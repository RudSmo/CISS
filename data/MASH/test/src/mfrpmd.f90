module mfrpmd

   ! Module for running Mean-field RPMD

   use types
   implicit double precision (a-h,o-z)

   integer :: nb ! number of beads
   integer :: nr ! number of replicas
   real(dp) :: hbar = 1.d0
   real(dp), allocatable :: betaN(:)
   real(dp), allocatable :: omegaN(:)

   real(dp), allocatable :: urp(:,:,:), vrp(:,:)
   character :: thermtyp ! thermostat type

contains

   subroutine init(beta,betamin,nb_,nr_,thermtyp_)
      use maths, only : symevp
      character :: thermtyp_
      ! integer, allocatable :: seed(:)
!
!     Initialize module variables
!
      nb = nb_
      nr = nr_
      thermtyp = thermtyp_
      allocate(betaN(nr),omegaN(nr))
      betaN(1) = beta/nb
      do m=2,nr
         betam = beta + (betamin-beta)*(m-1.0)/(nr-1.0)
         ! betam = beta * exp(log(betamin/beta)*(m-1.0)/(nr-1.0))
         betaN(m) = betam/nb
      end do
      print*, betaN
      omegaN = 1.d0/(betaN*hbar)
!
!     Initialize RP normal-mode transformation
!
      allocate(urp(nb,nb,nr),vrp(nb,nr))
      do m=1,nr
         rk = omegaN(m)**2
         urp = 0.d0
         do j = 1,nb
            urp(j,j,m) = 2*rk
         enddo
         do j = 2,nb
            urp(j-1,j,m) = -rk
            urp(j,j-1,m) = -rk
         enddo
         urp(1,nb,m) = urp(1,nb,m)-rk
         urp(nb,1,m) = urp(nb,1,m)-rk
         call symevp (urp(:,:,m),nb,nb,vrp(:,m),ierr)
         if (ierr .ne. 0) stop 'normal 2'
         vrp(1,m) = 0.d0! centroid mode has exactly zero frequency
      end do
      ! Initialize random seed
      ! call random_seed(size=nseed)
      ! allocate(seed(nseed))
      ! do n=1,nseed
      !    seed(n) = n
      ! end do
      ! call random_seed(put=seed)
      ! deallocate(seed)
   end subroutine


   subroutine runtrj(dt,nt,nth,Vt,Tt,Pt)
      use pes, only : nf, ns
      real(dp), intent(out) :: Vt(nt+1,nr), Tt(nt+1,nr), Pt(nt+1,ns)
!
!     Run thermostatted MF-RPMD trajectory, starting from GS equilibrium configuration.
!
      real(dp), allocatable :: q(:,:,:), p(:,:,:), force(:,:,:), Veff(:), pop(:)
      logical :: accepted
      allocate(q(nf,nb,nr),p(nf,nb,nr),force(nf,nb,nr),Veff(nr),pop(ns))
      q = 0.d0
      call samplep(p)
      call calcV(q,Veff,force,pop)
      call epot(Veff,q,Vt(1,:))
      call ekin(p,Tt(1,:))
      Pt(1,:) = pop(:)
      nacc = 0
      nrej = 0
      do it=1,nt
         if (thermtyp.eq.'a' .and. mod(it,nth).eq.0) call samplep(p) ! Andersen thermostat
         call step(p,q,Veff,force,pop,dt)
         if (nr.gt.1 .and. mod(it,10).eq.0) then
            call exchange(p,q,Veff,accepted)
            if (accepted) then
               nacc = nacc + 1
            else
               nrej = nrej + 1
            end if
         end if
         call epot(Veff,q,Vt(it+1,:))
         call ekin(p,Tt(it+1,:))
         Pt(it+1,:) = pop(:)
      end do
      if (nr.gt.1) print*, 'Exchange acceptance probability', float(nacc)/(nacc+nrej)
      deallocate(q,p,force,Veff,pop)
   end subroutine


   subroutine step(p,q,Veff,force,pop,dt)
      use pes, only : nf, mass, omega
      real(dp), intent(inout) :: p(:,:,:), q(:,:,:), force(:,:,:), Veff(:), pop(:)
! 
!     Evolve ring polymer through timestep dt
!
      dt2 = dt/2
      if (thermtyp.eq.'p') call pile(p,dt2)
      p = p + force * dt2
!
!     Cayley evolution of RP normal modes
!
      ! if (.not.cayley) stop 'mfrpmd step only implemented for cayley'
      if (nb.gt.1) then
         do m=1,nr
            do i=1,nf
               q(i,:,m) = matmul(q(i,:,m), Urp(:,:,m))
               p(i,:,m) = matmul(p(i,:,m), Urp(:,:,m))
            end do
         end do
      end if
      do m=1,nr
         ! q(:,1,m) = q(:,1,m) + p(:,1,m)/mass * dt
         do k=1,nb
            do i=1,nf
               em = mass(i)
               dtm = dt/em
               w2 = vrp(k,m) + omega(i)**2 ! don't include nuclear frequences
               x = w2*dt**2/4
               pp = (1-x)/(1+x)
               pq = -em*w2*dt/(1+x)
               qp = dtm/(1+x)
               qq = (1-x)/(1+x)
               pnew     = pq*q(i,k,m) + pp*p(i,k,m)  
               q(i,k,m) = qq*q(i,k,m) + qp*p(i,k,m)
               p(i,k,m) = pnew
            end do
         end do
      end do
      if (nb.gt.1) then
         do m=1,nr
            do i=1,nf
               q(i,:,m) = matmul(Urp(:,:,m), q(i,:,m))
               p(i,:,m) = matmul(Urp(:,:,m), p(i,:,m))
            end do
         end do
      end if
!
      call calcV(q,Veff,force,pop)
      p = p + force * dt2

      if (thermtyp.eq.'p') call pile(p,dt2)
   end subroutine

   subroutine calcV(q,Veff,force,pop)
      use pes, only : pot, grad, nf, ns
      use maths, only : trace,expm
      real(dp), intent(inout) :: q(:,:,:)
      real(dp), intent(inout) :: Veff(:)
      real(dp), intent(inout) :: force(:,:,:)
      real(dp), intent(inout) :: pop(:)
!
!     Calculate potential and force
!
      real(dp), allocatable :: V(:,:,:)    ! Diabatic potential for all beads
      real(dp), allocatable :: dVdq(:,:,:) ! Diabatic gradient for a single bead
      real(dp), allocatable :: M(:,:,:), H(:,:,:), D(:,:,:), Mk(:,:), Hk(:,:), Dki(:,:), Dk(:,:,:)
      allocate(V(ns,ns,nb),dVdq(nf,ns,ns),M(ns,ns,nb),H(ns,ns,nb),D(nf,ns,ns),Mk(ns,ns),Hk(ns,ns),Dki(ns,ns),Dk(nf,ns,ns))
      force = 0.d0
      do im=1,nr
         ttot = 0
         do k=1,nb
            call pot(q(:,k,im),V(:,:,k))
            ! Remove trace to avoid dividing zero by zero in force calculation.
            v0 = trace(V(:,:,k))/ns 
            do n=1,ns
               V(n,n,k) = V(n,n,k) - v0
            end do
            ttot = ttot + v0

            M(:,:,k) = -betaN(im)*V(:,:,k)
            call expm(M(:,:,k),ns)
         end do
         ! call calcM(V,M,im)
         call calcZ(M,Z)
         Veff(im) = -1.d0/betaN(im) * log(Z)  
         ! Keep trace when computing the potential
         Veff(im) = Veff(im) + ttot
         call calcH(M,H)
         ! call cpu_time(t1)
         do k=1,nb
            call grad(q(:,k,im),dVdq(:,:,:))
            ! Mk = M(:,:,k)
            ! Hk = H(:,:,k)
            ! ! do i=1,nf
            ! !    ! Trapezian rule of endpoints
            !    ! Dki = 0.5*(matmul(dVdq(i,:,:),Mk(:,:)) + matmul(Mk(:,:),dVdq(i,:,:)))
            ! !    force(i,k) = - trace(matmul(Dki,Hk))
            ! ! end do
            ! Dk = 0.d0
            ! do n3=1,ns
            !    do n2=1,ns
            !       do n1=1,ns
            !          Dk(:,n1,n3) = Dk(:,n1,n3) + dVdq(:,n1,n2)*Mk(n2,n3) + Mk(n1,n2)*dVdq(:,n2,n3)
            !          ! Dk(:,n1,n3) = Dk(:,n1,n3) + dVdq(:,n1,n2)*M12 + M23*dVdq(:,n2,n1)
            !          ! Dk(:,n1,n3) = Dk(:,n1,n3) + Mk(n1,n2)*dVdq(:,n2,n3)
            !          ! Dk(:,n1,n3) = Dk(:,n1,n3) + matmul(dVdq(:,n1,:),Mk(:,n3)) + matmul(dVdq(:,:,n3),Mk(n1,:))
            !       end do
            !    end do
            ! end do
            ! ! call cpu_time(tmid)
            ! Dk = Dk/2
            ! do n=1,ns
            !    do l=1,ns
            !       force(:,k) = force(:,k) - Dk(:,l,n)*Hk(n,l)
            !    end do
            ! end do

            ! OLD CODE - conserves energy for nb=1 and fairly OK for 8
            force(:,k,im) = 0.d0
            do i=1,ns
               force(:,k,im) = force(:,k,im) - matmul(dVdq(:,:,i),matmul(M(i,:,k),H(:,:,k)))
            end do
         end do
         ! call cpu_time(t2)
         force(:,:,im) = force(:,:,im) / Z
         if (im.eq.1) then
            ! Calculate populations for lowest replica
            call calcP(M,Z,pop)
         end if
      end do
      !!! TEST DERIVATIVE
      ! dq = 1.d-7
      ! q(1,1) = q(1,1) + dq
      ! do k=1,nb
      !    call pot(q(:,k),V(:,:,k))
      ! end do
      ! call calcM(V,M)
      ! call calcZ(M,Z)
      ! Veff2 = -1.d0/betaN * log(Z)
      ! print*, 'num', -(Veff2-Veff)/dq
      ! print*, 'ana', force(1,1)
      deallocate(V,dVdq,M,H,D,Mk,Hk,Dki,Dk)
   end subroutine

   subroutine exchange(p,q,Veff,accepted)
      use pes, only : nf
      real(dp), intent(inout) :: p(:,:,:), q(:,:,:), Veff(:)
      real(dp), allocatable :: qtmp(:,:), ptmp(:,:)
      logical :: accepted
      ! Randomly select pair of neighbouring temperatures
      call random_number(u)
      im = 1 + floor(u*(nr-1))
      jm = im + 1
      ! if (im.eq.1) then
      !    jm = 2
      ! else if (im.eq.nr) then
      !    jm = nr-1
      ! else
      !    call random_number(u)
      !    if (u.gt.0.5d0) then
      !       jm = im + 1
      !    else 
      !       jm = im - 1
      !    end if
      ! end if
      Delta = (betaN(im)-betaN(jm))*(Veff(im)-Veff(jm))
      call random_number(u)
      ! Accept with probability exp(-Delta), else reject. Note that if Delta<0, we always accept
      ! print*, (betaN(im)-betaN(jm)),(Veff(im)-Veff(jm))
      ! print*, it, im, jm, Delta, u.lt.exp(Delta)
      if (u.gt.exp(Delta)) then 
         accepted = .false.
         return 
      end if
      accepted = .true.
      allocate(qtmp(nf,nb),ptmp(nf,nb))
      qtmp = q(:,:,im)
      ptmp = p(:,:,im)
      q(:,:,im) = q(:,:,jm)
      p(:,:,im) = p(:,:,jm)
      q(:,:,jm) = qtmp
      p(:,:,jm) = ptmp
      deallocate(qtmp,ptmp)
      ! Rescale momenta (Sugita and Okamoto)
      fac = sqrt(betaN(im)/betaN(jm))
      p(:,:,jm) = p(:,:,jm)*fac
      p(:,:,im) = p(:,:,im)/fac
      ! Update Veff
      tmp = Veff(im)
      Veff(im) = Veff(jm)
      Veff(jm) = tmp
   end subroutine

   ! subroutine exchangewell(p,q,Veff,accepted)
   !    use pes, only : nf,ns
   !    use frexc, only : nf_bath
   !    real(dp), intent(inout) :: p(:,:,:), q(:,:,:), Veff(:)
   !    real(dp), allocatable :: qtmp(:,:), qtmp2(:,:)
   !    logical :: accepted
   !    ! Randomly select pair of states
   !    call random_number(u)
   !    im = 1 + floor(u*(ns-1))
   !    jm = im
   !    do while (jm.eq.im)
   !       call random_number(u)
   !       jm = 1 + floor(u*(ns-1))
   !    end do
   !    ! Swap states
   !    allocate(qtmp(nf_bath,ns),qtmp2(nf_bath,ns))
   !    qtmp = reshape(q(:,1,1),(/nf_bath,ns/))
   !    qtmp2 = qtmp
   !    qtmp2(:,im) =
      
   !    ttot = 0
   !    do k=1,nb
   !       call pot(q(:,k,im),V(:,:,k))
   !       ! Remove trace to avoid dividing zero by zero in force calculation.
   !       v0 = trace(V(:,:,k))/ns 
   !       do n=1,ns
   !          V(n,n,k) = V(n,n,k) - v0
   !       end do
   !       ttot = ttot + v0

   !       M(:,:,k) = -betaN(im)*V(:,:,k)
   !       call expm(M(:,:,k),ns)
   !    end do
   !    ! call calcM(V,M,im)
   !    call calcZ(M,Z)
   !    Veff(im) = -1.d0/betaN(im) * log(Z)  
   !    ! Keep trace when computing the potential
   !    Veff(im) = Veff(im) + ttot
   !    Delta = betaN(1)*(Veff(im)-Veff(jm))
   !    call random_number(u)
   !    ! Accept with probability exp(-Delta), else reject. Note that if Delta<0, we always accept
   !    ! print*, (betaN(im)-betaN(jm)),(Veff(im)-Veff(jm))
   !    ! print*, it, im, jm, Delta, u.lt.exp(Delta)
   !    if (u.gt.exp(Delta)) then 
   !       accepted = .false.
   !       return 
   !    end if
   !    accepted = .true.
   !    allocate(qtmp(nf,nb),ptmp(nf,nb))
   !    qtmp = q(:,:,im)
   !    ptmp = p(:,:,im)
   !    q(:,:,im) = q(:,:,jm)
   !    p(:,:,im) = p(:,:,jm)
   !    q(:,:,jm) = qtmp
   !    p(:,:,jm) = ptmp
   !    deallocate(qtmp,ptmp)
   !    ! Rescale momenta (Sugita and Okamoto)
   !    fac = sqrt(betaN(im)/betaN(jm))
   !    p(:,:,jm) = p(:,:,jm)*fac
   !    p(:,:,im) = p(:,:,im)/fac
   !    ! Update Veff
   !    tmp = Veff(im)
   !    Veff(im) = Veff(jm)
   !    Veff(jm) = tmp
   ! end subroutine

   ! subroutine calcM(V,M,im)
   !    use pes, only : ns
   !    use maths, only : expm
   !    real(dp), intent(in) :: V(:,:,:)
   !    real(dp), intent(inout) :: M(:,:,:)
   !    do k=1,nb
   !       M(:,:,k) = -betaN(im)*V(:,:,k)
   !       call expm(M(:,:,k),ns)
   !    end do
   ! end subroutine

   subroutine calcZ(M,Z)
      use pes, only : ns
      use maths, only : trace
      real(dp), intent(in) :: M(:,:,:)
!
!     Calculate local partition function for a given M(q)
!
      real(dp), allocatable :: work(:,:)
      allocate(work(ns,ns))
      work = M(:,:,1)
      do k=2,nb
         work = matmul(work,M(:,:,k))
      end do
      ! Take trace
      Z = trace(work)
      ! Z = 0.d0
      ! do n=1,ns 
      !    Z = Z + work(n,n)
      ! end do
      deallocate(work)
   end subroutine

   subroutine calcH(M,H)
      use pes, only : ns
      real(dp), intent(in) :: M(:,:,:)
      real(dp), intent(inout) :: H(:,:,:)
!
!     Use Bell's algorithm to calculate the 'hole' matrix H_k from
!     H_k = G_{k+1}F_{k-1}, which is built on the notation
!     G_k = M_k ... M_nb (right product) and
!     F_k = M_1 ... M_k  (left product) 
!
      real(dp), allocatable :: F(:,:,:), G(:,:,:)
      H = 0.d0
      if (nb.eq.1) then
         do n=1,ns
            H(n,n,1) = 1.d0
         end do
         return
      end if
      allocate(F(ns,ns,nb),G(ns,ns,nb))
      F(:,:,1) = M(:,:,1)
      G(:,:,nb) = M(:,:,nb)
      do k=2,nb-1
         F(:,:,k) = matmul(F(:,:,k-1),M(:,:,k))
         G(:,:,nb+1-k) = matmul(M(:,:,nb+1-k),G(:,:,nb+2-k))
      end do
      H(:,:,1) = G(:,:,2)
      H(:,:,nb) = F(:,:,nb-1)
      do k=2,nb-1
         H(:,:,k) = matmul(G(:,:,k+1),F(:,:,k-1))
      end do
      deallocate(F,G)
   end subroutine


   subroutine calcP(M,Z,pop)
      use pes, only : ns
      use maths, only : trace
      real(dp), intent(in) :: M(:,:,:)
      real(dp), intent(inout) :: pop(:)
! 
!     Calculate populations
!
      real(dp), allocatable :: work(:,:)
      allocate(work(ns,ns))
      ! pop = 0.d0
      do n=1,ns
         work = 0.d0
         work(n,n) = 1.d0
         do k=1,nb
            work = matmul(work,M(:,:,k))
         end do
         pop(n) = trace(work)
      end do
      pop = pop / Z
      deallocate(work)
   end subroutine

   subroutine epot(Veff,q,Ep)
      use pes, only : mass
      real(dp), intent(in) :: Veff(:), q(:,:,:)
      real(dp), intent(out) :: Ep(nr)
!
!     Potential energy of each replica
!      
      do m=1,nr
         Ep(m) = Veff(m)
         do k=1,nb
            next = k+1
            if (next.gt.nb) next = 1
            Ep(m) = Ep(m) + 0.5d0*omegaN(m)**2*sum(mass*(q(:,k,m)-q(:,next,m))**2)
         end do
      end do
   end subroutine

   subroutine ekin(p,T)
      use pes, only : mass
      real(dp), intent(in) :: p(:,:,:)
      real(dp), intent(out) :: T(nr)

!
!     Kinetic energy of each replica
!      
      T = 0.d0
      do m=1,nr
         do k=1,nb
            T(m) = T(m) + 0.5d0*sum(p(:,k,m)**2/mass)
         end do
      end do
   end subroutine


   subroutine samplep(p)
      use pes, only : nf, mass
      use maths, only : gasdev
      real(dp), intent(inout) :: p(:,:,:)
!
!     Resample momenta from thermal equilibrium (Anderson thermostat)
!
      real(dp), allocatable :: xi(:)
      allocate(xi(nf))
      do m=1,nr
         do k=1,nb
            call gasdev(xi,nf)
            p(:,k,m) = sqrt(mass/betaN(m)) * xi
         end do
      end do
      deallocate(xi)
   end subroutine

   subroutine pile(p,dt)
      use pes, only : nf, mass, omega
      use maths, only : gasdev
      real(dp), intent(inout) :: p(:,:,:)
!
!     Friction step for a PILE thermostat
!
      real(dp), allocatable :: xi(:)
      allocate(xi(nf))
!
      do m=1,nr
         if (nb.gt.1) then
            do i=1,nf
               p(i,:,m) = matmul(p(i,:,m), Urp(:,:,m))
            end do
         end if
   !     

         do k=1,nb
            call gasdev(xi,nf)
            do i=1,nf
               sigp = sqrt(mass(i)/betaN(m))
               gamma = sqrt(vrp(k,m)+omega(i)**2)
               c1 = exp(-dt*gamma)
               c2 = sqrt(1.d0-c1**2)
               p(i,k,m) = c1*p(i,k,m) + sigp*c2*xi(i)
            end do
         end do

   !
         if (nb.gt.1) then
            do i=1,nf
               p(i,:,m) = matmul(Urp(:,:,m), p(i,:,m))
            end do
         end if
      end do
      deallocate(xi)
   end subroutine

end module