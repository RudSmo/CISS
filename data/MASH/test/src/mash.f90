module mash

   ! Module for running mapping approach to surface hopping

   use types
   implicit double precision (a-h,o-z)

   real(dp) :: alpha,an  ! ns-dependent constants used in observables
   real(dp) :: beta      ! Reciprocal temperature
   logical :: useqtb = .false. ! Toggle quantum thermal bath


contains

! =============== Initialization =============
   subroutine init(beta_,useqtb_)
      use pes, only : ns
      real(dp) :: beta_
      logical :: useqtb_
!
!    Initialize module
!
      Hn = 0.d0
      do n=1,ns
         Hn = Hn + 1.d0/n
      end do
      alpha = (ns-1.d0)/(Hn-1.d0)
      an = (ns+1.d0)/alpha

      beta = beta_
      useqtb = useqtb_
   end subroutine


! ============= Potential, Hamiltonian =============

   subroutine mash_pot(q, qe, pe, vtot)
      use pes, only : nf, ns, potad
      real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
      real(dp), intent(out) :: vtot
!
!     Compute MASH potential at (q,qe,pe)
!
      integer :: a
      real(dp), allocatable :: U(:,:), Vad(:)
      allocate(U(ns,ns),Vad(ns))
      call potad(q,Vad,U)
      call cstate(qe,pe,U,a)
      vtot = Vad(a)
      deallocate(U,Vad)
   end subroutine


   real(dp) function ham(q, p, qe, pe)
      use pes, only : nf, ns, mass
      real(dp), intent(in) :: q(nf), p(nf), qe(ns), pe(ns)
!
!     Function to compute Hamiltonian at a phase-space point
!
      call mash_pot(q, qe, pe, vtot)
      ham = vtot + 0.5d0*sum(p**2/mass)
   end function

   
   real(dp) function ham_a(q, p, a)
      use pes, only : ns, mass, potad
      real(dp), intent(in) :: q(:), p(:)
      integer :: a
!
!     Function to compute Hamiltonian at a phase-space point
!
      real(dp), allocatable :: U(:,:), Vad(:)
      allocate(U(ns,ns),Vad(ns))
      call potad(q,Vad,U)
      ham_a = Vad(a) + 0.5d0*sum(p**2/mass)
      deallocate(U,Vad)
   end function

! =============== Observable-related subroutines =============

   subroutine cstate(qe, pe, U, a)
      use pes, only : ns
      real(dp), intent(in) :: qe(:), pe(:), U(:,:)
      integer, intent(out) :: a
!
!     Get current adiabatic state with the highest population |c_a|^2
!
      real(dp), allocatable :: qa(:), pa(:)
      complex(dpc), allocatable :: ca(:)
      allocate(qa(ns),pa(ns),ca(ns))
      ! Convert to adiabatic rep.
      qa = matmul(qe,U)
      pa = matmul(pe,U)
      ca = dcmplx(qa,pa)
      call cstate_ad(ca,a)
      deallocate(qa,pa,ca)
   end subroutine

   subroutine cstate_ad(ca, a)
      use pes, only : ns
      complex(dpc), intent(in) :: ca(:)
      integer, intent(out) :: a
!
!     Get current adiabatic state with the highest population |c_a|^2
!
      integer :: b
      real(dp), allocatable :: pop(:)
      allocate(pop(ns))
      ! Convert to adiabatic rep.
      pop = abs(ca)**2
      a = 1
      do b = 2,ns
         if (pop(b)>pop(a)) a=b
      end do
      deallocate(pop)
   end subroutine

   subroutine cstate2_ad(ca, a, b)
      use pes, only : ns
      complex(dpc), intent(in) :: ca(:)
      integer :: a, b
! 
!     Get state b with second-highest population excluding a
!
      real(dp), allocatable :: pop(:)
      allocate(pop(ns))
      pop = abs(ca)**2
      pop(a) = 0.d0
      b = 1
      do n = 2,ns
         if (pop(n)>pop(b)) b=n
      end do
      deallocate(pop)
   end subroutine

   subroutine obsbls_mash(c, obs, typ, poponly)
      use pes, only : ns
      ! use maths, only : iu
      complex(dpc), intent(inout) :: c(:)
      complex(dpc), intent(inout) :: obs(:,:)
      integer :: typ
      logical :: poponly
!
!  Calculate binned observables in given representation. 
!     obs(k,k) -- population
!     obs(k,l) -- coherence 
!
!     typ==1 -- binned (Theta function)
!     typ==2 -- weighted in the Phin way
!     typ==3 -- weighted in the (P-Pnext) way
!     typ==4 -- weighted in the rhon way
!
      ! complex(dpc), allocatable :: c2(:)
      overn = 1.d0/ns
      call cstate_ad(c,l) !highest populated state
      obs = 0.d0
      if (typ.eq.1) then
         obs(l,l) = 1.d0
      else if (typ.eq.2) then
         do k=1,ns
            obs(k,k) = overn + alpha*(abs(c(k))**2 - overn)
         end do
      else if (typ.eq.3) then
         call cstate2_ad(c,l,m)
         obs(l,l) = ns*(abs(c(l)**2)-abs(c(m)**2))
      else if (typ.eq.4) then
         do k=1,ns
            obs(k,k) = overn + an*(abs(c(k))**2 - overn)
         end do
      end if
      if (poponly) then
         return
      end if

      ! Coherences
      ort2 = 1.d0/sqrt(2.d0)
      do k=1,ns
         do l=k+1,ns
            if (typ.eq.1 .or. typ.eq.3) then
               ! TEST
               obs(k,l) = conjg(c(k))*c(l) * alpha
               obs(l,k) = conjg(obs(k,l))
               cycle
               ! END TEST
               ! allocate(c2(ns))
               ! c2 = c
               ! ! Sigma_x
               ! c2(l) = ort2*(c(l) + c(k))
               ! c2(k) = ort2*(c(l) - c(k))
               ! call cstate_ad(c2,l2)
               ! cohx = 0.d0
               ! w = 1.d0
               ! if (typ.eq.3) then
               !    call cstate2_ad(c2,l2,m2)
               !    w = ns*(abs(c2(l2)**2) - abs(c2(m2)**2))
               ! end if
               ! if (l2.eq.l) cohx = w
               ! if (l2.eq.k) cohx = -w
               ! ! Sigma_y
               ! c2(l) = ort2*(c(l) + iu*c(k))
               ! c2(k) = ort2*(c(l) - iu*c(k))
               ! call cstate_ad(c2,l2)
               ! cohy = 0.d0
               ! if (typ.eq.3) then
               !    call cstate2_ad(c2,l2,m2)
               !    w = ns*(abs(c2(l2)**2) - abs(c2(m2)**2))
               ! end if
               ! if (l2.eq.l) cohy = w
               ! if (l2.eq.k) cohy = -w
               ! ! Combine
               ! obs(k,l) = 0.5d0*(cohx + iu*cohy)
               ! obs(l,k) = 0.5d0*(cohx - iu*cohy)
               ! deallocate(c2)
            else if (typ.eq.2) then
               obs(k,l) = alpha*conjg(c(k))*c(l)
               obs(l,k) = conjg(obs(k,l))
            else if (typ.eq.4) then
               obs(k,l) = an*conjg(c(k))*c(l)
               obs(l,k) = conjg(obs(k,l))
            end if
         end do
      end do
   end subroutine

   subroutine obsbls(q, qe, pe, obs, rep, typ, poponly)
      use pes, only : ns, potad
      real(dp), intent(in) :: q(:), qe(:), pe(:)
      complex(dpc), intent(out) :: obs(:,:)
      character, intent(in) :: rep ! representation ('d' or 'a')
      integer, intent(in) :: typ ! projection type theta_n (1) or Phi_n (2)
      logical, intent(in) :: poponly ! Only compute populations
!
!  Observables in diabatic or adiabatic representation
!
      real(dp), allocatable :: Vad(:), U(:,:)
      complex(dpc), allocatable :: c(:), obsa(:,:)
      integer :: a
      allocate(c(ns))
      if (rep.eq.'a') then
         ! Adiabatic observables
         allocate(Vad(ns),U(ns,ns))
         call potad(q,Vad,U)
         c = dcmplx(matmul(qe,U),matmul(pe,U))
         call obsbls_mash(c, obs, typ, poponly)
         deallocate(Vad,U)
      else if (rep.eq.'d') then
         ! Diabatic (site) observables
         c = dcmplx(qe,pe)
         call obsbls_mash(c, obs, typ, poponly)
      else if (rep.eq.'e') then
         ! Diabatic (exciton) observables
         allocate(Vad(ns),U(ns,ns))
         call potad(q*0.d0,Vad,U)
         c = dcmplx(matmul(qe,U),matmul(pe,U))
         call obsbls_mash(c, obs, typ, poponly)
         deallocate(Vad,U)
      else if (rep.eq.'r') then
         ! Measure diabatic observables by rotating to an adiabatic basis
         allocate(Vad(ns),U(ns,ns),obsa(ns,ns))
         call potad(q,Vad,U)
         c = dcmplx(matmul(qe,U),matmul(pe,U))
         call obsbls_mash(c, obsa, typ, .false.)
         obs = 0.d0
         do n=1,ns
            do m=1,ns
               ! O_dia = U.O_ad.U^T
               do a=1,ns
                  obs(n,m) = obs(n,m) + sum(U(n,:)*obsa(:,a)*U(m,a))
               end do
            end do
         end do
         deallocate(Vad,U,obsa)
      else
         stop 'obsbls: Undefined representation'
      end if
      deallocate(c)
   end subroutine

   subroutine pops(q, qe, pe, pop, rep, typ)
      use pes, only : ns
      real(dp), intent(in) :: q(:), qe(:), pe(:)
      real(dp), intent(inout) :: pop(:)
      character, intent(in) :: rep
      integer, intent(in) :: typ
!
!  All populations |n><n|
!
      complex(dpc), allocatable :: obs(:,:)
      allocate(obs(ns,ns))
      call obsbls(q, qe, pe, obs, rep, typ, .true.)
      do n=1,ns
         pop(n) = real(obs(n,n))
      end do
      deallocate(obs)
   end subroutine

   subroutine obs_ad(c,obs,weight)
      use pes, only : ns
      complex(dpc), intent(in)  :: c(:)
      complex(dpc), intent(out) :: obs(:,:)
   !
   !  Observables in the adiabatic representation
   !
      integer :: a,aprime
      call cstate_ad(c,a) ! find index of highest populated state
      call cstate2_ad(c,a,aprime) ! find index of second highest populated state
      ! Compute populations with Theta (storing weight separately)
      obs = 0.d0
      obs(a,a) = 1.d0
      weight = ns*(abs(c(a)**2)-abs(c(aprime)**2))
      ! print*, 'weight', weight, abs(c(a)**2), abs(c(aprime)**2)
      ! Compute coherences (without scalings)
      do n=1,ns
         do m=1,ns
            if (n.eq.m) cycle
            obs(n,m) = conjg(c(n))*c(m)
         end do
      end do
   end subroutine
   
   subroutine corrfunc(qt, qet, pet, CAB, rep, init, nt)
      use pes, only : ns, potad
      real(dp), intent(in)  :: qt(:,:), qet(:,:), pet(:,:)
      ! complex(dpc), intent(out) :: CAB(nt+1,ns,ns,ns,ns)
      complex(dpc), intent(out) :: CAB(nt+1,ns)
      character :: rep
   !
   !  Correlation functions by going via the adiabatic representation
   !
      real(dp), allocatable :: Vad(:), U0(:,:), U(:,:), Uexc(:,:)
      complex(dpc), allocatable :: c(:), obs0(:,:), obst(:,:), CABad(:,:,:,:)
      ! complex(dpc) :: tmp
      allocate(Vad(ns),U0(ns,ns),U(ns,ns),Uexc(ns,ns))
      allocate(c(ns),obs0(ns,ns),obst(ns,ns),CABad(ns,ns,ns,ns))
      CAB = 0.d0
      ! Exciton transformation
      call potad(0.d0*qt(1,:),Vad,Uexc)
      ! Zero-time adiabatic observables
      call potad(qt(1,:),Vad,U0)
      c = dcmplx(matmul(qet(1,:),U0),matmul(pet(1,:),U0))
      call obs_ad(c,obs0,weight)
      if (rep.eq.'e') U0 = matmul(transpose(Uexc),U0)
      ! print*, 'obs_ad0',obs0
      do it=1,nt+1
         ! Time-dependent adiabatic observables
         call cpu_time(t1)
         call potad(qt(it,:),Vad,U)
         c = dcmplx(matmul(qet(it,:),U),matmul(pet(it,:),U))
         if (rep.eq.'e') U = matmul(transpose(Uexc),U)
         call obs_ad(c,obst,wtmp)
         ! print*, 'obs_adt',obst
         ! Put the adiabatic pieces together
         call cpu_time(t2)
         CABad = 0.d0
         do n=1,ns
            do m=1,ns
               do k=1,ns
                  do l=1,ns
                     if (n.eq.m .and. k.eq.l) then ! PP
                        CABad(n,n,k,k) = weight * obs0(n,n)*obst(k,k)
                     else if (n.eq.m .and. k.ne.l) then ! PC
                        CABad(n,n,k,l) = alpha  * obs0(n,n)*obst(k,l)
                     else if (n.ne.m .and. k.eq.l) then ! CP
                        CABad(n,m,k,k) = alpha  * obs0(n,m)*obst(k,k)
                     else ! CC
                        CABad(n,m,k,l) =(ns+1.d0)*obs0(n,m)*obst(k,l)
                     end if
                  end do
               end do
            end do
         end do
         call cpu_time(t3)
         ! print*, 'C11'
         ! print*, real(CABad(1,1,1,1)), real(CABad(1,1,1,2)), real(CABad(1,2,2,1))
         if (rep.eq.'a') then
            ! CAB(it,:,:,:,:) = CABad
            do k=1,ns
               CAB(it,k) = CABad(init,init,k,k)
            end do
         else
            ! Transform back to diabatic representation
            ! TODO: multiply U by exciton transform if requested
            ! do l=1,ns
            !    do k=1,ns
            !       do m=1,ns
            !          do n=1,ns
            !             tmp = 0.d0
            !             do ll=1,ns
            !                do kk=1,ns
            !                   do mm=1,ns
            !                      do nn=1,ns
            !                         tmp = tmp + U0(n,nn)*U0(m,mm)*U(k,kk)*U(l,ll)*CABad(nn,mm,kk,ll)
            !                      end do
            !                   end do
            !                end do
            !             end do
            !             CAB(it,n,m,k,l) = tmp
            !          end do
            !       end do
            !    end do
            ! end do

            do nn=1,ns
               do mm=1,ns
                  do k=1,ns
                     CAB(it,k) = CAB(it,k) + U(k,nn)*U(k,mm)* &
                        dot_product(U0(init,:),matmul(CABad(:,:,nn,mm),U0(init,:)))
                        ! dot_product(U(k,nn),matmul(CABad(:,:,nn,mm),U(k,mm)))
                  end do
               end do
            end do
            ! tmp = CAB(it,:,:)

   
            ! do n=1,ns
            !    do m=1,ns
            !       do nn=1,ns
            !          do mm=1,ns
            !             CAB(it,n,m,:,:) = CAB(it,n,m,:,:) + U0(n,nn)*U0(m,mm)*CABad(nn,mm,:,:)                     
            !          end do
            !       end do
            !    end do
            ! end do
            ! ! call cpu_time(t4)
            ! CABad = CAB(it,:,:,:,:) ! Store temporary result
            ! CAB(it,:,:,:,:) = 0.d0
            ! do n=1,ns
            !    do m=1,ns
            !       do nn=1,ns
            !          do mm=1,ns
            !             CAB(it,:,:,n,m) = CAB(it,:,:,n,m) + U(n,nn)*U(m,mm)*CABad(:,:,nn,mm)
            !          end do
            !       end do
            !    end do
            ! end do
         end if
         call cpu_time(t4)
         ! print*, t2-t1,t3-t2,t4-t3
         ! sum = 0.d0
         ! do n=1,ns
         !    sum = sum + real(CAB(1,1,n,n))
         ! end do
         ! print*, it, sum
         end do
      deallocate(Vad,U0,U,Uexc,c,obs0,obst,CABad)
   end subroutine



! =============== Dynamics-related subroutines ===============
   subroutine evolve(q, p, qe, pe, Vad, U, dvdq, a, dtbase, ierr)
      use pes, only : nf, ns, grad_a
      real(dp), intent(inout) :: q(:), p(:), qe(:), pe(:), Vad(:), &
                                 U(:,:), dvdq(:)
      integer, intent(inout) :: a, ierr
      real(dp), intent(in) :: dtbase
      real(dp), allocatable :: q0(:),p0(:),qe0(:),pe0(:), &
                               Vad0(:), U0(:,:), dvdq0(:)
      complex(dpc), allocatable :: ca0(:), ca1(:)
      real(dp), allocatable :: da(:,:), pa(:), pop(:)
      integer :: b
      logical :: accepted
!
!     Perform a time step
!
      allocate(q0(nf),p0(nf),qe0(ns),pe0(ns),Vad0(ns),U0(ns,ns),dvdq0(nf))
      allocate(ca0(ns),ca1(ns))
      allocate(da(nf,ns),pa(ns),pop(ns))

      if (useqtb) call step_T(p,dtbase/2)

      ierr = 0
      maxhop = 10
      dt = dtbase
      do ihop=1,maxhop ! Limit number of hops to look for in a timestep dt
         ! Store initial values
         call savetmp(q,p,qe,pe,Vad,U,dvdq,q0,p0,qe0,pe0,Vad0,U0,dvdq0)
      
         ! Store initial adiabatic wavefunction
         ca0 = dcmplx(matmul(qe,U),matmul(pe,U))
         
         ! Attempt full step
         call verlet(q,p,qe,pe,Vad,U,dvdq,a,dt)

         ! Calculate new active state
         ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
         call cstate_ad(ca1,b)

         if (b.eq.a) then
            ! Stayed on state - we're done
            exit
         else
            ! States have changed - find crossing time with bisection root search
            tl = 0.0
            tr = dt
            fl = deltaP(ca0,a)
            fr = deltaP(ca1,a)
            do iter=1,10
               tm = (tl+tr)/2 ! Mid point
               call savetmp(q0,p0,qe0,pe0,Vad0,U0,dvdq0,q,p,qe,pe,Vad,U,dvdq)
               call verlet(q,p,qe,pe,Vad,U,dvdq,a,tm)
               ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
               fm = deltaP(ca1,a)
               if (fm.gt.0) then
                  tl = tm
               else
                  tr = tm
               end if
               ! print*, iter, tm, dt, fm
            end do
            call cstate2_ad(ca1,a,b)
            call cross(q, p, ca1, a, b, Vad, U, accepted)
            if (accepted) then
               tx = tr
            else
               tx = tl
            end if
            call savetmp(q0,p0,qe0,pe0,Vad0,U0,dvdq0,q,p,qe,pe,Vad,U,dvdq)
            call verlet(q,p,qe,pe,Vad,U,dvdq,a,tx)
            ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
            call cross(q, p, ca1, a, b, Vad, U, accepted)
            if (accepted) then
               call grad_a(q,U,b,dvdq)
               a = b
            end if
            dt = dt - tx
         end if
         if (ihop.eq.maxhop) then
            ! Too many hops -- just finish it off
            ! ierr = 1
            call verlet(q,p,qe,pe,Vad,U,dvdq,a,dt)
            ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
            call cstate_ad(ca1,b)
            if (b.ne.a) then
               call cross(q, p, ca1, a, b, Vad, U, accepted)
               if (accepted) then
                  call grad_a(q,U,b,dvdq)
                  a = b
               end if
            end if                 
         end if
      end do

      if (useqtb) call step_T(p,dtbase/2)

      deallocate(da,pop,pa)
      deallocate(ca0,ca1)
      deallocate(q0,p0,qe0,pe0,Vad0,U0)
   end subroutine

   real(dp) function deltaP(ca,a)
      complex(dpc) :: ca(:)
      integer :: a,b
!
!     Difference in population between state a and highest populated state apart from a
!
      call cstate2_ad(ca,a,b)
      deltaP = abs(ca(a))**2 - abs(ca(b))**2
   end function

   subroutine savetmp(q,p,qe,pe,Vad,U,dvdq,q0,p0,qe0,pe0,Vad0,U0,dvdq0)
      real(dp), intent(inout) :: q(:),p(:),qe(:),pe(:),Vad(:),U(:,:),dvdq(:),&
            q0(:),p0(:),qe0(:),pe0(:), Vad0(:),U0(:,:),dvdq0(:)   
!
!     Store temporary variables and potential information
!
      q0 = q
      p0 = p
      qe0 = qe
      pe0 = pe
      Vad0 = Vad
      U0 = U
      dvdq0 = dvdq
   end subroutine


   subroutine verlet(q, p, qe, pe, Vad, U, dvdq, a, dt)
      use pes, only : potad, grad_a
      real(dp), intent(inout) :: q(:), p(:), qe(:), pe(:), Vad(:), &
                                 U(:,:), dvdq(:)
      integer, intent(in) :: a
      real(dp), intent(in) :: dt
!
!     Perform a time step on state n
!
      dt2 = dt/2
      call step_e(qe,pe,Vad,U,dt2)
      call step_p(p,dvdq,dt2)
      call step_q(q,p,dt)
      call potad(q, Vad, U)
      call grad_a(q,U,a,dvdq)
      call step_p(p,dvdq,dt2)
      call step_e(qe,pe,Vad,U,dt2)
   end subroutine


   subroutine step_p(p,dvdq,dt)
      real(dp), intent(inout) :: p(:)
      real(dp), intent(in) :: dvdq(:), dt
!
!    Evolve nuclear momenta
!
      p = p - dvdq * dt
   end subroutine

!    subroutine step_q(q, p, dt)
!       use pes, only : mass
!       real(dp), intent(inout) :: q(:), p(:)
!       real(dp), intent(in) :: dt
! !
! !     Evolve q and check for switch of adiabatic potential
! !
!       q = q + dt*p/mass
!    end subroutine

   subroutine step_q(q, p, dt)
      use pes, only : nf,mass,omega,cayley
      real(dp), intent(inout) :: q(nf), p(nf)
      real(dp), intent(in) :: dt
!
!     Cayley evolution through a time interval dt.
!
      if (cayley) then
         do i = 1,nf
            em = mass(i)
            dtm = dt/em
            w2 = omega(i)**2
            x = w2*dt**2/4
            pp = (1-x)/(1+x)
            pq = -em*w2*dt/(1+x)
            qp = dtm/(1+x)
            qq = (1-x)/(1+x)
            
            pnew = pq*q(i) + pp*p(i)  
            q(i) = qq*q(i) + qp*p(i)
            p(i) = pnew
         end do
      else
         q = q + dt*p/mass
      end if
   end subroutine

   subroutine step_e(qe, pe, Vad, U, dt)
      use pes, only : ns
      use maths, only : iu, symevp
      real(dp), intent(inout) :: qe(:), pe(:)
      real(dp), intent(in) :: Vad(:), U(:,:)
      real(dp), intent(in) :: dt
!
!     Evolve electronic coefficients c = qe + i*pe
!
      complex(dpc), allocatable :: c(:)
      allocate(c(ns))
      ! U transforms from adiabatic to diabatic basis
      qe = matmul(qe, U) ! Same as qe = transpose(U).qe
      pe = matmul(pe, U)
      c = dcmplx(qe, pe)
      c = exp(- iu * dt * Vad) * c
      qe = matmul(U, real(c)) ! qe = U.qe
      pe = matmul(U, aimag(c))
      deallocate(c)
   end subroutine

   subroutine step_T(p, dt)
      use pes, only : nf, mass, omega
      use maths, only : gasdev
      real(dp), intent(inout) :: p(:)
      real(dp), intent(in) :: dt
!
!     Friction step for a Langevin thermostat with strength gamma
!
      real(dp), allocatable :: xi(:)
      real(dp) :: sigp, betaeff
      allocate(xi(nf))
      call gasdev(xi,nf)
      do i=1,nf
         betaeff = tanh(beta*omega(i)/2)/(omega(i)/2.d0)
         sigp = sqrt(mass(i)/betaeff)
         gamma = 0.1*omega(i)
         c1 = exp(-dt*gamma)
         c2 = sqrt(1.d0-c1**2)
         p(i) = c1*p(i) + sigp*c2*xi(i)
      end do
      deallocate(xi)
   end subroutine

   subroutine cross(q,p,cad,n,m,Vad,U,accepted)
      use pes, only : nf, ns, mass, nacdir, omega
      real(dp), intent(inout) :: p(:)
      complex(dpc), intent(inout) :: cad(:)
      real(dp), intent(in) :: q(:), Vad(:), U(:,:)
      logical, intent(out) :: accepted
!
!     If energetically allowed, hop from n to m and rescale momentum
!     If not, reverse propagation direction along given component of the NAC vector
!
      real(dp), allocatable :: d(:,:,:), dj(:), pnac(:), porth(:)
      allocate(d(nf,ns,ns),dj(nf),pnac(nf),porth(nf))
      call nacdir(q,cad,Vad,U,n,m,dj)
      dj = dj/sqrt(mass)
      
      ! Use mass-scaled momenta
      p = p/sqrt(mass)
      if (nf.eq.1) then
         pnac = p
      else
         pnac = dot_product(p,dj)/dot_product(dj,dj)*dj
      end if
      porth = p - pnac
      ! TEST: rescale "classically adjusted" momenta
      if (useqtb) then
         do i=1,nf
            betaeff = tanh(beta*omega(i)/2)/(omega(i)/2.d0)
            pnac(i) = pnac(i) * sqrt(betaeff/beta)
         end do
      end if
      ! END TEST
      Ekin = 0.5d0*sum(pnac**2)
      Vdiff = Vad(m)-Vad(n)
      if ((Ekin-Vdiff).gt.0.d0) then
         ! Rescale momentum 
         pnac = sqrt(2.d0*(Ekin-Vdiff)) & ! (no masses since p is mass-scaled)
             * pnac/sqrt(dot_product(pnac,pnac))
         accepted = .true.
      else
         ! Reverse momentum along NAC vector
         pnac = -pnac
         accepted = .false.
      end if      
      ! TEST: scale back "classically adjusted" momenta
      if (useqtb) then
         do i=1,nf
            betaeff = tanh(beta*omega(i)/2)/(omega(i)/2.d0)
            pnac(i) = pnac(i) * sqrt(beta/betaeff)
         end do
      end if
      ! END TEST
      p = porth + pnac
      p = p*sqrt(mass)
      deallocate(d,dj,pnac,porth)
   end subroutine

   subroutine store(q, p, qe, pe, a, it, qt, pt, qet, pet, at)
      use types
      real(dp), intent(in) :: q(:), p(:), qe(:), pe(:)
      integer, intent(in) :: a, it
      real(dp), intent(inout) :: qt(:,:), pt(:,:), qet(:,:), pet(:,:)
      integer, intent(inout) :: at(:)
!
!  Store state at a given timestep
!
      qt(it,:) = q
      pt(it,:) = p
      qet(it,:) = qe
      pet(it,:) = pe
      at(it) = a
   end subroutine 


end module
