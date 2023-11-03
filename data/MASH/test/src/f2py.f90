!f2py

! =========== Potential-specific initializations =============
subroutine init_linvib(mass,omega,Vconst,Vlin,nf,ns)
   use types
   use linvib, only : init
   integer  :: nf,ns
   real(dp) :: mass(nf), omega(nf)
   real(dp) :: Vconst(ns,ns), Vlin(nf,ns,ns)
!
!  Initialize linear vibronic potential
!
   call init(nf,ns,mass,omega,Vconst,Vlin)
end subroutine

subroutine init_bilinvib(mass,omega,Vconst,Vlin,Vbilin,nf,ns)
   use types
   use bilinvib, only : init
   integer  :: nf,ns
   real(dp) :: mass(nf), omega(nf)
   real(dp) :: Vconst(ns,ns), Vlin(nf,ns,ns), Vbilin(nf,nf,ns,ns)
!
!  Initialize bilinear vibronic potential
!
   call init(nf,ns,mass,omega,Vconst,Vlin,Vbilin)
end subroutine

subroutine init_frexc(mass,omega,Vconst,kappa,nf,nf_bath,ns)
   use types
   use frexc, only : init
   integer  :: nf,ns
   real(dp) :: mass(nf), omega(nf)
   real(dp) :: Vconst(ns,ns), kappa(nf_bath)
!
!  Initialize Frenkel-exciton potential
!
   call init(nf,ns,mass,omega,Vconst,kappa)
end subroutine

subroutine init_tully(model,mass)
   use types
   use tully, only : init
   integer :: model
   real(dp) :: mass(1)
!
!  Initialize Tully model
!
   call init(model,mass)
end subroutine

! ================ General initializations =====================

subroutine init_mash(beta,useqtb)
   use types
   use mash, only : init
   real(dp) :: beta
   logical :: useqtb
!
!  Initialize MASH module.
!
   call init(beta,useqtb)
end subroutine

subroutine init_mfrpmd(beta,betamin,nb,nr,thermtyp)
   use types
   use mfrpmd, only : init
   real(dp) :: beta,betamin
   integer :: nb,nr
   character :: thermtyp
!
!  Initialize MF-RPMD module
!
   call init(beta,betamin,nb,nr,thermtyp)
end subroutine

! ================ Useful wrappers for potentials etc. =================

subroutine get_vad(q,nf,ns,Vad,U)
   use types
   use pes, only : potad
   integer, intent(in) :: nf,ns
   real(dp), intent(in) :: q(nf)
   real(dp), intent(out) :: Vad(ns), U(ns,ns)
!
!  Wrapper for adiabatic potential
!
   call potad(q,Vad,U)
end subroutine

subroutine pot_mat(q,nf,ns,V,dVdq)
   use types
   use pes, only : pot, grad
   integer, intent(in)  :: nf,ns
   real(dp), intent(in) :: q(nf)
   real(dp), intent(out):: v(ns,ns), dVdq(nf,ns,ns)
!
!  Wrapper for diabatic potential matrix
!
   call pot(q,v)
   call grad(q,dVdq)
end subroutine

subroutine mashpot(q,qe,pe,nf,ns,v)
   use types
   use mash, only : mash_pot
   integer, intent(in)  :: nf
   real(dp), intent(in) :: q(nf),qe(ns),pe(ns)
   real(dp), intent(out):: v
!
!  Wrapper for mash potential
!
   call mash_pot(q,qe,pe,v)
end subroutine


subroutine mashgrad(q,qe,pe,nf,ns,dvdq)
   use types
   use pes, only : potad, grad_a
   use mash, only : cstate
   integer, intent(in) :: nf, ns
   real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
   real(dp), intent(out) :: dvdq(nf)
!
!  Wrapper for gradient
!
   real(dp), allocatable :: Vad(:), U(:,:)
   integer :: a
   allocate(Vad(ns),U(ns,ns))
   call potad(q, Vad, U)
   call cstate(qe,pe,U,a)
   call grad_a(q, U, a, dvdq)
   deallocate(Vad, U)
end subroutine

subroutine nac(q,d,nf,ns)
   use types
   use pes, only : potad, mynac => nac
   integer, intent(in) :: nf, ns
   real(dp), intent(in) :: q(nf)
   real(dp), intent(out) :: d(nf,ns,ns)
!
!   Wrapper for nonadiabatic coupling vector
!
   real(dp), allocatable :: Vad(:), U(:,:)
   allocate(Vad(ns),U(ns,ns))
   call potad(q,Vad,U)
   call mynac(q,Vad,U,d)
   deallocate(Vad,U)
end subroutine

subroutine mashnacdir(q,cad,Vad,U,n,m,d,nf,ns)
   use types
   use pes, only : nacdir
   integer, intent(in) :: nf, ns
   real(dp), intent(in) :: q(nf), Vad(ns), U(ns,ns)
   complex(dpc), intent(in) :: cad(ns)
   integer, intent(in) :: n,m
   real(dp), intent(out) :: d(nf)
!
!   Wrapper for direction of momentum rescaling/reversal
!
   call nacdir(q,cad,Vad,U,n,m,d)
end subroutine

subroutine dia2ad(q,qe,pe,qa,pa,nf,ns)
   use types
   use pes, only : potad
   real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
   real(dp), intent(out) :: qa(ns), pa(ns)
   integer, intent(in) :: nf, ns
!
!  Convert diabatic amplitudes to adiabatic
!
   real(dp), allocatable :: Vad(:), U(:,:), eye(:,:)
   allocate(Vad(ns),U(ns,ns),eye(ns,ns))
   call potad(q,Vad,U)
   eye = 0.d0
   do n=1,ns
      eye(n,n) = 1.d0
   end do
   do n=1,ns
      if (dot_product(U(:,n),eye(:,n)).lt.0.d0) then
         U(:,n) = -U(:,n)
      end if
   end do
   qa(:) = matmul(qe, U)
   pa(:) = matmul(pe, U)
   deallocate(Vad,U,eye)
end subroutine

subroutine mash_pops(q, qe, pe, pop, rep, typ, nf, ns)
   use types
   use mash, only : pops
   real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
   integer, intent(in) :: nf, ns, typ
   character, intent(in) :: rep
   real(dp), intent(out) :: pop(ns)
!
!  Wrapper for mash.pops
!
   call pops(q, qe, pe, pop, rep, typ)
end subroutine


! =============== Main functions for running a trajectory ===============
subroutine runtrj(q, p, qe, pe, qt, pt, qet, pet, at, Et, &
   dt, ierr, nt, nf, ns)
   use types
   use pes, only : potad, grad_a
   use mash, only : store, evolve, ham, cstate, ham_a
   integer :: nt, nf, ns
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf), p(nf), qe(ns), pe(ns)
   real(dp), intent(out) :: qt(nt+1,nf), pt(nt+1,nf), &
      qet(nt+1,ns), pet(nt+1,ns), Et(nt+1)
   ! real(dp), intent(in) :: tol
   integer, intent(out) :: ierr, at(nt+1)
!
!  Run a trajectory
!
   real(dp), allocatable :: Vad(:), U(:,:), dvdq(:)
   ! real(dp) :: oldh, newh, refh
   logical :: fromfile = .false. ! Debugging option: 
   integer :: a , atmp
   ! refh = nf/beta

   allocate(Vad(ns),U(ns,ns),dvdq(nf))

   call potad(q,Vad,U)
   call cstate(qe,pe,U,a)

   if (fromfile) then
      print*, 'Loading from state.out'
      open (unit=4,file='state.out')
      read (4,*) q
      read (4,*) p
      read (4,*) qe
      read (4,*) pe
      read (4,*) a
      close (unit=4)
      call potad(q,Vad,U)
   end if

   ! oldh = ham(q,p,qe,pe)
   call grad_a(q,U,a,dvdq)
   do it = 1, nt
      call store(q, p, qe, pe, a, it, qt, pt, qet, pet, at)
      ! Et(it) = oldh
      Et(it) = ham_a(q,p,a)
      ! Et(it) = ham(q,p,qe,pe)
      atmp = a

      call evolve(q, p, qe, pe, Vad, U, dvdq, a, dt, ierr)

      ! newh = ham(q,p,qe,pe)
      ! if (abs((newh-oldh)/(refh)).gt.tol) then
      !    print*, 'energy not conserved'
      !    ! ierr = 1
      ! end if

      if (ierr.gt.0 .and. .not.fromfile) then
         open (unit=4,file='state.out')
         write (4,*) qt(it,:)
         write (4,*) pt(it,:)
         write (4,*) qet(it,:)
         write (4,*) pet(it,:)
         write (4,*) atmp
         close (unit=4)
      end if

      if (ierr.gt.0) exit ! Don't waste time if trajectory will be discarded

      if (q(1).ne.q(1)) then
         ! NaN
         ierr = 1
         ! stop 'Trajectory gives Not A Number!'
      end if
      ! oldh = newh
   end do
   if (ierr.eq.0) then
      call store(q, p, qe, pe, a, nt+1, qt, pt, qet, pet, at)
      ! Et(nt+1) = oldh
      Et(nt+1) = ham_a(q,p,a)
      ! Et(nt+1) = ham(q,p,qe,pe)
   end if

   deallocate(Vad,U,dvdq)
end subroutine

subroutine run_mfrp(dt,nt,nth,Vt,Tt,Pt,ns,nr)
   use types
   use mfrpmd, only : runtrj
   real(dp), intent(in) :: dt
   integer, intent(in) :: nt, nth, ns
   real(dp), intent(out) :: Vt(nt+1,nr), Tt(nt+1,nr), Pt(nt+1,ns)
!
!  Run thermostatted MF-RPMD trajectory and compute equilibrium averages
!
   call runtrj(dt,nt,nth,Vt,Tt,Pt)
end subroutine

! ============= Time-correlation functions ==============
subroutine tcf(q, p, qe, pe, cdt, & 
      Et, dt, ierr, rep, init, cftyp, nt, nf, ns)
   use types
   use mash, only : corrfunc, pops
   integer :: nt, nf, ns
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf), p(nf), qe(ns), pe(ns)
   real(dp), intent(out) :: cdt(nt+1,ns,ns)
   real(dp), intent(out) :: Et(nt+1)
   integer, intent(out) :: ierr
   character :: rep ! 'd','a','r' or 'e'
   character :: cftyp ! 'c' for Cleopatra, 'w' for WAB
   real(dp), allocatable :: qt(:,:), pt(:,:), qet(:,:), pet(:,:)
   integer, allocatable :: at(:)
   ! complex(dpc), allocatable :: CAB(:,:,:,:,:)
   complex(dpc), allocatable :: CAB(:,:)
   real(dp), allocatable :: popdtA(:), popdtB(:,:)

   allocate(qt(nt+1,nf), pt(nt+1,nf), qet(nt+1,ns), pet(nt+1,ns), at(nt+1))
   ! allocate(CAB(nt+1,ns,ns,ns,ns))

   call runtrj(q, p, qe, pe, qt, pt, qet, pet, at, Et, &
      dt, ierr, nt, nf, ns)

   if (cftyp.eq.'w') then
      ! ! MR generalization
      allocate(CAB(nt+1,ns))
      CAB = 0.d0
      call corrfunc(qt, qet, pet, CAB, rep, init, nt)

      do m=1,ns
         ! cdt(:,n,m) = real(CAB(:,n,n,m,m))
         cdt(:,init,m) = real(CAB(:,m))
      end do
      deallocate(CAB)
   else if (cftyp.eq.'c') then
      allocate(popdtA(ns),popdtB(nt+1,ns))
      ! call pops(qt(1,:), qet(1,:), pet(1,:), popdtA(:), 'r', 1)
      call pops(qt(1,:), qet(1,:), pet(1,:), popdtA(:), rep, 4)
      ! call pops(qt(1,:), qet(1,:), pet(1,:), popdtA(:), 'r', 3)
      do it=1,nt+1
         ! call pops(qt(it,:), qet(it,:), pet(it,:), popdtB(it,:), 'r', 1)
         call pops(qt(it,:), qet(it,:), pet(it,:), popdtB(it,:), rep, 2)
         ! call pops(qt(it,:), qet(it,:), pet(it,:), popdtB(it,:), 'r', 1)
      end do

      forall(n=1:ns,m=1:ns,it=1:nt+1)
         ! A(0)B(t) = rho(0) Phi(t) 
         cdt(it,n,m)=popdtA(n)*popdtB(it,m) 
      end forall
      deallocate(popdtA,popdtB)
   else if (cftyp.eq.'t') then
      allocate(popdtA(ns),popdtB(nt+1,ns))
      call pops(qt(1,:), qet(1,:), pet(1,:), popdtA(:), rep, 1)
      do it=1,nt+1
         call pops(qt(it,:), qet(it,:), pet(it,:), popdtB(it,:), rep, 2)
      end do
      forall(n=1:ns,m=1:ns,it=1:nt+1)
         cdt(it,n,m)=popdtA(n)*popdtB(it,m) 
      end forall
      deallocate(popdtA,popdtB)

   end if
   cdt = cdt * ns

   deallocate(qt,pt,qet,pet,at)
end subroutine


! ============= Parallelized functions ==============

subroutine runpar_poponly(q, p, qe, pe, bt, Et, ierr, &
   dt, nt, nf, ns, np)
   use types
   use mash, only : pops
   use omp_lib
   integer :: nt, nf, ns, np
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf,np), p(nf,np), qe(ns,np), pe(ns,np)
   real(dp), intent(out) :: bt(nt+1,ns)
   real(dp),allocatable :: dbt(:,:,:)
   real(dp),allocatable :: qt(:,:,:), pt(:,:,:), qet(:,:,:), pet(:,:,:)
   integer,allocatable :: at(:,:)
   real(dp), intent(out) :: Et(nt+1)
   integer, intent(out) :: ierr(np)
   real(dp), allocatable :: dEt(:,:)
!
!  Run a batch of trajectories in parallel
!  np = number of parallel threads
!
   allocate(dbt(nt+1,ns,np),qt(nt+1,nf,np),pt(nt+1,nf,np), &
            qet(nt+1,ns,np),pet(nt+1,ns,np),dEt(nt+1,np),at(nt+1,np))
!$omp parallel do default(shared) private(j)
   do j=1,np
      call runtrj(q(:,j), p(:,j), qe(:,j), pe(:,j), qt(:,:,j), &
         pt(:,:,j), qet(:,:,j), pet(:,:,j), at(:,j), dEt(:,j),  &
         dt, ierr(j), nt, nf, ns)

      do it=1,nt+1
         call pops(qt(it,:,j), qet(it,:,j), pet(it,:,j), dbt(it,:,j), 'd', 2)
      end do

      if (ierr(j).ne.0) then
         ! Error -> discard trajectory
         dbt(:,:,j) = 0.d0
      end if
   end do
!$omp end parallel do
   bt   = sum(dbt,dim=3)
   Et   = sum(dEt,dim=2)
   deallocate(dbt,dEt,qt,pt,qet,pet,at)
end subroutine

subroutine runpar_all(q, p, qe, pe, bdt, Et, ierr, &
   dt, nt, nf, ns, np)
   use types
   use mash, only : obsbls
   use omp_lib
   integer :: nt, nf, ns, np
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf,np), p(nf,np), qe(ns,np), pe(ns,np)
   complex(dpc), intent(out) :: bdt(nt+1,ns,ns)
   real(dp), intent(out) :: Et(nt+1)
   integer, intent(out) :: ierr(np)
   complex(dpc),allocatable :: dbt(:,:,:,:)
   real(dp), allocatable :: qt(:,:,:), pt(:,:,:), qet(:,:,:), pet(:,:,:)
   integer,allocatable :: at(:,:)
   real(dp), allocatable :: dEt(:,:)
!
!  Run a batch of trajectories in parallel
!  np = number of parallel threads
!
   allocate(dbt(nt+1,ns,ns,np),qt(nt+1,nf,np),pt(nt+1,nf,np), &
   qet(nt+1,ns,np),pet(nt+1,ns,np),dEt(nt+1,np))
   allocate(dEt(nt+1,np),at(nt+1,np))
!$omp parallel do default(shared) private(j)
   do j=1,np
      call runtrj(q(:,j), p(:,j), qe(:,j), pe(:,j), qt(:,:,j), &
         pt(:,:,j), qet(:,:,j), pet(:,:,j), at(:,j), dEt(:,j),  &
         dt, ierr(j), nt, nf, ns)
      do it=1,nt+1
         call obsbls(qt(it,:,j), qet(it,:,j), pet(it,:,j), & 
            dbt(it,:,:,j), 'd', 2, .false.)
      end do
      if (ierr(j).ne.0) then
         ! Error -> discard trajectory
         dbt(:,:,:,j) = 0.d0
      end if
   end do
!$omp end parallel do
   bdt = sum(dbt,dim=4)
   Et  = sum(dEt,dim=2)
   deallocate(dbt,qt,pt,qet,pet,dEt,at)
end subroutine

subroutine runpar_tcf(q, p, qe, pe, ct, Et, ierr, rep, init, cftyp, &
   dt, nt, nf, ns, np)
   use types
   use mash, only : pops
   use omp_lib
   integer :: nt, nf, ns, np
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf,np), p(nf,np), qe(ns,np), pe(ns,np)
   real(dp), intent(out) :: ct(nt+1,ns,ns)
   real(dp), intent(out) :: Et(nt+1)
   integer, intent(out) :: ierr(np)
   character :: rep, cftyp
   real(dp), allocatable :: dct(:,:,:,:), dEt(:,:)
!
!  Run a batch of trajectories in parallel
!  np = number of parallel threads
!
   allocate(dct(nt+1,ns,ns,np),dEt(nt+1,np))
!$omp parallel do default(shared) private(j)
   do j=1,np
      call tcf(q(:,j), p(:,j), qe(:,j), pe(:,j), &
         dct(:,:,:,j), dEt(:,j), dt, ierr(j), rep, init, cftyp, nt, nf, ns)

      if (ierr(j).ne.0) then
         ! Error -> discard trajectory
         dct(:,:,:,j) = 0.d0
      end if
   end do
!$omp end parallel do
   ct = sum(dct,dim=4)
   Et  = sum(dEt,dim=2)
   deallocate(dct,dEt)
end subroutine
