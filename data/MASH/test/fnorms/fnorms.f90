! Program to check if \int dc/norm rho_n Phi_m is equal to delta_nm
! where rho_n = a_N P_n + b_N and Phi_m = alpha_N P_m + beta_N
! and integral is over simplex


program fnorms
   implicit double precision (a-h,o-z)
   allocatable :: cpop(:,:), Pn(:)
   do ns=2,5
      hn = 0.d0
      do n=1,ns
         hn = hn + 1.d0/n
      end do
      alpha = (ns-1.d0)/(hn-1.d0)
      beta = (1.d0-alpha)/ns
      a = (ns+1.d0)/alpha
      b = -(1.d0+a*beta)

      ! ns = 4   ! number of states
      nsamp = 10000000

      allocate(cpop(ns,ns),Pn(ns))
      cpop = 0.d0
      do i=1,nsamp
         call sample(Pn,ns)
         do n=1,ns
            do m=1,ns
               cpop(n,m) = cpop(n,m) + (a*Pn(n)+b)*(alpha*Pn(m)+beta)
            end do
         end do
      end do
      cpop = cpop / nsamp * ns

      print*, 'N', ns
      do n=1,ns
         print*, cpop(n,:)
      end do
      deallocate(cpop,Pn)
   end do
end program

subroutine find_pop(Pn,n1,ns)
   implicit double precision (a-h,o-z)
   dimension Pn(ns)
!
!  Find index of highest populated state
!
   n1 = 1
   ! if (omit.eq.1) n1 = 2
   do i=1,ns
      if (Pn(i).gt.Pn(n1)) then
         n1 = i
      end if
   end do
   ! if (n1.eq.1) then
   !    n2 = 2
   ! else 
   !    n2 = 1
   ! end if
   ! do i=1,ns
   !    if ((Pn(i).gt.Pn(n2)).and.(i.ne.n1)) then
   !       n2 = i
   !    end if
   ! end do
   ! print*, 'n2', n2
end subroutine

subroutine find_coh(Pn,n1,n2,omit1,omit2,ns)
   implicit double precision (a-h,o-z)
   dimension Pn(ns)
!
!  Find indices of dominating coherence
!
   n1 = 1
   n2 = 2
   do i=1,ns
      do j=i+1,ns
         if (i.eq.j) cycle
         if (coh(Pn,i,j,ns).gt.coh(Pn,n1,n2,ns)) then
            if ((n1.eq.omit1).and.(n2.eq.omit2)) cycle
            n1 = i
            n2 = j
         end if
      end do
   end do
end subroutine

function coh(Pn,k,l,ns)
   implicit double precision (a-h,o-z)
   dimension Pn(ns)
!
!  Calculate population factor to coherence (k,l)
!  (excluding phase factor)
!
   coh = 2.d0*sqrt(Pn(k)*Pn(l))
end function

subroutine sample(Pn,ns)
   implicit double precision (a-h,o-z)
   dimension Pn(ns)
!
!  Sample Pn from simplex
!
   allocatable :: x(:), y(:)
   allocate(x(ns),y(ns))
   call gasdev(x,ns)
   call gasdev(y,ns)
   Pn = x**2+y**2
   Pn = Pn/sum(Pn) ! Pn are uniformly distributed and with sum 1
   deallocate(x,y)
end subroutine

subroutine gasdev(g,n)
   implicit double precision (a-h,o-z)
!
!  Generates an array of n normal deviates.
!
   dimension g(n)
   twopi = 2*dacos(-1.d0)
   do j = 1,n,2
      x = 1.d0-rand()! rand() is in [0,1) and can be exactly 0
      y = rand()
      r = dsqrt(-2*dlog(x))
      phi = twopi*y
      g(j) = r*dcos(phi)
      if (j .eq. n) exit
      g(j+1) = r*dsin(phi)
   enddo
end subroutine