module lproblem_qd
use file_OUT, only : LOUT
use precision, only : prec
use precision_qd
use memory
implicit none

private
public DecompositionData_qd
public simple_linsolver_qd
public init_linsolver_qd,free_linsolver_qd,use_linsolver_qd

type DecompositionData_qd
character(8) :: method
integer :: n
type(qd_real),allocatable :: A(:,:)
integer,allocatable :: perm(:)
end type DecompositionData_qd

contains

subroutine simple_linsolver_qd(LUStype,method,n,prec_A,prec_b,prec_x)
implicit none
character(1),intent(in) :: LUStype
character(*),intent(in) :: method
integer,intent(in) :: n
real(prec),intent(in) :: prec_A(:,:),prec_b(:)
real(prec),intent(out) :: prec_x(:)
type(qd_real),allocatable :: A(:,:),b(:),x(:)
integer,allocatable :: perm(:)
integer :: i

if(n<0) then

   write(LOUT,'(a)') 'Negative matrix dimension in linear equations solver!'
   stop

elseif(n==0) then

   write(LOUT,'(a)') 'WARNING!!! &
        &Zero matrix dimension in linear equations solver!'
   return

elseif(n==1) then

   prec_x(1) = prec_b(1)/prec_A(1,1)

else

   call mem_alloc(A,n,n)
   call mat_setU(LUStype,n,A,prec_A)

   select case(trim(method))
   case('CHOLESKY')

      do i=1,n
         x(i) = prec_b(i)
      enddo

      call cholesky_factor(n,A)
      call cholesky_solve(n,A,x)

   case('LDL1')

      call mem_alloc(perm,n)
      call mem_alloc(b,n)

      do i=1,n
         x(i) = prec_b(i)
      enddo

      call LDL1_factor(n,A,perm)
      call LDL1_solve(n,A,x,perm,b)

      call mem_dealloc(b)
      call mem_dealloc(perm)

   case('LDL2')

      call mem_alloc(perm,n)
      call mem_alloc(b,n)

      do i=1,n
         x(i) = prec_b(i)
      enddo

      call LDL2_factor(n,A,perm)
      call LDL2_solve(n,A,x,perm,b)

      call mem_dealloc(b)
      call mem_dealloc(perm)

   case('GAUSS')

      call mem_alloc(b,n)
      call mem_alloc(perm,n)

      call mat_setL(n,A)
      do i=1,n
         b(i) = prec_b(i)
      enddo

      call gauss(n,A,b,x,perm)

      call mem_dealloc(perm)
      call mem_dealloc(b)

   case default

      write(LOUT,'(2a)') 'Unknown method of solving linear equations: ',&
           trim(method)
      stop

   end select

   do i=1,n
      prec_x(i) = x(i)
   enddo

   call mem_dealloc(A)

endif

end subroutine simple_linsolver_qd

subroutine init_linsolver_qd(LHS,LUStype,method,n,prec_A)
implicit none
type(DecompositionData_qd) :: LHS
character(1),intent(in) :: LUStype
character(*),intent(in) :: method
integer,intent(in) :: n
real(prec),intent(in) :: prec_A(:,:)

if(n<0) then

   write(LOUT,'(a)') 'Negative matrix dimension in linear equations solver!'
   stop

elseif(n==0) then

   LHS%method = ''
   LHS%n = 0

   write(LOUT,'(a)') 'WARNING!!! &
        &Zero matrix dimension in linear equations solver!'
   return

elseif(n==1) then

   LHS%method = ''
   LHS%n = 1

   call mem_alloc(LHS%A,1,1)
   LHS%A(1,1) = 1._prec/prec_A(1,1)

else

   LHS%method = method
   LHS%n = n

   call mem_alloc(LHS%A,LHS%n,LHS%n)
   call mat_setU(LUStype,LHS%n,LHS%A,prec_A)

   select case(trim(LHS%method))
   case('CHOLESKY')

      call cholesky_factor(LHS%n,LHS%A)

   case('LDL1')

      call mem_alloc(LHS%perm,LHS%n)
      call LDL1_factor(LHS%n,LHS%A,LHS%perm)

   case('LDL2')

      call mem_alloc(LHS%perm,LHS%n)
      call LDL2_factor(LHS%n,LHS%A,LHS%perm)

   case('GAUSS')

      write(LOUT,'(a)') 'Gaussian elimination may be performed only by &
           &the simple_linsolver function!'
      stop

   case default

      write(LOUT,'(2a)') 'Unknown method of solving linear equations: ',&
           trim(LHS%method)
      stop

   end select

endif

end subroutine init_linsolver_qd

subroutine free_linsolver_qd(LHS)
implicit none
type(DecompositionData_qd) :: LHS

if(LHS%n>0) then

   if(allocated(LHS%perm)) call mem_dealloc(LHS%perm)
   call mem_dealloc(LHS%A)

endif

end subroutine free_linsolver_qd

subroutine use_linsolver_qd(LHS,n,prec_b,prec_x)
implicit none
type(DecompositionData_qd) :: LHS
integer,intent(in) :: n
real(prec),intent(in) :: prec_b(:)
real(prec),intent(out) :: prec_x(:)
type(qd_real),allocatable :: x(:),tmp(:)
integer :: i

if(LHS%n/=n) then

   write(LOUT,'(a)') 'RHS vector has improper size in use_linsolver!'
   stop

endif

if(LHS%n==1) then

   prec_x(1) = LHS%A(1,1)
   prec_x(1) = prec_x(1)*prec_b(1)

elseif(LHS%n>1) then

   call mem_alloc(x,n)

   do i=1,LHS%n
      x(i) = prec_b(i)
   enddo

   select case(trim(LHS%method))
   case('CHOLESKY')

      call cholesky_solve(LHS%n,LHS%A,x)

   case('LDL1')

      call mem_alloc(tmp,LHS%n)
      call LDL1_solve(LHS%n,LHS%A,x,LHS%perm,tmp)
      call mem_dealloc(tmp)

   case('LDL2')

      call mem_alloc(tmp,LHS%n)
      call LDL2_solve(LHS%n,LHS%A,x,LHS%perm,tmp)
      call mem_dealloc(tmp)

   case default

      write(LOUT,'(2a)') 'Unknown method of solving linear equations: ',&
           trim(LHS%method)
      stop

   end select

   do i=1,LHS%n
      prec_x(i) = x(i)
   enddo

   call mem_dealloc(x)

endif

end subroutine use_linsolver_qd

subroutine mat_setU(LUStype,n,A,prec_A)
implicit none
character(1),intent(in) :: LUStype
integer,intent(in) :: n
type(qd_real),intent(out) :: A(:,:)
real(prec),intent(in) :: prec_A(:,:)
integer :: i,j

select case(LUStype)
case('L','l')

   do j=1,n
      do i=1,j
         A(i,j) = prec_A(j,i)
      enddo
      do i=j+1,n
         A(i,j) = 0
      enddo
   enddo

case('U','u')

   do j=1,n
      do i=1,j
         A(i,j) = prec_A(i,j)
      enddo
      do i=j+1,n
         A(i,j) = 0
      enddo
   enddo

case('S','s')

   do j=1,n
      do i=1,j-1
         A(i,j) = 0.5_prec*(prec_A(i,j) + prec_A(j,i))
      enddo
      A(j,j) = prec_A(j,j)
      do i=j+1,n
         A(i,j) = 0
      enddo
   enddo

case default

   write(LOUT,'(2a)') 'Unknown part of matrix in linear equations solver: ',&
        trim(LUStype)
   stop

end select

end subroutine mat_setU

subroutine mat_setL(n,A)
implicit none
integer,intent(in) :: n
type(qd_real),intent(inout) :: A(:,:)
integer :: i,j

do j=1,n-1
   do i=j+1,n
      A(i,j) = A(j,i)
   enddo
enddo

end subroutine mat_setL

subroutine cholesky_factor(n,A)
implicit none
integer,intent(in) :: n
type(qd_real),intent(inout) :: A(n,n)
integer :: i,j,k
type(qd_real) :: rtmp

do j=1,n
   do i=1,j-1
      rtmp = A(i,j)
      do k=1,i-1
         rtmp = rtmp - A(k,i)*A(k,j)
      enddo
      A(i,j) = rtmp*A(i,i)
   enddo
   rtmp = A(j,j)
   do k=1,j-1
      rtmp = rtmp - A(k,j)**2
   enddo
   if(rtmp<0) then
      write(LOUT,'(a)') &
           'Matrix in linear equations solver is not positive definite!'
      stop
   endif
   A(j,j) = 1/sqrt(rtmp)
enddo

end subroutine cholesky_factor

subroutine cholesky_solve(n,A,x)
implicit none
integer,intent(in) :: n
type(qd_real),intent(in) :: A(:,:)
type(qd_real),intent(inout) :: x(:)
integer :: i,j
type(qd_real) :: rtmp

do j=1,n
   rtmp = x(j)
   do i=1,j-1
      rtmp = rtmp - A(i,j)*x(i)
   enddo
   x(j) = rtmp*A(j,j)
enddo

do j=n,1,-1
   rtmp = x(j)*A(j,j)
   do i=1,j-1
      x(i) = x(i) - rtmp*A(i,j)
   enddo
   x(j) = rtmp
enddo

end subroutine cholesky_solve

subroutine LDL1_factor(n,A,perm)
implicit none
integer,intent(in) :: n
type(qd_real),intent(inout) :: A(:,:)
integer,intent(out) :: perm(:)
integer :: i,j,k,k_max
type(qd_real) :: D,rtmp
logical :: do_warn

do i=1,n
   perm(i) = i
enddo

do_warn = .true.
do k=n,2,-1

   rtmp  = abs(A(k,k))
   k_max = k
   do i=k-1,1,-1
      if(abs(A(i,i))>rtmp) then
         rtmp  = abs(A(i,i))
         k_max = i
      endif
   enddo
   if(A(k_max,k_max)<=0.and.do_warn) then
      write(LOUT,'(a)') 'WARNING!!! &
           &Matrix in linear equations solver may be not positive definite!'
      do_warn = .false.
   endif

   call exchange(n,A,perm,k_max,k)

   D = 1/A(k,k)

   do j=1,k-1
      rtmp   = A(j,k)
      A(j,k) = A(j,k)*D
      do i=j,1,-1
         A(i,j) = A(i,j) - A(i,k)*rtmp
      enddo
   enddo

enddo

end subroutine LDL1_factor

subroutine LDL1_solve(n,A,x,perm,y)
implicit none
integer,intent(in) :: n
type(qd_real),intent(in) :: A(:,:)
type(qd_real),intent(inout) :: x(:)
integer,intent(in) :: perm(:)
type(qd_real) :: y(:)
integer :: i,j

do i=1,n
   y(i) = x(perm(i))
enddo

do j=n,2,-1
   do i=1,j-1
      y(i) = y(i) - A(i,j)*y(j)
   enddo
enddo

do i=1,n
   y(i) = y(i)/A(i,i)
enddo

do j=2,n
   do i=1,j-1
      y(j) = y(j) - y(i)*A(i,j)
   enddo
enddo

do i=1,n
   x(perm(i)) = y(i)
enddo

end subroutine LDL1_solve

subroutine LDL2_factor(n,A,perm)
implicit none
integer,intent(in) :: n
type(qd_real),intent(inout) :: A(:,:)
integer,intent(out) :: perm(:)
type(qd_real) :: alpha
integer :: i,j,k,i_max,j_max,k_max
type(qd_real) :: mu1,mu2,det,D11,D12,D22
type(qd_real) :: rtmp1,rtmp2

alpha = (1 + sqrt(17.d0))/8

do i=1,n
   perm(i) = i
enddo

k = n
do while(k>0)

   mu1   = -1
   k_max = 0
   mu2   = -1
   i_max = 0
   j_max = 0
   do j=k,1,-1
      if(abs(A(j,j))>mu1) then
         mu1   = abs(A(j,j))
         k_max = j
      endif
      do i=j-1,1,-1
         if(abs(A(i,j))>mu2) then
            mu2   = abs(A(i,j))
            i_max = i
            j_max = j
         endif
      enddo
   enddo
   if(mu2<=0) then
      if(mu1<=0) then
         write(LOUT,'(a)') 'Zero matrix in linear equations solver!'
         stop
      else
         mu2 = mu1*epsilon(mu1)
      endif
   endif

   if(mu1/mu2>alpha) then

      call exchange(n,A,perm,k_max,k)

      D11 = 1/A(k,k)

      do j=1,k-1
         rtmp1  = A(j,k)

         A(j,k) = rtmp1*D11
         do i=j,1,-1
            A(i,j) = A(i,j) - A(i,k)*rtmp1
         enddo

      enddo

      A(k,k) = D11

      k = k - 1

   else

      call exchange(n,A,perm,j_max,k)
      call exchange(n,A,perm,i_max,k-1)

      perm(k)   = -perm(k)
      perm(k-1) = -perm(k-1)

      D11 =  A(k,k)
      D12 = -A(k-1,k)
      D22 =  A(k-1,k-1)
      det = D11*D22 - D12**2
      D11 = D11/det
      D12 = D12/det
      D22 = D22/det

      do j=1,k-2
         rtmp1 = A(j,k-1)
         rtmp2 = A(j,k)

         A(j,k)   = rtmp1*D12 + rtmp2*D22
         do i=j,1,-1
            A(i,j) = A(i,j) - A(i,k)  *rtmp2
         enddo

         A(j,k-1) = rtmp1*D11 + rtmp2*D12
         do i=j,1,-1
            A(i,j) = A(i,j) - A(i,k-1)*rtmp1
         enddo

      enddo

      A(k-1,k-1) = D11
      A(k-1,k)   = D12
      A(k,k)     = D22

      k = k - 2

   endif
   
enddo

end subroutine LDL2_factor

subroutine LDL2_solve(n,A,x,perm,y)
implicit none
integer,intent(in) :: n
type(qd_real),intent(in) :: A(:,:)
type(qd_real),intent(inout) :: x(:)
integer,intent(in) :: perm(:)
type(qd_real) :: y(:)
integer :: i,j
type(qd_real) :: rtmp1,rtmp2

do i=1,n
   y(i) = x(abs(perm(i)))
enddo

j = n
do while(j>0)
   if(perm(j)>0) then

      rtmp1 = y(j)

      do i=1,j-1
         y(i) = y(i) - A(i,j)*rtmp1
      enddo

      y(j) = rtmp1*A(j,j)

      j = j - 1

   else

      rtmp1 = y(j-1)
      rtmp2 = y(j)

      do i=1,j-2
         y(i) = y(i) - A(i,j)  *rtmp2
      enddo
      do i=1,j-2
         y(i) = y(i) - A(i,j-1)*rtmp1
      enddo

      y(j-1) = A(j-1,j-1)*rtmp1 + A(j-1,j)*rtmp2
      y(j)   = A(j-1,j)  *rtmp1 + A(j,j)  *rtmp2

      j = j - 2

   endif
enddo

j = 0
do while(j<n)

   j = j + 1

   do i=1,j-1
      y(j) = y(j) - y(i)*A(i,j)
   enddo

   if(perm(j)<0) then

      j = j + 1

      do i=1,j-2
         y(j) = y(j) - y(i)*A(i,j)
      enddo

   endif

enddo

do i=1,n
   x(abs(perm(i))) = y(i)
enddo

end subroutine LDL2_solve

subroutine exchange(n,A,perm,j1_IN,j2_IN)
implicit none
integer,intent(in) :: n
type(qd_real),intent(inout) :: A(:,:)
integer,intent(inout) :: perm(:)
integer,intent(in) :: j1_IN,j2_IN
integer :: j1,j2
integer :: i
type(qd_real) :: rtmp
integer :: itmp

if(j1_IN/=j2_IN) then

   j1 = min(j1_IN,j2_IN)
   j2 = max(j1_IN,j2_IN)

   rtmp     = A(j1,j1)
   A(j1,j1) = A(j2,j2)
   A(j2,j2) = rtmp

   do i=1,j1-1
      rtmp    = A(i,j1)
      A(i,j1) = A(i,j2)
      A(i,j2) = rtmp
   enddo
   do i=j1+1,j2-1
      rtmp    = A(j1,i)
      A(j1,i) = A(i,j2)
      A(i,j2) = rtmp
   enddo
   do i=j2+1,n
      rtmp    = A(j1,i)
      A(j1,i) = A(j2,i)
      A(j2,i) = rtmp
   enddo

   itmp     = perm(j1)
   perm(j1) = perm(j2)
   perm(j2) = itmp

endif

end subroutine exchange

subroutine gauss(n,A,b,x,perm)
implicit none
integer,intent(in) :: n
type(qd_real),intent(inout) :: A(:,:),b(:)
type(qd_real),intent(out) :: x(:)
integer :: perm(:)
integer :: i,j,k,i_max,j_max
type(qd_real) :: rtmp
integer :: itmp

do i=1,n
   perm(i) = i
enddo

do k=n,2,-1

   rtmp  = abs(A(k,k))
   i_max = k
   j_max = k
   do j=1,k
      do i=1,k
         if(abs(A(i,j))>rtmp) then
            rtmp  = abs(A(i,j))
            i_max = i
            j_max = j
         endif
      enddo
   enddo

   if(i_max<k) then
      do j=1,n
         rtmp       = A(k,j)
         A(k,j)     = A(i_max,j)
         A(i_max,j) = rtmp
      enddo
      itmp        = perm(k)
      perm(k)     = perm(i_max)
      perm(i_max) = itmp
   endif

   if(j_max<k) then
      do i=1,k
         rtmp       = A(i,k)
         A(i,k)     = A(i,j_max)
         A(i,j_max) = rtmp
      enddo
      rtmp     = b(k)
      b(k)     = b(j_max)
      b(j_max) = rtmp
   endif

   do j=1,k-1
      rtmp = A(k,j)/A(k,k)
      do i=1,k-1
         A(i,j) = A(i,j) - rtmp*A(i,k)
      enddo
      A(k,j) = 0
      b(j) = b(j) - rtmp*b(k)
   enddo

enddo

do j=1,n
   rtmp = b(j)
   do i=1,j-1
      rtmp = rtmp - A(i,j)*b(i)
   enddo
   b(j) = rtmp/A(j,j)
enddo

do i=1,n
   x(perm(i)) = b(i)
enddo

end subroutine gauss

end module lproblem_qd
