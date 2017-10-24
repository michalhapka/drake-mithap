module eproblem
use file_OUT, only : LOUT
use precision, only : prec
use memory
implicit none

private
public symU_diagonalize
public test_diagonalize

integer,parameter :: max_QL_iterSingle = 30
real(prec),parameter :: thr_ZERO     = tiny(0._prec)
real(prec),parameter :: thr_EPSILON  = epsilon(0._prec)
real(prec),parameter :: thr_ROTATION = sqrt(1/epsilon(0._prec))

contains

subroutine symU_diagonalize(method,n,eval,evec,origH,S)
implicit none
character(*),intent(in) :: method
integer,intent(in) :: n
real(prec),intent(inout) :: eval(:)
real(prec),intent(inout) :: evec(:,:)
real(prec),intent(in) :: origH(:,:)
real(prec),intent(in),optional :: S(:,:)
logical :: generalized
real(prec),allocatable :: H(:,:),V(:,:),diagH(:),superH(:),U(:,:)
integer,allocatable :: perm(:)
integer :: i,j,pj

generalized = present(S)

if(n<0) then

   write(LOUT,'(a)') 'Negative matrix dimension in diagonalization!'
   stop

elseif(n==0) then

   write(LOUT,'(a)') 'WARNING!!! Zero matrix dimension in diagonalization!'
   return

elseif(n==1) then

   if(generalized) then
      if(S(1,1)<0._prec) then
         write(LOUT,'(a)') &
              'Overlap matrix in diagonalization is not positive definite!'
         stop
      endif
      eval(1) = origH(1,1)/S(1,1)
      evec(1,1) = 1._prec/sqrt(S(1,1))
   else
      eval(1) = origH(1,1)
      evec(1,1) = 1._prec
   endif

else

   call mem_alloc(H,n,n)
   call mem_alloc(V,n,n)
   call mem_alloc(diagH,n)
   call mem_alloc(superH,n)
   call mem_alloc(perm,n)

   H = 0._prec
   do j=1,n
      do i=1,j
         H(i,j) = origH(i,j)
      enddo
   enddo

   if(generalized) then
      call mem_alloc(U,n,n)
      U = 0._prec
      do j=1,n
         do i=1,j
            U(i,j) = S(i,j)
         enddo
      enddo
      call LDL_factor(n,U,perm)
      call LDL_transform(n,H,U,perm,V)
   endif

   V = 0._prec
   do i=1,n
      V(i,i) = 1._prec
   enddo

   select case(trim(method))
   case('QL')

      call tridiagonal_upper(n,H,V)
      do i=1,n
         diagH(i) = H(i,i)
      enddo
      do i=1,n-1
         superH(i) = H(i,i+1)
      enddo
      superH(n) = 0._prec
      call implicitQL(n,diagH,superH,V)

   case('JACOBI')

      call jacobi_upper(n,H,V)
      do i=1,n
         diagH(i) = H(i,i)
      enddo

   case default

      write(LOUT,'(2a)') 'Unknown method of diagonalization: ',trim(method)
      stop

   end select

   if(generalized) then
      call LDL_solve(n,V,U,perm,superH)
      call mem_dealloc(U)
   endif

   call CombSort(n,diagH,perm)
   do j=1,n
      eval(j) = diagH(j)
      pj = perm(j)
      do i=1,n
         evec(i,j) = V(i,pj)
      enddo
   enddo

   call mem_dealloc(perm)
   call mem_dealloc(superH)
   call mem_dealloc(diagH)
   call mem_dealloc(V)
   call mem_dealloc(H)

endif

end subroutine symU_diagonalize

subroutine LDL_factor(n,A,perm)
implicit none
integer,intent(in) :: n
real(prec),intent(inout) :: A(:,:)
integer,intent(out) :: perm(:)
integer :: i,j,k,k_max
real(prec) :: rtmp
integer :: itmp

do i=1,n
   perm(i) = i
enddo

do k=n,2,-1

   rtmp  = abs(A(k,k))
   k_max = k
   do i=1,k-1
      if(abs(A(i,i))>rtmp) then
         rtmp  = abs(A(i,i))
         k_max = i
      endif
   enddo

   if(k_max<k) then
      rtmp           = A(k,k)
      A(k,k)         = A(k_max,k_max)
      A(k_max,k_max) = rtmp
      do i=1,k_max-1
         rtmp       = A(i,k)
         A(i,k)     = A(i,k_max)
         A(i,k_max) = rtmp
      enddo
      do i=k_max+1,k-1
         rtmp       = A(i,k)
         A(i,k)     = A(k_max,i)
         A(k_max,i) = rtmp
      enddo
      do i=k+1,n
         rtmp       = A(k,i)
         A(k,i)     = A(k_max,i)
         A(k_max,i) = rtmp
      enddo
      itmp        = perm(k)
      perm(k)     = perm(k_max)
      perm(k_max) = itmp
   endif

   do j=1,k-1
      rtmp   = A(j,k)
      A(j,k) = A(j,k)/A(k,k)
      do i=j,1,-1
         A(i,j) = A(i,j) - A(i,k)*rtmp
      enddo
   enddo

enddo

do k=1,n
   if(A(k,k)<0._prec) then
      write(LOUT,'(a)') &
           'Overlap matrix in diagonalization is not positive definite!'
      stop
   endif
   A(k,k) = 1._prec/sqrt(A(k,k))
enddo

end subroutine LDL_factor

subroutine LDL_transform(n,A,U,perm,B)
implicit none
integer,intent(in) :: n
real(prec),intent(inout) :: A(:,:)
real(prec),intent(in) :: U(:,:)
integer,intent(in) :: perm(:)
real(prec) :: B(:,:)
integer :: i,j,k
real(prec) :: rtmp

do i=1,n
   k = perm(i)
   do j=1,k
      B(i,j) = A(j,k)
   enddo
   do j=k+1,n
      B(i,j) = A(k,j)
   enddo
enddo

do k=1,n
   do j=n,1,-1
      rtmp = B(j,k)
      do i=1,j-1
         B(i,k) = B(i,k) - U(i,j)*rtmp
      enddo
      B(j,k) = U(j,j)*rtmp
   enddo
enddo

do j=1,n
   k = perm(j)
   do i=1,j
      A(i,j) = B(i,k)
   enddo
enddo

do j=n,1,-1
   do i=1,j-1
      rtmp = U(i,j)
      do k=1,i
         A(k,i) = A(k,i) - A(k,j)*rtmp
      enddo
   enddo
   rtmp = U(j,j)
   do k=1,j
      A(k,j) = A(k,j)*rtmp
   enddo
enddo

end subroutine LDL_transform

subroutine LDL_solve(n,A,U,perm,b)
implicit none
integer,intent(in) :: n
real(prec),intent(inout) :: A(:,:)
real(prec),intent(in) :: U(:,:)
integer,intent(in) :: perm(:)
real(prec) :: b(:)
integer :: i,j,k
real(prec) :: rtmp

do k=1,n
   do j=1,n
      rtmp = A(j,k)*U(j,j)
      do i=j-1,1,-1
         rtmp = rtmp - b(i)*U(i,j)
      enddo
      b(j) = rtmp
   enddo
   do j=1,n
      A(perm(j),k) = b(j)
   enddo
enddo

end subroutine LDL_solve

subroutine tridiagonal_upper(n,A,V)
implicit none
integer,intent(in) :: n
real(prec),intent(inout) :: A(:,:),V(:,:)
integer :: i,j,k,j_1
real(prec) :: sumoff,scale,xnorm,xnorm2,kfac
real(prec) :: rtmp,rtmp1,rtmp2

do j=n,3,-1
   j_1 = j - 1

   sumoff = 0._prec
   do i=1,j-2
      sumoff = sumoff + abs(A(i,j))
   enddo
   if(sumoff<thr_ZERO) cycle

   scale = 1._prec/(sumoff + abs(A(j_1,j)))
   do i=1,j_1
      A(i,j) = A(i,j)*scale
   enddo

   xnorm2 = 0._prec
   do i=1,j_1
      xnorm2 = xnorm2 + A(i,j)**2
   enddo
   xnorm = sqrt(xnorm2); if(A(j_1,j)<0._prec) xnorm = -xnorm
   rtmp = 1._prec/sqrt(xnorm2 + xnorm*A(j_1,j))
   do i=1,j-2
      V(i,j) = A(i,j)*rtmp
   enddo
   V(j_1,j) = (A(j_1,j) + xnorm)*rtmp

   do i=1,j_1
      rtmp1 = 0._prec
      rtmp2 = V(i,j)
      do k=1,i-1
         rtmp1  = rtmp1  + A(k,i)*V(k,j)
         A(k,j) = A(k,j) + A(k,i)*rtmp2
      enddo
      A(i,j) = rtmp1 + A(i,i)*rtmp2
   enddo

   rtmp = 0._prec
   do i=1,j_1
      rtmp = rtmp + A(i,j)*V(i,j)
   enddo
   kfac = 0.5_prec*rtmp
   do i=1,j_1
      A(i,j) = A(i,j) - kfac*V(i,j)
   enddo

   do i=1,j_1
      rtmp1 = V(i,j)
      rtmp2 = A(i,j)
      do k=1,i
         A(k,i) = A(k,i) - rtmp1*A(k,j) - rtmp2*V(k,j)
      enddo
   enddo

   do i=1,j-2
      A(i,j) = 0._prec
   enddo
   A(j_1,j) = -xnorm/scale

enddo

do j=3,n
   j_1 = j - 1

   if(abs(V(j_1,j))<thr_ZERO) cycle

   do i=1,j_1
      rtmp = 0._prec
      do k=1,j_1
         rtmp = rtmp + V(k,i)*V(k,j)
      enddo
      do k=1,j_1
         V(k,i) = V(k,i) - rtmp*V(k,j)
      enddo
   enddo

   do i=1,j_1
      V(i,j) = 0._prec
   enddo

enddo

end subroutine tridiagonal_upper

subroutine implicitQL(n,A,B,V)
implicit none
integer,intent(in) :: n
real(prec),intent(inout) :: A(:),B(:),V(:,:)
integer :: jstart,jend,j,j_1,i
integer :: iter,iterSingle,Nrot
real(prec) :: theta,t,c,s,shift,rot
real(prec) :: d1,d2,e,f1,f2,af1,af2,pythag,l12

iter = 0
iterSingle = 0
Nrot = 0

jstart = 1
do while(jstart<n)

   if(iterSingle<max_QL_iterSingle) then
      jend = jstart
      do while(jend<n)
         if(abs(B(jend))<thr_EPSILON*(abs(A(jend))+abs(A(jend+1)))) exit
         jend = jend + 1
      enddo
   else
      write(LOUT,'(a,i6)') &
           'WARNING!!! Too many iterations to converge eigenvalue: ',jstart
      jend = jstart
   endif

   if(jend==jstart) then

      jstart = jstart + 1
      iterSingle = 0

   else

      iter = iter + 1
      iterSingle = iterSingle + 1

      d1 = A(jstart)
      d2 = A(jstart+1)
      e  = B(jstart)
      theta = (d2 - d1)/(2._prec*e)
      if(abs(theta)<thr_ROTATION) then
         if(theta>0._prec) then
            t = 1._prec/(theta + sqrt(theta**2 + 1._prec))
         else
            t = 1._prec/(theta - sqrt(theta**2 + 1._prec))
         endif
      else
         t = 1._prec/(2._prec*theta)
      endif
      shift = d1 - t*e

      c = 1._prec
      s = 1._prec
      d1 = A(jend)
      f2 = d1 - shift
      do j=jend,jstart+1,-1
         j_1 = j - 1

         d2 = d1
         d1 = A(j_1)
         e  = c*B(j_1)
         f1 = s*B(j_1)

         af1 = abs(f1)
         af2 = abs(f2)
         if(af1>af2) then
            pythag = af1*sqrt(1._prec + (af2/af1)**2)
         else
            pythag = af2*sqrt(1._prec + (af1/af2)**2)
         endif
         c = f2/pythag
         s = f1/pythag

         l12 = s*(d1 - d2) + 2._prec*c*e
         d1   = d1 - s*l12
         A(j) = d2 + s*l12
         B(j) = pythag
         f2   = c*l12 - e

         Nrot = Nrot + 1
         do i=1,n
            rot = V(i,j_1)
            V(i,j_1) = c*rot - s*V(i,j)
            V(i,j)   = s*rot + c*V(i,j)
         enddo
      enddo
      A(jstart) = d1
      B(jstart) = f2
      B(jend) = 0._prec

   endif

enddo

end subroutine implicitQL

subroutine jacobi_upper(n,A,V)
implicit none
integer,intent(in) :: n
real(prec),intent(inout) :: A(:,:),V(:,:)
integer :: iter
integer :: i,j,k,imax,jmax
real(prec) :: rtmp,sumoff2,maxoff
real(prec) :: theta,t,c,s,tau
real(prec) :: ival,jval

iter = 0
do

   imax = 0
   jmax = 0
   sumoff2 = 0._prec
   maxoff  = 0._prec
   do j=2,n
      do i=1,j-1
         rtmp = abs(A(i,j))
         sumoff2 = sumoff2 + rtmp**2
         if(rtmp>maxoff) then
            maxoff = rtmp
            imax = i
            jmax = j
         endif
      enddo
   enddo
   sumoff2 = sqrt(sumoff2)

   if(sumoff2<thr_EPSILON) exit
   iter = iter + 1

   theta = (A(jmax,jmax) - A(imax,imax))/(2._prec*A(imax,jmax))
   if(abs(theta)<thr_ROTATION) then
      if(theta>0._prec) then
         t = 1._prec/(theta + sqrt(theta**2 + 1._prec))
      else
         t = 1._prec/(theta - sqrt(theta**2 + 1._prec))
      endif
   else
      t = 1._prec/(2._prec*theta)
   endif
   c = 1._prec/sqrt(t**2 + 1._prec)
   s = t*c
   tau = s/(1._prec + c)

   A(imax,imax) = A(imax,imax) - t*A(imax,jmax)
   A(jmax,jmax) = A(jmax,jmax) + t*A(imax,jmax)
   A(imax,jmax) = 0._prec
   do k=1,imax-1
      ival = A(k,imax)
      jval = A(k,jmax)
      A(k,imax) = ival - s*(jval + tau*ival)
      A(k,jmax) = jval + s*(ival - tau*jval)
   enddo
   do k=imax+1,jmax-1
      ival = A(imax,k)
      jval = A(k,jmax)
      A(imax,k) = ival - s*(jval + tau*ival)
      A(k,jmax) = jval + s*(ival - tau*jval)
   enddo
   do k=jmax+1,n
      ival = A(imax,k)
      jval = A(jmax,k)
      A(imax,k) = ival - s*(jval + tau*ival)
      A(jmax,k) = jval + s*(ival - tau*jval)
   enddo

   do k=1,n
      ival = V(k,imax)
      jval = V(k,jmax)
      V(k,imax) = ival - s*(jval + tau*ival)
      V(k,jmax) = jval + s*(ival - tau*jval)
   enddo

enddo

end subroutine jacobi_upper

subroutine CombSort(n,val,perm)
implicit none
integer,intent(in) :: n
real(prec),intent(inout) :: val(:)
integer,intent(out) :: perm(:)
integer :: i,gap,itmp
real(prec) :: rtmp
logical :: swapped

do i=1,n
   perm(i) = i
enddo

gap = n
swapped = .true.
do while(gap>1.or.swapped)
   gap = max(1,(gap*10)/13)
   if(gap==9.or.gap==10) gap = 11
   swapped = .false.

   do i=1,n-gap
      if (val(i+gap)<val(i)) then

         rtmp = val(i)
         val(i) = val(i+gap)
         val(i+gap) = rtmp

         itmp = perm(i)
         perm(i) = perm(i+gap)
         perm(i+gap) = itmp

         swapped = .true.
      endif
   enddo

enddo

end subroutine CombSort

subroutine test_diagonalize(n,eval,evec,H,S)
implicit none
integer,intent(in) :: n
real(prec),intent(in) :: eval(:)
real(prec),intent(in) :: evec(:,:)
real(prec),intent(in) :: H(:,:)
real(prec),intent(in),optional :: S(:,:)
logical :: generalized
real(prec),allocatable :: HV(:,:),SV(:,:)
real(prec) :: max_resid,mean_resid,max_dia,mean_dia,max_off,mean_off
real(prec) :: rtmp
integer :: i,j

generalized = present(S)

call mem_alloc(HV,n,n)
call mem_alloc(SV,n,n)

call multiply_upper(n,H,evec,HV)
if(generalized) then
   call multiply_upper(n,S,evec,SV)
else
   SV(:,:) = evec(1:n,1:n)
endif

max_resid  = 0._prec
mean_resid = 0._prec
do i=1,n
   rtmp = sqrt(sum((HV(:,i) - eval(i)*SV(:,i))**2))
   max_resid  = max(max_resid,rtmp)
   mean_resid = mean_resid + rtmp
enddo
mean_resid = mean_resid/n

max_dia  = 0._prec
mean_dia = 0._prec
do i=1,n
   rtmp = sqrt(sum(evec(1:n,i)*SV(:,i))) - 1._prec
   max_dia  = max(max_dia,abs(rtmp))
   mean_dia = mean_dia + rtmp
enddo
mean_dia = mean_dia/n

max_off  = 0._prec
mean_off = 0._prec
if(generalized) then
   do j=1,n
      do i=1,n
         if(i==j) cycle
         rtmp = sum(evec(1:n,i)*SV(:,j))
         max_off  = max(max_off,abs(rtmp))
         mean_off = mean_off + rtmp
      enddo
   enddo
   mean_off = mean_off/max(n*(n-1),1)
else
   do j=2,n
      do i=1,j-1
         rtmp = sum(evec(1:n,i)*SV(:,j))
         max_off  = max(max_off,abs(rtmp))
         mean_off = mean_off + rtmp
      enddo
   enddo
   mean_off = mean_off/max((n*(n-1))/2,1)
endif

call mem_dealloc(SV)
call mem_dealloc(HV)

write(LOUT,'(52x,a)') 'maximal     mean'
if(generalized) then
   write(LOUT,'(a,2es11.3)') &
        'Norm of residual vector:          ||HV - e SV|| =',max_resid,mean_resid
else
   write(LOUT,'(a,2es11.3)') &
        'Norm of residual vector:          ||HV -  e V|| =',max_resid,mean_resid
endif
write(LOUT,'(a,2es11.3)') &
     'Norm of eigenvector:               ||i|| - 1 =   ',max_dia,mean_dia
write(LOUT,'(a,2es11.3)') &
     'Product of different eigenvectors:    <i|j> =    ',max_off,mean_off

end subroutine test_diagonalize

subroutine multiply_upper(n,A,B,C)
implicit none
integer,intent(in) :: n
real(prec),intent(in) :: A(:,:),B(:,:)
real(prec),intent(out) :: C(:,:)
integer :: i,j,k
real(prec) :: rtmp1,rtmp2

do j=1,n
   do i=1,n
      rtmp1 = 0._prec
      rtmp2 = B(i,j)
      do k=1,i-1
         rtmp1  = rtmp1  + A(k,i)*B(k,j)
         C(k,j) = C(k,j) + A(k,i)*rtmp2
      enddo
      C(i,j) = rtmp1 + A(i,i)*rtmp2
   enddo
enddo

end subroutine multiply_upper

end module eproblem
