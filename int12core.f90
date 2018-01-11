module int12core
use file_OUT, only : LOUT
use precision, only : prec
use precision_qd
implicit none

private
public int1_slater,int2_slater

contains

function int1_slater(a,prec_alpha) result(prec_intval)
implicit none
real(prec) :: prec_intval
integer,intent(in) :: a
real(prec),intent(in) :: prec_alpha
type(qd_real) :: alpha,intval

if(a<-2) then
   write(LOUT,'(a,i5)') 'Wrong argument in int1_slater: ',a
   stop
endif

alpha = prec_alpha

intval = regfactorial(a+2,alpha)

prec_intval = intval

end function int1_slater

function int2_slater(a,b,c,prec_alpha,prec_beta) result(prec_intval)
implicit none
real(prec) :: prec_intval
integer,intent(in) :: a,b,c
real(prec),intent(in) :: prec_alpha,prec_beta
type(qd_real) :: alpha,beta,intval
type(qd_real) :: ratio,tmp
integer :: i

if(c<-1.or.a+b+c<-5) then
   write(LOUT,'(a,i5,a,i5,a,i5)') 'Wrong arguments in int2_slater: ',&
        a,', ',b,', ',c
   stop
endif

alpha = prec_alpha
beta  = prec_beta

ratio = c + 2
ratio = 2/ratio

tmp = 0
do i=0,(c+1)/2
   tmp = tmp + binomial(c+2,2*i+1) &
        *(int2_F(a+2*i+2,b+c-2*i+2,alpha,beta) &
        + int2_F(b+2*i+2,a+c-2*i+2,beta,alpha))
enddo
tmp = ratio*tmp

intval = tmp

prec_intval = intval

end function int2_slater

function int2_F(p,q,alpha,beta) result(Fpq)
implicit none
type(qd_real) :: Fpq
integer,intent(in) :: p,q
type(qd_real),intent(in) :: alpha,beta
type(qd_real) :: regfac_p,regfac_q,regfac_pq,inv_p,ratio,tmp
integer :: j

if(p>=0) then

   if(q>=0) then

      regfac_p = regfactorial(p,alpha + beta)
      regfac_q = regfactorial(q,beta)
      ratio = beta/(alpha + beta)

      tmp = 1
      do j=1,q
         tmp = tmp + binomial(p+j,p)*ratio**j
      enddo

      Fpq = regfac_p*regfac_q*tmp

   elseif(q==-1.and.alpha>beta) then

      regfac_p = regfactorial(p,alpha)

      ratio = beta/(alpha + beta)
      tmp = -log(ratio)

      ratio = alpha/(alpha + beta)
      do j=1,p
         tmp = tmp - ratio**j/j
      enddo

      Fpq = regfac_p*tmp

   else

      regfac_pq = regfactorial(p+q+1,alpha + beta)
      inv_p = p + 1
      inv_p = 1/inv_p
      ratio = alpha/(alpha + beta)

      tmp = confluent_2F1(1,p+q+2,p+2,ratio)

      Fpq = regfac_pq*inv_p*tmp

   endif

else
   write(LOUT,'(a,i5,a,i5)') 'Case not covered in int2_F: ',p,', ',q
   stop
endif

end function int2_F

function confluent_2F1(a,b,c,x) result(val)
implicit none
integer,parameter :: CONFRAC_maxit  = 10000
integer,parameter :: CONFRAC_thrmlt = 16
type(qd_real) :: val
integer,intent(in) :: a,b,c
type(qd_real),intent(in) :: x
type(qd_real) :: qd_eps,qd_tiny
type(qd_real) :: qd_a,qd_b,qd_c,ratio
type(qd_real) :: ai,bi,num,num_1,num_2,den,den_2,inv_den,prev,res
integer :: i

qd_eps  = CONFRAC_thrmlt*epsilon(val)
qd_tiny = tiny(val)

qd_a = a
qd_b = b
qd_c = c

num_2 = 1
num_1 = 1
den_2 = 0

prev = num_1

do i=1,CONFRAC_maxit

   ratio = ((qd_a + i)/(1 + i)) * ((qd_b + i)/(qd_c + i)) * x
   ai = -ratio
   bi = 1 + ratio

   num = bi*num_1 + ai*num_2
   den = bi + ai*den_2
   inv_den = 1/den

   res = num*inv_den
   if(abs(prev-res)<max(res*qd_eps,qd_tiny)) exit
   prev = res

   num_2 = num_1*inv_den
   num_1 = num*inv_den
   den_2 = inv_den

enddo
if(i>CONFRAC_maxit) then
   write(LOUT,'(a)') 'WARNING: too many iterations in continued_fraction!'
endif

val = 1 + (qd_a*(qd_b/qd_c) * x)/res

end function confluent_2F1

function regfactorial(n,alpha) result(regfac)
implicit none
type(qd_real) :: regfac
integer,intent(in) :: n
type(qd_real),intent(in) :: alpha
integer :: i

if(n<0) then
   write(LOUT,'(a,i5)') 'Wrong argument in regfactorial: ',n
   stop
endif

regfac = 1/alpha
do i=1,n
   regfac = regfac/alpha
   regfac = regfac*i
enddo

end function regfactorial

function binomial(n,k) result(bin)
implicit none
type(qd_real) :: bin
integer,intent(in) :: n,k
integer :: i

if(n<0.or.k<0.or.k>n) then
   write(LOUT,'(a,i5,a,i5)') 'Wrong arguments in binomial: ',n,', ',k
   stop
endif

bin = 1
do i=1,min(k,n-k)
   bin = bin*(n-i+1)
   bin = bin/i
enddo

end function binomial

end module int12core
