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
type(qd_real) :: regfac_p,regfac_q,ratio,tmp
integer :: j

if(p<0.or.q<0) then
   write(LOUT,'(a,i5,a,i5)') 'Case not covered in int2_F: ',p,', ',q
   stop
endif

regfac_p = regfactorial(p,alpha + beta)
regfac_q = regfactorial(q,beta)
ratio = beta/(alpha + beta)

tmp = 1
do j=1,q
   tmp = tmp + binomial(p+j,p)*ratio**j
enddo

Fpq = regfac_p*regfac_q*tmp

end function int2_F

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
