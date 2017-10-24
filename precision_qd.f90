module precision_qd
use precision, only : prec
use qdmodule
implicit none

integer,parameter :: qd_size = (4*storage_size(0.d0))/8

interface assignment(=)
module procedure qd_from_prec
module procedure prec_from_qd
end interface assignment(=)

contains

subroutine qd_from_prec(qd_a,a)
implicit none
type(qd_real),intent(out) :: qd_a
real(prec),intent(in) :: a
real(prec) :: tmp
integer :: i

tmp = a
do i=1,3
   qd_a%re(i) = tmp
   tmp = tmp - qd_a%re(i)
enddo

qd_a%re(4) = 0.d0

end subroutine qd_from_prec

subroutine prec_from_qd(a,qd_a)
implicit none
real(prec),intent(out) :: a
type(qd_real),intent(in) :: qd_a

a = sum(real(qd_a%re(1:3),prec))

end subroutine prec_from_qd

end module precision_qd
