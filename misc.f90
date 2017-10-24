module misc
use precision, only : prec
implicit none

private
public toUPPER,swap,pack_vector,unpack_vector

interface swap
module procedure swap_int,swap_log,swap_real
end interface swap

contains

subroutine toUPPER(string)
implicit none
integer,parameter :: codea = iachar('a'), codez = iachar('z')
integer,parameter :: offset = iachar('A') - iachar('a')
character(*),intent(inout) :: string
integer :: i,code

do i=1,len_trim(string)
   code = iachar(string(i:i))
   if(code>=codea .and. code<=codez) string(i:i) = achar(code + offset)
enddo

end subroutine toUPPER

elemental subroutine swap_int(v1,v2)
implicit none
integer,intent(inout) :: v1,v2
integer :: tmp

tmp = v1
v1  = v2
v2  = tmp

end subroutine swap_int

elemental subroutine swap_log(v1,v2)
implicit none
logical,intent(inout) :: v1,v2
logical :: tmp

tmp = v1
v1  = v2
v2  = tmp

end subroutine swap_log

elemental subroutine swap_real(v1,v2)
implicit none
real(prec),intent(inout) :: v1,v2
real(prec) :: tmp

tmp = v1
v1  = v2
v2  = tmp

end subroutine swap_real

subroutine pack_vector(vecS,n,vecL,mask)
implicit none
real(prec),intent(out) :: vecS(:)
integer,intent(in) :: n
real(prec),intent(in) :: vecL(:)
logical,intent(in) :: mask(:)
integer :: iS,iL

vecS = 0._prec

iS = 0
do iL=1,n
   if(mask(iL)) then
      iS = iS + 1
      vecS(iS) = vecL(iL)
   endif
enddo

end subroutine pack_vector

subroutine unpack_vector(vecS,n,vecL,mask)
implicit none
real(prec),intent(in) :: vecS(:)
integer,intent(in) :: n
real(prec),intent(inout) :: vecL(:)
logical,intent(in) :: mask(:)
integer :: iS,iL

iS = 0
do iL=1,n
   if(mask(iL)) then
      iS = iS + 1
      vecL(iL) = vecS(iS)
   endif
enddo

end subroutine unpack_vector

end module misc
