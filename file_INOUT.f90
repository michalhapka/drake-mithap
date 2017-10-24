module file_OUT
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

private
public LOUT, LERR
public open_OUT, close_OUT

character(*),parameter :: OUT_name = 'DRAKE.OUT'
integer,save :: OUT_opened = 0

integer,save,protected :: LOUT = output_unit
integer,parameter :: LERR = error_unit

contains

subroutine open_OUT
implicit none

select case(OUT_opened)
case(1:)

   write(LERR,'(a)') &
        '!!!ERROR!!! Trying to open an already opened output file!'
   stop

case(0)

   open(newunit=LOUT,file=OUT_name,&
        action='write',status='replace')
   OUT_opened = 1

case(:-1)

   open(newunit=LOUT,file=OUT_name,&
        action='write',status='old',position='append')
   OUT_opened = -OUT_opened + 1

end select

end subroutine open_OUT

subroutine close_OUT
implicit none

select case(OUT_opened)
case(1:)

   close(LOUT,status='keep')
   LOUT = output_unit
   OUT_opened = -OUT_opened

case(0)

   ! case ignored

case(:-1)

   write(LERR,'(a)') &
        '!!!ERROR!!! Trying to close an output file that has not been opened!'
   stop

end select

end subroutine close_OUT

end module file_OUT

!-------------------------------------------------------------------------------

module file_IN
use file_OUT, only : LOUT, LERR
implicit none

private
public LIN
public open_IN, close_IN

character(*),parameter :: IN_name = 'DRAKE.IN'
logical,save :: IN_opened = .false.

integer,save,protected :: LIN = huge(0)

contains

subroutine open_IN
implicit none

if(IN_opened) then

   write(LERR,'(a)') &
        '!!!ERROR!!! Trying to open an already opened input file!'
   stop

else

   write(LOUT,'(2a)') 'Reading data from the input file ',IN_name

   open(newunit=LIN,file=IN_name,&
        action='read',status='old',position='rewind')
   IN_opened = .true.

endif

end subroutine open_IN

subroutine close_IN
implicit none

if(IN_opened) then

   close(LIN,status='keep')
   LIN = huge(0)
   IN_opened = .false.

else

   write(LERR,'(a)') &
        '!!!ERROR!!! Trying to close an input file that has not been opened!'
   stop

endif

end subroutine close_IN

end module file_IN
