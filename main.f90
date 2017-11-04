subroutine f_main
use file_OUT
use precision, only : dble
use memory
use time
use commontypes
use inputread
use systemdef
use calcdriver
implicit none
type(InputData) :: Input
type(SystemData) :: System
type(ControlData) :: Control
real(dble) :: Tcpu, Twall

!call open_OUT
call timer('START',Tcpu,Twall)

call read_Input(Input)
call fill_Control(Control,Input)
!call create_System(System,Input,Control)
call free_Input(Input)

!call calculate(System,Control)

!call free_System(System)

call mem_report
call timer('Total program',Tcpu,Twall)
call close_OUT

end subroutine f_main
