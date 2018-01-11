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
call create_System(System,Input,Control)
call free_Input(Input)
call calculate(System,Control)

call free_System(System)
call mem_report
call timer('Total program',Tcpu,Twall)
call close_OUT

end subroutine f_main

!!$subroutine f_main
!!$use precision, only : prec
!!$use int12core
!!$implicit none
!!$integer :: a,b,c
!!$real(prec) :: alpha,beta
!!$real(prec) :: f0,f1,f2
!!$
!!$a = 0
!!$b = 0
!!$c = 5
!!$alpha = 3._prec
!!$beta  = 13._prec
!!$
!!$f0 = int2_slater(a,b,c,alpha,beta)
!!$write(*,*) f0
!!$
!!$f1 = 0.5_prec*(&
!!$     int2_slater(a+1,b-1,c,alpha,beta) + &
!!$     int2_slater(a-1,b+1,c,alpha,beta) - &
!!$     int2_slater(a-1,b-1,c+2,alpha,beta) &
!!$     )
!!$write(*,*) f1
!!$
!!$f2 = 1.5_prec/(c+2._prec)*(&
!!$     int2_slater(a,b-2,c+2,alpha,beta) + &
!!$     int2_slater(a-2,b,c+2,alpha,beta) - &
!!$     int2_slater(a-2,b-2,c+4,alpha,beta) &
!!$     ) + &
!!$     int2_slater(a,b,c,alpha,beta)
!!$write(*,*) f2
!!$
!!$end subroutine f_main
