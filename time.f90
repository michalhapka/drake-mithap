module time
use file_OUT, only : LOUT
use precision, only : prec => dble
implicit none

private
public timer

integer,save :: refYear = 0

real(prec),parameter :: second   = 1._prec
real(prec),parameter :: millisec = second * 1.e-3_prec
real(prec),parameter :: minute   = second * 60._prec
real(prec),parameter :: hour     = minute * 60._prec
real(prec),parameter :: day      = hour   * 24._prec

integer,parameter :: string_out_LENGTH = 20

contains

subroutine timer(action,Tcpu,Twall)
implicit none
character(*),intent(in) :: action
real(prec),intent(inout) :: Tcpu,Twall
real(prec) :: Tcpu_local,Twall_local
integer :: DateTime(8)
character(string_out_LENGTH) :: string_out

call date_and_time(values=DateTime)
if(refYear==0) refYear = DateTime(1)

call cpu_time(Tcpu_local)
call wall_time(Twall_local,DateTime)

if(adjustl(trim(action))/='START') then

   string_out = action(1:min(string_out_LENGTH,len(action)))

   write(LOUT,'(3a)',advance='no') 'CPU  time in ',string_out,' is '
   call print_time(Tcpu_local-Tcpu)
   write(LOUT,'(3a)',advance='no') 'Wall time in ',string_out,' is '
   call print_time(Twall_local-Twall)

endif

Tcpu  = Tcpu_local
Twall = Twall_local

end subroutine timer

subroutine wall_time(time,DateTime)
implicit none
real(prec) :: time
integer :: DateTime(8)
integer :: imonth,iyear

time = millisec * DateTime(8) &
     + second   * DateTime(7) &
     + minute   * DateTime(6) &
     + hour     * DateTime(5) &
     + day      *(DateTime(3) - 1)

do imonth=1,DateTime(2)-1
   select case(imonth)
   case(1,3,5,7,8,10,12)
      time = time + day * 31._prec
   case(4,6,9,11)
      time = time + day * 30._prec
   case(2)
      if(leapYear(DateTime(1))) then
         time = time + day * 29._prec
      else
         time = time + day * 28._prec
      endif
   end select
enddo

do iyear=refYear,DateTime(1)-1
   if(leapYear(iyear)) then
      time = time + day * 366._prec
   else
      time = time + day * 365._prec
   endif
enddo

end subroutine wall_time

function leapYear(year)
implicit none
logical :: leapYear
integer :: year

if(mod(year,4)/=0) then
   leapYear = .false.
elseif(mod(year,100)/=0) then
   leapYear = .true.
elseif(mod(year,400)/=0) then
   leapYear = .false.
else
   leapYear = .true.
endif

end function leapYear

subroutine print_time(full)
implicit none
real(prec) :: full
integer :: dec_hour,dec_minute
real(prec) :: dec_second

dec_hour   = int(full/hour)
dec_minute = int((full - dec_hour*hour)/minute)
dec_second = full - dec_hour*hour - dec_minute*minute

if(dec_hour/=0) then
   write(LOUT,'(i6,a2, i3,a2, f7.3,a2)') &
        dec_hour,  ' h', &
        dec_minute,' m', &
        dec_second,' s'
elseif(dec_minute/=0) then
   write(LOUT,'(8x, i3,a2, f7.3,a2)') &
        dec_minute,' m', &
        dec_second,' s'
else
   write(LOUT,'(8x, 5x, f7.3,a2)') &
        dec_second,' s'
endif

end subroutine print_time

end module time
