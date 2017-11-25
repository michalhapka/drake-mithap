module optimizer
use file_OUT, only : LOUT
use precision, only : prec
use memory
use commontypes

private
public powell

real(prec),parameter :: tauL = 0.5_prec*(sqrt(5._prec) - 1._prec)
real(prec),parameter :: tauS = 1._prec - tauL

interface
function energyInterface(vector,System,Control) result(energy)
import :: prec,SystemData,ControlData
real(prec) :: energy
real(prec) :: vector(:)
type(SystemData) :: System
type(ControlData) :: Control
end function energyInterface
end interface

contains

subroutine powell(n,xIO,energy,System,Control)
implicit none
integer,intent(in) :: n
real(prec) :: xIO(:)
procedure(energyInterface) :: energy
type(SystemData) :: System
type(ControlData) :: Control
integer :: i,iter
integer :: calls_counter
real(prec) :: f0,fi,delta_f
real(prec),allocatable :: x0(:),xi(:),xTMP(:),di(:),dALL(:,:)

write(LOUT,'()')
write(LOUT,'()')
write(LOUT,'(5x,a)') '--- Powell algorithm ---'

call mem_alloc(x0,n)
call mem_alloc(xi,n)
call mem_alloc(xTMP,n)
call mem_alloc(di,n)
call mem_alloc(dALL,n,n)

iter = 0

x0(:) = xIO
f0    = energy(x0,System,Control)
calls_counter = 1


write(LOUT,'()')
write(LOUT,'(1x,a,t11,12x,a,12x,a,4x,a)') &
     'iter','energy','energy change','function calls'
write(LOUT,'(i5,t11,f25.18,t61,i5)') iter,f0,calls_counter
flush(LOUT)

iter = 1

do
   if(mod(iter-1,n)==0) then
      dALL = 0._prec
      do i=1,n
         dALL(i,i) = x0(i)*Control%OPTSTEPMLT
      enddo
      write(LOUT,'()')
   endif

   xi(:) = x0
   fi    = f0
   do i=1,n
      di(:) = dALL(:,i)
      if(i>1) dALL(:,i-1) = di
      call linmin(xi,fi,di,xTMP,calls_counter,energy,System,Control)
      write(LOUT,'(2i5,f25.18,t61,i5)') iter,i,fi,calls_counter
      flush(LOUT)
   enddo

   di(:) = x0 - xi
   dALL(:,n) = di
   if(sqrt(dot_product(di,di))<Control%OPTNORMTHR) then
      calls_counter = 0
      delta_f = 0._prec
   else
      call linmin(xi,fi,di,xTMP,calls_counter,energy,System,Control)
      delta_f = abs(f0 - fi)
   endif
   write(LOUT,'(i5,t11,f25.18,es18.6,t61,i5)') iter,fi,delta_f,calls_counter
   flush(LOUT)

   x0(:) = xi
   f0    = fi

   iter = iter + 1

   if(delta_f<Control%OPTTHR) exit
enddo

xIO(:) = x0

call mem_dealloc(dALL)
call mem_dealloc(di)
call mem_dealloc(xTMP)
call mem_dealloc(xi)
call mem_dealloc(x0)

write(LOUT,'()')

end subroutine powell

subroutine linmin(xIO,fIO,d,x,calls_counter,energy,System,Control)
implicit none
real(prec),intent(inout) :: xIO(:),fIO
real(prec),intent(in) :: d(:)
real(prec) :: x(:)
integer,intent(out) :: calls_counter
procedure(energyInterface) :: energy
type(SystemData) :: System
type(ControlData) :: Control
real(prec) :: f_l,f_r,f,f_prev
real(prec) :: a,b,l,r,m
logical :: go_left

calls_counter = 0

x   = xIO - d
f_l = energy(x,System,Control)
x   = xIO + d

f_r = energy(x,System,Control)
calls_counter = calls_counter + 2

if(f_l>fIO.and.f_r>fIO) then
   a = -1._prec
   b =  1._prec
elseif(f_l<fIO.and.f_r<fIO) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Moving in both directions minimizes energy in linmin!'
   stop
elseif(abs(f_l-fIO)>epsilon(0._prec).or.abs(f_r-fIO)>epsilon(0._prec)) then
   go_left = (f_l<fIO)
   if(go_left) then
      f = f_l
      m = -1._prec
   else
      f = f_r
      m = +1._prec
   endif
   do
      f_prev = f
      m = m*Control%OPTLINMLT
      x = xIO + m*d
      f = energy(x,System,Control)
      calls_counter = calls_counter + 1
      if(f>f_prev) exit
   enddo
   if(go_left) then
      a = m
      b = m/Control%OPTLINMLT**2
      if(b>-0.5_prec*(1._prec + 1._prec/Control%OPTLINMLT)) b = 0._prec
   else
      a = m/Control%OPTLINMLT**2
      if(a<+0.5_prec*(1._prec + 1._prec/Control%OPTLINMLT)) a = 0._prec
      b = m
   endif
else
   return
endif

m = b - a
l = a + tauS*m
r = a + tauL*m

x   = xIO + l*d
f_l = energy(x,System,Control)
x   = xIO + r*d
f_r = energy(x,System,Control)
calls_counter = calls_counter + 2

do while(abs(f_l-f_r)>Control%OPTTHR)
   if(f_l>f_r) then
      a = l
      l = r
      r = a + tauL*(b - a)
      f_l = f_r
      x   = xIO + r*d
      f_r = energy(x,System,Control)
   else
      b = r
      r = l
      l = a + tauS*(b - a)
      f_r = f_l
      x   = xIO + l*d
      f_l = energy(x,System,Control)
   endif
   calls_counter = calls_counter + 1
enddo

m = 0.5_prec*(l + r)
x = xIO + m*d
f = energy(x,System,Control)
calls_counter = calls_counter + 1

xIO = x
fIO = f

end subroutine linmin

end module optimizer
