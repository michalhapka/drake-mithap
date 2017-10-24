module memory
use file_OUT, only : LOUT
use precision
use precision_qd
implicit none

private
public mem_report
public mem_alloc,mem_dealloc

interface mem_alloc
module procedure alloc_log1  ,alloc_log2  ,alloc_log3
module procedure alloc_int1  ,alloc_int2  ,alloc_int3
module procedure alloc_real1 ,alloc_real2 ,alloc_real3 ,alloc_real4
module procedure alloc_qd1   ,alloc_qd2
end interface mem_alloc

interface mem_dealloc
module procedure dealloc_log1  ,dealloc_log2  ,dealloc_log3
module procedure dealloc_int1  ,dealloc_int2  ,dealloc_int3
module procedure dealloc_real1 ,dealloc_real2 ,dealloc_real3 ,dealloc_real4
module procedure dealloc_qd1   ,dealloc_qd2
end interface mem_dealloc

integer(longint),save :: max_allocated_log  = 0_longint
integer(longint),save :: max_allocated_int  = 0_longint
integer(longint),save :: max_allocated_real = 0_longint
integer(longint),save :: max_allocated_qd   = 0_longint
integer(longint),save :: allocated_log  = 0_longint
integer(longint),save :: allocated_int  = 0_longint
integer(longint),save :: allocated_real = 0_longint
integer(longint),save :: allocated_qd   = 0_longint

integer(longint),parameter :: kilobyte = ibset(0_longint,10)
integer(longint),parameter :: megabyte = ibset(0_longint,20)
integer(longint),parameter :: gigabyte = ibset(0_longint,30)
integer(longint),parameter :: kilobyte_check = 10_longint ** 3
integer(longint),parameter :: megabyte_check = 10_longint ** 6
integer(longint),parameter :: gigabyte_check = 10_longint ** 9

integer,parameter :: message_LENGTH = 80
character(message_LENGTH),parameter :: &
     message_negative = 'Negative amount of allocated memory!'
character(message_LENGTH),parameter :: &
     message_alloc    = 'Trying to allocate an already allocated object!'
character(message_LENGTH),parameter :: &
     message_dealloc  = 'Trying to deallocate an unallocated object!'

contains

subroutine mem_report
implicit none

write(LOUT,'()')

write(LOUT,'(a)',advance='no') 'Max. memory allocated for LOGICAL:  '
call mem_report_number(max_allocated_log)
write(LOUT,'(a)',advance='no') 'Max. memory allocated for INTEGER:  '
call mem_report_number(max_allocated_int)
write(LOUT,'(a)',advance='no') 'Max. memory allocated for REAL:     '
call mem_report_number(max_allocated_real)
write(LOUT,'(a)',advance='no') 'Max. memory allocated for QD:       '
call mem_report_number(max_allocated_qd)

if(allocated_log/=0_longint) then
   write(LOUT,'(a)',advance='no') &
        '!!!MEMORY IS STILL ALLOCATED!!! for LOGICAL:  '
   call mem_report_number(allocated_log)
endif
if(allocated_int/=0_longint) then
   write(LOUT,'(a)',advance='no') &
        '!!!MEMORY IS STILL ALLOCATED!!! for INTEGER:  '
   call mem_report_number(allocated_int)
endif
if(allocated_real/=0_longint) then
   write(LOUT,'(a)',advance='no') &
        '!!!MEMORY IS STILL ALLOCATED!!! for REAL:     '
   call mem_report_number(allocated_real)
endif
if(allocated_qd/=0_longint) then
   write(LOUT,'(a)',advance='no') &
        '!!!MEMORY IS STILL ALLOCATED!!! for QD:       '
   call mem_report_number(allocated_qd)
endif

end subroutine mem_report

subroutine mem_report_number(number)
implicit none
integer(longint) :: number
character(9) :: string

if(number>gigabyte_check) then
   call prepare_string(gigabyte)
   write(LOUT,'(a9,a)') string,' Gb'
elseif(number>megabyte_check) then
   call prepare_string(megabyte)
   write(LOUT,'(a9,a)') string,' Mb'
elseif(number>kilobyte_check) then
   call prepare_string(kilobyte)
   write(LOUT,'(a9,a)') string,' kb'
else
   write(LOUT,'(i9,a)') number,'  b'
endif

contains

  subroutine prepare_string(base)
  implicit none
  integer(longint) :: base

  write(string,'(i5,a,i3.3)') number/base,'.', &
       (1000_longint*mod(number,base) + base/2_longint)/base

  end subroutine prepare_string

end subroutine mem_report_number

subroutine mem_modify_log(xsize,add)
implicit none
integer(longint) :: xsize
logical :: add

if(add) then
   allocated_log = allocated_log + xsize
   max_allocated_log = max(max_allocated_log,allocated_log)
else
   allocated_log = allocated_log - xsize
   if(allocated_log<0_longint) call write_message_and_stop(message_negative)
endif

end subroutine mem_modify_log

subroutine mem_modify_int(xsize,add)
implicit none
integer(longint) :: xsize
logical :: add

if(add) then
   allocated_int = allocated_int + xsize
   max_allocated_int = max(max_allocated_int,allocated_int)
else
   allocated_int = allocated_int - xsize
   if(allocated_int<0_longint) call write_message_and_stop(message_negative)
endif

end subroutine mem_modify_int

subroutine mem_modify_real(xsize,add)
implicit none
integer(longint) :: xsize
logical :: add

if(add) then
   allocated_real = allocated_real + xsize
   max_allocated_real = max(max_allocated_real,allocated_real)
else
   allocated_real = allocated_real - xsize
   if(allocated_real<0_longint) call write_message_and_stop(message_negative)
endif

end subroutine mem_modify_real

subroutine mem_modify_qd(xsize,add)
implicit none
integer(longint) :: xsize
logical :: add

if(add) then
   allocated_qd = allocated_qd + xsize
   max_allocated_qd = max(max_allocated_qd,allocated_qd)
else
   allocated_qd = allocated_qd - xsize
   if(allocated_qd<0_longint) call write_message_and_stop(message_negative)
endif

end subroutine mem_modify_qd

subroutine alloc_log1(x,n1)
implicit none
logical,allocatable :: x(:)
integer :: n1
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a)') n1,' logical'
      call write_message_and_stop(message,object)
   endif
   xsize = log_size*size(x,kind=longint)
   call mem_modify_log(xsize,.true.)
endif

end subroutine alloc_log1

subroutine alloc_log2(x,n1,n2)
implicit none
logical,allocatable :: x(:,:)
integer :: n1,n2
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a)') n1,' x',n2,' logical'
      call write_message_and_stop(message,object)
   endif
   xsize = log_size*size(x,kind=longint)
   call mem_modify_log(xsize,.true.)
endif

end subroutine alloc_log2

subroutine alloc_log3(x,n1,n2,n3)
implicit none
logical,allocatable :: x(:,:,:)
integer :: n1,n2,n3
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2,n3),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a,i10,a)') n1,' x',n2,' x',n3,' logical'
      call write_message_and_stop(message,object)
   endif
   xsize = log_size*size(x,kind=longint)
   call mem_modify_log(xsize,.true.)
endif

end subroutine alloc_log3

subroutine dealloc_log1(x)
implicit none
logical,allocatable :: x(:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = log_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_log(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_log1

subroutine dealloc_log2(x)
implicit none
logical,allocatable :: x(:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = log_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_log(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_log2

subroutine dealloc_log3(x)
implicit none
logical,allocatable :: x(:,:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = log_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_log(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_log3

subroutine alloc_int1(x,n1)
implicit none
integer,allocatable :: x(:)
integer :: n1
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a)') n1,' integer'
      call write_message_and_stop(message,object)
   endif
   xsize = int_size*size(x,kind=longint)
   call mem_modify_int(xsize,.true.)
endif

end subroutine alloc_int1

subroutine alloc_int2(x,n1,n2)
implicit none
integer,allocatable :: x(:,:)
integer :: n1,n2
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a)') n1,' x',n2,' integer'
      call write_message_and_stop(message,object)
   endif
   xsize = int_size*size(x,kind=longint)
   call mem_modify_int(xsize,.true.)
endif

end subroutine alloc_int2

subroutine alloc_int3(x,n1,n2,n3)
implicit none
integer,allocatable :: x(:,:,:)
integer :: n1,n2,n3
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2,n3),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a,i10,a)') n1,' x',n2,' x',n3,' integer'
      call write_message_and_stop(message,object)
   endif
   xsize = int_size*size(x,kind=longint)
   call mem_modify_int(xsize,.true.)
endif

end subroutine alloc_int3

subroutine dealloc_int1(x)
implicit none
integer,allocatable :: x(:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = int_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_int(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_int1

subroutine dealloc_int2(x)
implicit none
integer,allocatable :: x(:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = int_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_int(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_int2

subroutine dealloc_int3(x)
implicit none
integer,allocatable :: x(:,:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = int_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_int(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_int3

subroutine alloc_real1(x,n1)
implicit none
real(prec),allocatable :: x(:)
integer :: n1
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a)') n1,' real'
      call write_message_and_stop(message,object)
   endif
   xsize = real_size*size(x,kind=longint)
   call mem_modify_real(xsize,.true.)
endif

end subroutine alloc_real1

subroutine alloc_real2(x,n1,n2)
implicit none
real(prec),allocatable :: x(:,:)
integer :: n1,n2
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a)') n1,' x',n2,' real'
      call write_message_and_stop(message,object)
   endif
   xsize = real_size*size(x,kind=longint)
   call mem_modify_real(xsize,.true.)
endif

end subroutine alloc_real2

subroutine alloc_real3(x,n1,n2,n3)
implicit none
real(prec),allocatable :: x(:,:,:)
integer :: n1,n2,n3
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2,n3),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a,i10,a)') n1,' x',n2,' x',n3,' real'
      call write_message_and_stop(message,object)
   endif
   xsize = real_size*size(x,kind=longint)
   call mem_modify_real(xsize,.true.)
endif

end subroutine alloc_real3

subroutine alloc_real4(x,n1,n2,n3,n4)
implicit none
real(prec),allocatable :: x(:,:,:,:)
integer :: n1,n2,n3,n4
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2,n3,n4),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a,i10,a,i10,a)') &
           n1,' x',n2,' x',n3,' x',n4,' real'
      call write_message_and_stop(message,object)
   endif
   xsize = real_size*size(x,kind=longint)
   call mem_modify_real(xsize,.true.)
endif

end subroutine alloc_real4

subroutine dealloc_real1(x)
implicit none
real(prec),allocatable :: x(:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = real_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_real(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_real1

subroutine dealloc_real2(x)
implicit none
real(prec),allocatable :: x(:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = real_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_real(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_real2

subroutine dealloc_real3(x)
implicit none
real(prec),allocatable :: x(:,:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = real_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_real(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_real3

subroutine dealloc_real4(x)
implicit none
real(prec),allocatable :: x(:,:,:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = real_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_real(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_real4

subroutine alloc_qd1(x,n1)
implicit none
type(qd_real),allocatable :: x(:)
integer :: n1
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a)') n1,' qd'
      call write_message_and_stop(message,object)
   endif
   xsize = qd_size*size(x,kind=longint)
   call mem_modify_qd(xsize,.true.)
endif

end subroutine alloc_qd1

subroutine alloc_qd2(x,n1,n2)
implicit none
type(qd_real),allocatable :: x(:,:)
integer :: n1,n2
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message,object

if(allocated(x)) then
   call write_message_and_stop(message_alloc)
else
   allocate(x(n1,n2),stat=ierr,errmsg=message)
   if(ierr/=0) then
      write(object,'(i10,a,i10,a)') n1,' x',n2,' qd'
      call write_message_and_stop(message,object)
   endif
   xsize = qd_size*size(x,kind=longint)
   call mem_modify_qd(xsize,.true.)
endif

end subroutine alloc_qd2

subroutine dealloc_qd1(x)
implicit none
type(qd_real),allocatable :: x(:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = qd_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_qd(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_qd1

subroutine dealloc_qd2(x)
implicit none
type(qd_real),allocatable :: x(:,:)
integer(longint) :: xsize
integer :: ierr
character(message_LENGTH) :: message

if(allocated(x)) then
   xsize = qd_size*size(x,kind=longint)
   deallocate(x,stat=ierr,errmsg=message)
   if(ierr/=0) call write_message_and_stop(message)
   call mem_modify_qd(xsize,.false.)
else
   call write_message_and_stop(message_dealloc)
endif

end subroutine dealloc_qd2

subroutine write_message_and_stop(message1,message2)
implicit none
character(message_LENGTH) :: message1
character(message_LENGTH),optional :: message2

write(LOUT,'(a)') message1
if(present(message2)) write(LOUT,'(a)') message2

call mem_report

stop

end subroutine write_message_and_stop

end module memory
