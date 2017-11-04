module inputread
use file_OUT, only : LOUT
use precision, only : prec
use file_IN
use memory
use commontypes
implicit none

private
public read_Input

integer,parameter :: slength = 80
character(*),parameter :: sfmt = '(a80)'
character(*),parameter :: set_comment = '!#'
character(*),parameter :: set_int     = '+-0123456789'
character(*),parameter :: set_real    = set_int//'.ED'
character(*),parameter :: end_flag    = '$'

interface read_value_and_advance_line
module procedure read_int_and_advance_line
module procedure read_real_and_advance_line
module procedure read_string_and_advance_line
end interface read_value_and_advance_line

interface read_value_random
module procedure read_log_random
module procedure read_int_random
module procedure read_real_random
end interface read_value_random

contains

subroutine read_Input(Input)
implicit none
type(InputData) :: Input

call init_Input(Input)

call open_IN
call read_InputFile(Input)
call close_IN

call print_Input(Input)

call print_EndSection

end subroutine read_Input

subroutine read_InputFile(Input)
use misc, only : swap
use angmom
implicit none
type(InputData) :: Input
character(slength) :: line,keyword
integer :: icnt,iorb,ipair1,ipair2,imult
integer :: i,j
logical :: swapped
character(8) :: smult, lorb_tmp
logical,allocatable :: isOrb(:),isPair(:,:)

call read_line(line)
call read_keyword_and_advance_line(line,'Z','=')
call read_value_and_advance_line(line,Input%nucZ,end_flag)

call read_line(line)
call read_keyword_and_advance_line(line,'MAXL','=')
call read_value_and_advance_line(line,Input%maxl,end_flag)

call read_line(line)
call read_keyword_and_advance_line(line,'NORB','=')
!call read_value_and_advance_line(line,Input%norb,end_flag)
call mem_alloc(Input%norb,Input%maxl+1)
!read(line, *) (Input%norb(i), i=1,Input%maxl+1)
call read_norbs(line,Input%maxl+1,Input%norb,icnt)

if(icnt.lt.Input%maxl+1) then
   write(LOUT,'(a, i5, a)') 'Expected', Input%maxl+1, ' orbitals!'
   call incorrect_data
endif

if(minval(Input%norb)<1) then
   write(LOUT,'(a)') 'Incorrect number of orbitals!'
   call incorrect_data
endif

call read_line(line)
call read_keyword_and_advance_line(line,'CALC','=')
call read_value_and_advance_line(line,Input%calc_type,end_flag)
if(all(trim(Input%calc_type)/=possible_calc_type)) then
   write(LOUT,'(a)') 'Unrecognizable type of calculations! Choose one of:'
   write(LOUT,'(a)') possible_calc_type
   call incorrect_data
endif

call read_line(line)
if(check_keyword(line,'ETA','=')) then

   call read_keyword_and_advance_line(line,'ETA','=')
   call read_value_and_advance_line(line,Input%Neta,end_flag)
   if(Input%Neta<1) then
      write(LOUT,'(a)') 'Incorrect number of eta parameters!'
      call incorrect_data
   endif
   call mem_alloc(Input%eta,Input%Neta)
   call read_value_random(Input%Neta,Input%eta)
   icnt = 1
   swapped = .true.
   do while(icnt<Input%Neta.and.swapped)
      swapped = .false.
      i = Input%Neta
      do while(icnt<i)
         if(Input%eta(i-1)<Input%eta(i)) then
            call swap(Input%eta(i-1),Input%eta(i))
            swapped = .true.
         elseif(Input%eta(i-1)==Input%eta(i)) then
            do j=i,Input%Neta-1
               Input%eta(j) = Input%eta(j+1)
            enddo
            Input%Neta = Input%Neta - 1
         endif
         i = i - 1
      enddo
      icnt = icnt + 1
   enddo

else

   if(index(Input%calc_type,'SCF')==0) then
      write(LOUT,'(a)') &
           'Eta parameters must be specified for post SCF calculations!'
      call incorrect_data
   endif
   
   call cancel_line
   Input%Neta = 0
   call mem_alloc(Input%eta,Input%Neta)

endif


call read_line(line)
call read_keyword_and_advance_line(line,'EXPONENTS','=')
call read_value_and_advance_line(line,Input%Nexponents,end_flag)
if(Input%Nexponents<1) then
   write(LOUT,'(a)') 'Incorrect number of exponents!'
   call incorrect_data
endif
call mem_alloc(Input%exponents,Input%Nexponents)
call mem_alloc(Input%optimized,Input%Nexponents)

call read_value_random(Input%Nexponents,Input%exponents)

call read_line(line)
if(check_keyword(line,'OPTIMIZE',end_flag)) then
   call read_keyword_and_advance_line(line,'OPTIMIZE',end_flag)
   call read_value_random(Input%Nexponents,Input%optimized)
else
   call cancel_line
   Input%optimized = .true.
endif
if(index(Input%calc_type,'OPT')==0) Input%optimized = .false.

call read_line(line)
call read_keyword_and_advance_line(line,'ORBS','=')
call read_value_and_advance_line(line,Input%orbs_basis,end_flag)
if(all(trim(Input%orbs_basis)/=possible_basis)) then
   write(LOUT,'(a)') 'Unrecognizable type of orbital basis set! Choose one of:'
   write(LOUT,'(a)') possible_basis
   call incorrect_data
endif

select case(trim(Input%orbs_basis))
case('COMMON')

   allocate(Input%OrbInput_L(Input%maxl+1))

do i=1,size(Input%OrbInput_L)
   allocate(Input%OrbInput_L(i)%OrbInput(1))

   call read_line(line)
   call read_keyword_and_advance_line(line,'LORB','=')
   call read_value_and_advance_line(line,lorb_tmp,end_flag)
   Input%OrbInput_L(i)%lorb =  get_angmom_l(lorb_tmp)

   if(Input%OrbInput_L(i)%lorb.gt.Input%maxl) then
   write(LOUT,'(2a)') 'Expecting orbitals up to ', get_angmom_name(Input%maxl)
   call incorrect_data 
   endif

   call read_line(line)
   call read_OrbInput(Input%OrbInput_L(i)%OrbInput(1))

enddo

! hapka:  old_drake
!   allocate(Input%OrbInput(1))
!
!   call read_line(line)
!   call read_OrbInput(Input%OrbInput(1))

case('SEPARATE')

   allocate(Input%OrbInput_L(Input%maxl+1))
do i=1,size(Input%OrbInput_L)
   allocate(Input%OrbInput_L(i)%OrbInput(Input%norb(i)))
   call mem_alloc(isOrb,Input%norb(i))
   isOrb = .false.
   icnt = 0
   do
      call read_line(line)

      call read_keyword_and_advance_line(line,'ORB','=')
      call read_value_and_advance_line(line,iorb,',',.false.)
      if(iorb<1.or.iorb>Input%norb(i)) then
         write(LOUT,'(a)') 'Incorrect orbital number!'
         call incorrect_data
      endif

      if(isOrb(iorb)) then
         write(LOUT,'(a,i3,a)') 'Orbital ',iorb,' has already been defined!'
         call incorrect_data
      else
         call read_OrbInput(Input%OrbInput_L(i)%OrbInput(iorb))
         isOrb(iorb) = .true.
         icnt = icnt + 1
      endif

      if(icnt==Input%norb(i)) then
         if(.not.all(isOrb)) then
            write(LOUT,'(a)') 'ERROR!!! Not enough orbitals read from input!'
            stop
         endif
         exit
      elseif(icnt>Input%norb(i)) then
         write(LOUT,'(a)') 'ERROR!!! Too many orbitals read from input!'
         stop
      endif
   enddo
   call mem_dealloc(isOrb)

enddo
! hapka: old_drake
!   allocate(Input%OrbInput(Input%norb))
!
!   call mem_alloc(isOrb,Input%norb)
!   isOrb = .false.
!   icnt = 0
!   do
!      call read_line(line)
!
!      call read_keyword_and_advance_line(line,'ORB','=')
!      call read_value_and_advance_line(line,iorb,',',.false.)
!      if(iorb<1.or.iorb>Input%norb) then
!         write(LOUT,'(a)') 'Incorrect orbital number!'
!         call incorrect_data
!      endif
!
!      if(isOrb(iorb)) then
!         write(LOUT,'(a,i3,a)') 'Orbital ',iorb,' has already been defined!'
!         call incorrect_data
!      else
!         call read_OrbInput(Input%OrbInput(iorb))
!         isOrb(iorb) = .true.
!         icnt = icnt + 1
!      endif
!
!      if(icnt==Input%norb) then
!         if(.not.all(isOrb)) then
!            write(LOUT,'(a)') 'ERROR!!! Not enough orbitals read from input!'
!            stop
!         endif
!         exit
!      elseif(icnt>Input%norb) then
!         write(LOUT,'(a)') 'ERROR!!! Too many orbitals read from input!'
!         stop
!      endif
!   enddo
!   call mem_dealloc(isOrb)

case default

   write(LOUT,'(a)') 'Unrecognizable type of orbital basis set: &
        &reading inputfile!'
   stop

end select

call read_line(line)
! hapka: old_drake
!call read_keyword_and_advance_line(line,'END',end_flag)
call get_keyword_and_advance_line(line,keyword,'=')
if(trim(keyword).eq.'LORB') then
   write(LOUT,'(a)') 'Too many orbital blocks!'
   call incorrect_data  
endif

! hapka: post HF part (to be done later)
!if(index(Input%calc_type,'SCF')==0) then
   
 !  call read_line(line)
 !  call read_keyword_and_advance_line(line,'PAIRS','=')
 !  call read_value_and_advance_line(line,Input%pairs_basis,end_flag)
 !  if(all(trim(Input%pairs_basis)/=possible_basis)) then
 !     write(LOUT,'(a)') 'Unrecognizable type of pairs basis set! Choose one of:'
 !     write(LOUT,'(a)') possible_basis
 !     call incorrect_data
 !  endif

!   select case(trim(Input%pairs_basis))
!   case('COMMON')
!
!      allocate(Input%PairInput(1,1))
!
!      Input%PairInput(1,1)%mult = 0
!      call read_line(line)
!      call read_PairInput(Input%PairInput(1,1))
!
!   case('SEPARATE')
!
!      allocate(Input%PairInput(Input%norb,Input%norb))
!
!      call mem_alloc(isPair,Input%norb,Input%norb)
!      isPair = .false.
!      icnt = 0
!      do
!         call read_line(line)
!
!         call read_keyword_and_advance_line(line,'PAIR','=')
!         call read_value_and_advance_line(line,ipair1,',',.false.)
!         call read_value_and_advance_line(line,ipair2,',',.false.)
!         if((ipair1<1.or.ipair1>Input%norb).or. &
!              (ipair2<1.or.ipair2>Input%norb)) then
!            write(LOUT,'(a)') 'Incorrect pair number!'
!            call incorrect_data
!         endif
!         if(ipair1>ipair2) call swap(ipair1,ipair2)
!
!         if(check_keyword(line,'MULT','=')) then
!
!            call read_keyword_and_advance_line(line,'MULT','=')
!            call read_value_and_advance_line(line,smult,',',.false.)
!            do imult=1,3
!               if(trim(smult)==possible_mult(imult)) exit
!            enddo
!            if(imult>3) then
!               write(LOUT,'(a)') 'Unrecognizable multiplicity! Choose one of:'
!               write(LOUT,'(a)') possible_mult(1:3)
!               call incorrect_data
!            endif
!
!            if(ipair1==ipair2.and.imult/=1) then
!               write(LOUT,'(2a)') 'This pair can only be ',possible_mult(1)
!               call incorrect_data
!            endif
!
!         else
!
!            if(ipair1/=ipair2) then
!               write(LOUT,'(a)') &
!                    'Multiplicity must be specified for this pair! &
!                    &Choose one of:'
!               write(LOUT,'(a)') possible_mult(1:3)
!               call incorrect_data
!            endif
!
!            imult = 1
!
!         endif
!
!         select case(imult)
!         case(1)
!
!            if(isPair(ipair2,ipair1)) then
!               write(LOUT,'(a,2i3,3a)') &
!                    'Pair ',ipair1,ipair2,', ',trim(possible_mult(1)),&
!                    ' has already been defined!'
!               call incorrect_data
!            else
!               Input%PairInput(ipair2,ipair1)%mult = 1
!               call read_PairInput(Input%PairInput(ipair2,ipair1))
!               isPair(ipair2,ipair1) = .true.
!               icnt = icnt + 1
!            endif
!
!         case(2)
!
!            if(isPair(ipair2,ipair1).or.isPair(ipair1,ipair2)) then
!               if(isPair(ipair2,ipair1)) then
!                  write(LOUT,'(a,2i3,3a)') &
!                       'Pair ',ipair1,ipair2,', ',trim(possible_mult(1)),&
!                       ' has already been defined!'
!               endif
!               if(isPair(ipair1,ipair2)) then
!                  write(LOUT,'(a,2i3,3a)') &
!                       'Pair ',ipair1,ipair2,', ',trim(possible_mult(3)),&
!                       ' has already been defined!'
!               endif
!               call incorrect_data
!            else
!               Input%PairInput(ipair2,ipair1)%mult = 2
!               Input%PairInput(ipair1,ipair2)%mult = -1
!               call read_PairInput(Input%PairInput(ipair2,ipair1))
!               call fake_PairInput(Input%PairInput(ipair1,ipair2))
!               isPair(ipair2,ipair1) = .true.
!               isPair(ipair1,ipair2) = .true.
!               icnt = icnt + 2
!            endif
!               
!         case(3)
!
!            if(isPair(ipair1,ipair2)) then
!               write(LOUT,'(a,2i3,3a)') &
!                    'Pair ',ipair1,ipair2,', ',trim(possible_mult(3)),&
!                    ' has already been defined!'
!               call incorrect_data
!            else
!               Input%PairInput(ipair1,ipair2)%mult = 3
!               call read_PairInput(Input%PairInput(ipair1,ipair2))
!               isPair(ipair1,ipair2) = .true.
!               icnt = icnt + 1
!            endif
!
!         end select
!
!         if(icnt==Input%norb**2) then
!            if(.not.all(isPair)) then
!               write(LOUT,'(a)') 'ERROR!!! Not enough pairs read from input!'
!               stop
!            endif
!            exit
!         elseif(icnt>Input%norb**2) then
!            write(LOUT,'(a)') 'ERROR!!! Too many pairs read from input!'
!            stop
!         endif
!      enddo
!      call mem_dealloc(isPair)
!
!   case default
!
!      write(LOUT,'(a)') 'Unrecognizable type of pairs basis set: &
!           &reading inputfile!'
!      stop
!
!   end select
!
!   call read_line(line)
!   call read_keyword_and_advance_line(line,'END',end_flag)
!
!else
!
!   call read_line(line)
!   if(check_keyword(line,'PAIRS','=')) then
!
!      do
!         call read_line(line)
!         if(check_keyword(line,'END',end_flag)) exit
!      enddo
!
!   else
!
!      call cancel_line
!
!   endif
!
!endif

do
   call read_line(line)
   call get_keyword_and_advance_line(line,keyword,'=')
   select case(trim(keyword))
   case('END')
      exit
   case('INT3')
      if(Input%set_INT3) then
         write(LOUT,'(a)') 'INT3 has already been set!'
         call incorrect_data
      else
         Input%INT3 = line(1:1)
         if(all(Input%INT3/=possible_int3)) then
            write(LOUT,'(a)') 'Unrecognizable source of 3-el integrals! &
                 &Choose one of:'
            write(LOUT,'(a)') possible_int3
            call incorrect_data
         endif
         Input%set_INT3 = .true.
      endif
   case('PRINT')
      if(Input%set_LPRINT) then
         write(LOUT,'(a)') 'PRINT has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%LPRINT,end_flag)
         Input%set_LPRINT = .true.
      endif
   case('SIMTHR')
      if(Input%set_SIMTHR) then
         write(LOUT,'(a)') 'SIMTHR has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%SIMTHR,end_flag)
         Input%set_SIMTHR = .true.
      endif
   case('REDTHR')
      if(Input%set_REDTHR) then
         write(LOUT,'(a)') 'REDTHR has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%REDTHR,end_flag)
         Input%set_REDTHR = .true.
      endif
   case('SCFTHR')
      if(Input%set_SCFTHR) then
         write(LOUT,'(a)') 'SCFTHR has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%SCFTHR,end_flag)
         Input%set_SCFTHR = .true.
      endif
   case('SCFCNT')
      if(Input%set_SCFCNT) then
         write(LOUT,'(a)') 'SCFCNT has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%SCFCNT,end_flag)
         Input%set_SCFCNT = .true.
      endif
   case('SCFMAXIT')
      if(Input%set_SCFMAXIT) then
         write(LOUT,'(a)') 'SCFMAXIT has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%SCFMAXIT,end_flag)
         Input%set_SCFMAXIT = .true.
      endif
   case('SCFSHIFTVAL')
      if(Input%set_SCFSHIFTVAL) then
         write(LOUT,'(a)') 'SCFSHIFTVAL has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%SCFSHIFTVAL,end_flag)
         Input%set_SCFSHIFTVAL = .true.
      endif
   case('SCFSHIFTBRK')
      if(Input%set_SCFSHIFTBRK) then
         write(LOUT,'(a)') 'SCFSHIFTBRK has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%SCFSHIFTBRK,end_flag)
         Input%set_SCFSHIFTBRK = .true.
      endif
   case('ENERGTHR')
      if(Input%set_ENERGTHR) then
         write(LOUT,'(a)') 'ENERGTHR has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%ENERGTHR,end_flag)
         Input%set_ENERGTHR = .true.
      endif
   case('ENERGMAXIT')
      if(Input%set_ENERGMAXIT) then
         write(LOUT,'(a)') 'ENERGMAXIT has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%ENERGMAXIT,end_flag)
         Input%set_ENERGMAXIT = .true.
      endif
   case('DIVERTHR')
      if(Input%set_DIVERTHR) then
         write(LOUT,'(a)') 'DIVERTHR has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%DIVERTHR,end_flag)
         Input%set_DIVERTHR = .true.
      endif
   case('OPTTHRMLT')
      if(Input%set_OPTTHRMLT) then
         write(LOUT,'(a)') 'OPTTHRMLT has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%OPTTHRMLT,end_flag)
         Input%set_OPTTHRMLT = .true.
      endif
   case('OPTSTEPMLT')
      if(Input%set_OPTSTEPMLT) then
         write(LOUT,'(a)') 'OPTSTEPMLT has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%OPTSTEPMLT,end_flag)
         Input%set_OPTSTEPMLT = .true.
      endif
   case('OPTLINMLT')
      if(Input%set_OPTLINMLT) then
         write(LOUT,'(a)') 'OPTLINMLT has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%OPTLINMLT,end_flag)
         Input%set_OPTLINMLT = .true.
      endif
   case('OPTNORMTHR')
      if(Input%set_OPTNORMTHR) then
         write(LOUT,'(a)') 'OPTNORMTHR has already been set!'
         call incorrect_data
      else
         call read_value_and_advance_line(line,Input%OPTNORMTHR,end_flag)
         Input%set_OPTNORMTHR = .true.
      endif
   case default
      write(LOUT,'(a)') 'Unknown keyword!'
      call incorrect_data
   end select
enddo

contains

  subroutine read_OrbInput(OrbInput)
  implicit none
  type(OrbInputData) :: OrbInput
  integer :: i,iexp,irange

  call read_keyword_and_advance_line(line,'PRIMITIVES','=')
  call read_value_and_advance_line(line,OrbInput%Nprimitives,end_flag)
  if(OrbInput%Nprimitives<0) then
     write(LOUT,'(a)') 'Incorrect number of orbital primitives!'
     call incorrect_data
  endif

  call mem_alloc(OrbInput%orb_control,2,OrbInput%Nprimitives)
  OrbInput%orb_control = 0

  do i=1,OrbInput%Nprimitives
     call read_line(line)

     call read_value_and_advance_line(line,iexp,'|')
     if(iexp<1.or.iexp>Input%Nexponents) then
        write(LOUT,'(a)') 'The exponent does not exist!'
        call incorrect_data
     endif

     call read_value_and_advance_line(line,irange,end_flag)
     if(irange<0) then
        write(LOUT,'(a)') 'Incorrect r-range in the definition of orbital!'
        call incorrect_data
     endif

     OrbInput%orb_control(1,i) = iexp
     OrbInput%orb_control(2,i) = irange
  enddo

  end subroutine read_OrbInput

!  subroutine fake_PairInput(PairInput)
!  implicit none
!  type(PairInputData) :: PairInput
!
!  PairInput%Nprimitives = 0
!
!  call mem_alloc(PairInput%pair_control,6,PairInput%Nprimitives)
!
!  end subroutine fake_PairInput
!
!  subroutine read_PairInput(PairInput)
!  implicit none
!  type(PairInputData) :: PairInput
!  integer :: iexp(2),itype,irange(3)
!  integer :: i,j
!
!  call read_keyword_and_advance_line(line,'PRIMITIVES','=')
!  call read_value_and_advance_line(line,PairInput%Nprimitives,end_flag)
!  if(PairInput%Nprimitives<0) then
!     write(LOUT,'(a)') 'Incorrect number of pair primitives!'
!     call incorrect_data
!  endif
!
!  call mem_alloc(PairInput%pair_control,6,PairInput%Nprimitives)
!  PairInput%pair_control = 0
!
!  do i=1,PairInput%Nprimitives
!     call read_line(line)
!
!     call read_value_and_advance_line(line,iexp(1),',',.false.)
!     call read_value_and_advance_line(line,iexp(2),'|')
!     if(any(iexp<1).or.any(iexp>Input%Nexponents)) then
!        write(LOUT,'(a)') 'At least one of the exponents does not exist!'
!        call incorrect_data
!     endif
!
!     call read_value_and_advance_line(line,itype,':')
!     if(itype<1.or.itype>3) then
!        write(LOUT,'(a)') &
!             'Incorrect range type specification in the definition of pair!'
!        call incorrect_data
!     endif
!
!     do j=1,itype-1
!        call read_value_and_advance_line(line,irange(j),',',.false.)
!     enddo
!     call read_value_and_advance_line(line,irange(itype),end_flag)
!     if(any(irange(1:itype)<0)) then
!        write(LOUT,'(a)') 'Incorrect r-range in the definition of pair!'
!        call incorrect_data
!     endif
!
!     if(iexp(1)>iexp(2)) then
!        call swap(iexp(1),iexp(2))
!        if(itype>2) call swap(irange(1),irange(2))
!     endif
!
!     PairInput%pair_control(1:2,i) = iexp
!     PairInput%pair_control(3,i) = itype
!     PairInput%pair_control(4:3+itype,i) = irange(1:itype)
!  enddo
!
!  end subroutine read_PairInput

end subroutine read_InputFile

subroutine read_line(line,doUP)
use misc, only : toUPPER
implicit none
character(slength) :: line
logical,optional :: doUP
logical :: transform
integer :: ios,comment_pos

if(present(doUP)) then
   transform = doUP
else
   transform = .true.
endif

do
   read(LIN,sfmt,iostat=ios) line
   if(ios<0) line = 'END'

   line = adjustl(line)
   comment_pos = scan(line,set_comment)
   if(comment_pos/=0) line = line(1:comment_pos-1)

   if(len_trim(line)/=0) exit
enddo

if(transform) call toUPPER(line)

end subroutine read_line

subroutine cancel_line
implicit none

backspace(LIN)

end subroutine cancel_line

function check_keyword(line,keyword,separator) result(is_present)
implicit none
logical :: is_present
character(slength) :: line
character(*),intent(in) :: keyword
character(1),intent(in) :: separator
integer :: pos,ios

pos = scan(line,' '//separator)
ios = len(keyword)

is_present = (line(1:ios)==keyword.and.pos==ios+1)

end function check_keyword

subroutine read_keyword_and_advance_line(line,keyword,separator,mandatory)
implicit none
character(slength) :: line
character(*),intent(in) :: keyword
character(1),intent(in) :: separator
logical,intent(in),optional :: mandatory
integer :: pos,ios

pos = scan(line,' '//separator)
ios = len(keyword)

if(.not.(line(1:ios)==keyword.and.pos==ios+1)) then
   write(LOUT,'(3a)') &
        'Mandatory keyword ',keyword,' not found in the inputfile!'
   call incorrect_data
endif

line = line(pos:)
call advance_line(line,separator,mandatory)

end subroutine read_keyword_and_advance_line

subroutine get_keyword_and_advance_line(line,keyword,separator,mandatory)
implicit none
character(slength) :: line
character(*),intent(out) :: keyword
character(1),intent(in) :: separator
logical,intent(in),optional :: mandatory
integer :: pos

pos = scan(line,' '//separator)

keyword = line(1:pos-1)

line = line(pos:)
if(trim(keyword)=='END') then
   call advance_line(line,end_flag,mandatory)
else
   call advance_line(line,separator,mandatory)
endif

end subroutine get_keyword_and_advance_line

subroutine read_int_and_advance_line(line,val,separator,mandatory)
implicit none
character(slength) :: line
integer,intent(out) :: val
character(1),intent(in) :: separator
logical,intent(in),optional :: mandatory
integer :: pos,ios

pos = scan(line,' '//separator)
ios = verify(line,set_int)
if(ios/=pos) call incorrect_data

read(line(1:pos-1),*,iostat=ios) val
if(ios/=0) call incorrect_data

line = line(pos:)
call advance_line(line,separator,mandatory)

end subroutine read_int_and_advance_line

subroutine read_real_and_advance_line(line,val,separator,mandatory)
implicit none
character(slength) :: line
real(prec),intent(out) :: val
character(1),intent(in) :: separator
logical,intent(in),optional :: mandatory
integer :: pos,ios

pos = scan(line,' '//separator)
ios = verify(line,set_real)
if(ios/=pos) call incorrect_data

read(line(1:pos-1),*,iostat=ios) val
if(ios/=0) call incorrect_data

line = line(pos:)
call advance_line(line,separator,mandatory)

end subroutine read_real_and_advance_line

subroutine read_string_and_advance_line(line,val,separator,mandatory)
implicit none
character(slength) :: line
character(*),intent(out) :: val
character(1),intent(in) :: separator
logical,intent(in),optional :: mandatory
integer :: pos

pos = scan(line,' '//separator)

val = line(1:min(pos-1,len(val)))
pos = len_trim(val) + 1

line = line(pos:)
call advance_line(line,separator,mandatory)

end subroutine read_string_and_advance_line

subroutine advance_line(line,separator,mandatory_option)
implicit none
character(slength) :: line
character(1),intent(in) :: separator
logical,intent(in),optional :: mandatory_option
logical :: mandatory
integer :: pos
integer :: i,cnt

if(separator==end_flag) then

   if(len_trim(line)/=0) call incorrect_data

else

   if(present(mandatory_option)) then
      mandatory = mandatory_option
   else
      mandatory = .true.
   endif

   pos = verify(line,' '//separator)
   if(pos==0.or.pos==1) call incorrect_data

   cnt = 0
   do i=1,pos-1
      if(line(i:i)==separator) cnt = cnt + 1
   enddo

   if(mandatory) then
      if(cnt/=1) call incorrect_data
   else
      if(cnt>1) call incorrect_data
   endif

   line = line(pos:)

endif

end subroutine advance_line

subroutine read_log_random(n,val)
implicit none
integer,intent(in) :: n
logical,intent(out) :: val(:)
character(slength) :: line
integer :: icnt,pos

icnt = 0
do
   call read_line(line)
   call prepare_line(line)

   do while(len_trim(line)>0)
      icnt = icnt + 1
      if(icnt>n) then
         write(LOUT,'(a,i5)') 'Too many values given! Expected: ',n
         call incorrect_data
      endif
      read(line,*) val(icnt)
      pos = scan(line,' ')
      line = adjustl(line(pos:))
   enddo

   if(icnt==n) exit
enddo

end subroutine read_log_random

subroutine read_int_random(n,val)
implicit none
integer,intent(in) :: n
integer,intent(out) :: val(:)
character(slength) :: line
integer :: icnt,pos

icnt = 0
do
   call read_line(line)
   call prepare_line(line)
   do while(len_trim(line)>0)
      icnt = icnt + 1
      if(icnt>n) then
         write(LOUT,'(a,i5)') 'Too many values given! Expected: ',n
         call incorrect_data
      endif
      read(line,*) val(icnt)
      pos = verify(line,set_int)
      line = adjustl(line(pos:))
   enddo

   if(icnt==n) exit
enddo

end subroutine read_int_random

subroutine read_real_random(n,val)
implicit none
integer,intent(in) :: n
real(prec),intent(out) :: val(:)
character(slength) :: line
integer :: icnt,pos

icnt = 0
do
   call read_line(line)
   call prepare_line(line)

   do while(len_trim(line)>0)
      icnt = icnt + 1
      if(icnt>n) then
         write(LOUT,'(a,i5)') 'Too many values given! Expected: ',n
         call incorrect_data
      endif
      read(line,*) val(icnt)
      pos = verify(line,set_real)
      line = adjustl(line(pos:))
   enddo

   if(icnt==n) exit
enddo

end subroutine read_real_random

subroutine prepare_line(line)
implicit none
character(slength) :: line
integer :: pos

pos = scan(line,':')
if(pos/=0) line = adjustl(line(pos+1:))
do
   pos = scan(line,',')
   if(pos==0) exit
   line(pos:pos) = ' '
enddo

end subroutine prepare_line

subroutine read_norbs(line,n,val,icnt)
implicit none
integer,intent(in) :: n
integer            :: val(:)
character(slength) :: line
integer,intent(out):: icnt
integer :: pos

icnt = 0
do
   do while(len_trim(line)>0)
      icnt = icnt + 1
      if(icnt>n) then
         write(LOUT,'(a,i5)') 'Too many orbitals given! Expected: ',n
         call incorrect_data
      endif
      read(line,*) val(icnt)
      pos = verify(line,set_int)
      line = adjustl(line(pos:))
   enddo

   if(icnt==n.or.icnt.le.n) exit
enddo

end subroutine read_norbs

subroutine incorrect_data
implicit none
character(slength) :: line

backspace(LIN)
read(LIN,sfmt) line

write(LOUT,'(a)') 'ERROR!!! &
     &Incorrect format or inconsistent data in the input file.'
write(LOUT,'()')
write(LOUT,'(a)') trim(line)

stop

end subroutine incorrect_data

end module inputread
