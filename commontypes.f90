module commontypes
use file_OUT, only : LOUT
use precision, only : prec
use memory
implicit none

private
public possible_calc_type,possible_basis,possible_mult,possible_int3
public print_EndSection
public InputData,OrbInputData,PairInputData,init_Input,free_Input,print_Input
public print_orb_control,print_pair_control,compare_prim_orbs,compare_prim_pairs
public ControlData,fill_Control
public SystemData,init_System,free_System,print_System,print_SystemExponents
public OrbSystemData,OrbSpecData,OrbReducedData
public PairSystemData,PairSpecData,PairReducedData
public check_PairSystem,check_consistency_PairReduced
public increment_PairSpec
public maxOmega_PairSpec,maxt_PairSpec,maxuv_PairSpec
public maxOmega_PairReduced,maxt_PairReduced,maxuv_PairReduced
public TripletData,init_Triplet,free_Triplet,shrink_Triplet,expand_Triplet
public SCForbitalsData,init_SCForbitals,free_SCForbitals,print_SCForbitals

interface maxuv_PairReduced
module procedure maxuv_PairReduced_not
module procedure maxuv_PairReduced_t
end interface maxuv_PairReduced

interface shrink_Triplet
module procedure shrink_Triplet_mat
module procedure shrink_Triplet_mat_dim
module procedure shrink_Triplet_vec
end interface shrink_Triplet

interface expand_Triplet
module procedure expand_Triplet_mat
module procedure expand_Triplet_vec
end interface expand_Triplet

character(*),parameter :: possible_calc_type(*) = &
     [character(8) :: &
     'SCFOPT', 'SCF', &
     'MP2OPT', 'MP2', 'IEPAOPT','IEPA', &
     'CID'  , 'LCCD'  , 'FCCD'  , &
     'CID-L', 'LCCD-L', 'FCCD-L']

character(*),parameter :: possible_basis(2) = &
     [character(8) :: 'COMMON', 'SEPARATE']

character(*),parameter :: possible_mult(-1:3) = &
     [character(8) :: 'UNDEF', 'ALL', 'SINGLET', 'BOTH', 'TRIPLET']

character,parameter :: possible_int3(*) = ['L','P','Q']

type OrbInputData
integer :: Nprimitives
integer,allocatable :: orb_control(:,:)
end type OrbInputData

type PairInputData
integer :: mult
integer :: Nprimitives
integer,allocatable :: pair_control(:,:)
end type PairInputData

type InputData
real(prec) :: nucZ
integer :: norb
character(8) :: calc_type
integer :: Neta
integer,allocatable :: eta(:)
integer :: Nexponents
real(prec),allocatable :: exponents(:)
logical,allocatable :: optimized(:)
character(8) :: orbs_basis
type(OrbInputData),allocatable :: OrbInput(:)
character(8) :: pairs_basis
type(PairInputData),allocatable :: PairInput(:,:)
character  :: INT3;        logical :: set_INT3
integer    :: LPRINT;      logical :: set_LPRINT
real(prec) :: SIMTHR;      logical :: set_SIMTHR
real(prec) :: REDTHR;      logical :: set_REDTHR
real(prec) :: SCFTHR;      logical :: set_SCFTHR
integer    :: SCFCNT;      logical :: set_SCFCNT
integer    :: SCFMAXIT;    logical :: set_SCFMAXIT
real(prec) :: SCFSHIFTVAL; logical :: set_SCFSHIFTVAL
integer    :: SCFSHIFTBRK; logical :: set_SCFSHIFTBRK
real(prec) :: ENERGTHR;    logical :: set_ENERGTHR
integer    :: ENERGMAXIT;  logical :: set_ENERGMAXIT
real(prec) :: DIVERTHR;    logical :: set_DIVERTHR
real(prec) :: OPTTHRMLT;   logical :: set_OPTTHRMLT
real(prec) :: OPTSTEPMLT;  logical :: set_OPTSTEPMLT
real(prec) :: OPTLINMLT;   logical :: set_OPTLINMLT
real(prec) :: OPTNORMTHR;  logical :: set_OPTNORMTHR
end type InputData

type ControlData
character  :: INT3
integer    :: LPRINT
real(prec) :: SIMTHR
real(prec) :: REDTHR
real(prec) :: SCFTHR
integer    :: SCFCNT
integer    :: SCFMAXIT
real(prec) :: SCFSHIFTVAL
integer    :: SCFSHIFTBRK
real(prec) :: ENERGTHR
integer    :: ENERGMAXIT
real(prec) :: DIVERTHR
real(prec) :: OPTTHR
real(prec) :: OPTSTEPMLT
real(prec) :: OPTLINMLT
real(prec) :: OPTNORMTHR
end type ControlData

type OrbSpecData
integer :: iexp
integer :: irange
integer :: nbas
integer :: offset
end type OrbSpecData

type OrbSystemData
integer :: nbas
integer :: n_prim
type(OrbSpecData),allocatable :: OrbSpec(:)
end type OrbSystemData

type OrbReducedData
logical :: isUsed
integer :: iexp
integer :: maxrange
end type OrbReducedData

type PairSpecData
integer :: iexp1,iexp2,iexpSQ
integer :: itype
logical :: restrict
integer :: n_range
integer :: irange1
integer,allocatable :: irange23(:,:)
integer :: nbas
integer :: offset
end type PairSpecData

type PairSystemData
integer :: mult
integer :: nbas
integer :: n_prim
type(PairSpecData),allocatable :: PairSpec(:)
end type PairSystemData

type PairReducedData
logical :: isUsed
integer :: iexp1,iexp2,iexpSQ
logical :: istype1,istype2,istype3
integer :: maxrange1
integer :: n_range2
integer,allocatable :: maxrange2(:,:)
integer :: n_range3
integer,allocatable :: maxrange3(:,:)
end type PairReducedData

type SystemData
real(prec) :: nucZ
integer :: norb
character(8) :: calc_type
logical :: post_SCF,optimize
logical :: common_orbs,common_pairs
integer :: n_orbs,n_pairs
integer :: Neta
integer,allocatable :: eta(:)
integer :: Nexp,NexpSQ
real(prec),allocatable :: exponents(:)
logical,allocatable :: isUsed_orbs(:),isUsed_pairs(:)
logical,allocatable :: isOpt_orbs(:) ,isOpt_pairs(:)
type(OrbSystemData),allocatable :: OrbSystem(:)
type(OrbReducedData),allocatable :: OrbReduced(:)
type(PairSystemData),allocatable :: PairSystem(:,:)
type(PairReducedData),allocatable :: PairReduced(:),OrbPairReduced(:)
end type SystemData

type TripletData
logical :: restrict
integer :: nbas,nbas_orig
integer,allocatable :: idxS(:),idxE(:)
end type TripletData

type SCForbitalsData
integer :: norb,nbas
real(prec) :: energy
real(prec) :: HOMO,LUMO
real(prec),allocatable :: orb_energy(:)
real(prec),allocatable :: orb_vector(:,:)
real(prec),allocatable :: pair_energy(:,:,:)
end type SCForbitalsData

contains

subroutine print_EndSection
implicit none
integer :: i

write(LOUT,'()')
write(LOUT,'(8a10)') ('----------',i=1,8)
write(LOUT,'(8a10)') ('**********',i=1,8)
write(LOUT,'(8a10)') ('----------',i=1,8)

end subroutine print_EndSection

subroutine init_Input(Input)
implicit none
type(InputData) :: Input

Input%nucZ        = 0._prec
Input%norb        = 0
Input%calc_type   = ''
Input%Neta        = 0
Input%Nexponents  = 0
Input%orbs_basis  = ''
Input%pairs_basis = ''

Input%INT3        = 'Q';           Input%set_INT3        = .false.
Input%LPRINT      = 0;             Input%set_LPRINT      = .false.
Input%SIMTHR      = 1.e-6_prec;    Input%set_SIMTHR      = .false.
Input%REDTHR      = 1.e-1_prec;    Input%set_REDTHR      = .false.
Input%SCFTHR      = 1.e-12_prec;   Input%set_SCFTHR      = .false.
Input%SCFCNT      = 6;             Input%set_SCFCNT      = .false.
Input%SCFMAXIT    = 200;           Input%set_SCFMAXIT    = .false.
Input%SCFSHIFTVAL = -1._prec;      Input%set_SCFSHIFTVAL = .false.
Input%SCFSHIFTBRK = -1;            Input%set_SCFSHIFTBRK = .false.
Input%ENERGTHR    = 1.e-10_prec;   Input%set_ENERGTHR    = .false.
Input%ENERGMAXIT  = 40;            Input%set_ENERGMAXIT  = .false.
Input%DIVERTHR    = 1.e1_prec;     Input%set_DIVERTHR    = .false.
Input%OPTTHRMLT   = 1.e-1_prec;    Input%set_OPTTHRMLT   = .false.
Input%OPTSTEPMLT  = 1.e-2_prec;    Input%set_OPTSTEPMLT  = .false.
Input%OPTLINMLT   = 2._prec;       Input%set_OPTLINMLT   = .false.
Input%OPTNORMTHR  = 1.e-7_prec;    Input%set_OPTNORMTHR  = .false.

end subroutine init_Input

subroutine free_Input(Input)
implicit none
type(InputData) :: Input
integer :: i,j

if(allocated(Input%PairInput)) then

   do j=1,size(Input%PairInput,2)
      do i=1,size(Input%PairInput,1)
         call mem_dealloc(Input%PairInput(i,j)%pair_control)
      enddo
   enddo
   deallocate(Input%PairInput)

endif

if(allocated(Input%OrbInput)) then

   do i=1,size(Input%OrbInput)
      call mem_dealloc(Input%OrbInput(i)%orb_control)
   enddo
   deallocate(Input%OrbInput)

endif

call mem_dealloc(Input%optimized)
call mem_dealloc(Input%exponents)
call mem_dealloc(Input%eta)

end subroutine free_Input

subroutine print_Input(Input)
use misc, only : swap
implicit none
type(InputData) :: Input
integer :: i,j,idx1,idx2

write(LOUT,'()')
write(LOUT,'(a,f7.3)') 'Z    = ',Input%nucZ
write(LOUT,'(a,i3)') 'NORB = ',Input%norb
write(LOUT,'()')
write(LOUT,'(2a)') 'CALC = ',Input%calc_type
if(Input%Neta>0) then
   write(LOUT,'()')
   write(LOUT,'(a,i3)') 'ETA = ',Input%Neta
   do i=1,Input%Neta
      write(LOUT,'(i3,a,i5)') i,' : ',Input%eta(i)
   enddo
endif
write(LOUT,'()')
write(LOUT,'(a,i3)') 'EXPONENTS = ',Input%Nexponents
do i=1,Input%Nexponents
   write(LOUT,'(i3,a,f16.8)') i,' : ',Input%exponents(i)
enddo
if(index(Input%calc_type,'OPT')/=0) then
   write(LOUT,'(a)') 'OPTIMIZE'
   do i=1,Input%Nexponents
      write(LOUT,'(i3,a,l3)') i,' : ',Input%optimized(i)
   enddo
endif

write(LOUT,'()')
write(LOUT,'(2a)') 'ORBS = ',Input%orbs_basis
select case(trim(Input%orbs_basis))
case('COMMON')

   associate(OrbInput => Input%OrbInput(1))
     call print_orb_control(OrbInput%Nprimitives,OrbInput%orb_control)
   end associate

case('SEPARATE')

   do i=1,Input%norb

      associate(OrbInput => Input%OrbInput(i))
        write(LOUT,'(a,i3,a)',advance='no') 'ORB = ',i,', '
        call print_orb_control(OrbInput%Nprimitives,OrbInput%orb_control)
      end associate

   enddo

case default

   write(LOUT,'(a)') 'ERROR!!! &
        &Unrecognizable type of orbital basis set: printing read inputfile!'
   stop

end select

write(LOUT,'(a)') 'END'

if(index(Input%calc_type,'SCF')==0) then

   write(LOUT,'()')
   write(LOUT,'(2a)') 'PAIRS = ',Input%pairs_basis
   select case(trim(Input%pairs_basis))
   case('COMMON')

      associate(PairInput => Input%PairInput(1,1))
        call print_pair_control(PairInput%Nprimitives,&
             PairInput%pair_control)
      end associate

   case('SEPARATE')

      do j=1,Input%norb
         do i=1,j

            idx1 = j
            idx2 = i
            do
               associate(PairInput => Input%PairInput(idx1,idx2))
                 if(PairInput%mult>0) then
                    write(LOUT,'(2(a,i3),3a)',advance='no') &
                         'PAIR = ',i,', ',j,&
                         ', MULT = ',trim(possible_mult(PairInput%mult)),', '
                    call print_pair_control(PairInput%Nprimitives,&
                         PairInput%pair_control)
                 elseif(PairInput%mult==0) then
                    write(LOUT,'(a)') 'ERROR!!! &
                         &Wrong multiplicity for SEPARATE pair basis set!'
                    stop
                 endif
               end associate

               if(idx1<=idx2) exit
               call swap(idx1,idx2)
            enddo

         enddo
      enddo

   case default

      write(LOUT,'(a)') 'ERROR!!! &
           &Unrecognizable type of pair basis set: printing read inputfile!'
      stop

   end select

   write(LOUT,'(a)') 'END'

endif

write(LOUT,'()')
if(Input%set_INT3) then
   write(LOUT,'(a,a8)') 'INT3       = ',Input%INT3
else
   write(LOUT,'(a,a8,t35,a)') 'INT3       = ',Input%INT3,'! default'
endif
if(Input%set_LPRINT) then
   write(LOUT,'(a,i8)') 'PRINT      = ',Input%LPRINT
else
   write(LOUT,'(a,i8,t35,a)') 'PRINT      = ',Input%LPRINT,'! default'
endif
if(Input%set_SIMTHR) then
   write(LOUT,'(a,es16.3)') 'SIMTHR     = ',Input%SIMTHR
else
   write(LOUT,'(a,es16.3,t35,a)') 'SIMTHR     = ',Input%SIMTHR,'! default'
endif
if(Input%set_REDTHR) then
   write(LOUT,'(a,es16.3)') 'REDTHR     = ',Input%REDTHR
else
   write(LOUT,'(a,es16.3,t35,a)') 'REDTHR     = ',Input%REDTHR,'! default'
endif
if(Input%set_SCFTHR) then
   write(LOUT,'(a,es16.3)') 'SCFTHR     = ',Input%SCFTHR
else
   write(LOUT,'(a,es16.3,t35,a)') 'SCFTHR     = ',Input%SCFTHR,'! default'
endif
if(Input%set_SCFCNT) then
   write(LOUT,'(a,i8)') 'SCFCNT     = ',Input%SCFCNT
else
   write(LOUT,'(a,i8,t35,a)') 'SCFCNT     = ',Input%SCFCNT,'! default'
endif
if(Input%set_SCFMAXIT) then
   write(LOUT,'(a,i8)') 'SCFMAXIT   = ',Input%SCFMAXIT
else
   write(LOUT,'(a,i8,t35,a)') 'SCFMAXIT   = ',Input%SCFMAXIT,'! default'
endif
if(Input%set_SCFSHIFTVAL) then
   write(LOUT,'(a,es16.3)') 'SCFSHIFTVAL= ',Input%SCFSHIFTVAL
else
   write(LOUT,'(a,es16.3,t35,a)') 'SCFSHIFTVAL= ',Input%SCFSHIFTVAL,'! default'
endif
if(Input%set_SCFSHIFTBRK) then
   write(LOUT,'(a,i8)') 'SCFSHIFTBRK= ',Input%SCFSHIFTBRK
else
   write(LOUT,'(a,i8,t35,a)') 'SCFSHIFTBRK= ',Input%SCFSHIFTBRK,'! default'
endif
if(Input%set_ENERGTHR) then
   write(LOUT,'(a,es16.3)') 'ENERGTHR   = ',Input%ENERGTHR
else
   write(LOUT,'(a,es16.3,t35,a)') 'ENERGTHR   = ',Input%ENERGTHR,'! default'
endif
if(Input%set_ENERGMAXIT) then
   write(LOUT,'(a,i8)') 'ENERGMAXIT = ',Input%ENERGMAXIT
else
   write(LOUT,'(a,i8,t35,a)') 'ENERGMAXIT = ',Input%ENERGMAXIT,'! default'
endif
if(Input%set_DIVERTHR) then
   write(LOUT,'(a,es16.3)') 'DIVERTHR   = ',Input%DIVERTHR
else
   write(LOUT,'(a,es16.3,t35,a)') 'DIVERTHR   = ',Input%DIVERTHR,'! default'
endif
if(Input%set_OPTTHRMLT) then
   write(LOUT,'(a,es16.3)') 'OPTTHRMLT  = ',Input%OPTTHRMLT
else
   write(LOUT,'(a,es16.3,t35,a)') 'OPTTHRMLT  = ',Input%OPTTHRMLT,'! default'
endif
if(Input%set_OPTSTEPMLT) then
   write(LOUT,'(a,es16.3)') 'OPTSTEPMLT = ',Input%OPTSTEPMLT
else
   write(LOUT,'(a,es16.3,t35,a)') 'OPTSTEPMLT = ',Input%OPTSTEPMLT,'! default'
endif
if(Input%set_OPTLINMLT) then
   write(LOUT,'(a,es16.3)') 'OPTLINMLT  = ',Input%OPTLINMLT
else
   write(LOUT,'(a,es16.3,t35,a)') 'OPTLINMLT  = ',Input%OPTLINMLT,'! default'
endif
if(Input%set_OPTNORMTHR) then
   write(LOUT,'(a,es16.3)') 'OPTNORMTHR = ',Input%OPTNORMTHR
else
   write(LOUT,'(a,es16.3,t35,a)') 'OPTNORMTHR = ',Input%OPTNORMTHR,'! default'
endif

end subroutine print_Input

subroutine print_orb_control(Nprim,orb_control)
implicit none
integer :: Nprim
integer :: orb_control(:,:)
integer :: i

write(LOUT,'(a,i3)') 'PRIMITIVES = ',Nprim
do i=1,Nprim
   write(LOUT,'(i3,a,i3)') orb_control(1,i),' |',orb_control(2,i)
enddo

end subroutine print_orb_control

subroutine print_pair_control(Nprim,pair_control)
implicit none
integer :: Nprim
integer :: pair_control(:,:)
integer :: itype
integer :: i,j

write(LOUT,'(a,i3)') 'PRIMITIVES = ',Nprim
do i=1,Nprim
   write(LOUT,'(3(i3,a))',advance='no') &
        pair_control(1,i),', ', pair_control(2,i),' |', pair_control(3,i),' :'
   itype = pair_control(3,i)
   do j=1,itype-1
      write(LOUT,'(i3,a)',advance='no') pair_control(3+j,i),', '
   enddo
   write(LOUT,'(i3)') pair_control(3+itype,i)
enddo

end subroutine print_pair_control

function compare_prim_orbs(orb1,orb2) result(diff_type)
implicit none
integer :: diff_type
integer,intent(in) :: orb1(:),orb2(:)
integer :: key1,key2

key1 = orb1(1)
key2 = orb2(1)
if(key1/=key2) then
   diff_type = sign(1,key2-key1)
   return
endif

key1 = orb1(2)
key2 = orb2(2)
if(key1/=key2) then
   diff_type = sign(22,key2-key1)
   return
endif

diff_type = 0

end function compare_prim_orbs

function compare_prim_pairs(pair1,pair2) result(diff_type)
implicit none
integer :: diff_type
integer,intent(in) :: pair1(:),pair2(:)
integer :: key1,key2
integer :: itype

if(pair1(1)>pair1(2).or.pair2(1)>pair2(2)) then
   write(LOUT,'(a)') 'ERROR!!! Incorrect exponents order in compare_prim_pairs!'
   stop
endif
key1 = pair1(1) + (pair1(2)-1)*pair1(2)/2
key2 = pair2(1) + (pair2(2)-1)*pair2(2)/2
if(key1/=key2) then
   diff_type = sign(1,key2-key1)
   return
endif

key1 = pair1(3)
key2 = pair2(3)
if(key1/=key2) then
   diff_type = sign(6,key2-key1)
   return
endif

itype = key1
select case(itype)
case(1)

   key1 = pair1(4)
   key2 = pair2(4)
   if(key1/=key2) then
      diff_type = sign(44,key2-key1)
      return
   endif

case(2)

   key1 = pair1(5)
   key2 = pair2(5)
   if(key1/=key2) then
      diff_type = sign(2,key2-key1)
      return
   endif

   key1 = pair1(4)
   key2 = pair2(4)
   if(key1/=key2) then
      diff_type = sign(44,key2-key1)
      return
   endif

case(3)

   key1 = pair1(6)
   key2 = pair2(6)
   if(key1/=key2) then
      diff_type = sign(2,key2-key1)
      return
   endif

   key1 = pair1(5)
   key2 = pair2(5)
   if(key1/=key2) then
      diff_type = sign(45,key2-key1)
      return
   endif

   key1 = pair1(4)
   key2 = pair2(4)
   if(key1/=key2) then
      diff_type = sign(45,key2-key1)
      return
   endif

case default

   write(LOUT,'(a)') &
        'ERROR!!! Incorrect range type specification in compare_prim_pairs!'
   stop

end select

diff_type = 0

end function compare_prim_pairs

subroutine fill_Control(Control,Input)
implicit none
type(ControlData) :: Control
type(InputData) :: Input

Control%INT3        = Input%INT3
Control%LPRINT      = abs(Input%LPRINT)
Control%SIMTHR      = abs(Input%SIMTHR)
Control%REDTHR      = abs(Input%REDTHR)
Control%SCFTHR      = abs(Input%SCFTHR)
Control%SCFCNT      = abs(Input%SCFCNT)
Control%SCFMAXIT    = abs(Input%SCFMAXIT)
Control%SCFSHIFTVAL = Input%SCFSHIFTVAL
Control%SCFSHIFTBRK = Input%SCFSHIFTBRK
Control%ENERGTHR    = abs(Input%ENERGTHR)
Control%ENERGMAXIT  = abs(Input%ENERGMAXIT)
Control%DIVERTHR    = abs(Input%DIVERTHR)
if(index(Input%calc_type,'OPT')/=0) then
   if(index(Input%calc_type,'SCF')==0) then
      Control%OPTTHR      = Control%ENERGTHR
      Control%ENERGTHR    = Control%ENERGTHR * abs(Input%OPTTHRMLT)
   else
      Control%OPTTHR      = Control%SCFTHR
      Control%SCFTHR      = Control%SCFTHR   * abs(Input%OPTTHRMLT)
   endif
else
   Control%OPTTHR      = huge(0._prec)
endif
Control%OPTSTEPMLT  = abs(Input%OPTSTEPMLT)
Control%OPTLINMLT   = abs(Input%OPTLINMLT)
Control%OPTNORMTHR  = abs(Input%OPTNORMTHR)

end subroutine fill_Control

subroutine init_System(System)
implicit none
type(SystemData) :: System
integer :: i,j

call mem_alloc(System%eta         ,System%Neta)
call mem_alloc(System%exponents   ,System%Nexp)
call mem_alloc(System%isUsed_orbs ,System%Nexp)
call mem_alloc(System%isUsed_pairs,System%Nexp)
call mem_alloc(System%isOpt_orbs  ,System%Nexp)
call mem_alloc(System%isOpt_pairs ,System%Nexp)

allocate(System%OrbSystem(System%n_orbs))
allocate(System%OrbReduced(System%Nexp))

do i=1,System%n_orbs
   call init_OrbSystem(System%OrbSystem(i))
enddo

do i=1,System%Nexp
   call init_OrbReduced(System%OrbReduced(i))
enddo

if(System%post_SCF) then

   allocate(System%PairSystem(System%n_pairs,System%n_pairs))
   allocate(System%PairReduced(System%NexpSQ))
   allocate(System%OrbPairReduced(System%NexpSQ))

   do j=1,System%n_pairs
      do i=1,System%n_pairs
         call init_PairSystem(System%PairSystem(i,j))
      enddo
   enddo

   do i=1,System%NexpSQ
      call init_PairReduced(System%PairReduced(i))
   enddo
   do i=1,System%NexpSQ
      call init_PairReduced(System%OrbPairReduced(i))
   enddo

endif

end subroutine init_System

subroutine free_System(System)
implicit none
type(SystemData) :: System
integer :: i,j

if(System%post_SCF) then

   do i=1,System%NexpSQ
      call free_PairReduced(System%OrbPairReduced(i))
   enddo
   do i=1,System%NexpSQ
      call free_PairReduced(System%PairReduced(i))
   enddo

   do j=1,System%n_pairs
      do i=1,System%n_pairs
         call free_PairSystem(System%PairSystem(i,j))
      enddo
   enddo

   deallocate(System%OrbPairReduced)
   deallocate(System%PairReduced)
   deallocate(System%PairSystem)

endif

do i=1,System%n_orbs
   call free_OrbSystem(System%OrbSystem(i))
enddo

deallocate(System%OrbReduced)
deallocate(System%OrbSystem)

call mem_dealloc(System%isOpt_pairs)
call mem_dealloc(System%isOpt_orbs)
call mem_dealloc(System%isUsed_pairs)
call mem_dealloc(System%isUsed_orbs)
call mem_dealloc(System%exponents)
call mem_dealloc(System%eta)

end subroutine free_System

subroutine init_OrbSystem(OrbSystem)
implicit none
type(OrbSystemData) :: OrbSystem

OrbSystem%nbas    = 0
OrbSystem%n_prim  = 0

end subroutine init_OrbSystem

subroutine free_OrbSystem(OrbSystem)
implicit none
type(OrbSystemData) :: OrbSystem

deallocate(OrbSystem%OrbSpec)

end subroutine free_OrbSystem

subroutine init_OrbReduced(OrbReduced)
implicit none
type(OrbReducedData) :: OrbReduced

OrbReduced%isUsed   = .false.
OrbReduced%iexp     = 0
OrbReduced%maxrange = -1

end subroutine init_OrbReduced

subroutine init_PairSystem(PairSystem)
implicit none
type(PairSystemData) :: PairSystem

PairSystem%mult   = -1
PairSystem%nbas   = 0
PairSystem%n_prim = 0

end subroutine init_PairSystem

subroutine free_PairSystem(PairSystem)
implicit none
type(PairSystemData) :: PairSystem
integer :: i

if(PairSystem%mult<0) return

do i=1,PairSystem%n_prim
   associate(PairSpec => PairSystem%PairSpec(i))
     if(PairSpec%itype==2.or.&
          PairSpec%itype==3) call mem_dealloc(PairSpec%irange23)
   end associate
enddo
deallocate(PairSystem%PairSpec)

end subroutine free_PairSystem

subroutine init_PairReduced(PairReduced)
implicit none
type(PairReducedData) :: PairReduced

PairReduced%isUsed    = .false.
PairReduced%iexp1     = 0
PairReduced%iexp2     = 0
PairReduced%iexpSQ    = 0
PairReduced%istype1   = .false.
PairReduced%istype2   = .false.
PairReduced%istype3   = .false.
PairReduced%maxrange1 = -1
PairReduced%n_range2  = 0
PairReduced%n_range3  = 0

end subroutine init_PairReduced

subroutine free_PairReduced(PairReduced)
implicit none
type(PairReducedData) :: PairReduced

if(PairReduced%istype2) call mem_dealloc(PairReduced%maxrange2)
if(PairReduced%istype3) call mem_dealloc(PairReduced%maxrange3)

end subroutine free_PairReduced

subroutine check_PairSystem(PairSystem)
implicit none
type(PairSystemData) :: PairSystem
integer :: i_prim
integer :: i_exact,j_exact,k_exact
integer :: ibas,ibasTOT,iend,ik
integer :: i,j,k,flag
logical :: correct

if(PairSystem%mult<0) return

ibasTOT = 0

do i_prim=1,PairSystem%n_prim
   associate(PairSpec => PairSystem%PairSpec(i_prim))

     if(ibasTOT/=PairSpec%offset) then
        write(LOUT,'(a)') 'ERROR!!! &
             &Offset for basis functions predicted incorrectly!'
        stop
     endif

     select case(PairSpec%itype)
     case(1)

        flag = -1

        ibas = 0
        do k_exact=0,PairSpec%irange1
           do j_exact=0,PairSpec%irange1-k_exact
              iend = PairSpec%irange1-k_exact-j_exact
              if(PairSpec%restrict) iend = min(iend,j_exact)
              do i_exact=0,iend
                 ibas    = ibas    + 1
                 ibasTOT = ibasTOT + 1

                 correct = increment_PairSpec(i,j,k,flag,PairSpec)
                 if(i/=i_exact.or.j/=j_exact.or.k/=k_exact.or..not.correct) then
                    write(LOUT,'(a)') 'ERROR!!! &
                         &Incorrect indices for range type 1!'
                    write(LOUT,'(a,i6,4x,a,3i5,4x,a,3i5,4x,a,l2)') &
                         'ibas=',ibas,&
                         'expected:',i_exact,j_exact,k_exact,&
                         'found:',i,j,k,&
                         'correct:',correct
                    stop
                 endif

              enddo
           enddo
        enddo
        if(ibas/=PairSpec%nbas) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Number of basis functions for range type 1 &
                &predicted incorrectly!'
           stop
        endif

        correct = increment_PairSpec(i,j,k,flag,PairSpec)
        if(correct) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Increment procedure not finished for range type 1!'
           stop
        endif

     case(2)

        flag = -1

        ibas = 0
        do ik=1,PairSpec%n_range
           k_exact = PairSpec%irange23(2,ik)
           do j_exact=0,PairSpec%irange23(1,ik)
              iend = PairSpec%irange23(1,ik)-j_exact
              if(PairSpec%restrict) iend = min(iend,j_exact)
              do i_exact=0,iend
                 ibas    = ibas    + 1
                 ibasTOT = ibasTOT + 1

                 correct = increment_PairSpec(i,j,k,flag,PairSpec)
                 if(i/=i_exact.or.j/=j_exact.or.k/=k_exact.or..not.correct) then
                    write(LOUT,'(a)') 'ERROR!!! &
                         &Incorrect indices for range type 2!'
                    write(LOUT,'(a,i6,4x,a,3i5,4x,a,3i5,4x,a,l2)') &
                         'ibas=',ibas,&
                         'expected:',i_exact,j_exact,k_exact,&
                         'found:',i,j,k,&
                         'correct:',correct
                    stop
                 endif

              enddo
           enddo
        enddo
        if(ibas/=PairSpec%nbas) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Number of basis functions for range type 2 &
                &predicted incorrectly!'
           stop
        endif

        correct = increment_PairSpec(i,j,k,flag,PairSpec)
        if(correct) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Increment procedure not finished for range type 2!'
           stop
        endif

     case(3)

        flag = -1

        ibas = 0
        do ik=1,PairSpec%n_range
           k_exact = PairSpec%irange23(3,ik)
           do j_exact=0,PairSpec%irange23(2,ik)
              iend = PairSpec%irange23(1,ik)
              if(PairSpec%restrict) iend = j_exact
              do i_exact=0,iend
                 ibas    = ibas    + 1
                 ibasTOT = ibasTOT + 1

                 correct = increment_PairSpec(i,j,k,flag,PairSpec)
                 if(i/=i_exact.or.j/=j_exact.or.k/=k_exact.or..not.correct) then
                    write(LOUT,'(a)') 'ERROR!!! &
                         &Incorrect indices for range type 3!'
                    write(LOUT,'(a,i6,4x,a,3i5,4x,a,3i5,4x,a,l2)') &
                         'ibas=',ibas,&
                         'expected:',i_exact,j_exact,k_exact,&
                         'found:',i,j,k,&
                         'correct:',correct
                    stop
                 endif

              enddo
           enddo
        enddo
        if(ibas/=PairSpec%nbas) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Number of basis functions for range type 3 &
                &predicted incorrectly!'
           stop
        endif

        correct = increment_PairSpec(i,j,k,flag,PairSpec)
        if(correct) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Increment procedure not finished for range type 3!'
           stop
        endif

     case default

        write(LOUT,'(a)') 'ERROR!!! &
             &Incorrect range type specification in check_PairSystem!'
        stop

     end select

   end associate
enddo

if(ibasTOT/=PairSystem%nbas) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Total number of basis functions predicted incorrectly!'
   stop
endif

end subroutine check_PairSystem

function increment_PairSpec(u,v,t,flag,PairSpec) result(correct)
implicit none
logical :: correct
integer,intent(inout) :: u,v,t
integer,intent(inout) :: flag
type(PairSpecData),intent(in) :: PairSpec

correct = .true.

select case(PairSpec%itype)
case(1)

   if(flag<0) then

      u = 0
      v = 0
      t = 0
      flag = 0
      return

   else

      u = u + 1
      if(PairSpec%restrict) then
         if(u+v+t<=PairSpec%irange1.and.u<=v) return
      else
         if(u+v+t<=PairSpec%irange1) return
      endif
      u = 0
      v = v + 1
      if(v+t<=PairSpec%irange1) return
      v = 0
      t = t + 1
      if(t<=PairSpec%irange1) return
      t = 0
      flag = -1

   endif

case(2)

   if(flag<0) then

      u = 0
      v = 0
      flag = 1
      t = PairSpec%irange23(2,flag)
      return

   else

      u = u + 1
      if(PairSpec%restrict) then
         if(u+v<=PairSpec%irange23(1,flag).and.u<=v) return
      else
         if(u+v<=PairSpec%irange23(1,flag)) return
      endif
      u = 0
      v = v + 1
      if(v<=PairSpec%irange23(1,flag)) return
      v = 0
      flag = flag + 1
      if(flag<=PairSpec%n_range) then
         t = PairSpec%irange23(2,flag)
         return
      endif
      t = 0
      flag = -1

   endif

case(3)

   if(flag<0) then

      u = 0
      v = 0
      flag = 1
      t = PairSpec%irange23(3,flag)
      return

   else

      u = u + 1
      if(PairSpec%restrict) then
         if(u<=v) return
      else
         if(u<=PairSpec%irange23(1,flag)) return
      endif
      u = 0
      v = v + 1
      if(v<=PairSpec%irange23(2,flag)) return
      v = 0
      flag = flag + 1
      if(flag<=PairSpec%n_range) then
         t = PairSpec%irange23(3,flag)
         return
      endif
      t = 0
      flag = -1

   endif

case default

   write(LOUT,'(a)') 'ERROR!!! &
        &Incorrect range type specification in increment_PairSystem!'
   stop

end select

correct = .false.

end function increment_PairSpec

function maxOmega_PairSpec(PairSpec) result(Omega)
implicit none
integer :: Omega
type(PairSpecData),intent(in) :: PairSpec
integer :: i

Omega = 0

select case(PairSpec%itype)
case(1)

   Omega = PairSpec%irange1

case(2)

   do i=1,PairSpec%n_range
      Omega = max(Omega,sum(PairSpec%irange23(1:2,i)))
   enddo

case(3)

   do i=1,PairSpec%n_range
      Omega = max(Omega,sum(PairSpec%irange23(1:3,i)))
   enddo

case default

   write(LOUT,'(a)') 'ERROR!!! &
        &Incorrect range type specification in maxOmega_PairSpec!'
   stop

end select

end function maxOmega_PairSpec

function maxt_PairSpec(PairSpec) result(maxt)
implicit none
integer :: maxt
type(PairSpecData),intent(in) :: PairSpec

select case(PairSpec%itype)
case(1)

   maxt = PairSpec%irange1

case(2)

   maxt = PairSpec%irange23(2,PairSpec%n_range)

case(3)

   maxt = PairSpec%irange23(3,PairSpec%n_range)

case default

   write(LOUT,'(a)') 'ERROR!!! &
        &Incorrect range type specification in maxt_PairSpec!'
   stop

end select

end function maxt_PairSpec

function maxuv_PairSpec(PairSpec) result(maxuv)
implicit none
integer :: maxuv
type(PairSpecData),intent(in) :: PairSpec
integer :: i

maxuv = 0

select case(PairSpec%itype)
case(1)

   maxuv = PairSpec%irange1

case(2)

   do i=1,PairSpec%n_range
      maxuv = max(maxuv,PairSpec%irange23(1,i))
   enddo

case(3)

   do i=1,PairSpec%n_range
      maxuv = max(maxuv,PairSpec%irange23(1,i),PairSpec%irange23(2,i))
   enddo

case default

   write(LOUT,'(a)') 'ERROR!!! &
        &Incorrect range type specification in maxuv_PairSpec!'
   stop

end select

end function maxuv_PairSpec

subroutine check_consistency_PairReduced(NexpSQ,OrbPairReduced,PairReduced)
implicit none
integer :: NexpSQ
type(PairReducedData) :: OrbPairReduced(:),PairReduced(:)
integer :: i,sOrbPair,sPair

do i=1,NexpSQ
   associate(OrbPairSpec => OrbPairReduced(i), PairSpec => PairReduced(i))

     if(.not.OrbPairSpec%isUsed.and.PairSpec%isUsed) then
        write(LOUT,'(a)') 'ERROR!!! &
             &Test 0a in check_consistency_PairReduced!'
        stop
     endif
     if(PairSpec%isUsed) then
        if(OrbPairSpec%iexpSQ/=PairSpec%iexpSQ.or.&
             OrbPairSpec%iexp1/=PairSpec%iexp1.or.&
             OrbPairSpec%iexp2/=PairSpec%iexp2) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 0b in check_consistency_PairReduced!'
           stop
        endif
     endif

     if(OrbPairSpec%istype1.neqv.PairSpec%istype1) then
        write(LOUT,'(a)') 'ERROR!!! &
             &Test 1a in check_consistency_PairReduced!'
        stop
     endif
     if(PairSpec%istype1) then
        if(OrbPairSpec%maxrange1/=PairSpec%maxrange1) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 1b in check_consistency_PairReduced!'
           stop
        endif
     endif

     if(OrbPairSpec%istype2.neqv.PairSpec%istype2) then
        write(LOUT,'(a)') 'ERROR!!! &
             &Test 2a in check_consistency_PairReduced!'
        stop
     endif
     if(PairSpec%istype2) then
        if(OrbPairSpec%n_range2/=PairSpec%n_range2) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 2b in check_consistency_PairReduced!'
           stop
        endif
        if(any(OrbPairSpec%maxrange2/=PairSpec%maxrange2)) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 2c in check_consistency_PairReduced!'
           stop
        endif
     endif

     if(.not.OrbPairSpec%istype3.and.PairSpec%istype3) then
        write(LOUT,'(a)') 'ERROR!!! &
             &Test 3a in check_consistency_PairReduced!'
        stop
     endif
     if(PairSpec%istype3) then
        sOrbPair = merge(2,1,OrbPairSpec%maxrange3(3,1)==0)
        sPair    = merge(2,1,   PairSpec%maxrange3(3,1)==0)
        if((OrbPairSpec%n_range3-sOrbPair+1)/=(PairSpec%n_range3-sPair+1)) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 3b in check_consistency_PairReduced!'
           stop
        endif
        if(any(OrbPairSpec%maxrange3(:,sOrbPair:)/= &
             PairSpec%maxrange3(:,sPair:))) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 3c in check_consistency_PairReduced!'
           stop
        endif
     elseif(OrbPairSpec%istype3) then
        if(OrbPairSpec%n_range3/=1) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 3b in check_consistency_PairReduced!'
           stop
        endif
        if(OrbPairSpec%maxrange3(3,1)/=0) then
           write(LOUT,'(a)') 'ERROR!!! &
                &Test 3c in check_consistency_PairReduced!'
           stop
        endif
     endif

   end associate
enddo

end subroutine check_consistency_PairReduced

function maxOmega_PairReduced(PairReduced) result(Omega)
implicit none
integer :: Omega
type(PairReducedData),intent(in) :: PairReduced
integer :: i

Omega = 0

if(PairReduced%istype1) &
     Omega = max(Omega,PairReduced%maxrange1)

if(PairReduced%istype2) then
   do i=1,PairReduced%n_range2
      Omega = max(Omega,sum(PairReduced%maxrange2(1:2,i)))
   enddo
endif

if(PairReduced%istype3) then
   do i=1,PairReduced%n_range3
      Omega = max(Omega,sum(PairReduced%maxrange3(1:3,i)))
   enddo
endif

end function maxOmega_PairReduced

function maxt_PairReduced(PairReduced,smallest) result(maxt)
implicit none
integer :: maxt
type(PairReducedData),intent(in) :: PairReduced
integer,intent(in),optional :: smallest

if(present(smallest)) then
   maxt = smallest - 1
else
   maxt = 0
endif

if(PairReduced%istype1) &
     maxt = max(maxt,PairReduced%maxrange1)

if(PairReduced%istype2) &
     maxt = max(maxt,PairReduced%maxrange2(2,PairReduced%n_range2))

if(PairReduced%istype3) &
     maxt = max(maxt,PairReduced%maxrange3(3,PairReduced%n_range3))

end function maxt_PairReduced

function maxuv_PairReduced_not(PairReduced,smallest) result(maxuv)
implicit none
integer :: maxuv
type(PairReducedData),intent(in) :: PairReduced
integer,intent(in),optional :: smallest
integer :: i

if(present(smallest)) then
   maxuv = smallest - 1
else
   maxuv = 0
endif

if(PairReduced%istype1) &
     maxuv = max(maxuv,PairReduced%maxrange1)

if(PairReduced%istype2) then
   do i=1,PairReduced%n_range2
      maxuv = max(maxuv,PairReduced%maxrange2(1,i))
   enddo
endif

if(PairReduced%istype3) then
   do i=1,PairReduced%n_range3
      maxuv = max(maxuv,PairReduced%maxrange3(1,i),PairReduced%maxrange3(2,i))
   enddo
endif

end function maxuv_PairReduced_not

function maxuv_PairReduced_t(t,PairReduced,smallest) result(maxuv)
implicit none
integer :: maxuv(2)
integer,intent(in) :: t
type(PairReducedData),intent(in) :: PairReduced
integer,intent(in) :: smallest
integer :: i

maxuv = smallest - 1

if(PairReduced%istype1) then
   maxuv(1) = max(maxuv(1),PairReduced%maxrange1-t)
   maxuv(2) = max(maxuv(2),PairReduced%maxrange1-t)
endif

if(PairReduced%istype2) then
   do i=1,PairReduced%n_range2
      if(PairReduced%maxrange2(2,i)==t) then
         maxuv(1) = max(maxuv(1),PairReduced%maxrange2(1,i))
         maxuv(2) = max(maxuv(2),PairReduced%maxrange2(1,i))
      end if
   enddo
endif

if(PairReduced%istype3) then
   do i=1,PairReduced%n_range3
      if(PairReduced%maxrange3(3,i)==t) then
         maxuv(1) = max(maxuv(1),PairReduced%maxrange3(1,i))
         maxuv(2) = max(maxuv(2),PairReduced%maxrange3(2,i))
      end if
   enddo
endif

if((maxuv(1)<smallest).neqv.(maxuv(2)<smallest)) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Ranges of u and v cannot be determined in maxuv_PairReduced!'
   stop
endif

end function maxuv_PairReduced_t

subroutine print_System(System,LPRINT)
use misc, only : swap
implicit none
type(SystemData) :: System
integer :: LPRINT
integer :: i,j,idx1,idx2

write(LOUT,'()')
write(LOUT,'(5x,a)') '--- Calculation specifications ---'

write(LOUT,'()')
write(LOUT,'(a,f7.3)') 'Nuclear charge:     ',System%nucZ
write(LOUT,'(a,i3)') 'Number of orbitals: ',System%norb
write(LOUT,'(a,3x,a,t35,a)',advance='no') &
     'Calculation type:   ',trim(System%calc_type),'('
if(System%post_SCF) then
   if(System%optimize) then
      write(LOUT,'(a)',advance='no') 'pair optimization'
   else
      write(LOUT,'(a)',advance='no') 'correlation energy at a given level'
   endif
else
   if(System%optimize) then
      write(LOUT,'(a)',advance='no') 'orbital optimization'
   else
      write(LOUT,'(a)',advance='no') 'only SCF energy'
   endif
endif
write(LOUT,'(a)') ')'

if(System%Neta>0) then
   write(LOUT,'()')
   write(LOUT,'(a)') 'Symmetry forcing parameter eta:'
   do i=1,System%Neta
      write(LOUT,'(i3,a,i5)') i,' :',System%eta(i)
   enddo
endif

call print_SystemExponents(System)

write(LOUT,'()')
write(LOUT,'()')
write(LOUT,'(5x,a)',advance='no') 'Orbital specification - basis set is '
if(System%common_orbs) then
   write(LOUT,'(a)') 'COMMON'
else
   write(LOUT,'(a)') 'SEPARATE'
endif
do i=1,System%n_orbs
   associate(OrbSystem => System%OrbSystem(i))
     write(LOUT,'()')
     if(.not.System%common_orbs) write(LOUT,'(a,i3)') 'orbital no. ',i
     call print_OrbSystem(OrbSystem,LPRINT)
   end associate
enddo

if(LPRINT>=5) then

   write(LOUT,'()')
   write(LOUT,'(a)') 'reduced obital information'
   write(LOUT,'()')
   write(LOUT,'(11x,a,2x,a,1x,a)') 'no.','exp','maxrange'
   do i=1,System%Nexp
      associate(OrbReduced => System%OrbReduced(i))
        write(LOUT,'(10x,i3,a)',advance='no') i,' : '
        if(OrbReduced%isUsed) then
           write(LOUT,'(i3,2x,i3)') OrbReduced%iexp,OrbReduced%maxrange
        else
           write(LOUT,'(a)') 'NOT USED'
        endif
      end associate
   enddo

endif

if(System%post_SCF) then

   write(LOUT,'()')
   write(LOUT,'()')
   write(LOUT,'(5x,a)',advance='no') 'Pair specification - basis set is '
   if(System%common_pairs) then
      write(LOUT,'(a)') 'COMMON'
   else
      write(LOUT,'(a)') 'SEPARATE'
   endif
   do j=1,System%n_pairs
      do i=1,j

         idx1 = j
         idx2 = i
         do
            associate(PairSystem => System%PairSystem(idx1,idx2))
              if(PairSystem%mult>=0) then
                 write(LOUT,'()')
                 if(.not.System%common_pairs) &
                      write(LOUT,'(a,2(i3,a),a)') &
                      'pair no. ',i,', ',j,', ',possible_mult(PairSystem%mult)
                 call print_PairSystem(PairSystem,LPRINT)
              endif
            end associate

            if(idx1<=idx2) exit
            call swap(idx1,idx2)
         enddo

      enddo
   enddo

   if(LPRINT>=5) then

      write(LOUT,'()')
      write(LOUT,'(a)') 'reduced pair information'
      write(LOUT,'()')
      write(LOUT,'(11x,a,4x,a,1x,a,t40,a,11x,a,4x,a,1x,a)') &
           'no.','exps','(expSQ)','||','no.','exps','(expSQ)'
      do i=1,System%NexpSQ
         call print_two_PairReduced(i,&
              System%OrbPairReduced(i),System%PairReduced(i))
         if(i/=System%NexpSQ) write(LOUT,'(t40,a)') '||'
      enddo

   endif

endif

end subroutine print_System

subroutine print_SystemExponents(System)
implicit none
type(SystemData) :: System
integer :: i

write(LOUT,'()')
write(LOUT,'(a)') 'Exponents'
write(LOUT,'(1x,a,7x,a,10x,a,7x,a)') &
     'no.','value','orbs','pairs'
do i=1,System%Nexp
   write(LOUT,'(i3,a,f16.8,2(l8,1x,a))') i,' :',&
        System%exponents(i),&
        System%isUsed_orbs(i) ,merge('OPT','   ',System%isOpt_orbs(i)),&
        System%isUsed_pairs(i),merge('OPT','   ',System%isOpt_pairs(i))
enddo

end subroutine print_SystemExponents

subroutine print_OrbSystem(OrbSystem,LPRINT)
implicit none
type(OrbSystemData) :: OrbSystem
integer :: LPRINT
integer :: i

write(LOUT,'(a,i5)') 'Number of basis functions: ',OrbSystem%nbas
write(LOUT,'(a,i5)') 'Number of primitives:      ',OrbSystem%n_prim

if(LPRINT>=2) then

   write(LOUT,'()')
   write(LOUT,'(11x,a,2x,a,1x,a,5x,a,2x,a)') &
        'no.','exp','range','nbas','functions'
   do i=1,OrbSystem%n_prim
      associate(OrbSpec => OrbSystem%OrbSpec(i))
        write(LOUT,'(10x,i3,a,i3,2x,i3,5x,i4,1x,i4,a,i4)') &
             i,' : ',OrbSpec%iexp,OrbSpec%irange,&
             OrbSpec%nbas,&
             OrbSpec%offset+1,' - ',OrbSpec%offset+OrbSpec%nbas
      end associate
   enddo

endif

end subroutine print_OrbSystem

subroutine print_PairSystem(PairSystem,LPRINT)
implicit none
type(PairSystemData) :: PairSystem
integer :: LPRINT
integer :: i,j

write(LOUT,'(a,i5)') 'Number of basis functions: ',PairSystem%nbas
write(LOUT,'(a,i5)') 'Number of primitives:      ',PairSystem%n_prim

if(LPRINT>=2) then

   write(LOUT,'()')
   write(LOUT,'(11x,a,4x,a,1x,a,26x,a,5x,a)') &
        'no.','exps','(expSQ)','nbas','functions'
   do j=1,PairSystem%n_prim
      associate(PairSpec => PairSystem%PairSpec(j))
        write(LOUT,'(10x,i3,a,2i3,a,i5,a,8x,a,i1,10x,i5,2x,i5,a,i5)') &
             j,' : ',PairSpec%iexp1,PairSpec%iexp2,' (',PairSpec%iexpSQ,') ',&
             'TYPE',PairSpec%itype,&
             PairSpec%nbas,&
             PairSpec%offset+1,' - ',PairSpec%offset+PairSpec%nbas
        select case(PairSpec%itype)
        case(1)
           write(LOUT,'(31x,a,i3)',advance='no') &
                't+u+v<= ',PairSpec%irange1
           if(PairSpec%restrict) then
              write(LOUT,'(a)') ', u<=v'
           else
              write(LOUT,'()')
           endif
        case(2)
           do i=1,PairSpec%n_range
              write(LOUT,'(31x,2(a,i3))',advance='no') &
                   't= ',PairSpec%irange23(2,i),&
                   ', u+v<= ',PairSpec%irange23(1,i)
              if(PairSpec%restrict) then
                 write(LOUT,'(a)') ', u<=v'
              else
                 write(LOUT,'()')
              endif
           enddo
        case(3)
           do i=1,PairSpec%n_range
              write(LOUT,'(31x,3(a,i3))',advance='no') &
                   't= ',PairSpec%irange23(3,i),&
                   ', u<= ',PairSpec%irange23(1,i),&
                   ', v<= ',PairSpec%irange23(2,i)
              if(PairSpec%restrict) then
                 write(LOUT,'(a)') ', u<=v'
              else
                 write(LOUT,'()')
              endif
           enddo
        case default
           write(LOUT,'(a)') 'ERROR!!! &
                &Incorrect range type specification in print_PairSystem!'
           stop
        end select
      end associate
      if(j/=PairSystem%n_prim) write(LOUT,'()')
   enddo

endif

end subroutine print_PairSystem

subroutine print_two_PairReduced(inum,OrbPairSpec,PairSpec)
implicit none
integer :: inum
type(PairReducedData) :: OrbPairSpec,PairSpec
character(39) :: line1,line2
integer :: i

call header(line1,OrbPairSpec)
call header(line2,PairSpec)
write(LOUT,'(3a)') line1,'||',line2

if(.not.OrbPairSpec%isUsed) return

if(PairSpec%istype1) then
   call type1(line1,OrbPairSpec)
   call type1(line2,PairSpec)
   write(LOUT,'(3a)') line1,'||',line2
endif

if(PairSpec%istype2) then
   do i=1,PairSpec%n_range2
      call type2(line1,i,OrbPairSpec)
      call type2(line2,i,PairSpec)
      write(LOUT,'(3a)') line1,'||',line2
   enddo
endif

if(PairSpec%istype3) then
   if(OrbPairSpec%n_range3==PairSpec%n_range3) then
      do i=1,PairSpec%n_range3
         call type3(line1,i,OrbPairSpec)
         call type3(line2,i,PairSpec)
         write(LOUT,'(3a)') line1,'||',line2
      enddo
   else
      call type3(line1,1,OrbPairSpec)
      call type3header(line2)
      write(LOUT,'(3a)') line1,'||',line2
      do i=1,PairSpec%n_range3
         call type3(line1,i+1,OrbPairSpec)
         call type3(line2,i  ,PairSpec,.false.)
         write(LOUT,'(3a)') line1,'||',line2
      enddo
   endif
elseif(OrbPairSpec%istype3) then
   do i=1,OrbPairSpec%n_range3
      call type3(line1,i,OrbPairSpec)
      write(LOUT,'(2a)') line1,'||'
   enddo
endif

contains

  subroutine header(line,PairSpec)
  implicit none
  character(*) :: line
  type(PairReducedData) :: PairSpec
  integer :: ipos

  write(line,'(10x,i3,a)') inum,' :'
  ipos = len_trim(line)
  if(PairSpec%isUsed) then
     write(line(ipos+1:),'(1x,2i3,a,i5,a)') &
          PairSpec%iexp1,PairSpec%iexp2,' (',PairSpec%iexpSQ,') '
  else
     write(line(ipos+1:),'(1x,a)') 'NOT USED'
  endif

  end subroutine header

  subroutine type1(line,PairSpec)
  implicit none
  character(*) :: line
  type(PairReducedData) :: PairSpec

  write(line,'(14x,a,2x,i3)') 'TYPE1',PairSpec%maxrange1

  end subroutine type1

  subroutine type2(line,irange,PairSpec)
  implicit none
  character(*) :: line
  integer :: irange
  type(PairReducedData) :: PairSpec

  write(line,'(14x,a,2x,2i3)') &
       merge('TYPE2','     ',irange==1),PairSpec%maxrange2(:,irange)

  end subroutine type2

  subroutine type3header(line)
  implicit none
  character(*) :: line

  write(line,'(14x,a,2x,3(2x,"-"))') 'TYPE3'

  end subroutine type3header

  subroutine type3(line,irange,PairSpec,do_modify)
  implicit none
  character(*) :: line
  integer :: irange
  type(PairReducedData) :: PairSpec
  logical,optional :: do_modify
  logical :: modifier

  if(present(do_modify)) then
     modifier = do_modify
  else
     modifier = .true.
  endif

  write(line,'(14x,a,2x,3i3)') &
       merge('TYPE3','     ',irange==1.and.modifier),&
       PairSpec%maxrange3(:,irange)

  end subroutine type3

end subroutine print_two_PairReduced

subroutine init_Triplet(Triplet,PairSystem)
implicit none
type(TripletData) :: Triplet
type(PairSystemData),intent(in) :: PairSystem
integer :: i_prim
integer :: ibas_orig,ibas
integer :: u,v,t,flag

Triplet%restrict  = .false.

Triplet%nbas      = PairSystem%nbas
Triplet%nbas_orig = PairSystem%nbas

do i_prim=1,PairSystem%n_prim
   associate(PairSpec => PairSystem%PairSpec(i_prim))

     Triplet%restrict = Triplet%restrict.or.PairSpec%restrict

   end associate
enddo

if(Triplet%restrict) then

   call mem_alloc(Triplet%idxS,Triplet%nbas_orig)
   call mem_alloc(Triplet%idxE,Triplet%nbas_orig)

   Triplet%idxS = 0
   Triplet%idxE = 0

   ibas_orig = 0
   ibas      = 0
   do i_prim=1,PairSystem%n_prim
      associate(PairSpec => PairSystem%PairSpec(i_prim))

        if(PairSpec%restrict) then

           flag = -1
           do while(increment_PairSpec(u,v,t,flag,PairSpec))
              ibas_orig = ibas_orig + 1
              if(u/=v) then
                 ibas      = ibas + 1
                 Triplet%idxS(ibas)      = ibas_orig
                 Triplet%idxE(ibas_orig) = ibas
              endif
           enddo

        else

           do flag=1,PairSpec%nbas
              ibas_orig = ibas_orig + 1
              ibas      = ibas + 1
              Triplet%idxS(ibas)      = ibas_orig
              Triplet%idxE(ibas_orig) = ibas
           enddo

        endif

      end associate
   enddo
   if(ibas_orig/=PairSystem%nbas) then
      write(LOUT,'(a)') 'ERROR!!! &
           &Incorrect original number of basis functions in init_Triplet!'
      stop
   endif

   Triplet%nbas = ibas

endif

end subroutine init_Triplet

subroutine free_Triplet(Triplet)
implicit none
type(TripletData) :: Triplet

if(Triplet%restrict) then

   call mem_dealloc(Triplet%idxE)
   call mem_dealloc(Triplet%idxS)

endif

end subroutine free_Triplet

subroutine shrink_Triplet_mat(mat,Triplet)
implicit none
real(prec) :: mat(:,:)
type(TripletData),intent(in) :: Triplet
integer :: ibas,jbas
integer :: j_idx

if(Triplet%restrict) then

   do jbas=1,Triplet%nbas

      j_idx = Triplet%idxS(jbas)
      do ibas=1,Triplet%nbas
         mat(ibas,jbas) = mat(Triplet%idxS(ibas),j_idx)
      enddo
      mat(Triplet%nbas+1:,jbas) = 0._prec

   enddo
   mat(:,Triplet%nbas+1:) = 0._prec

endif

end subroutine shrink_Triplet_mat

subroutine shrink_Triplet_mat_dim(mat,Triplet,dim)
implicit none
real(prec) :: mat(:,:)
type(TripletData),intent(in) :: Triplet
integer,intent(in) :: dim
integer :: ibas,jbas
integer :: j_idx

select case(dim)
case(1)

   if(Triplet%restrict) then

      do jbas=1,Triplet%nbas_orig
         do ibas=1,Triplet%nbas
            mat(ibas,jbas) = mat(Triplet%idxS(ibas),jbas)
         enddo
         mat(Triplet%nbas+1:,jbas) = 0._prec
      enddo

   endif

case(2)

   if(Triplet%restrict) then

      do jbas=1,Triplet%nbas
         j_idx = Triplet%idxS(jbas)
         do ibas=1,Triplet%nbas_orig
            mat(ibas,jbas) = mat(ibas,j_idx)
         enddo
      enddo
      mat(:,Triplet%nbas+1:) = 0._prec

   endif

case default

   write(LOUT,'(a)') 'ERROR!!! Incorrect dimension in shrink_Triplet!'
   stop

end select

end subroutine shrink_Triplet_mat_dim

subroutine shrink_Triplet_vec(vec,Triplet)
implicit none
real(prec) :: vec(:)
type(TripletData),intent(in) :: Triplet
integer :: ibas

if(Triplet%restrict) then

   do ibas=1,Triplet%nbas
      vec(ibas) = vec(Triplet%idxS(ibas))
   enddo
   vec(Triplet%nbas+1:) = 0._prec

endif

end subroutine shrink_Triplet_vec

subroutine expand_Triplet_mat(mat,Triplet)
implicit none
real(prec) :: mat(:,:)
type(TripletData),intent(in) :: Triplet
integer :: ibas,jbas
integer :: j_idx

if(Triplet%restrict) then

   do jbas=Triplet%nbas_orig,1,-1
      if(Triplet%idxE(jbas)>0) then

         j_idx = Triplet%idxE(jbas)
         do ibas=Triplet%nbas_orig,1,-1
            if(Triplet%idxE(ibas)>0) then
               mat(ibas,jbas) = mat(Triplet%idxE(ibas),j_idx)
            else
               mat(ibas,jbas) = 0._prec
            endif
         enddo

      else
         mat(:,jbas) = 0._prec
      endif
   enddo
      
endif

end subroutine expand_Triplet_mat

subroutine expand_Triplet_vec(vec,Triplet)
implicit none
real(prec) :: vec(:)
type(TripletData),intent(in) :: Triplet
integer :: ibas

if(Triplet%restrict) then

   do ibas=Triplet%nbas_orig,1,-1
      if(Triplet%idxE(ibas)>0) then
         vec(ibas) = vec(Triplet%idxE(ibas))
      else
         vec(ibas) = 0._prec
      endif
   enddo

endif

end subroutine expand_Triplet_vec

subroutine init_SCForbitals(SCForbitals,System)
implicit none
type(SCForbitalsData) :: SCForbitals
type(SystemData) :: System
integer :: i

SCForbitals%norb = System%norb
SCForbitals%nbas = 0
do i=1,System%n_orbs
   associate(OrbSystem => System%OrbSystem(i))
     SCForbitals%nbas = max(SCForbitals%nbas,OrbSystem%nbas)
   end associate
enddo

call mem_alloc(SCForbitals%orb_energy,SCForbitals%norb)
call mem_alloc(SCForbitals%orb_vector,SCForbitals%nbas,SCForbitals%norb)
call mem_alloc(SCForbitals%pair_energy,&
     SCForbitals%norb,SCForbitals%norb,&
     SCForbitals%norb*(SCForbitals%norb + 1)/2)

SCForbitals%energy = 0._prec

SCForbitals%orb_energy  = 0._prec
SCForbitals%orb_vector  = 0._prec
SCForbitals%pair_energy = 0._prec

end subroutine init_SCForbitals

subroutine free_SCForbitals(SCForbitals)
implicit none
type(SCForbitalsData) :: SCForbitals

call mem_dealloc(SCForbitals%pair_energy)
call mem_dealloc(SCForbitals%orb_vector)
call mem_dealloc(SCForbitals%orb_energy)

end subroutine free_SCForbitals

subroutine print_SCForbitals(SCForbitals,System)
implicit none
type(SCForbitalsData) :: SCForbitals
type(SystemData) :: System
integer :: ibas,iorb
integer :: i1,j1,i2,j2,ij2
integer :: i,j
real(prec) :: energy
integer,allocatable :: nbas(:)

write(LOUT,'()')
write(LOUT,'(5x,a)') '--- SCF orbitals and energies ---'

write(LOUT,'()')
write(LOUT,'(a,i5)') 'Number of orbitals: ',SCForbitals%norb
if(System%common_orbs) then
   write(LOUT,'(a,i5)') 'Basis set size:     ',SCForbitals%nbas
else
   write(LOUT,'(a,i5)') 'Max basis set size: ',SCForbitals%nbas
endif

write(LOUT,'()')
write(LOUT,'(13x,a,i2,9(:15x,a,i2))') ('orb no. ',iorb,iorb=1,SCForbitals%norb)

write(LOUT,'()')
write(LOUT,'(a)') 'AO coefficients'
if(System%common_orbs) then
   do ibas=1,SCForbitals%nbas
      write(LOUT,'(i5,10es25.14)') ibas,&
           (SCForbitals%orb_vector(ibas,iorb),iorb=1,SCForbitals%norb)
   enddo
else
   call mem_alloc(nbas,SCForbitals%norb)
   do iorb=1,SCForbitals%norb
      nbas(iorb) = System%OrbSystem(iorb)%nbas
   enddo
   do ibas=1,SCForbitals%nbas
      write(LOUT,'(i5)',advance='no') ibas
      do iorb=1,SCForbitals%norb
         if(ibas<=nbas(iorb)) then
            write(LOUT,'(es25.14)',advance='no') &
                 SCForbitals%orb_vector(ibas,iorb)
         else
            write(LOUT,'(a25)',advance='no') ''
         endif
      enddo
      write(LOUT,'()')
   enddo
   call mem_dealloc(nbas)
endif

write(LOUT,'()')
write(LOUT,'(a)') 'orbital energies'
write(LOUT,'(5x,10f25.18)') &
     (SCForbitals%orb_energy(iorb),iorb=1,SCForbitals%norb)

write(LOUT,'()') 
write(LOUT,'(a,f25.18)') 'LUMO: ',SCForbitals%LUMO
write(LOUT,'(a,f25.18)') 'HOMO: ',SCForbitals%HOMO

write(LOUT,'()')
write(LOUT,'(a)') 'pair energies'
ij2 = 0
do j2=1,SCForbitals%norb
   do i2=1,j2
      ij2 = ij2 + 1
      do j1=1,SCForbitals%norb
         do i1=1,SCForbitals%norb
            write(LOUT,'(5x,a,2i2,a,2i2,a,f25.18)') &
                 '<',i1,j1,' |',i2,j2,' >',SCForbitals%pair_energy(i1,j1,ij2)
         enddo
      enddo
      write(LOUT,'()')
   enddo
enddo

energy = 0._prec
do i=1,SCForbitals%norb
   energy = energy &
        + 2._prec*(SCForbitals%orb_energy(i) &
        - 0.5_prec*SCForbitals%pair_energy(i,i,i*(i+1)/2))
enddo
do j=1,SCForbitals%norb
   do i=1,j-1
      energy = energy - 2._prec*( &
           2*SCForbitals%pair_energy(i,j,i+j*(j-1)/2) &
           - SCForbitals%pair_energy(j,i,i+j*(j-1)/2))
   enddo
enddo

write(LOUT,'(a)') 'SCF energy'
write(LOUT,'(a,f25.18)') 'from orbital and pair energies: ',energy
write(LOUT,'(a,f25.18)') 'from the SCF procedure:         ',SCForbitals%energy

end subroutine print_SCForbitals

end module commontypes
