module calcdriver
use file_OUT, only : LOUT
use precision, only : prec
use memory
use commontypes
use optimizer
use SCFdriver
use CCdriver
implicit none

private
public calculate

contains

 subroutine calculate(System,Control)
 use misc, only : pack_vector,unpack_vector
 type(SystemData) :: System
 type(ControlData) :: Control
 type(SCForbitalsData) :: SCForbitals
 integer :: Noptimize
 real(prec),allocatable :: optimize(:)
 
 if(System%post_SCF) then
 
    call init_SCForbitals(SCForbitals,System)
 
    write(LOUT,'()')
    write(LOUT,'(5x,a)') '--- SCF calculation ---'
    call calculateSCF(System,Control,SCForbitals)
    if(Control%LPRINT>=5) call print_SCForbitals(SCForbitals,System)
 
    call print_EndSection
 
    call CC_fix_SCForbitals(SCForbitals)
    if(System%optimize) then
 
       write(LOUT,'()')
       write(LOUT,'(a)') 'NON AVAILABLE!'
!       write(LOUT,'(5x,3a)') &
!            '--- Initial ',&
!            trim(System%calc_type(1:index(System%calc_type,'OPT')-1)),&
!            ' calculation ---'
!       call calculateCC(System,Control)
! 
!       Noptimize = count(System%isOpt_pairs)
!       call mem_alloc(optimize,Noptimize)
!       call pack_vector(optimize,&
!            System%Nexp,System%exponents,System%isOpt_pairs)
! 
!       call powell(Noptimize,optimize,energyCC,System,Control)
! 
!       call unpack_vector(optimize,&
!            System%Nexp,System%exponents,System%isOpt_pairs)
!       call mem_dealloc(optimize)
! 
!       write(LOUT,'()')
!       write(LOUT,'(5x,3a)') &
!            '--- Final ',&
!            trim(System%calc_type(1:index(System%calc_type,'OPT')-1)),&
!            ' calculation ---'
!       call calculateCC(System,Control)
! 
!       call print_SystemExponents(System)
! 
    else
 
       write(LOUT,'()')
       write(LOUT,'(5x,3a)') '--- ',trim(System%calc_type),' calculation ---'
       call calculateCC(System,Control)
 
    endif
    call CC_release_SCForbitals
 
    call free_SCForbitals(SCForbitals,System%maxl)
 
 else
 
    if(System%optimize) then
 
       write(LOUT,'()')
       write(LOUT,'(5x,a)') '--- Initial SCF calculation ---'
       call calculateSCF(System,Control)
 
       Noptimize = count(System%isOpt_orbs)
       call mem_alloc(optimize,Noptimize)
       call pack_vector(optimize,&
            System%Nexp,System%exponents,System%isOpt_orbs)
 
       call powell(Noptimize,optimize,energySCF,System,Control)
 
       call unpack_vector(optimize,&
            System%Nexp,System%exponents,System%isOpt_orbs)
       call mem_dealloc(optimize)
 
       write(LOUT,'()')
       write(LOUT,'(5x,a)') '--- Final SCF calculation ---'
       call calculateSCF(System,Control)
 
       call print_SystemExponents(System)
 
    else
 
       write(LOUT,'()')
       write(LOUT,'(5x,a)') '--- SCF calculation ---'
       call calculateSCF(System,Control)
! 
    endif
 
 endif
 
 end subroutine calculate
 
end module calcdriver
