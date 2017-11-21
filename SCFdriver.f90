module SCFdriver
use file_OUT, only : LOUT
use precision, only : prec
use memory
use eproblem
use eproblem_qd
use commontypes
use SCFint
implicit none

private
public calculateSCF,energySCF

integer,parameter :: number_OVERLAP = 3

real(prec) :: nucZ
integer :: Nexp
real(prec),allocatable :: exponents(:)

real(prec) :: RESULT_energy
type(SCForbitalsData),pointer :: RESULT_extend

type Problem_l_Data
real(prec), allocatable :: H(:,:)
real(prec), allocatable :: S(:,:)
real(prec), allocatable :: D(:,:)
end type Problem_l_Data 

type Problem_ll_Data
real(prec), allocatable :: J(:,:)
real(prec), allocatable :: K(:,:)
end type Problem_ll_Data

type RefLamData
integer :: nlambda
integer,allocatable :: lambda(:)
real(prec),allocatable :: prefac(:)
end type RefLamData

type(Problem_l_Data),allocatable :: problem_l(:)
type(Problem_ll_Data),allocatable :: problem_ll(:,:)

type(RefLamData),allocatable :: refJ(:,:), refK(:,:)

contains

subroutine calculateSCF(System,Control,SCForbitals)
implicit none
type(SystemData) :: System
type(ControlData) :: Control
type(SCForbitalsData),target,optional :: SCForbitals

nucZ = System%nucZ

Nexp = System%Nexp
call mem_alloc(exponents,Nexp)
exponents(:) = System%exponents

RESULT_energy = 0._prec
nullify(RESULT_extend)
if(present(SCForbitals)) RESULT_extend => SCForbitals

call chooseSCF(System,Control,.true.)

if(associated(RESULT_extend)) nullify(RESULT_extend)

call mem_dealloc(exponents)

end subroutine calculateSCF

function energySCF(MODexponents,System,Control) result(energy)
use misc, only : unpack_vector
implicit none
real(prec) :: energy
real(prec) :: MODexponents(:)
type(SystemData) :: System
type(ControlData) :: Control

nucZ = System%nucZ

Nexp = System%Nexp
call mem_alloc(exponents,Nexp)
exponents(:) = System%exponents
call unpack_vector(MODexponents,Nexp,exponents,System%isOpt_orbs)

RESULT_energy = 0._prec
nullify(RESULT_extend)

call chooseSCF(System,Control,.false.)

energy = RESULT_energy

call mem_dealloc(exponents)

end function energySCF

subroutine chooseSCF(System,Control,fullPRINT)
implicit none
type(SystemData) :: System
type(ControlData) :: Control
logical :: fullPRINT
integer :: i

 if(System%common_orbs) then
    do i=1,System%maxl+1
       if(System%n_orbs(i)/=1) then
          write(LOUT,'(a)') 'ERROR!!! Inconsistent data in chooseSCF!'
          stop
       endif
    enddo

    call SCF_energy_common(System%maxl,System%norb,System%OrbSystem_L,System%OrbReduced,&
         Control,fullPRINT)

  else
    do i=1,System%maxl+1
       if(System%n_orbs(i)/=System%norb(i)) then
          write(LOUT,'(a)') 'ERROR!!! Inconsistent data in chooseSCF!'
          stop
       endif
    enddo
 
    do i=1,System%maxl+1
       if(System%n_orbs(i)/=2) then
          write(LOUT,'(a)') 'ERROR!!! Separate SCF orbitals only for norb = 2 !'
          stop
       endif
    enddo

!    call SCF_energy_separate(System%OrbSystem,System%OrbReduced,&
!         Control,fullPRINT)

  endif

! hapka: old_drake
! if(System%common_orbs) then
!    if(System%n_orbs/=1) then
!       write(LOUT,'(a)') 'ERROR!!! Inconsistent data in chooseSCF!'
!       stop
!    endif
! 
!    call SCF_energy_common(System%norb,System%OrbSystem(1),System%OrbReduced,&
!         Control,fullPRINT)
! 
! else
!    if(System%n_orbs/=System%norb) then
!       write(LOUT,'(a)') 'ERROR!!! Inconsistent data in chooseSCF!'
!       stop
!    
!    if(System%norb/=2) then
!       write(LOUT,'(a)') 'ERROR!!! Separate SCF orbitals only for norb = 2 !'
!       stop
!    endif
! 
!    call SCF_energy_separate(System%OrbSystem,System%OrbReduced,&
!         Control,fullPRINT)
! 
! endif

end subroutine chooseSCF

subroutine SCF_energy_common(maxl,norb,OrbSystem_L,OrbReduced,Control,fullPRINT)
implicit none
integer :: maxl
!integer,intent(in) :: norb
integer,intent(in) :: norb(:)
!type(OrbSystemData),intent(in) :: OrbSystem
type(OrbSystemData_L),intent(in) :: OrbSystem_L(:)
type(OrbReducedData),intent(in) :: OrbReduced(:)
type(ControlData),intent(in) :: Control
logical,intent(in) :: fullPRINT
integer :: LPRINT_mod
!integer :: nbas
integer,allocatable :: nbas(:)
logical :: converged
real(prec) :: energy
real(prec),allocatable :: matH(:,:),matS(:,:),twoel(:,:)
real(prec),allocatable :: eval(:),evec(:,:)
real(prec),allocatable :: matF(:,:),matG(:,:),dens(:,:)
real(prec) :: prefac
integer :: lambda
integer :: i, j, k
integer :: i1,j1,i2,j2,ij2

LPRINT_mod = merge(Control%LPRINT,0,fullPRINT)

!nbas = OrbSystem%nbas

allocate(nbas(maxl+1))

do i=1,maxl+1
 nbas(i) = OrbSystem_L(i)%OrbSystem(1)%nbas
enddo

allocate(problem_l(maxl+1))
allocate(problem_ll(maxl+1,maxl+1))
do i=1,maxl+1
   call mem_alloc(problem_l(i)%H,nbas(i),nbas(i))
   call mem_alloc(problem_l(i)%S,nbas(i),nbas(i))
! add twoel!
enddo

do i=1,maxl+1
   do j=1,maxl+1
      call mem_alloc(problem_ll(i,j)%J,nbas(i)**2,nbas(j)**2)
      call mem_alloc(problem_ll(i,j)%K,nbas(i)**2,nbas(j)**2)
   enddo
enddo

! hapka: old_drake
!call mem_alloc(matH,nbas,nbas)
!call mem_alloc(matS,nbas,nbas)
!call mem_alloc(twoel,nbas**2,nbas**2)

call create_SCFint(Nexp,exponents,OrbReduced,LPRINT_mod)
call create_refJK(maxl)

do i=1,maxl+1
  associate(&
       matH => problem_l(i)%H, &
       matS => problem_l(i)%S, &
       OrbSystem => OrbSystem_L(i)%OrbSystem(1) )

       call SCFint_matS(matS,OrbSystem,OrbSystem)
       call SCFint_matH(matH,OrbSystem,OrbSystem,nucZ)

!       do j=1,nbas(i)
!          write(*,*) matH(j,j)
!       enddo
!       write(*,*) ''

  end associate
enddo
  
! hapka: old_drake
!call SCFint_matH(matH,OrbSystem,OrbSystem,nucZ)
!call SCFint_matS(matS,OrbSystem,OrbSystem)

do i=0,maxl
   do j=0,maxl
   associate(&
      refJ => refJ(i+1,j+1), &
      refK => refK(i+1,j+1), &
      matJ => problem_ll(i+1,j+1)%J, &
      matK => problem_ll(i+1,j+1)%K, &
      ijOrbSystem => OrbSystem_L(i+1)%OrbSystem(1), &
      klOrbSystem => OrbSystem_L(j+1)%OrbSystem(1)  )
 
      matJ = 0._prec
      do k=1,refJ%nlambda 
         prefac = refJ%prefac(k)
         lambda = refJ%lambda(k)
         call SCFint_matJ(matJ,ijOrbSystem,ijOrbSystem, &
                       &  klOrbSystem,klOrbSystem,lambda,prefac)
      enddo

      matK = 0._prec
      do k=1,refK%nlambda
         prefac = refK%prefac(k)
         lambda = refK%lambda(k)
         call SCFint_matK(matK,ijOrbSystem,ijOrbSystem, &
                    &  klOrbSystem,klOrbSystem,lambda,prefac)
      enddo

   end associate
   enddo
enddo

!call SCFint_matJ(twoel,OrbSystem,OrbSystem,OrbSystem,OrbSystem)

call free_refJK(maxl)
call free_SCFint

!call mem_alloc(eval,nbas)
!call mem_alloc(evec,nbas,nbas)
!call mem_alloc(matF,nbas,nbas)
!call mem_alloc(matG,nbas,nbas)
!call mem_alloc(dens,nbas,nbas)

!call Roothaan_iterations(norb,nbas,converged,energy,eval,evec,&
!     matH,matS,twoel,matF,matG,dens,Control,LPRINT_mod)

!if(fullPRINT) then
!   write(LOUT,'()')
!   write(LOUT,'(a)') 'Orbital energies'
!   do i=1,norb
!      write(LOUT,'(i5,f25.18)') i,eval(i)
!   enddo
!   write(LOUT,'()')
!   if(converged) then
!      write(LOUT,'(a,f25.18)') '    Converged SCF energy: ',energy
!   else
!      write(LOUT,'(a,f25.18)') 'NOT converged SCF energy: ',energy
!   endif
!endif
!
!RESULT_energy = energy
!
!if(associated(RESULT_extend)) then
!   RESULT_extend%energy = energy
!   RESULT_extend%HOMO = eval(norb)
!   RESULT_extend%LUMO = merge(eval(norb+1),0._prec,nbas>norb)
!   RESULT_extend%orb_energy(1:norb) = eval(1:norb)
!   RESULT_extend%orb_vector(1:nbas,1:norb) = evec(1:nbas,1:norb)
!   ij2 = 0
!   do j2=1,norb
!      do i2=1,j2
!         ij2 = ij2 + 1
!         do j1=1,norb
!            do i1=1,norb
!               RESULT_extend%pair_energy(i1,j1,ij2) = make_contract(&
!                    nbas,RESULT_extend%orb_vector(:,i1),&
!                    nbas,RESULT_extend%orb_vector(:,i2),&
!                    nbas,RESULT_extend%orb_vector(:,j1),&
!                    nbas,RESULT_extend%orb_vector(:,j2),&
!                    twoel,dens)
!            enddo
!         enddo
!      enddo
!   enddo
!endif
!
!call mem_dealloc(dens)
!call mem_dealloc(matG)
!call mem_dealloc(matF)
!call mem_dealloc(evec)
!call mem_dealloc(eval)

!call mem_dealloc(twoel)
!call mem_dealloc(matS)
!call mem_dealloc(matH)

do i=1,maxl+1
   call mem_dealloc(problem_l(i)%H)
   call mem_dealloc(problem_l(i)%S)
enddo

do i=1,maxl+1
   do j=1,maxl+1
      call mem_dealloc(problem_ll(i,j)%J)
      call mem_dealloc(problem_ll(i,j)%K)
   enddo
enddo


deallocate(problem_l)
deallocate(problem_ll)
deallocate(nbas)

end subroutine SCF_energy_common

subroutine SCF_energy_separate(OrbSystem,OrbReduced,Control,fullPRINT)
implicit none
type(OrbSystemData),intent(in) :: OrbSystem(:)
type(OrbReducedData),intent(in) :: OrbReduced(:)
type(ControlData),intent(in) :: Control
logical,intent(in) :: fullPRINT
integer :: LPRINT_mod
integer :: nbas1,nbas2,nbas_max
integer :: iter,cnt_min
logical :: converged
real(prec) :: energy1,energy2,energy,energy_prev,energy_delta,energy_delta_min
real(prec) :: s11,s12,s21,s22,norms
real(prec),allocatable :: matH1(:,:),matS1(:,:),matJ1(:,:)
real(prec),allocatable :: matH2(:,:),matS2(:,:),matJ2(:,:)
real(prec),allocatable :: matS12(:,:),matJ12(:,:),matK12(:,:)
real(prec),allocatable :: vector1(:),eval1(:),evec1(:,:)
real(prec),allocatable :: vector2(:),eval2(:),evec2(:,:)
real(prec),allocatable :: matF(:,:),matG(:,:),dens(:,:)
real(prec),allocatable :: vecP(:),matFmod(:,:),matSmod(:,:)
real(prec),allocatable :: twoel(:,:)
integer :: i

LPRINT_mod = merge(Control%LPRINT,0,fullPRINT)

nbas1 = OrbSystem(1)%nbas
nbas2 = OrbSystem(2)%nbas
nbas_max = max(nbas1,nbas2)

call mem_alloc(matH1,nbas1,nbas1)
call mem_alloc(matS1,nbas1,nbas1)
call mem_alloc(matJ1,nbas1**2,nbas1**2)

call mem_alloc(matH2,nbas2,nbas2)
call mem_alloc(matS2,nbas2,nbas2)
call mem_alloc(matJ2,nbas2**2,nbas2**2)

call mem_alloc(matS12,nbas1,nbas2)
call mem_alloc(matJ12,nbas1**2,nbas2**2)
call mem_alloc(matK12,nbas1**2,nbas2**2)

! call create_SCFint(Nexp,exponents,OrbReduced,LPRINT_mod)
! 
! call SCFint_matH(matH1,OrbSystem(1),OrbSystem(1),nucZ)
! call SCFint_matS(matS1,OrbSystem(1),OrbSystem(1))
! call SCFint_matJ(matJ1,OrbSystem(1),OrbSystem(1),OrbSystem(1),OrbSystem(1))
! 
! call SCFint_matH(matH2,OrbSystem(2),OrbSystem(2),nucZ)
! call SCFint_matS(matS2,OrbSystem(2),OrbSystem(2))
! call SCFint_matJ(matJ2,OrbSystem(2),OrbSystem(2),OrbSystem(2),OrbSystem(2))
! 
! call SCFint_matS(matS12,OrbSystem(1),OrbSystem(2))
! call SCFint_matJ(matJ12,OrbSystem(1),OrbSystem(1),OrbSystem(2),OrbSystem(2))
! call SCFint_matK(matK12,OrbSystem(1),OrbSystem(1),OrbSystem(2),OrbSystem(2))
! 
! call free_SCFint

call mem_alloc(vector1,nbas1)
call mem_alloc(eval1,nbas1)
call mem_alloc(evec1,nbas1,nbas1)
call mem_alloc(vector2,nbas2)
call mem_alloc(eval2,nbas2)
call mem_alloc(evec2,nbas2,nbas2)
call mem_alloc(matF,nbas_max,nbas_max)
call mem_alloc(matG,nbas_max,nbas_max)
call mem_alloc(dens,nbas_max,nbas_max)
call mem_alloc(vecP,nbas_max)
call mem_alloc(matFmod,nbas_max,nbas_max)
call mem_alloc(matSmod,nbas_max,nbas_max)

if(LPRINT_mod>=1) then
   write(LOUT,'()')
   write(LOUT,'()')
   write(LOUT,'(10x,a)') 'SCF for primitives of the FIRST  orbital'
   write(LOUT,'()')
endif

call Roothaan_iterations(2,nbas1,converged,energy1,eval1,evec1,&
     matH1,matS1,matJ1,matF,matG,dens,Control,LPRINT_mod)

if(LPRINT_mod>=1) then
   write(LOUT,'()')
   write(LOUT,'(a)') 'Orbital energies'
   do i=1,2
      write(LOUT,'(i5,f25.18)') i,eval1(i)
   enddo
   write(LOUT,'()')
   if(converged) then
      write(LOUT,'(a,f25.18)') '    Converged SCF energy: ',energy1
   else
      write(LOUT,'(a,f25.18)') 'NOT converged SCF energy: ',energy1
   endif
endif

if(LPRINT_mod>=1) then
   write(LOUT,'()')
   write(LOUT,'()')
   write(LOUT,'(10x,a)') 'SCF for primitives of the SECOND orbital'
   write(LOUT,'()')
endif

call Roothaan_iterations(2,nbas2,converged,energy2,eval2,evec2,&
     matH2,matS2,matJ2,matF,matG,dens,Control,LPRINT_mod)

if(LPRINT_mod>=1) then
   write(LOUT,'()')
   write(LOUT,'(a)') 'Orbital energies'
   do i=1,2
      write(LOUT,'(i5,f25.18)') i,eval2(i)
   enddo
   write(LOUT,'()')
   if(converged) then
      write(LOUT,'(a,f25.18)') '    Converged SCF energy: ',energy2
   else
      write(LOUT,'(a,f25.18)') 'NOT converged SCF energy: ',energy2
   endif
endif

call make_dens(nbas1,nbas2,evec1(:,1),evec2(:,1),dens)
s11 = make_trace(nbas1,nbas2,matS12,dens)
call make_dens(nbas1,nbas2,evec1(:,1),evec2(:,2),dens)
s12 = make_trace(nbas1,nbas2,matS12,dens)

norms = hypot(s11,s12)
vector1(:) = evec1(:,1)
vector2(:) = (s11/norms)*evec2(:,2) + (-s12/norms)*evec2(:,1)

if(LPRINT_mod>=1) then
   call make_dens(nbas1,nbas2,evec1(:,2),evec2(:,1),dens)
   s21 = make_trace(nbas1,nbas2,matS12,dens)
   call make_dens(nbas1,nbas2,evec1(:,2),evec2(:,2),dens)
   s22 = make_trace(nbas1,nbas2,matS12,dens)

   write(LOUT,'()')
   write(LOUT,'()')
   write(LOUT,'(10x,a)') 'Overlap between raw FIRST and SECOND orbitals'
   write(LOUT,'()')
   write(LOUT,'(10x,5x,2(6x,i5,7x))') 1,2
   write(LOUT,'(10x,i5,2es18.6)') 1,s11,s12
   write(LOUT,'(10x,i5,2es18.6)') 2,s21,s22
endif

if(LPRINT_mod>=1) then
   write(LOUT,'()')
   write(LOUT,'()')
   write(LOUT,'(10x,a)') 'SCF: different basis sets for different orbitals'
   write(LOUT,'()')
   write(LOUT,'(1x,a,12x,a,12x,a,2x,a)') &
        'iter','energy','energy change','min-change counter'
endif

iter = 0

cnt_min = 0
energy_delta_min = huge(0._prec)
converged = .false.
do
   energy = 0._prec
   call make_dens(nbas2,nbas2,vector2,vector2,dens)
   call make_matG_2( 1._prec,nbas2,nbas2,matG,1,matJ2,dens)
   energy = energy &
        + 2._prec*(make_trace(nbas2,nbas2,matH2,dens) &
        + 0.5_prec*make_trace(nbas2,nbas2,matG,dens))
   call make_matG_2( 2._prec,nbas1,nbas2,matG,1,matJ12,dens,.true.)
   call make_matG_2(-1._prec,nbas1,nbas2,matG,1,matK12,dens,.false.)
   call make_dens(nbas1,nbas1,vector1,vector1,dens)
   energy = energy &
        + 2._prec*make_trace(nbas1,nbas1,matG,dens)
   call make_matG_2( 1._prec,nbas1,nbas1,matG,1,matJ1,dens)
   energy = energy &
        + 2._prec*(make_trace(nbas1,nbas1,matH1,dens) &
        + 0.5_prec*make_trace(nbas1,nbas1,matG,dens))

   if(iter>0) then

      energy_delta = abs(energy-energy_prev)
      if(energy_delta<energy_delta_min) then
         cnt_min = 0
         energy_delta_min = energy_delta
      else
         cnt_min = cnt_min + 1
      endif

      if(LPRINT_mod>=1) write(LOUT,'(i5,f25.18,es18.6,6x,i5)') &
           iter,energy,energy_delta,cnt_min

      if(energy_delta<Control%SCFTHR) converged = .true.

   else

      if(LPRINT_mod>=1) write(LOUT,'(i5,f25.18)') iter,energy

   endif
   energy_prev = energy

   iter = iter + 1

   if(iter>Control%SCFMAXIT.and..not.converged) then
      write(LOUT,'(a,i5,a,es10.3,a)') 'WARNING!!! &
           &SCF did not converge in SCFMAXIT =', &
           Control%SCFMAXIT,' iterations!     (',energy_delta,' )'
      exit
   endif
   if(cnt_min>=Control%SCFCNT.and..not.converged) then
      write(LOUT,'(a,i5,a,es10.3,a)') 'WARNING!!! &
           &No significant change in last SCFCNT =', &
           Control%SCFCNT,' iterations! (',energy_delta,' )'
      exit
   endif

   call make_dens(nbas2,nbas2,vector2,vector2,dens)
   call make_matG_2( 1._prec,nbas2,nbas2,matG,2,matJ2,dens)
   call make_dens(nbas1,nbas1,vector1,vector1,dens)
   call make_matG_2( 2._prec,nbas1,nbas2,matG,2,matJ12,dens,.false.)
   call make_matG_2(-1._prec,nbas1,nbas2,matG,2,matK12,dens,.false.)

   matF(1:nbas2,1:nbas2) = matH2(1:nbas2,1:nbas2) + matG(1:nbas2,1:nbas2)
   call make_over(nbas1,nbas2,vecP,2,matS12,vector1)
   call transform(nbas2,vecP,matF,matFmod)
   call transform(nbas2,vecP,matS2,matSmod)
   call symU_diagonalize_qd('QL',nbas2-1,eval2,evec2,matFmod,matSmod)
   energy2 = eval2(1)
   if(.not.converged) call backtransform(nbas2,vecP,evec2(:,1),vector2)

   call make_dens(nbas1,nbas1,vector1,vector1,dens)
   call make_matG_2( 1._prec,nbas1,nbas1,matG,1,matJ1,dens)
   call make_dens(nbas2,nbas2,vector2,vector2,dens)
   call make_matG_2( 2._prec,nbas1,nbas2,matG,1,matJ12,dens,.false.)
   call make_matG_2(-1._prec,nbas1,nbas2,matG,1,matK12,dens,.false.)

   matF(1:nbas1,1:nbas1) = matH1(1:nbas1,1:nbas1) + matG(1:nbas1,1:nbas1)
   call make_over(nbas1,nbas2,vecP,1,matS12,vector2)
   call transform(nbas1,vecP,matF,matFmod)
   call transform(nbas1,vecP,matS1,matSmod)
   call symU_diagonalize_qd('QL',nbas1-1,eval1,evec1,matFmod,matSmod)
   energy1 = eval1(1)
   if(.not.converged) call backtransform(nbas1,vecP,evec1(:,1),vector1)

   if(converged) exit
enddo

if(fullPRINT) then
   call make_dens(nbas1,nbas2,vector1,vector2,dens)
   s12 = make_trace(nbas1,nbas2,matS12,dens)

   write(LOUT,'()')
   write(LOUT,'(a)') 'Orbital energies'
   write(LOUT,'(i5,f25.18)') 1,energy1
   write(LOUT,'(i5,f25.18)') 2,energy2
   write(LOUT,'()')
   write(LOUT,'(a,es18.6)') 'Overlap check: ',s12
   write(LOUT,'()')
   if(converged) then
      write(LOUT,'(a,f25.18)') '    Converged SCF energy: ',energy
   else
      write(LOUT,'(a,f25.18)') 'NOT converged SCF energy: ',energy
   endif
endif

RESULT_energy = energy

if(associated(RESULT_extend)) then
   RESULT_extend%energy = energy
   RESULT_extend%HOMO = max(energy1,energy2)
   RESULT_extend%LUMO = min(&
        merge(eval1(2),0._prec,nbas1>2),&
        merge(eval2(2),0._prec,nbas2>2))
   RESULT_extend%orb_energy(1) = energy1
   RESULT_extend%orb_energy(2) = energy2
   RESULT_extend%orb_vector(1:nbas1,1) = vector1
   RESULT_extend%orb_vector(1:nbas2,2) = vector2
endif

call mem_dealloc(matSmod)
call mem_dealloc(matFmod)
call mem_dealloc(vecP)
call mem_dealloc(dens)
call mem_dealloc(matG)
call mem_dealloc(matF)
call mem_dealloc(evec2)
call mem_dealloc(eval2)
call mem_dealloc(vector2)
call mem_dealloc(evec1)
call mem_dealloc(eval1)
call mem_dealloc(vector1)

call mem_dealloc(matK12)
call mem_dealloc(matJ12)
call mem_dealloc(matS12)

call mem_dealloc(matJ2)
call mem_dealloc(matS2)
call mem_dealloc(matH2)

call mem_dealloc(matJ1)
call mem_dealloc(matS1)
call mem_dealloc(matH1)

if(associated(RESULT_extend)) then
   call mem_alloc(dens,nbas_max,nbas_max)
   call mem_alloc(twoel,nbas_max**2,nbas_max**2)
   call create_SCFint(Nexp,exponents,OrbReduced,0)

!   call SCFint_matJ(twoel,OrbSystem(1),OrbSystem(1),OrbSystem(1),OrbSystem(1))
   energy = make_contract(&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas1,RESULT_extend%orb_vector(:,1),&
        twoel,dens)
   RESULT_extend%pair_energy(1,1,1) = energy

!   call SCFint_matJ(twoel,OrbSystem(1),OrbSystem(1),OrbSystem(1),OrbSystem(2))
   energy = make_contract(&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        twoel,dens)
   RESULT_extend%pair_energy(2,1,1) = energy
   RESULT_extend%pair_energy(1,2,1) = energy
   RESULT_extend%pair_energy(1,1,2) = energy

!   call SCFint_matJ(twoel,OrbSystem(1),OrbSystem(2),OrbSystem(1),OrbSystem(2))
   energy = make_contract(&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        twoel,dens)
   RESULT_extend%pair_energy(2,2,1) = energy
   RESULT_extend%pair_energy(2,1,2) = energy
   RESULT_extend%pair_energy(1,1,3) = energy

!   call SCFint_matJ(twoel,OrbSystem(1),OrbSystem(1),OrbSystem(2),OrbSystem(2))
   energy = make_contract(&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        twoel,dens)
   RESULT_extend%pair_energy(1,2,2) = energy

!   call SCFint_matJ(twoel,OrbSystem(1),OrbSystem(2),OrbSystem(2),OrbSystem(2))
   energy = make_contract(&
        nbas1,RESULT_extend%orb_vector(:,1),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        twoel,dens)
   RESULT_extend%pair_energy(2,2,2) = energy
   RESULT_extend%pair_energy(2,1,3) = energy
   RESULT_extend%pair_energy(1,2,3) = energy

!   call SCFint_matJ(twoel,OrbSystem(2),OrbSystem(2),OrbSystem(2),OrbSystem(2))
   energy = make_contract(&
        nbas2,RESULT_extend%orb_vector(:,2),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        nbas2,RESULT_extend%orb_vector(:,2),&
        twoel,dens)
   RESULT_extend%pair_energy(2,2,3) = energy

   call free_SCFint
   call mem_dealloc(twoel)
   call mem_dealloc(dens)
endif

end subroutine SCF_energy_separate

subroutine Roothaan_iterations(norb,nbas,converged,energy,eval,evec,&
     matH,matS,twoel,matF,matG,dens,Control,LPRINT)
implicit none
integer,intent(in) :: norb,nbas
logical,intent(out) :: converged
real(prec),intent(out) :: energy,eval(:),evec(:,:)
real(prec),intent(in) :: matH(:,:),matS(:,:),twoel(:,:)
real(prec) :: matF(:,:),matG(:,:),dens(:,:)
type(ControlData),intent(in) :: Control
integer,intent(in) :: LPRINT
logical :: level_shifting
integer :: iter,cnt_min
real(prec) :: energy_prev,energy_delta,energy_delta_min
integer :: i

level_shifting = (Control%SCFSHIFTVAL>0._prec.and.Control%SCFSHIFTBRK>0)
if(level_shifting.and.LPRINT>=1) then
   write(LOUT,'()')
   write(LOUT,'(a)') 'Level shifting procedure enabled'
   write(LOUT,'(a,es12.3)') 'VAL = ',Control%SCFSHIFTVAL
   write(LOUT,'(a,i4)') 'BRK = ',Control%SCFSHIFTBRK
endif

if(LPRINT>=2) then
   write(LOUT,'()')
   write(LOUT,'(5x,a)') 'Checking the overlap matrix'

   call symU_diagonalize_qd('QL',nbas,eval,evec,matS)
   call test_diagonalize_qd(nbas,eval,evec,matS)

   write(LOUT,'()')
   if(nbas>2*number_OVERLAP+1) then
      write(LOUT,'(a)') 'A few smallest and largest eigenvalues &
           &of the overlap matrix'
      do i=1,number_OVERLAP
         write(LOUT,'(i5,es18.6)') i,eval(i)
      enddo
      write(LOUT,'(2x,a)') '---'
      do i=nbas-number_OVERLAP+1,nbas
         write(LOUT,'(i5,es18.6)') i,eval(i)
      enddo
   else
      write(LOUT,'(a)') 'Eigenvalues of the overlap matrix'
      do i=1,nbas
         write(LOUT,'(i5,es18.6)') i,eval(i)
      enddo
   endif
endif

iter = 0

if(LPRINT>=2) then
   write(LOUT,'()')
   write(LOUT,'(5x,a)') 'Checking zeroth SCF iteration'
endif

call symU_diagonalize_qd('QL',nbas,eval,evec,matH,matS)
if(LPRINT>=2) call test_diagonalize_qd(nbas,eval,evec,matH,matS)
do i=1,norb
   call make_dens(nbas,nbas,evec(:,i),evec(:,i),dens,i==1)
enddo
energy = 2._prec*make_trace(nbas,nbas,matH,dens)

if(LPRINT>=2) then
   write(LOUT,'()')
   write(LOUT,'(a)') 'Orbital energies'
   do i=1,norb
      write(LOUT,'(i5,f25.18)') i,eval(i)
   enddo
endif

if(LPRINT>=1) then
   write(LOUT,'()')
   write(LOUT,'(1x,a,12x,a,12x,a,2x,a)') &
        'iter','energy','energy change','min-change counter'
   write(LOUT,'(i5,f25.18)') iter,energy
endif

iter = 1

cnt_min = 0
energy_delta_min = huge(0._prec)
converged = .false.
do
   energy_prev = energy

   call make_matG_1(2._prec,-1._prec,nbas,matG,twoel,dens)

   if(level_shifting.and.iter<=Control%SCFSHIFTBRK) then
      call make_matmul(nbas,matS,evec,matF)
      do i=norb+1,nbas
         call make_dens(nbas,nbas,matF(:,i),matF(:,i),dens,i==(norb+1))
      enddo
      matF(1:nbas,1:nbas) = matH(1:nbas,1:nbas) + matG(1:nbas,1:nbas) &
           + Control%SCFSHIFTVAL*dens(1:nbas,1:nbas)
   else
      matF(1:nbas,1:nbas) = matH(1:nbas,1:nbas) + matG(1:nbas,1:nbas)
   endif
   call symU_diagonalize_qd('QL',nbas,eval,evec,matF,matS)
   do i=1,norb
      call make_dens(nbas,nbas,evec(:,i),evec(:,i),dens,i==1)
   enddo

   energy = &
        + 2._prec*(make_trace(nbas,nbas,matH,dens) &
        + 0.5_prec*make_trace(nbas,nbas,matG,dens))

   energy_delta = abs(energy-energy_prev)
   if(energy_delta<energy_delta_min) then
      cnt_min = 0
      energy_delta_min = energy_delta
   else
      cnt_min = cnt_min + 1
   endif

   if(LPRINT>=1) write(LOUT,'(i5,f25.18,es18.6,6x,i5)') &
        iter,energy,energy_delta,cnt_min

   iter = iter + 1

   if(energy_delta<Control%SCFTHR) then
      converged = .true.
      exit
   endif
   if(iter>Control%SCFMAXIT) then
      write(LOUT,'(a,i5,a,es10.3,a)') 'WARNING!!! &
           &SCF did not converge in SCFMAXIT =', &
           Control%SCFMAXIT,' iterations!     (',energy_delta,' )'
      exit
   endif
   if(cnt_min>=Control%SCFCNT) then
      write(LOUT,'(a,i5,a,es10.3,a)') 'WARNING!!! &
           &No significant change in last SCFCNT =', &
           Control%SCFCNT,' iterations! (',energy_delta,' )'
      exit
   endif
enddo

if(LPRINT>=2) then
   write(LOUT,'()')
   write(LOUT,'(5x,a)') 'Checking last SCF iteration'

   call test_diagonalize_qd(nbas,eval,evec,matF,matS)
endif

end subroutine Roothaan_iterations

subroutine make_matmul(nbas,A,B,C)
implicit none
integer,intent(in) :: nbas
real(prec),intent(in) :: A(:,:),B(:,:)
real(prec),intent(out) :: C(:,:)
integer :: i,j,k
real(prec) :: rtmp

C = 0._prec

do j=1,nbas
   do k=1,nbas
      rtmp = B(k,j)
      do i=1,nbas
         C(i,j) = C(i,j) + A(i,k)*rtmp
      enddo
   enddo
enddo

end subroutine make_matmul

subroutine make_dens(nbas1,nbas2,vec1,vec2,dens,do_cleaning)
implicit none
integer,intent(in) :: nbas1,nbas2
real(prec),intent(in) :: vec1(:),vec2(:)
real(prec) :: dens(:,:)
logical,intent(in),optional :: do_cleaning
logical :: cleaning
integer :: i,j

if(present(do_cleaning)) then
   cleaning = do_cleaning
else
   cleaning = .true.
endif

if(cleaning) dens = 0._prec

do j=1,nbas2
   do i=1,nbas1
      dens(i,j) = dens(i,j) + vec1(i)*vec2(j)
   enddo
enddo

end subroutine make_dens

function make_trace(nbas1,nbas2,mat,dens) result(trace)
implicit none
real(prec) :: trace
integer,intent(in) :: nbas1,nbas2
real(prec),intent(in) :: mat(:,:)
real(prec),intent(in) :: dens(:,:)
integer :: i,j

trace = 0._prec

do j=1,nbas2
   do i=1,nbas1
      trace = trace + mat(i,j)*dens(i,j)
   enddo
enddo

end function make_trace

function make_contract(nbas_i,vec_i,nbas_j,vec_j,nbas_k,vec_k,nbas_l,vec_l,&
     twoel,dens) result(contract)
implicit none
real(prec) :: contract
integer,intent(in) :: nbas_i,nbas_j,nbas_k,nbas_l
real(prec),intent(in) :: vec_i(:),vec_j(:),vec_k(:),vec_l(:)
real(prec),intent(in) :: twoel(:,:)
real(prec) :: dens(:,:)
integer :: i,j,k,l,itwo,jtwo
real(prec) :: rtmp

contract = 0._prec

do j=1,nbas_j
   do i=1,nbas_i
      dens(i,j) = vec_i(i)*vec_j(j)
   enddo
enddo

jtwo = 0
do l=1,nbas_l
   do k=1,nbas_k
      jtwo = jtwo + 1
      rtmp = vec_k(k)*vec_l(l)

      itwo = 0
      do j=1,nbas_j
         do i=1,nbas_i
            itwo = itwo + 1

            contract = contract + dens(i,j)*twoel(itwo,jtwo)*rtmp

         enddo
      enddo

   enddo
enddo

end function make_contract

subroutine make_matG_1(coeffJ,coeffK,nbas,matG,twoel,dens)
implicit none
real(prec),intent(in) :: coeffJ,coeffK
integer,intent(in) :: nbas
real(prec) :: matG(:,:)
real(prec),intent(in) :: twoel(:,:),dens(:,:)
integer :: i,j,k,l,itwo,jtwo
real(prec) :: rtmp

matG = 0._prec

jtwo = 0
do l=1,nbas
   do k=1,nbas
      jtwo = jtwo + 1

      rtmp = coeffJ*dens(k,l)

      itwo = 0
      do j=1,nbas
         do i=1,nbas
            itwo = itwo + 1

            matG(i,j) = matG(i,j) + twoel(itwo,jtwo)*rtmp

         enddo
      enddo

   enddo
enddo

jtwo = 0
do l=1,nbas
   do j=1,nbas
      jtwo = jtwo + 1

      itwo = 0
      do k=1,nbas
         rtmp = coeffK*dens(k,l)
         do i=1,nbas
            itwo = itwo + 1

            matG(i,j) = matG(i,j) + twoel(itwo,jtwo)*rtmp

         enddo
      enddo

   enddo
enddo

end subroutine make_matG_1

subroutine make_matG_2(coeffX,nbas1,nbas2,matG,dimtype,twoel,dens,do_cleaning)
implicit none
real(prec),intent(in) :: coeffX
integer,intent(in) :: nbas1,nbas2
real(prec) :: matG(:,:)
integer,intent(in) :: dimtype
real(prec),intent(in) :: twoel(:,:),dens(:,:)
logical,intent(in),optional :: do_cleaning
logical :: cleaning
integer :: i,j,k,l,itwo,jtwo
real(prec) :: rtmp

if(present(do_cleaning)) then
   cleaning = do_cleaning
else
   cleaning = .true.
endif

if(cleaning) matG = 0._prec

select case(dimtype)
case(1)

   jtwo = 0
   do l=1,nbas2
      do k=1,nbas2
         jtwo = jtwo + 1

         rtmp = coeffX*dens(k,l)

         itwo = 0
         do j=1,nbas1
            do i=1,nbas1
               itwo = itwo + 1

               matG(i,j) = matG(i,j) + twoel(itwo,jtwo)*rtmp

            enddo
         enddo

      enddo
   enddo

case(2)

   jtwo = 0
   do l=1,nbas2
      do k=1,nbas2
         jtwo = jtwo + 1

         rtmp = 0._prec

         itwo = 0
         do j=1,nbas1
            do i=1,nbas1
               itwo = itwo + 1

               rtmp = rtmp + dens(i,j)*twoel(itwo,jtwo)

            enddo
         enddo

         matG(k,l) = matG(k,l) + coeffX*rtmp

      enddo
   enddo

case default

   write(LOUT,'(a,i5)') 'ERROR!!! &
        &Incorrect dimtype in make_matG_2: ',dimtype

end select

end subroutine make_matG_2

subroutine create_refJK(maxl)
! to be adapted to general case later
implicit none
integer,intent(in) :: maxl
integer :: i, j
real(prec) :: prefac

allocate(refJ(maxl+1,maxl+1), &
         refK(maxl+1,maxl+1))

do i=1,maxl+1
   do j=1,maxl+1
      prefac = 1._prec*j - 0.5_prec
      associate(refJ => refJ(i,j) )
      allocate(refJ%lambda(1),refJ%prefac(1))
      refJ%nlambda = 1
      refJ%lambda = 0
      refJ%prefac = prefac
      end associate
   enddo
enddo

do i=1,maxl+1
   do j=1,maxl+1
      associate(refK => refK(i,j) )
      allocate(refK%lambda((i+j)/2), &
               refK%prefac((i+j)/2) )
      refK%nlambda = (i+j)/2
      end associate
   enddo
enddo

refK(1,1)%lambda(1) = 0
refK(1,2)%lambda(1) = 1
refK(2,1)%lambda(1) = 1
refK(2,2)%lambda(1) = 0
refK(2,2)%lambda(2) = 2

refK(1,1)%prefac(1) = 0.5_prec
refK(1,2)%prefac(1) = 1.5_prec
refK(2,1)%prefac(1) = 0.5_prec
refK(2,2)%prefac(1) = 0.5_prec
refK(2,2)%prefac(2) = 1._prec

!deallocate(refJ,refK)

end subroutine create_refJK

subroutine free_refJK(maxl)
implicit none
integer,intent(in) :: maxl
integer :: i, j

! dellocate refJK
do i=1,maxl+1
   do j=1,maxl+1
      associate(refK => refK(i,j), &
                refJ => refJ(i,j) )
      deallocate(refK%lambda, refK%prefac, &
                 refJ%lambda, refJ%prefac )
      end associate
   enddo
enddo
deallocate(refJ,refK)

end subroutine free_refJK

subroutine make_over(nbas1,nbas2,p,dimtype,matS,vec)
implicit none
integer,intent(in) :: nbas1,nbas2
real(prec) :: p(:)
integer,intent(in) :: dimtype
real(prec),intent(in) :: matS(:,:),vec(:)
integer :: i,j

p = 0._prec

select case(dimtype)
case(1)

   do j=1,nbas2
      do i=1,nbas1
         p(i) = p(i) + matS(i,j)*vec(j)
      enddo
   enddo

case(2)

   do j=1,nbas2
      do i=1,nbas1
         p(j) = p(j) + vec(i)*matS(i,j)
      enddo
   enddo

case default

   write(LOUT,'(a,i5)') 'ERROR!!! &
        &Incorrect dimtype in make_over: ',dimtype

end select

end subroutine make_over

subroutine transform(nbas,p,mat,mat_mod)
implicit none
integer,intent(in) :: nbas
real(prec),intent(in) :: p(:),mat(:,:)
real(prec) :: mat_mod(:,:)
integer :: i,j

mat_mod = 0._prec

do j=1,nbas
   do i=1,nbas-1
      mat_mod(i,j) = p(i+1)*mat(i,j) - p(i)*mat(i+1,j)
   enddo
enddo

do j=1,nbas-1
   do i=1,nbas-1
      mat_mod(i,j) = mat_mod(i,j)*p(j+1) - mat_mod(i,j+1)*p(j)
   enddo
enddo
mat_mod(:,nbas) = 0._prec

end subroutine transform

subroutine backtransform(nbas,p,vec,vec_mod)
implicit none
integer,intent(in) :: nbas
real(prec),intent(in) :: p(:),vec(:)
real(prec) :: vec_mod(:)
integer :: i

vec_mod = 0._prec

vec_mod(1) = p(2)*vec(1)
do i=2,nbas-1
   vec_mod(i) = -p(i-1)*vec(i-1) + p(i+1)*vec(i) 
enddo
vec_mod(nbas) = -p(nbas-1)*vec(nbas-1)

end subroutine backtransform

end module SCFdriver
