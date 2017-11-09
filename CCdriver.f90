module CCdriver
use file_OUT, only : LOUT
use precision, only : prec
use memory
use lproblem
use commontypes
use CCint2
use CCint3
implicit none

private
! public CC_fix_SCForbitals,CC_release_SCForbitals
! public calculateCC,energyCC
! 
! type(SCForbitalsData),pointer :: SCF
! 
! integer :: Nexp
! real(prec),allocatable :: exponents(:)
! 
! real(prec) :: RESULT_energy
! 
! type CCpairData
! integer :: nbas
! real(prec) :: multE,multN
! real(prec),allocatable :: vecE(:),vecN(:)
! real(prec),allocatable :: tau(:),Qtau(:),prev(:)
! real(prec),allocatable :: energy(:,:,:),orthog(:,:,:)
! integer :: DIIS_size
! real(prec),allocatable :: DIIS_vec(:,:)
! end type CCpairData
! 
! type OSdata
! integer :: size
! integer :: start,end
! end type OSdata
! 
! integer,parameter :: power_LIMIT = -16
! 
! contains
! 
! subroutine CC_fix_SCForbitals(SCForbitals)
! implicit none
! type(SCForbitalsData),target :: SCForbitals
! 
! call checkSCF(.false.)
! 
! SCF => SCForbitals
! 
! end subroutine CC_fix_SCForbitals
! 
! subroutine CC_release_SCForbitals
! implicit none
! 
! call checkSCF(.true.)
! 
! nullify(SCF)
! 
! end subroutine CC_release_SCForbitals
! 
! subroutine calculateCC(System,Control)
! implicit none
! type(SystemData) :: System
! type(ControlData) :: Control
! 
! call checkSCF(.true.)
! 
! Nexp = System%Nexp
! call mem_alloc(exponents,Nexp)
! exponents(:) = System%exponents
! 
! RESULT_energy = 0._prec
! 
! call chooseCC(System,Control,.true.)
! 
! call mem_dealloc(exponents)
! 
! end subroutine calculateCC
! 
! function energyCC(MODexponents,System,Control) result(energy)
! use misc, only : unpack_vector
! implicit none
! real(prec) :: energy
! real(prec) :: MODexponents(:)
! type(SystemData) :: System
! type(ControlData) :: Control
! 
! call checkSCF(.true.)
! 
! Nexp = System%Nexp
! call mem_alloc(exponents,Nexp)
! exponents(:) = System%exponents
! call unpack_vector(MODexponents,Nexp,exponents,System%isOpt_pairs)
! 
! RESULT_energy = 0._prec
! 
! call chooseCC(System,Control,.false.)
! 
! energy = RESULT_energy
! 
! call mem_dealloc(exponents)
! 
! end function energyCC
! 
! subroutine checkSCF(expected)
! implicit none
! logical,intent(in) :: expected
! 
! if(associated(SCF).neqv.expected) then
! 
!    if(expected) then
!       write(LOUT,'(a)') 'ERROR!!! &
!            &SCF data for CC calculations has not been associated!'
!    else
!       write(LOUT,'(a)') 'ERROR!!! &
!            &SCF data for CC calculations has already been associated!'
!    endif
! 
!    stop
! 
! endif
! 
! end subroutine checkSCF
! 
! subroutine chooseCC(System,Control,fullPRINT)
! implicit none
! type(SystemData),intent(in) :: System
! type(ControlData),intent(in) :: Control
! logical,intent(in) :: fullPRINT
! 
! if(System%common_pairs) then
!    if(System%n_pairs/=1) then
!       write(LOUT,'(a)') 'ERROR!!! Inconsistent data in chooseCC!'
!       stop
!    endif
! 
!    call CC_energy_common(System,Control,fullPRINT)
! 
! else
! 
!    write(LOUT,'(a)') 'NOT YET!!!'
!    stop
! 
! endif
! 
! end subroutine chooseCC
! 
! subroutine CC_energy_common(System,Control,fullPRINT)
! use precision, only : dble
! use time
! implicit none
! type(SystemData),intent(in) :: System
! type(ControlData),intent(in) :: Control
! logical,intent(in) :: fullPRINT
! integer :: LPRINT
! integer :: stage
! integer :: MAXIT
! integer :: norb,nbas
! integer :: ieta,power,iter
! integer :: ipair,jpair,ijpair,iorb,jorb,korb
! real(prec) :: eta,factor
! real(prec) :: energy_this,energy_prev,edelta_this,edelta_prev
! real(prec),allocatable :: matS_S(:,:),matS_T(:,:)
! real(prec),allocatable :: matF_S(:,:),matF_T(:,:)
! real(prec),allocatable :: matP1_S(:,:),matP1_T(:,:)
! real(prec),allocatable :: matPe_S(:,:),matPe_T(:,:)
! real(prec),allocatable :: matP2_S(:,:),matP2_T(:,:)
! real(prec),allocatable :: matL1_S(:,:),matL1_T(:,:)
! real(prec),allocatable :: matJ_S(:,:,:,:),matJ_T(:,:,:,:)
! real(prec),allocatable :: matK_S(:,:,:,:),matK_T(:,:,:,:)
! real(prec),allocatable :: matM_S(:,:,:,:),matM_T(:,:,:,:)
! real(prec),allocatable :: vec0(:,:,:),vec1(:,:,:),vecP(:,:,:)
! real(prec),allocatable :: pair_energy(:,:,:,:)
! real(prec),allocatable :: TMP(:,:)
! type(CCpairData),allocatable :: CCpairs(:,:)
! real(prec),allocatable :: energy(:,:,:),orthog(:,:,:)
! integer,allocatable :: final_iter(:)
! type(TripletData) :: Triplet
! real(prec),allocatable :: LHS(:,:),RHS(:)
! type(DecompositionData) :: project_S,project_T
! type(DecompositionData),allocatable :: pairLHS(:,:)
! logical :: do_OS
! integer :: OS_size
! type(OSdata),allocatable :: OS(:,:)
! integer :: DIIS_start,DIIS_size
! logical :: DIIS
! integer :: DIIS_off,DIIS_n
! real(dble) :: Tcpu,Twall
! 
! LPRINT = merge(Control%LPRINT,0,fullPRINT)
! 
! select case(trim(System%calc_type))
! case('FCCD','FCCD-L')
!    stage = 3
! case('LCCD','LCCD-L')
!    stage = 2
! case('CID','CID-L')
!    stage = 1
! case default
!    stage = 0
! end select
! 
! do_OS = (index(System%calc_type,'-L')/=0)
! 
! MAXIT = merge(Control%ENERGMAXIT,1,stage>0)
! 
! if(do_OS) then
!    DIIS_start = 3
!    DIIS_size  = 6
! else
!    DIIS_start = -1
!    DIIS_size  = -1
! endif
! DIIS = (stage>0).and.(DIIS_start>1).and.(DIIS_size>1)
! if(.not.DIIS) DIIS_size = -1
! 
! norb = System%norb
! 
! associate(PairSystem => System%PairSystem(1,1))
! 
!   nbas = PairSystem%nbas
! 
!   call mem_alloc(matS_S,nbas,nbas);  call mem_alloc(matS_T,nbas,nbas)
!   call mem_alloc(matF_S,nbas,nbas);  call mem_alloc(matF_T,nbas,nbas)
!   call mem_alloc(matP1_S,nbas,nbas); call mem_alloc(matP1_T,nbas,nbas)
!   call mem_alloc(matPe_S,nbas,nbas); call mem_alloc(matPe_T,nbas,nbas)
!   call mem_alloc(matP2_S,nbas,nbas); call mem_alloc(matP2_T,nbas,nbas)
!   if(stage>0) then
!      call mem_alloc(matL1_S,nbas,nbas); call mem_alloc(matL1_T,nbas,nbas)
!      call mem_alloc(matJ_S,nbas,nbas,norb,norb)
!      call mem_alloc(matJ_T,nbas,nbas,norb,norb)
!      call mem_alloc(matK_S,nbas,nbas,norb,norb)
!      call mem_alloc(matK_T,nbas,nbas,norb,norb)
!      call mem_alloc(matM_S,nbas,nbas,norb,norb)
!      call mem_alloc(matM_T,nbas,nbas,norb,norb)
!   endif
! 
!   call mem_alloc(vec0,nbas,norb,norb)
!   call mem_alloc(vec1,nbas,norb,norb)
!   call mem_alloc(vecP,nbas,norb,norb)
! 
!   if(stage>2) call mem_alloc(pair_energy,norb,norb,norb,norb)
! 
!   call mem_alloc(TMP,nbas,nbas)
! 
!   if(fullPRINT) then
!      write(LOUT,'()')
!      call timer('START',Tcpu,Twall)
!   endif
! 
!   call create_CCint2(Nexp,exponents,System%OrbPairReduced,LPRINT)
!   if(fullPRINT) then
!      call timer('2-el integrals init',Tcpu,Twall)
!      flush(LOUT)
!   endif
! 
!   call CCint2_matS('A',0,TMP,PairSystem,PairSystem)
!   matS_S(:,:) = TMP
!   matS_T(:,:) = TMP
!   call CCint2_matS('B',0,TMP,PairSystem,PairSystem)
!   matS_S(:,:) = matS_S + TMP
!   matS_T(:,:) = matS_T - TMP
! 
!   if(stage>0) then
!      call CCint2_matS('A',-1,TMP,PairSystem,PairSystem)
!      matL1_S(:,:) = TMP
!      matL1_T(:,:) = TMP
!      call CCint2_matS('B',-1,TMP,PairSystem,PairSystem)
!      matL1_S(:,:) = matL1_S + TMP
!      matL1_T(:,:) = matL1_T - TMP
!   endif
! 
!   call CCint2_matH('A',TMP,PairSystem,PairSystem,System%nucZ)
!   matF_S(:,:) = TMP
!   matF_T(:,:) = TMP
!   call CCint2_matH('B',TMP,PairSystem,PairSystem,System%nucZ)
!   matF_S(:,:) = matF_S + TMP
!   matF_T(:,:) = matF_T - TMP
! 
!   do jorb=1,norb
!      do iorb=1,norb
!         associate(&
!              i_vect => SCF%orb_vector(:,iorb),&
!              j_vect => SCF%orb_vector(:,jorb),&
!              i_OrbSystem => System%OrbSystem(min(iorb,System%n_orbs)),&
!              j_OrbSystem => System%OrbSystem(min(jorb,System%n_orbs)))
! 
!           call CCint2_vec(0,vec0(:,iorb,jorb),PairSystem,&
!                i_vect,i_OrbSystem,j_vect,j_OrbSystem)
! 
!           call CCint2_vec(-1,vec1(:,iorb,jorb),PairSystem,&
!                i_vect,i_OrbSystem,j_vect,j_OrbSystem)
! 
!         end associate
!      enddo
!   enddo
! 
!   call free_CCint2
!   if(fullPRINT) then
!      call timer('2-el integrals use',Tcpu,Twall)
!      flush(LOUT)
!   endif
! 
!   call create_CCint3(Nexp,exponents,System%OrbReduced,System%PairReduced,&
!        Control%INT3,LPRINT)
!   if(fullPRINT) then
!      call timer('3-el integrals init',Tcpu,Twall)
!      flush(LOUT)
!   endif
! 
!   if(stage>0) then
! 
!      do jorb=1,norb
!         do iorb=1,norb
!            associate(&
!                 i_vect => SCF%orb_vector(:,iorb),&
!                 i_OrbSystem => System%OrbSystem(min(iorb,System%n_orbs)),&
!                 j_vect => SCF%orb_vector(:,jorb),&
!                 j_OrbSystem => System%OrbSystem(min(jorb,System%n_orbs)))
! 
!              call CCint3_matJ('A',TMP,PairSystem,PairSystem,&
!                   i_vect,i_OrbSystem,j_vect,j_OrbSystem)
!              matJ_S(:,:,iorb,jorb) = TMP
!              matJ_T(:,:,iorb,jorb) = TMP
!              if(iorb==jorb) then
!                 matF_S(:,:) = matF_S + TMP*2._prec
!                 matF_T(:,:) = matF_T + TMP*2._prec
!              endif
!              call CCint3_matJ('B',TMP,PairSystem,PairSystem,&
!                   i_vect,i_OrbSystem,j_vect,j_OrbSystem)
!              matJ_S(:,:,iorb,jorb) = matJ_S(:,:,iorb,jorb) + TMP
!              matJ_T(:,:,iorb,jorb) = matJ_T(:,:,iorb,jorb) - TMP
!              if(iorb==jorb) then
!                 matF_S(:,:) = matF_S + TMP*2._prec
!                 matF_T(:,:) = matF_T - TMP*2._prec
!              endif
! 
!              call CCint3_matK('A',TMP,PairSystem,PairSystem,&
!                   i_vect,i_OrbSystem,j_vect,j_OrbSystem)
!              matK_S(:,:,iorb,jorb) = TMP
!              matK_T(:,:,iorb,jorb) = TMP
!              if(iorb==jorb) then
!                 matF_S(:,:) = matF_S - TMP
!                 matF_T(:,:) = matF_T - TMP
!              endif
!              call CCint3_matK('B',TMP,PairSystem,PairSystem,&
!                   i_vect,i_OrbSystem,j_vect,j_OrbSystem)
!              matK_S(:,:,iorb,jorb) = matK_S(:,:,iorb,jorb) + TMP
!              matK_T(:,:,iorb,jorb) = matK_T(:,:,iorb,jorb) - TMP
!              if(iorb==jorb) then
!                 matF_S(:,:) = matF_S - TMP
!                 matF_T(:,:) = matF_T + TMP
!              endif
! 
!              call CCint3_matM('A',TMP,PairSystem,PairSystem,&
!                   i_vect,i_OrbSystem,j_vect,j_OrbSystem)
!              matM_S(:,:,iorb,jorb) = TMP
!              matM_T(:,:,iorb,jorb) = TMP
!              call CCint3_matM('B',TMP,PairSystem,PairSystem,&
!                   i_vect,i_OrbSystem,j_vect,j_OrbSystem)
!              matM_S(:,:,iorb,jorb) = matM_S(:,:,iorb,jorb) - TMP
!              matM_T(:,:,iorb,jorb) = matM_T(:,:,iorb,jorb) + TMP
! 
!            end associate
!         enddo
!      enddo
! 
!   else
! 
!      do korb=1,norb
!         associate(&
!              k_vect => SCF%orb_vector(:,korb),&
!              k_OrbSystem => System%OrbSystem(min(korb,System%n_orbs)))
! 
!           call CCint3_matJ('A',TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!           matF_S(:,:) = matF_S + TMP*2._prec
!           matF_T(:,:) = matF_T + TMP*2._prec
!           call CCint3_matJ('B',TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!           matF_S(:,:) = matF_S + TMP*2._prec
!           matF_T(:,:) = matF_T - TMP*2._prec
! 
!           call CCint3_matK('A',TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!           matF_S(:,:) = matF_S - TMP
!           matF_T(:,:) = matF_T - TMP
!           call CCint3_matK('B',TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!           matF_S(:,:) = matF_S - TMP
!           matF_T(:,:) = matF_T + TMP
! 
!         end associate
!      enddo
! 
!   endif
! 
!   matP1_S = 0._prec
!   matP1_T = 0._prec
!   matPe_S = 0._prec
!   matPe_T = 0._prec
!   do korb=1,norb
!      associate(&
!           k_ener => SCF%orb_energy(korb),&
!           k_vect => SCF%orb_vector(:,korb),&
!           k_OrbSystem => System%OrbSystem(min(korb,System%n_orbs)))
! 
!        call CCint3_matP('A',0,TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!        matP1_S(:,:) = matP1_S + TMP
!        matP1_T(:,:) = matP1_T + TMP
!        matPe_S(:,:) = matPe_S + TMP*k_ener
!        matPe_T(:,:) = matPe_T + TMP*k_ener
!        call CCint3_matP('B',0,TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!        matP1_S(:,:) = matP1_S + TMP
!        matP1_T(:,:) = matP1_T - TMP
!        matPe_S(:,:) = matPe_S + TMP*k_ener
!        matPe_T(:,:) = matPe_T - TMP*k_ener
! 
!        if(stage>0) then
!           call CCint3_matP('A',-1,TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!           matL1_S(:,:) = matL1_S - TMP
!           matL1_T(:,:) = matL1_T - TMP
!           call CCint3_matP('B',-1,TMP,PairSystem,PairSystem,k_vect,k_OrbSystem)
!           matL1_S(:,:) = matL1_S - TMP
!           matL1_T(:,:) = matL1_T + TMP
!        endif
! 
!      end associate
!   enddo
! 
!   vecP = 0._prec
!   do jorb=1,norb
!      do iorb=1,norb
!         associate(&
!              i_vect => SCF%orb_vector(:,iorb),&
!              j_vect => SCF%orb_vector(:,jorb),&
!              i_OrbSystem => System%OrbSystem(min(iorb,System%n_orbs)),&
!              j_OrbSystem => System%OrbSystem(min(jorb,System%n_orbs)))
!           do korb=1,norb
!              associate(&
!                   k_vect => SCF%orb_vector(:,korb),&
!                   k_OrbSystem => System%OrbSystem(min(korb,System%n_orbs)))
! 
!                call CCint3_vecP(TMP(:,1),PairSystem,&
!                     i_vect,i_OrbSystem,j_vect,j_OrbSystem,&
!                     k_vect,k_OrbSystem)
!                vecP(:,iorb,jorb) = vecP(:,iorb,jorb) + TMP(:,1)
! 
!              end associate
!           enddo
!         end associate
!      enddo
!   enddo
! 
!   call free_CCint3
!   if(fullPRINT) then
!      call timer('3-el integrals use',Tcpu,Twall)
!      flush(LOUT)
!   endif
! 
!   matP2_S = 0._prec
!   matP2_T = 0._prec
!   do jorb=1,norb
!      do iorb=1,norb
! 
!         call outerMat(nbas,1._prec,vec0(:,iorb,jorb),vec0(:,iorb,jorb),TMP)
!         matP2_S(:,:) = matP2_S + TMP
!         if(iorb/=jorb) matP2_T(:,:) = matP2_T + TMP
!         call outerMat(nbas,1._prec,vec0(:,iorb,jorb),vec0(:,jorb,iorb),TMP)
!         matP2_S(:,:) = matP2_S + TMP
!         if(iorb/=jorb) matP2_T(:,:) = matP2_T - TMP
! 
!         if(stage>0) then
!            call outerMat(nbas,1._prec,vec0(:,iorb,jorb),vec1(:,iorb,jorb),TMP)
!            matL1_S(:,:) = matL1_S + TMP
!            if(iorb/=jorb) matL1_T(:,:) = matL1_T + TMP
!            call outerMat(nbas,1._prec,vec0(:,iorb,jorb),vec1(:,jorb,iorb),TMP)
!            matL1_S(:,:) = matL1_S + TMP
!            if(iorb/=jorb) matL1_T(:,:) = matL1_T - TMP
!         endif
! 
!      enddo
!   enddo
! 
!   call mem_dealloc(TMP)
! 
!   allocate(CCpairs(norb,norb))
!   ijpair = 0
!   do jpair=1,norb
!      do ipair=1,jpair
!         ijpair = ijpair + 1
! 
!         associate(CCpair => CCpairs(jpair,ipair))
! 
!           call init_CCpair(CCpair,nbas,MAXIT,System%Neta,DIIS_size)
! 
!           CCpair%multE = merge(0.5_prec,1._prec,ipair==jpair)
!           CCpair%multN = -CCpair%multE
! 
!           CCpair%vecE(:) = vec1(:,ipair,jpair) + vec1(:,jpair,ipair)
!           CCpair%vecN(:) = &
!                - vec1(:,ipair,jpair) - vec1(:,jpair,ipair) &
!                + vecP(:,ipair,jpair) + vecP(:,jpair,ipair)
!           do jorb=1,norb
!              do iorb=1,norb
!                 CCpair%vecN(:) = CCpair%vecN - vec0(:,iorb,jorb)&
!                      * (SCF%pair_energy(iorb,jorb,ijpair) &
!                      +  SCF%pair_energy(jorb,iorb,ijpair))
!              enddo
!           enddo
! 
!         end associate
! 
!         if(ipair==jpair) cycle
! 
!         associate(CCpair => CCpairs(ipair,jpair))
! 
!           call init_CCpair(CCpair,nbas,MAXIT,System%Neta,DIIS_size)
! 
!           CCpair%multE = 3._prec
!           CCpair%multN = -CCpair%multE
! 
!           CCpair%vecE(:) = vec1(:,ipair,jpair) - vec1(:,jpair,ipair)
!           CCpair%vecN(:) = &
!                - vec1(:,ipair,jpair) + vec1(:,jpair,ipair) &
!                + vecP(:,ipair,jpair) - vecP(:,jpair,ipair)
!           do jorb=1,norb
!              do iorb=1,norb
!                 if(iorb==jorb) cycle
!                 CCpair%vecN(:) = CCpair%vecN - vec0(:,iorb,jorb)&
!                      * (SCF%pair_energy(iorb,jorb,ijpair) &
!                      -  SCF%pair_energy(jorb,iorb,ijpair))
!              enddo
!           enddo
! 
!         end associate
! 
!      enddo
!   enddo
! 
!   call mem_alloc(energy,3,MAXIT,System%Neta)
!   call mem_alloc(orthog,2,MAXIT,System%Neta)
!   call mem_alloc(final_iter,System%Neta)
! 
!   call init_Triplet(Triplet,PairSystem)
!   if(fullPRINT.and.norb>1) then
!      write(LOUT,'()')
!      write(LOUT,'(2(a,i5))') &
!           'Triplet reduction: ',Triplet%nbas_orig,' -> ',Triplet%nbas
!   endif
! 
!   call mem_alloc(LHS,nbas,nbas)
!   call mem_alloc(RHS,nbas)
! 
!   LHS(:,:) = matS_S
!   if(fullPRINT) call symmetry_check(0,0,nbas,LHS)
!   call init_linsolver(project_S,'U','LDL1',nbas,LHS)
!   LHS(:,:) = matS_T
!   call shrink_Triplet(LHS,Triplet)
!   if(fullPRINT) call symmetry_check(0,0,Triplet%nbas,LHS)
!   call init_linsolver(project_T,'U','LDL1',Triplet%nbas,LHS)
! 
!   if(fullPRINT) flush(LOUT)
! 
!   if(stage>1) allocate(pairLHS(norb,norb))
! 
!   do ieta=1,System%Neta
! 
!      power = System%eta(ieta)
!      if(power>power_LIMIT) then
!         eta = 10._prec**power
!      else
!         eta = 0._prec
!      endif
!      if(fullPRINT) call print_eta_full(power)
! 
!      iter = 1
! 
!      if(DIIS) then
!         DIIS_off = 0
!         DIIS_n   = 0
!      endif
! 
!      if(do_OS.and.ieta==1) then
! 
!         allocate(OS(norb,norb))
! 
!         OS_size = 0
!         do jpair=1,norb
!            do ipair=1,norb
!               associate(OS => OS(ipair,jpair))
! 
!                 if(ipair>=jpair) then
!                    OS%size = nbas
!                 else
!                    OS%size = Triplet%nbas
!                 endif
! 
!                 OS%start = OS_size + 1
!                 OS%end   = OS_size + OS%size
!                 OS_size  = OS_size + OS%size
! 
!               end associate
!            enddo
!         enddo
! 
!         call mem_alloc(TMP,OS_size,OS_size+2)
!         TMP = 0._prec
! 
!         ijpair = 0
!         do jpair=1,norb
!            do ipair=1,jpair
!               ijpair = ijpair + 1
! 
!               associate(OS1 => OS(jpair,ipair))
! 
!                 RHS(:) = CCpairs(jpair,ipair)%vecN
!                 TMP(OS1%start:OS1%end,OS_size+1) = RHS(1:nbas)
! 
!                 call MWO(LHS,jpair,ipair,eta,&
!                      matF_S,matS_S,matP1_S,matPe_S,matP2_S)
!                 LHS(:,:) = LHS + matL1_S
!                 TMP(OS1%start:OS1%end,OS1%start:OS1%end) = LHS(1:nbas,1:nbas)
! 
!                 LHS(:,:) = matS_S
!                 do jorb=1,norb
!                    do iorb=1,jorb
!                       factor = &
!                            SCF%pair_energy(iorb,jorb,ijpair) + &
!                            SCF%pair_energy(jorb,iorb,ijpair)
!                       if(iorb==jorb) factor = factor/2._prec
!                       associate(OS2 => OS(jorb,iorb))
!                         TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                              TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                              + factor*LHS(1:nbas,1:nbas)
!                       end associate
!                    enddo
!                 enddo
! 
!                 do korb=1,norb
!                    LHS(:,:) = matJ_S(:,:,korb,jpair) &
!                         - 0.5_prec*matK_S(:,:,korb,jpair)
!                    associate(OS2 => OS(max(ipair,korb),min(ipair,korb)))
!                      TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                           TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                           - LHS(1:nbas,1:nbas)
!                    end associate
!                    if(korb/=ipair) then
!                       factor =  1.5_prec*merge( 1._prec,-1._prec,ipair<korb)
!                       LHS(:,:) = matM_S(:,:,korb,jpair)
!                       call shrink_Triplet(LHS,Triplet,2)
!                       associate(OS2 => OS(min(ipair,korb),max(ipair,korb)))
!                         TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                              TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                              - factor*LHS(1:nbas,1:Triplet%nbas)
!                       end associate
!                    endif
!                    LHS(:,:) = matJ_S(:,:,korb,ipair) &
!                         - 0.5_prec*matK_S(:,:,korb,ipair)
!                    associate(OS2 => OS(max(korb,jpair),min(korb,jpair)))
!                      TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                           TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                           - LHS(1:nbas,1:nbas)
!                    end associate
!                    if(korb/=jpair) then
!                       factor = -1.5_prec*merge( 1._prec,-1._prec,korb<jpair)
!                       LHS(:,:) = matM_S(:,:,korb,ipair)
!                       call shrink_Triplet(LHS,Triplet,2)
!                       associate(OS2 => OS(min(korb,jpair),max(korb,jpair)))
!                         TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                              TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                              - factor*LHS(1:nbas,1:Triplet%nbas)
!                       end associate
!                    endif
!                 enddo
! 
!               end associate
! 
!               if(ipair==jpair) cycle
! 
!               associate(OS1 => OS(ipair,jpair))
! 
!                 RHS(:) = CCpairs(ipair,jpair)%vecN
!                 call shrink_Triplet(RHS,Triplet)
!                 TMP(OS1%start:OS1%end,OS_size+1) = RHS(1:Triplet%nbas)
! 
!                 call MWO(LHS,ipair,jpair,eta,&
!                      matF_T,matS_T,matP1_T,matPe_T,matP2_T)
!                 LHS(:,:) = LHS + matL1_T
!                 call shrink_Triplet(LHS,Triplet)
!                 TMP(OS1%start:OS1%end,OS1%start:OS1%end) = &
!                      LHS(1:Triplet%nbas,1:Triplet%nbas)
! 
!                 LHS(:,:) = matS_T
!                 call shrink_Triplet(LHS,Triplet)
!                 do jorb=2,norb
!                    do iorb=1,jorb-1
!                       factor = &
!                            SCF%pair_energy(iorb,jorb,ijpair) - &
!                            SCF%pair_energy(jorb,iorb,ijpair)
!                       associate(OS2 => OS(iorb,jorb))
!                         TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                              TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                              + factor*LHS(1:Triplet%nbas,1:Triplet%nbas)
!                       end associate
!                    enddo
!                 enddo
! 
!                 do korb=1,norb
!                    if(korb/=ipair) then
!                       factor = merge( 1._prec,-1._prec,ipair<korb)
!                       LHS(:,:) = matJ_T(:,:,korb,jpair) &
!                            - 1.5_prec*matK_T(:,:,korb,jpair)
!                       call shrink_Triplet(LHS,Triplet)
!                       associate(OS2 => OS(min(ipair,korb),max(ipair,korb)))
!                         TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                              TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                              - factor*LHS(1:Triplet%nbas,1:Triplet%nbas)
!                       end associate
!                    endif
!                    LHS(:,:) =  0.5_prec*matM_T(:,:,korb,jpair)
!                    call shrink_Triplet(LHS,Triplet,1)
!                    associate(OS2 => OS(max(ipair,korb),min(ipair,korb)))
!                      TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                           TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                           - LHS(1:Triplet%nbas,1:nbas)
!                    end associate
!                    if(korb/=jpair) then
!                       factor = merge( 1._prec,-1._prec,korb<jpair)
!                       LHS(:,:) = matJ_T(:,:,korb,ipair) &
!                            - 1.5_prec*matK_T(:,:,korb,ipair)
!                       call shrink_Triplet(LHS,Triplet)
!                       associate(OS2 => OS(min(korb,jpair),max(korb,jpair)))
!                         TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                              TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                              - factor*LHS(1:Triplet%nbas,1:Triplet%nbas)
!                       end associate
!                    endif
!                    LHS(:,:) = -0.5_prec*matM_T(:,:,korb,ipair)
!                    call shrink_Triplet(LHS,Triplet,1)
!                    associate(OS2 => OS(max(korb,jpair),min(korb,jpair)))
!                      TMP(OS1%start:OS1%end,OS2%start:OS2%end) = &
!                           TMP(OS1%start:OS1%end,OS2%start:OS2%end) &
!                           - LHS(1:Triplet%nbas,1:nbas)
!                    end associate
!                 enddo
! 
!               end associate
! 
!            enddo
!         enddo
! 
!         call simple_linsolver('G','GAUSS',OS_size,&
!              TMP(:,1:OS_size),TMP(:,OS_size+1),TMP(:,OS_size+2))
! 
!         do jpair=1,norb
!            do ipair=1,norb
!               associate(OS => OS(ipair,jpair),CCpair => CCpairs(ipair,jpair))
!                 if(ipair>=jpair) then
! 
!                    CCpair%tau(1:nbas) = TMP(OS%start:OS%end,OS_size+2)
! 
!                 else
! 
!                    CCpair%tau(1:Triplet%nbas) = TMP(OS%start:OS%end,OS_size+2)
!                    call expand_Triplet(CCpair%tau,Triplet)
! 
!                 endif
!               end associate
!            enddo
!         enddo
! 
!         call mem_dealloc(TMP)
! 
!         deallocate(OS)
! 
!      endif
!      if(do_OS) then
! 
!         call summarize_iteration
! 
!         energy_prev = energy(3,iter,ieta)
!         edelta_prev = abs(energy_prev)
!         iter = iter + 1
! 
!      endif
! 
!      if(stage>1) then
!         do jpair=1,norb
!            do ipair=1,jpair
! 
!               call MWO(LHS,jpair,ipair,eta,&
!                    matF_S,matS_S,matP1_S,matPe_S,matP2_S)
!               if(fullPRINT) call symmetry_check(jpair,ipair,nbas,LHS)
!               call init_linsolver(pairLHS(jpair,ipair),&
!                    'U','LDL2',nbas,LHS)
! 
!               if(ipair==jpair) cycle
! 
!               call MWO(LHS,ipair,jpair,eta,&
!                    matF_T,matS_T,matP1_T,matPe_T,matP2_T)
!               call shrink_Triplet(LHS,Triplet)
!               if(fullPRINT) call symmetry_check(ipair,jpair,Triplet%nbas,LHS)
!               call init_linsolver(pairLHS(ipair,jpair),&
!                    'U','LDL2',Triplet%nbas,LHS)
! 
!            enddo
!         enddo
!      endif
! 
!      if(fullPRINT) flush(LOUT)
! 
!      do
! 
!         if(stage>2.and.iter>1) then
!            do jpair=1,norb
!               do ipair=1,jpair
! 
!                  associate(CCpair => CCpairs(jpair,ipair),&
!                       elms => pair_energy(:,:,jpair,ipair))
!                    elms = 0._prec
!                    do jorb=1,norb
!                       do iorb=1,jorb
!                          factor = &
!                               inner(nbas,vec1(:,iorb,jorb),CCpair%prev) + &
!                               inner(nbas,vec1(:,jorb,iorb),CCpair%prev)
!                          elms(iorb,jorb) = factor
!                          elms(jorb,iorb) = factor
!                       enddo
!                    enddo
!                  end associate
! 
!                  if(ipair==jpair) cycle
! 
!                  associate(CCpair => CCpairs(ipair,jpair),&
!                       elms => pair_energy(:,:,ipair,jpair))
!                    elms = 0._prec
!                    do jorb=2,norb
!                       do iorb=1,jorb-1
!                          factor = &
!                               inner(nbas,vec1(:,iorb,jorb),CCpair%prev) - &
!                               inner(nbas,vec1(:,jorb,iorb),CCpair%prev)
!                          elms(iorb,jorb) =  factor
!                          elms(jorb,iorb) = -factor
!                       enddo
!                    enddo
!                  end associate
! 
!               enddo
!            enddo
!         endif
! 
!         ijpair = 0
!         do jpair=1,norb
!            do ipair=1,jpair
!               ijpair = ijpair + 1
! 
!               associate(CCpair => CCpairs(jpair,ipair))
! 
!                 RHS(:) = CCpair%vecN
! 
!                 if(stage>0.and.iter>1) then
! 
!                    call MatVec(nbas,-1._prec,matL1_S,CCpair%prev,RHS,.false.)
! 
!                    do jorb=1,norb
!                       do iorb=1,jorb
!                          factor = &
!                               SCF%pair_energy(iorb,jorb,ijpair) + &
!                               SCF%pair_energy(jorb,iorb,ijpair)
!                          if(iorb==jorb) factor = factor/2._prec
!                          call MatVec(nbas,-factor,matS_S,&
!                               CCpairs(jorb,iorb)%prev,RHS,.false.)
!                       enddo
!                    enddo
! 
!                    do korb=1,norb
!                       call MatVec(nbas, 1._prec,&
!                            matJ_S(:,:,korb,jpair),&
!                            CCpairs(max(ipair,korb),min(ipair,korb))%prev,&
!                            RHS,.false.)
!                       call MatVec(nbas,-0.5_prec,&
!                            matK_S(:,:,korb,jpair),&
!                            CCpairs(max(ipair,korb),min(ipair,korb))%prev,&
!                            RHS,.false.)
!                       if(korb/=ipair) then
!                          factor = merge( 1._prec,-1._prec,ipair<korb)
!                          call MatVec(nbas, 1.5_prec*factor,&
!                               matM_S(:,:,korb,jpair),&
!                               CCpairs(min(ipair,korb),max(ipair,korb))%prev,&
!                               RHS,.false.)
!                       endif
!                       call MatVec(nbas, 1._prec,&
!                            matJ_S(:,:,korb,ipair),&
!                            CCpairs(max(korb,jpair),min(korb,jpair))%prev,&
!                            RHS,.false.)
!                       call MatVec(nbas,-0.5_prec,&
!                            matK_S(:,:,korb,ipair),&
!                            CCpairs(max(korb,jpair),min(korb,jpair))%prev,&
!                            RHS,.false.)
!                       if(korb/=jpair) then
!                          factor = merge( 1._prec,-1._prec,korb<jpair)
!                          call MatVec(nbas,-1.5_prec*factor,&
!                               matM_S(:,:,korb,ipair),&
!                               CCpairs(min(korb,jpair),max(korb,jpair))%prev,&
!                               RHS,.false.)
!                       endif
!                    enddo
! 
!                    if(stage>2) then
! 
!                       do jorb=1,norb
!                          do iorb=1,jorb
!                             factor = pair_energy(iorb,jorb,jpair,ipair)
!                             if(iorb==jorb) factor = factor/2._prec
!                             call MatVec(nbas,-factor,matS_S,&
!                                  CCpairs(jorb,iorb)%prev,RHS,.false.)
!                          enddo
!                       enddo
! 
!                       do jorb=1,norb
!                          do iorb=1,norb
!                             factor = pair_energy(iorb,jorb,&
!                                  max(jpair,jorb),min(jpair,jorb))
!                             call MatVec(nbas,0.5_prec*factor,&
!                                  matS_S,CCpairs(&
!                                  max(ipair,iorb),min(ipair,iorb)&
!                                  )%prev,RHS,.false.)
!                             factor = pair_energy(iorb,jorb,&
!                                  max(ipair,jorb),min(ipair,jorb))
!                             call MatVec(nbas,0.5_prec*factor,&
!                                  matS_S,CCpairs(&
!                                  max(iorb,jpair),min(iorb,jpair)&
!                                  )%prev,RHS,.false.)
!                             if(iorb==jorb) cycle
!                             if(jpair/=jorb) then
!                                factor = pair_energy(iorb,jorb,&
!                                     min(jpair,jorb),max(jpair,jorb))
!                                factor = merge(factor,-factor,jpair<jorb)
!                                call MatVec(nbas,1.5_prec*factor,&
!                                     matS_S,CCpairs(&
!                                     max(ipair,iorb),min(ipair,iorb)&
!                                     )%prev,RHS,.false.)
!                             endif
!                             if(ipair/=jorb) then
!                                factor = pair_energy(iorb,jorb,&
!                                     min(ipair,jorb),max(ipair,jorb))
!                                factor = merge(factor,-factor,ipair<jorb)
!                                call MatVec(nbas,1.5_prec*factor,&
!                                     matS_S,CCpairs(&
!                                     max(iorb,jpair),min(iorb,jpair)&
!                                     )%prev,RHS,.false.)
!                             endif
!                          enddo
!                       enddo
! 
!                    endif
! 
!                 endif
! 
!                 if(stage>1) then
!                    call use_linsolver(pairLHS(jpair,ipair),nbas,RHS,CCpair%tau)
!                 else
!                    call MWO(LHS,jpair,ipair,eta,&
!                         matF_S,matS_S,matP1_S,matPe_S,matP2_S)
!                    if(stage==1.and.iter>1) &
!                         LHS(:,:) = LHS - energy(3,iter-1,ieta)*matS_S
!                    if(iter==1.and.fullPRINT) &
!                         call symmetry_check(jpair,ipair,nbas,LHS)
!                    call simple_linsolver('U','LDL2',nbas,LHS,RHS,CCpair%tau)
!                 endif
! 
!               end associate
! 
!               if(ipair==jpair) cycle
! 
!               associate(CCpair => CCpairs(ipair,jpair))
! 
!                 RHS(:) = CCpair%vecN
! 
!                 if(stage>0.and.iter>1) then
! 
!                    call MatVec(nbas,-1._prec,matL1_T,CCpair%prev,RHS,.false.)
! 
!                    do jorb=2,norb
!                       do iorb=1,jorb-1
!                          factor = &
!                               SCF%pair_energy(iorb,jorb,ijpair) - &
!                               SCF%pair_energy(jorb,iorb,ijpair)
!                          call MatVec(nbas,-factor,matS_T,&
!                               CCpairs(iorb,jorb)%prev,RHS,.false.)
!                       enddo
!                    enddo
! 
!                    do korb=1,norb
!                       if(korb/=ipair) then
!                          factor = merge( 1._prec,-1._prec,ipair<korb)
!                          call MatVec(nbas, 1._prec*factor,&
!                               matJ_T(:,:,korb,jpair),&
!                               CCpairs(min(ipair,korb),max(ipair,korb))%prev,&
!                               RHS,.false.)
!                          call MatVec(nbas,-1.5_prec*factor,&
!                               matK_T(:,:,korb,jpair),&
!                               CCpairs(min(ipair,korb),max(ipair,korb))%prev,&
!                               RHS,.false.)
!                       endif
!                       call MatVec(nbas, 0.5_prec,&
!                            matM_T(:,:,korb,jpair),&
!                            CCpairs(max(ipair,korb),min(ipair,korb))%prev,&
!                            RHS,.false.)
!                       if(korb/=jpair) then
!                          factor = merge( 1._prec,-1._prec,korb<jpair)
!                          call MatVec(nbas, 1._prec*factor,&
!                               matJ_T(:,:,korb,ipair),&
!                               CCpairs(min(korb,jpair),max(korb,jpair))%prev,&
!                               RHS,.false.)
!                          call MatVec(nbas,-1.5_prec*factor,&
!                               matK_T(:,:,korb,ipair),&
!                               CCpairs(min(korb,jpair),max(korb,jpair))%prev,&
!                               RHS,.false.)
!                       endif
!                       call MatVec(nbas,-0.5_prec,&
!                            matM_T(:,:,korb,ipair),&
!                            CCpairs(max(korb,jpair),min(korb,jpair))%prev,&
!                            RHS,.false.)
!                    enddo
! 
!                    if(stage>2) then
! 
!                       do jorb=2,norb
!                          do iorb=1,jorb-1
!                             factor = pair_energy(iorb,jorb,ipair,jpair)
!                             call MatVec(nbas,-factor,matS_T,&
!                                  CCpairs(iorb,jorb)%prev,RHS,.false.)
!                          enddo
!                       enddo
! 
!                       do jorb=1,norb
!                          do iorb=1,norb
!                             if(ipair/=iorb) then
!                                factor = pair_energy(iorb,jorb,&
!                                     max(jpair,jorb),min(jpair,jorb))
!                                factor = merge(factor,-factor,ipair<iorb)
!                                call MatVec(nbas,0.5_prec*factor,&
!                                     matS_T,CCpairs(&
!                                     min(ipair,iorb),max(ipair,iorb)&
!                                     )%prev,RHS,.false.)
!                             endif
!                             if(iorb/=jpair) then
!                                factor = pair_energy(iorb,jorb,&
!                                     max(ipair,jorb),min(ipair,jorb))
!                                factor = merge(factor,-factor,iorb<jpair)
!                                call MatVec(nbas,0.5_prec*factor,&
!                                     matS_T,CCpairs(&
!                                     min(iorb,jpair),max(iorb,jpair)&
!                                     )%prev,RHS,.false.)
!                             endif
!                             if(iorb==jorb) cycle
!                             if(ipair/=iorb.and.jpair/=jorb) then
!                                factor = pair_energy(iorb,jorb,&
!                                     min(jpair,jorb),max(jpair,jorb))
!                                factor = merge(factor,-factor,ipair<iorb)
!                                factor = merge(factor,-factor,jpair<jorb)
!                                call MatVec(nbas,1.5_prec*factor,&
!                                     matS_T,CCpairs(&
!                                     min(ipair,iorb),max(ipair,iorb)&
!                                     )%prev,RHS,.false.)
!                             endif
!                             if(iorb/=jpair.and.ipair/=jorb) then
!                                factor = pair_energy(iorb,jorb,&
!                                     min(ipair,jorb),max(ipair,jorb))
!                                factor = merge(factor,-factor,iorb<jpair)
!                                factor = merge(factor,-factor,ipair<jorb)
!                                call MatVec(nbas,1.5_prec*factor,&
!                                     matS_T,CCpairs(&
!                                     min(iorb,jpair),max(iorb,jpair)&
!                                     )%prev,RHS,.false.)
!                             endif
!                          enddo
!                       enddo
! 
!                    endif
! 
!                 endif
! 
!                 call shrink_Triplet(RHS,Triplet)
!                 if(stage>1) then
!                    call use_linsolver(pairLHS(ipair,jpair),&
!                         Triplet%nbas,RHS,CCpair%tau)
!                 else
!                    call MWO(LHS,ipair,jpair,eta,&
!                         matF_T,matS_T,matP1_T,matPe_T,matP2_T)
!                    if(stage==1.and.iter>1) &
!                         LHS(:,:) = LHS - energy(3,iter-1,ieta)*matS_T
!                    call shrink_Triplet(LHS,Triplet)
!                    if(iter==1.and.fullPRINT) &
!                         call symmetry_check(ipair,jpair,Triplet%nbas,LHS)
!                    call simple_linsolver('U','LDL2',&
!                         Triplet%nbas,LHS,RHS,CCpair%tau)
!                 endif
!                 call expand_Triplet(CCpair%tau,Triplet)
! 
!               end associate
! 
!            enddo
!         enddo
! 
!         call summarize_iteration
! 
!         if(stage>0) then
! 
!            energy_this = energy(3,iter,ieta)
!            if(iter>1) then
!               edelta_this = abs(energy_this - energy_prev)
!               if(edelta_this<Control%ENERGTHR) then
!                  if(fullPRINT) write(LOUT,'(a)') 'CONVERGENCE'
!                  exit
!               endif
!               if(edelta_this/edelta_prev>Control%DIVERTHR) then
!                  if(fullPRINT) write(LOUT,'(a)') 'DIVERGENCE'
!                  exit
!               endif
!            else
!               edelta_this = abs(energy_this)
!            endif
!            energy_prev = energy_this
!            edelta_prev = edelta_this
! 
!            if(iter>=MAXIT) then
!               if(fullPRINT) write(LOUT,'(a)') 'NO CONVERGENCE'
!               exit
!            endif
! 
!            iter = iter + 1
! 
!         else
! 
!            exit
! 
!         endif
! 
!      enddo
!      final_iter(ieta) = iter
! 
!      if(stage>1) then
!         do jpair=1,norb
!            do ipair=1,norb
!               call free_linsolver(pairLHS(ipair,jpair))
!            enddo
!         enddo
!      endif
! 
!      if(LPRINT>=1) then
!         write(LOUT,'()')
!         if(stage>0) then
!            if(do_OS) then
!               write(LOUT,'(a)') 'First iteration'
!            else
!               write(LOUT,'(a)') 'MP2'
!            endif
!         endif
!         do jpair=1,norb
!            do ipair=1,jpair
!               associate(CCpair => CCpairs(jpair,ipair))
!                 write(LOUT,'(2i1,1x,a,3f21.15,2f6.1)') ipair,jpair,'S',&
!                      CCpair%energy(:,1,ieta),&
!                      CCpair%orthog(:,1,ieta)
!               end associate
!               if(ipair==jpair) cycle
!               associate(CCpair => CCpairs(ipair,jpair))
!                 write(LOUT,'(2i1,1x,a,3f21.15,2f6.1)') ipair,jpair,'T',&
!                      CCpair%energy(:,1,ieta),&
!                      CCpair%orthog(:,1,ieta)
!               end associate
!            enddo
!         enddo
!         if(stage>0) then
!            write(LOUT,'()')
!            write(LOUT,'(a)') trim(System%calc_type)
!            do jpair=1,norb
!               do ipair=1,jpair
!                  associate(CCpair => CCpairs(jpair,ipair))
!                    write(LOUT,'(2i1,1x,a,3f21.15,2f6.1)') ipair,jpair,'S',&
!                         CCpair%energy(:,final_iter(ieta),ieta),&
!                         CCpair%orthog(:,final_iter(ieta),ieta)
!                  end associate
!                  if(ipair==jpair) cycle
!                  associate(CCpair => CCpairs(ipair,jpair))
!                    write(LOUT,'(2i1,1x,a,3f21.15,2f6.1)') ipair,jpair,'T',&
!                         CCpair%energy(:,final_iter(ieta),ieta),&
!                         CCpair%orthog(:,final_iter(ieta),ieta)
!                  end associate
!               enddo
!            enddo
!         endif
!      endif
! 
!      if(fullPRINT) then
!         write(LOUT,'()')
!         flush(LOUT)
!      endif
!   enddo
! 
!   if(stage>1) deallocate(pairLHS)
! 
!   call free_linsolver(project_T)
!   call free_linsolver(project_S)
! 
!   call mem_dealloc(RHS)
!   call mem_dealloc(LHS)
! 
!   call free_Triplet(Triplet)
! 
!   if(fullPRINT) then
!      write(LOUT,'()')
!      if(do_OS) then
!         write(LOUT,'(5x,a)') 'FINAL RESULTS : first iteration'
!      else
!         write(LOUT,'(5x,a)') 'FINAL RESULTS : MP2'
!      endif
!      call header_eta
!      do ieta=1,System%Neta
!         power = System%eta(ieta)
!         call print_eta_partial(power)
!         write(LOUT,'(3f21.15,2f6.1)') energy(:,1,ieta),orthog(:,1,ieta)
!      enddo
!      if(stage>0) then
!         write(LOUT,'()')
!         write(LOUT,'(5x,2a)') 'FINAL RESULTS : ',trim(System%calc_type)
!         call header_eta
!         do ieta=1,System%Neta
!            power = System%eta(ieta)
!            call print_eta_partial(power)
!            write(LOUT,'(3f21.15,2f6.1)') &
!                 energy(:,final_iter(ieta),ieta),&
!                 orthog(:,final_iter(ieta),ieta)
!         enddo
!      endif
!      flush(LOUT)
!   endif
! 
!   RESULT_energy = energy(3,final_iter(System%Neta),System%Neta)
! 
!   call mem_dealloc(final_iter)
!   call mem_dealloc(orthog)
!   call mem_dealloc(energy)
! 
!   do jpair=1,norb
!      do ipair=1,norb
!         associate(CCpair => CCpairs(ipair,jpair))
! 
!           call free_CCpair(CCpair)
! 
!         end associate
!      enddo
!   enddo
!   deallocate(CCpairs)
! 
!   if(stage>2) call mem_dealloc(pair_energy)
! 
!   call mem_dealloc(vecP)
!   call mem_dealloc(vec1)
!   call mem_dealloc(vec0)
! 
!   if(stage>0) then
!      call mem_dealloc(matM_T)
!      call mem_dealloc(matM_S)
!      call mem_dealloc(matK_T)
!      call mem_dealloc(matK_S)
!      call mem_dealloc(matJ_T)
!      call mem_dealloc(matJ_S)
!      call mem_dealloc(matL1_S); call mem_dealloc(matL1_T)
!   endif
!   call mem_dealloc(matP2_S); call mem_dealloc(matP2_T)
!   call mem_dealloc(matPe_S); call mem_dealloc(matPe_T)
!   call mem_dealloc(matP1_S); call mem_dealloc(matP1_T)
!   call mem_dealloc(matF_S);  call mem_dealloc(matF_T)
!   call mem_dealloc(matS_S);  call mem_dealloc(matS_T)
! 
! end associate
! 
! contains
! 
!   subroutine summarize_iteration
!   implicit none
! 
!   call collect_properties
! 
!   if(DIIS.and.iter>=DIIS_start-1) then
!      if(.false.) then
!         call common_DIIS
!      else
!         call separate_DIIS
!      endif
!   endif
! 
!   if(fullPRINT) then
!      if(iter==1) call header_iter
!      write(LOUT,'(i4,3f21.15,2f6.1)') &
!           iter,energy(:,iter,ieta),orthog(:,iter,ieta)
!      flush(LOUT)
!   endif
! 
!   end subroutine summarize_iteration
! 
!   subroutine collect_properties
!   implicit none
! 
!   energy(:,iter,ieta) = 0._prec
!   orthog(:,iter,ieta) = -huge(0._prec)
!   do jpair=1,norb
!      do ipair=1,norb
!         associate(CCpair => CCpairs(ipair,jpair))
! 
!           if(ipair>=jpair) then
! 
!              LHS(:,:) = matS_S - matP1_S + matP2_S
!              call MatVec(nbas,1._prec,LHS,CCpair%tau,RHS)
!              call use_linsolver(project_S,nbas,RHS,CCpair%Qtau)
! 
!              associate(energy => CCpair%energy(:,iter,ieta))
!                energy(1) = &
!                     CCpair%multE*inner(CCpair%nbas,CCpair%vecE,CCpair%tau)
!                energy(2) = &
!                     CCpair%multE*inner(CCpair%nbas,CCpair%vecE,CCpair%Qtau)
!                energy(3) = &
!                     CCpair%multN*inner(CCpair%nbas,CCpair%vecN,CCpair%tau)
!              end associate
!              associate(orthog => CCpair%orthog(:,iter,ieta))
!                orthog(1) = log10(abs(&
!                     innerMat(CCpair%nbas,CCpair%tau,CCpair%tau,matP1_S)/&
!                     innerMat(CCpair%nbas,CCpair%tau,CCpair%tau,matS_S)))
!                orthog(2) = log10(abs(&
!                     innerMat(CCpair%nbas,CCpair%Qtau,CCpair%Qtau,matP1_S)/&
!                     innerMat(CCpair%nbas,CCpair%Qtau,CCpair%Qtau,matS_S)))
!              end associate
! 
!           else
! 
!              LHS(:,:) = matS_T - matP1_T + matP2_T
!              call MatVec(nbas,1._prec,LHS,CCpair%tau,RHS)
!              call shrink_Triplet(RHS,Triplet)
!              call use_linsolver(project_T,Triplet%nbas,RHS,CCpair%Qtau)
!              call expand_Triplet(CCpair%Qtau,Triplet)
! 
!              associate(energy => CCpair%energy(:,iter,ieta))
!                energy(1) = &
!                     CCpair%multE*inner(CCpair%nbas,CCpair%vecE,CCpair%tau)
!                energy(2) = &
!                     CCpair%multE*inner(CCpair%nbas,CCpair%vecE,CCpair%Qtau)
!                energy(3) = &
!                     CCpair%multN*inner(CCpair%nbas,CCpair%vecN,CCpair%tau)
!              end associate
!              associate(orthog => CCpair%orthog(:,iter,ieta))
!                orthog(1) = log10(abs(&
!                     innerMat(CCpair%nbas,CCpair%tau,CCpair%tau,matP1_T)/&
!                     innerMat(CCpair%nbas,CCpair%tau,CCpair%tau,matS_T)))
!                orthog(2) = log10(abs(&
!                     innerMat(CCpair%nbas,CCpair%Qtau,CCpair%Qtau,matP1_T)/&
!                     innerMat(CCpair%nbas,CCpair%Qtau,CCpair%Qtau,matS_T)))
!              end associate
! 
!           endif
! 
!           energy(:,iter,ieta) = &
!                energy(:,iter,ieta) + CCpair%energy(:,iter,ieta)
!           orthog(:,iter,ieta) = &
!                max(orthog(:,iter,ieta),CCpair%orthog(:,iter,ieta))
! 
!           if(CCpair%orthog(1,iter,ieta)<CCpair%orthog(2,iter,ieta)) then
!              CCpair%prev(:) = CCpair%tau
!           else
!              CCpair%prev(:) = CCpair%Qtau
!           endif
! 
!         end associate
!      enddo
!   enddo
! 
!   end subroutine collect_properties
! 
!   subroutine common_DIIS
!   implicit none
!   real(prec),allocatable :: mat(:,:),b(:),x(:)
!   integer :: i,j
!   real(prec) :: rtmp
! 
!   DIIS_n = DIIS_n + 1
!   if(DIIS_n>DIIS_size+1) then
!      DIIS_n   = DIIS_size + 1
!      DIIS_off = idx(DIIS_off + 1)
!   endif
! 
!   do jpair=1,norb
!      do ipair=1,norb
!         associate(CCpair => CCpairs(ipair,jpair))
! 
!           CCpair%DIIS_vec(:,idx(DIIS_off + DIIS_n)) = CCpair%prev
! 
!         end associate
!      enddo
!   enddo
! 
!   if(DIIS_n<2) return
! 
!   call mem_alloc(mat,DIIS_n,DIIS_n)
!   call mem_alloc(b,DIIS_n)
!   call mem_alloc(x,DIIS_n)
! 
!   mat = 0._prec
!   do jpair=1,norb
!      do ipair=1,norb
!         associate(CCpair => CCpairs(ipair,jpair))
! 
!           do j=1,DIIS_n-1
!              RHS(:) =&
!                   - CCpair%DIIS_vec(:,idx(DIIS_off + j)) &
!                   + CCpair%DIIS_vec(:,idx(DIIS_off + j + 1))
!              rtmp = inner(nbas,RHS,CCpair%DIIS_vec(:,idx(DIIS_off + 1)))
!              do i=1,j
!                 mat(i,j) = mat(i,j) - rtmp
!                 rtmp = inner(nbas,RHS,CCpair%DIIS_vec(:,idx(DIIS_off + i + 1)))
!                 mat(i,j) = mat(i,j) + rtmp
!              enddo
!           enddo
! 
!         end associate
!      enddo
!   enddo
!   do i=1,DIIS_n-1
!      mat(i,DIIS_n) = -1._prec
!   enddo
! 
!   b = 0._prec
!   b(DIIS_n) = -1._prec
! 
!   call simple_linsolver('U','LDL2',DIIS_n,mat,b,x)
! 
! !!$  do i=1,DIIS_n
! !!$     write(*,'(100a10)',advance='no') ('',j=1,i-1)
! !!$     write(*,'(100es10.2)') (mat(i,j),j=i,DIIS_n)
! !!$  enddo
! !!$  write(*,'(100es10.2)') b
! !!$  write(*,'(100es10.2)') x
! !!$  write(*,*) sum(x(1:DIIS_n-1))
! 
!   do jpair=1,norb
!      do ipair=1,norb
!         associate(CCpair => CCpairs(ipair,jpair))
! 
!           RHS = 0._prec
!           do i=1,DIIS_n-1
!              RHS(:) = RHS + x(i)*CCpair%DIIS_vec(:,idx(DIIS_off + i + 1))
!           enddo
! 
!           CCpair%prev(:) = RHS
!           !CCpair%DIIS_vec(:,idx(DIIS_off + DIIS_n)) = RHS
! 
!         end associate
!      enddo
!   enddo
! 
!   call mem_dealloc(x)
!   call mem_dealloc(b)
!   call mem_dealloc(mat)
! 
!   end subroutine common_DIIS
! 
!   subroutine separate_DIIS
!   implicit none
!   real(prec),allocatable :: mat(:,:),b(:),x(:)
!   integer :: i,j
!   real(prec) :: rtmp
! 
!   DIIS_n = DIIS_n + 1
!   if(DIIS_n>DIIS_size+1) then
!      DIIS_n   = DIIS_size + 1
!      DIIS_off = idx(DIIS_off + 1)
!   endif
! 
!   do jpair=1,norb
!      do ipair=1,norb
!         associate(CCpair => CCpairs(ipair,jpair))
! 
!           CCpair%DIIS_vec(:,idx(DIIS_off + DIIS_n)) = CCpair%prev
! 
!           if(DIIS_n<2) cycle
! 
!           call mem_alloc(mat,DIIS_n,DIIS_n)
!           call mem_alloc(b,DIIS_n)
!           call mem_alloc(x,DIIS_n)
! 
!           mat = 0._prec
!           do j=1,DIIS_n-1
!              RHS(:) =&
!                   - CCpair%DIIS_vec(:,idx(DIIS_off + j)) &
!                   + CCpair%DIIS_vec(:,idx(DIIS_off + j + 1))
!              rtmp = inner(nbas,RHS,CCpair%DIIS_vec(:,idx(DIIS_off + 1)))
!              do i=1,j
!                 mat(i,j) = mat(i,j) - rtmp
!                 rtmp = inner(nbas,RHS,CCpair%DIIS_vec(:,idx(DIIS_off + i + 1)))
!                 mat(i,j) = mat(i,j) + rtmp
!              enddo
!           enddo
!           do i=1,DIIS_n-1
!              mat(i,DIIS_n) = -1._prec
!           enddo
! 
!           b = 0._prec
!           b(DIIS_n) = -1._prec
! 
!           call simple_linsolver('U','LDL2',DIIS_n,mat,b,x)
! 
! !!$          do i=1,DIIS_n
! !!$             write(*,'(100a10)',advance='no') ('',j=1,i-1)
! !!$             write(*,'(100es10.2)') (mat(i,j),j=i,DIIS_n)
! !!$          enddo
! !!$          write(*,'(100es10.2)') b
! !!$          write(*,'(100es10.2)') x
! !!$          write(*,*) sum(x(1:DIIS_n-1))
! 
!           RHS = 0._prec
!           do i=1,DIIS_n-1
!              RHS(:) = RHS + x(i)*CCpair%DIIS_vec(:,idx(DIIS_off + i + 1))
!           enddo
! 
!           CCpair%prev(:) = RHS
!           !CCpair%DIIS_vec(:,idx(DIIS_off + DIIS_n)) = RHS
! 
!           call mem_dealloc(x)
!           call mem_dealloc(b)
!           call mem_dealloc(mat)
! 
!         end associate
!      enddo
!   enddo
! 
!   end subroutine separate_DIIS
! 
!   function idx(i)
!   implicit none
!   integer :: idx
!   integer,intent(in) :: i
!   idx = mod(i-1,DIIS_size+1) + 1
!   if(idx==0) idx = DIIS_size + 1
!   end function idx
! 
! end subroutine CC_energy_common
! 
! subroutine header_iter
! implicit none
! write(LOUT,'()')
! write(LOUT,'(a,3(6x,a,11x),2(2x,a))') &
!      'iter','E(0)','E(1)','E(q)','S(1)','S(2)'
! end subroutine header_iter
! 
! subroutine header_eta
! implicit none
! write(LOUT,'()')
! write(LOUT,'(1x,a,3(6x,a,11x),2(2x,a))') &
!      'eta','E(0)','E(1)','E(q)','S(1)','S(2)'
! end subroutine header_eta
! 
! subroutine print_eta_full(power)
! implicit none
! integer :: power
! write(LOUT,'()')
! if(power>power_LIMIT) then
!    write(LOUT,'(5x,a,i5)') 'log10(ETA) = ',power
! else
!    write(LOUT,'(5x,a)') 'ETA = 0'
! endif
! write(LOUT,'()')
! end subroutine print_eta_full
! 
! subroutine print_eta_partial(power)
! implicit none
! integer :: power
! if(power>power_LIMIT) then
!    write(LOUT,'(i4)',advance='no') power
! else
!    write(LOUT,'(a)',advance='no') '-inf'
! endif
! end subroutine print_eta_partial
! 
! subroutine init_CCpair(CCpair,nbas,MAXIT,Neta,DIIS_size)
! implicit none
! type(CCpairData) :: CCpair
! integer :: nbas,MAXIT,Neta,DIIS_size
! 
! CCpair%nbas = nbas
! 
! call mem_alloc(CCpair%vecE,CCpair%nbas)
! call mem_alloc(CCpair%vecN,CCpair%nbas)
! call mem_alloc(CCpair%tau,CCpair%nbas)
! call mem_alloc(CCpair%Qtau,CCpair%nbas)
! call mem_alloc(CCpair%prev,CCpair%nbas)
! 
! call mem_alloc(CCpair%energy,3,MAXIT,Neta)
! call mem_alloc(CCpair%orthog,2,MAXIT,Neta)
! 
! CCpair%DIIS_size = DIIS_size
! if(CCpair%DIIS_size>0) &
!      call mem_alloc(CCpair%DIIS_vec,CCpair%nbas,CCpair%DIIS_size+1)
! 
! end subroutine init_CCpair
! 
! subroutine free_CCpair(CCpair)
! implicit none
! type(CCpairData) :: CCpair
! 
! if(CCpair%DIIS_size>0) call mem_dealloc(CCpair%DIIS_vec)
! 
! call mem_dealloc(CCpair%orthog)
! call mem_dealloc(CCpair%energy)
! 
! call mem_dealloc(CCpair%prev)
! call mem_dealloc(CCpair%Qtau)
! call mem_dealloc(CCpair%tau)
! call mem_dealloc(CCpair%vecN)
! call mem_dealloc(CCpair%vecE)
! 
! end subroutine free_CCpair
! 
! subroutine MWO(LHS,ipair,jpair,eta,matF,matS,matP1,matPe,matP2)
! implicit none
! real(prec) :: LHS(:,:)
! integer :: ipair,jpair
! real(prec) :: eta
! real(prec) :: matF(:,:),matS(:,:),matP1(:,:),matPe(:,:),matP2(:,:)
! real(prec) :: ei,ej,Delta1,Delta2,Delta3
! 
! ei = SCF%orb_energy(min(ipair,jpair))
! ej = SCF%orb_energy(max(ipair,jpair))
! 
! Delta1 = ei + ej - SCF%HOMO + eta
! Delta2 = -1._prec
! Delta3 = 2*SCF%HOMO - ei - ej
! 
! LHS = matF - (ei + ej)*matS + Delta1*matP1 + Delta2*matPe + Delta3*matP2
! 
! end subroutine MWO
! 
! subroutine symmetry_check(ipair,jpair,n,Mat)
! implicit none
! integer,intent(in) :: ipair,jpair
! integer,intent(in) :: n
! real(prec),intent(in) :: Mat(:,:)
! real(prec) :: test,val1,val2
! integer :: i,j
! 
! test = 0._prec
! 
! do j=2,n
!    do i=1,j-1
!       val1 = Mat(i,j)
!       val2 = Mat(j,i)
!       test = max(test,abs((val1 - val2)/(0.5_prec*(val1 + val2))))
!    enddo
! enddo
! 
! if(ipair>0.and.jpair>0) then
!    write(LOUT,'(a,2i1,a,a,es13.3)') 'Symmetry test for pair ',&
!         min(ipair,jpair),max(ipair,jpair),merge('S','T',ipair>=jpair),' :',test
! else
!    write(LOUT,'(a,es13.3)') 'Symmetry test : ',test
! endif
! 
! end subroutine symmetry_check
! 
! function inner(n,Vec1,Vec2) result(dot)
! implicit none
! real(prec) :: dot
! integer,intent(in) :: n
! real(prec),intent(in) :: Vec1(:),Vec2(:)
! integer :: i
! 
! dot = 0._prec
! 
! do i=1,n
!    dot = dot + Vec1(i)*Vec2(i)
! enddo
! 
! end function inner
! 
! function innerMat(n,Vec1,Vec2,Mat) result(dot)
! implicit none
! real(prec) :: dot
! integer,intent(in) :: n
! real(prec),intent(in) :: Vec1(:),Vec2(:),Mat(:,:)
! integer :: i,j
! 
! dot = 0._prec
! 
! do j=1,n
!    do i=1,n
!       dot = dot + Vec1(i)*Mat(i,j)*Vec2(j)
!    enddo
! enddo
! 
! end function innerMat
! 
! subroutine outerMat(n,alpha,Vec1,Vec2,Mat,do_cleaning)
! implicit none
! integer,intent(in) :: n
! real(prec),intent(in) :: alpha
! real(prec),intent(in) :: Vec1(:),Vec2(:)
! real(prec),intent(inout) :: Mat(:,:)
! logical,intent(in),optional :: do_cleaning
! logical :: cleaning
! integer :: i,j
! real(prec) :: rtmp
! 
! if(present(do_cleaning)) then
!    cleaning = do_cleaning
! else
!    cleaning = .true.
! endif
! 
! if(cleaning) Mat = 0._prec
! 
! do j=1,n
!    rtmp = alpha*Vec2(j)
!    do i=1,n
!       Mat(i,j) = Mat(i,j) + Vec1(i)*rtmp
!    enddo
! enddo
! 
! end subroutine outerMat
! 
! subroutine MatVec(n,alpha,Mat,Vec,Res,do_cleaning)
! implicit none
! integer,intent(in) :: n
! real(prec),intent(in) :: alpha
! real(prec),intent(in) :: Mat(:,:),Vec(:)
! real(prec),intent(inout) :: Res(:)
! logical,intent(in),optional :: do_cleaning
! logical :: cleaning
! integer :: i,j
! real(prec) :: rtmp
! 
! if(present(do_cleaning)) then
!    cleaning = do_cleaning
! else
!    cleaning = .true.
! endif
! 
! if(cleaning) Res = 0._prec
! 
! do j=1,n
!    rtmp = alpha*Vec(j)
!    do i=1,n
!       Res(i) = Res(i) + Mat(i,j)*rtmp
!    enddo
! enddo
! 
! end subroutine MatVec
! 
end module CCdriver
