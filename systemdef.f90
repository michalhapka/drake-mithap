module systemdef
use file_OUT, only : LOUT
use precision, only : prec
use memory
use commontypes
implicit none

private
public create_System

contains

subroutine create_System(System,Input,Control)
implicit none
type(SystemData) :: System
type(InputData),intent(in) :: Input
type(ControlData),intent(in) :: Control
logical :: post_SCF,optimize
logical :: common_orbs,common_pairs
integer :: n_orbs,n_pairs
integer :: maxprim_orbs,maxprim_pairs
integer :: sumprim_orbs,sumprim_pairs
integer :: nprim_orbs,nprim_pairs
integer :: Neta
integer :: Nexp_orig,Nexp
logical,allocatable :: isUsed_orbs(:),isUsed_pairs(:),isUsed(:)
logical,allocatable :: isOpt_orbs(:),isOpt_pairs(:)
integer,allocatable :: exp_reorder(:)
integer,allocatable :: prim(:,:),prim_all(:,:)
integer,allocatable :: diff(:),diff_all(:)
logical,allocatable :: expSQ_restrict(:)
character(2) :: sepspace
logical :: written
integer :: i,j,k

! post_SCF = (index(Input%calc_type,'SCF')==0)
! optimize = (index(Input%calc_type,'OPT')/=0)
! 
! select case(trim(Input%orbs_basis))
! case('COMMON')
!    common_orbs = .true.
!    n_orbs = 1
! case('SEPARATE')
!    common_orbs = .false.
!    n_orbs = Input%norb
!    if(n_orbs/=2) then
!       write(LOUT,'(a)') &
!            'ERROR!!! Separate SCF orbitals possible only for norb = 2 !'
!       stop
!    endif
! case default
!    write(LOUT,'(a)') 'Unrecognizable type of orbital basis set: &
!         &creating system data!'
!    stop
! end select
! 
! if(post_SCF) then
!    select case(trim(Input%pairs_basis))
!    case('COMMON')
!       common_pairs = .true.
!       n_pairs = 1
!    case('SEPARATE')
!       common_pairs = .false.
!       n_pairs = Input%norb
!    case default
!       write(LOUT,'(a)') 'Unrecognizable type of pairs basis set: &
!            &creating system data!'
!       stop
!    end select
! else
!    common_pairs = .true.
!    n_pairs = 0
! endif
! 
! if(post_SCF) then
!    if(optimize) then
!       if(Input%Neta>1) then
!          write(LOUT,'(a)') 'WARNING!!! &
!               &In post SCF optimizations &
!               &only the smallest eta value will be used!'
!       endif
!       Neta = 1
!    else
!       Neta = Input%Neta
!    endif
! else
!    Neta = 0
! endif
! 
! Nexp_orig = Input%Nexponents
! 
! call mem_alloc(exp_reorder,Nexp_orig)
! 
! call mem_alloc(isUsed_orbs ,Nexp_orig)
! call mem_alloc(isUsed_pairs,Nexp_orig)
! call mem_alloc(isUsed      ,Nexp_orig)
! call mem_alloc(isOpt_orbs  ,Nexp_orig)
! call mem_alloc(isOpt_pairs ,Nexp_orig)
! 
! maxprim_orbs = 0
! sumprim_orbs = 0
! isUsed_orbs  = .false.
! do j=1,n_orbs
!    associate(OrbInput => Input%OrbInput(j))
! 
!      nprim_orbs = OrbInput%Nprimitives
!      if(nprim_orbs<1) then
!         write(LOUT,'(a)') 'Every orbital should be specified!'
!         stop
!      endif
! 
!      maxprim_orbs = max(maxprim_orbs,nprim_orbs)
!      sumprim_orbs = sumprim_orbs + nprim_orbs
! 
!      do i=1,nprim_orbs
!         isUsed_orbs(OrbInput%orb_control(1,i)) = .true.
!      enddo
! 
!    end associate
! enddo
! 
! maxprim_pairs = 0
! sumprim_pairs = 0
! isUsed_pairs  = .false.
! if(post_SCF) then
!    do k=1,n_pairs
!       do j=1,n_pairs
!          associate(PairInput => Input%PairInput(j,k))
! 
!            if(PairInput%mult<0) cycle
! 
!            nprim_pairs = PairInput%Nprimitives
! 
!            maxprim_pairs = max(maxprim_pairs,nprim_pairs)
!            sumprim_pairs = sumprim_pairs + nprim_pairs
! 
!            do i=1,nprim_pairs
!               isUsed_pairs(PairInput%pair_control(1,i)) = .true.
!               isUsed_pairs(PairInput%pair_control(2,i)) = .true.
!            enddo
! 
!          end associate
!       enddo
!    enddo
!    if(maxprim_pairs<1) then
!       write(LOUT,'(a)') 'At least one pair should be specified!'
!       stop
!    endif
! endif
! 
! isOpt_orbs  = .false.
! isOpt_pairs = .false.
! if(optimize) then
!    if(post_SCF) then
! 
!       isOpt_pairs(:) = (isUsed_pairs.and.Input%optimized)
! 
!       if(.not.any(isOpt_pairs)) then
!          write(LOUT,'(a)') 'No pair exponents to optimize!'
!          stop
!       endif
! 
!       isUsed(:) = (isUsed_orbs.and.isOpt_pairs)
!       if(any(isUsed)) then
!          write(LOUT,'()')
!          write(LOUT,'(a)',advance='no') 'Some orbital and &
!               &optimized pair exponents overlap'
!          sepspace = ': '
!          do i=1,Nexp_orig
!             if(isUsed(i)) then
!                write(LOUT,'(a,i3)',advance='no') sepspace,i
!                sepspace = ', '
!             endif
!          enddo
!          write(LOUT,'()')
!          stop
!       endif
! 
!    else
! 
!       isOpt_orbs(:) = (isUsed_orbs.and.Input%optimized)
! 
!       if(.not.any(isOpt_orbs)) then
!          write(LOUT,'(a)') 'No orbital exponents to optimize!'
!          stop
!       endif
! 
!    endif
! endif
! 
! if(post_SCF) then
!    isUsed(:) = (isUsed_pairs.and..not.isOpt_pairs)
!    isUsed(:) = (isUsed_orbs.or.isUsed)
! else
!    isUsed(:) = (isUsed_orbs.and..not.isOpt_orbs)
! endif
! written = .false.
! do j=2,Nexp_orig
!    do i=1,j-1
!       if(isUsed(i).and.isUsed(j)) then
!          if(abs(Input%exponents(i) - Input%exponents(j))<Control%SIMTHR) then
!             if(.not.written) then
!                write(LOUT,'(a,es10.3,a)') 'ERROR!!! &
!                     &Nonopt exponents are too close to each other &
!                     &(SIMTHR = ',Control%SIMTHR,')'
!                written = .true.
!             endif
!             write(LOUT,'(2i3,5x,2f16.8)') &
!                  i,j,Input%exponents(i),Input%exponents(j)
!          endif
!       endif
!    enddo
! enddo
! if(written) stop
! 
! isUsed(:) = (isUsed_orbs.or.isUsed_pairs)
! Nexp = count(isUsed)
! 
! if(Nexp/=Nexp_orig) then
!    write(LOUT,'()')
!    write(LOUT,'(a)',advance='no') &
!         'Input exponents that are not used will be ignored'
!    sepspace = ': '
!    do i=1,Nexp_orig
!       if(.not.isUsed(i)) then
!          write(LOUT,'(a,i3)',advance='no') sepspace,i
!          sepspace = ', '
!       endif
!    enddo
!    write(LOUT,'()')
! endif
! 
! System%nucZ         = Input%nucZ
! System%maxl         = Input%maxl
! System%norb         = Input%norb
! System%calc_type    = Input%calc_type
! System%post_SCF     = post_SCF
! System%optimize     = optimize
! System%common_orbs  = common_orbs
! System%common_pairs = common_pairs
! System%n_orbs       = n_orbs
! System%n_pairs      = n_pairs
! System%Neta         = Neta
! System%Nexp         = Nexp
! System%NexpSQ       = System%Nexp*(System%Nexp + 1)/2
! 
! call init_System(System)
! 
! System%eta(:) = Input%eta(Input%Neta-System%Neta+1:Input%Neta)
! 
! i=0
! do j=1,Nexp_orig
!    if(isUsed(j)) then
!       i = i + 1
! 
!       System%exponents(i)    = Input%exponents(j)
!       System%isUsed_orbs(i)  = isUsed_orbs(j)
!       System%isUsed_pairs(i) = isUsed_pairs(j)
!       System%isOpt_orbs(i)   = isOpt_orbs(j)
!       System%isOpt_pairs(i)  = isOpt_pairs(i)
! 
!       exp_reorder(j) = i
!    else
!       exp_reorder(j) = -1
!    endif
! enddo
! if(i/=System%Nexp) then
!    write(LOUT,'(a)') 'ERROR!!! Wrong number of exponents in create_System!'
!    stop
! endif
! 
! call mem_dealloc(isOpt_pairs)
! call mem_dealloc(isOpt_orbs)
! call mem_dealloc(isUsed)
! call mem_dealloc(isUsed_pairs)
! call mem_dealloc(isUsed_orbs)
! 
! call mem_alloc(prim,2,maxprim_orbs)
! call mem_alloc(prim_all,2,sumprim_orbs)
! call mem_alloc(diff,maxprim_orbs)
! call mem_alloc(diff_all,sumprim_orbs)
! 
! sumprim_orbs = 0
! do j=1,System%n_orbs
!    associate(OrbInput => Input%OrbInput(j), &
!         OrbSystem => System%OrbSystem(j))
! 
!      nprim_orbs = OrbInput%Nprimitives
!      prim(:,1:nprim_orbs) = OrbInput%orb_control
!      do i=1,nprim_orbs
!         prim(1,i) = exp_reorder(prim(1,i))
!      enddo
! 
!      call sort_prim(nprim_orbs,prim,diff,compare_prim_orbs)
!      if(any(diff(1:nprim_orbs)>1)) then
!         write(LOUT,'(a)') &
!              'ERROR!!! Some of the orbital primitives are redefined!'
!         call print_orb_control(nprim_orbs,prim)
!         stop
!      endif
!      if(any(diff(1:nprim_orbs)==0)) then
!         write(LOUT,'(a)') &
!              'ERROR!!! Some of the orbital primitives appear more than once!'
!         call print_orb_control(nprim_orbs,prim)
!         stop
!      endif
! 
!      call create_OrbSystem(OrbSystem,nprim_orbs,prim)
! 
!      prim_all(:,sumprim_orbs+1:sumprim_orbs+nprim_orbs) = &
!           prim(:,1:nprim_orbs)
!      sumprim_orbs = sumprim_orbs + nprim_orbs
! 
!    end associate
! enddo
! 
! call sort_prim(sumprim_orbs,prim_all,diff_all,compare_prim_orbs)
! call reduce_prim(sumprim_orbs,prim_all,diff_all)
! call create_OrbReduced(System%OrbReduced,sumprim_orbs,prim_all,&
!      System%Nexp,System%isUsed_orbs)
! 
! call mem_dealloc(diff_all)
! call mem_dealloc(diff)
! call mem_dealloc(prim_all)
! call mem_dealloc(prim)
! 
! if(System%post_SCF) then
! 
!    sumprim_pairs = sumprim_pairs + sumprim_orbs*(sumprim_orbs + 1)/2
! 
!    call mem_alloc(prim,6,maxprim_pairs)
!    call mem_alloc(prim_all,6,sumprim_pairs)
!    call mem_alloc(diff,maxprim_pairs)
!    call mem_alloc(diff_all,sumprim_pairs)
! 
!    call mem_alloc(expSQ_restrict,System%NexpSQ)
!    i = 0
!    do k=1,System%Nexp
!       do j=1,k
!          i = i + 1
!          expSQ_restrict(i) = &
!               (abs(System%exponents(j) - System%exponents(k))<Control%REDTHR)
!       enddo
!    enddo
!    if(i/=System%NexpSQ) then
!       write(LOUT,'(a)') 'ERROR!!! &
!            &Wrong number of exponents(SQ) in create_System!'
!       stop
!    endif
! 
!    sumprim_pairs = 0
!    do k=1,System%n_pairs
!       do j=1,System%n_pairs
!          associate(PairInput => Input%PairInput(j,k), &
!               PairSystem => System%PairSystem(j,k))
! 
!            PairSystem%mult = PairInput%mult
! 
!            if(PairSystem%mult<0) cycle
! 
!            write(LOUT,'(a,2(i3,a),a)') 'Processing pair: ',&
!                 min(j,k),', ',max(j,k),', ',possible_mult(PairSystem%mult)
! 
!            nprim_pairs = PairInput%Nprimitives
!            prim(:,1:nprim_pairs) = PairInput%pair_control
!            do i=1,nprim_pairs
!               prim(1,i) = exp_reorder(prim(1,i))
!               prim(2,i) = exp_reorder(prim(2,i))
!            enddo
! 
!            call sort_prim(nprim_pairs,prim,diff,compare_prim_pairs)
!            if(any(diff(1:nprim_pairs)>2)) then
!               write(LOUT,'(a)') &
!                    'ERROR!!! Some of the pair primitives are redefined!'
!               call print_pair_control(nprim_pairs,prim)
!               stop
!            endif
!            if(any(diff(1:nprim_pairs)==0)) then
!               write(LOUT,'(a)') &
!                    'ERROR!!! Some of the pair primitives appear more than once!'
!               call print_pair_control(nprim_pairs,prim)
!               stop
!            endif
! 
!            call create_PairSystem(PairSystem,nprim_pairs,prim,&
!                 diff,expSQ_restrict)
!            call check_PairSystem(PairSystem)
! 
!            prim_all(:,sumprim_pairs+1:sumprim_pairs+nprim_pairs) = &
!                 prim(:,1:nprim_pairs)
!            sumprim_pairs = sumprim_pairs + nprim_pairs
! 
!          end associate
!       enddo
!    enddo
! 
!    call sort_prim(sumprim_pairs,prim_all,diff_all,compare_prim_pairs)
!    call reduce_prim(sumprim_pairs,prim_all,diff_all)
!    call create_PairReduced(System%PairReduced,&
!         sumprim_pairs,prim_all,diff_all,System%Nexp,System%isUsed_pairs)
! 
!    call add_OrbReduced(sumprim_pairs,prim_all,System%Nexp,System%OrbReduced)
!    call sort_prim(sumprim_pairs,prim_all,diff_all,compare_prim_pairs)
!    call reduce_prim(sumprim_pairs,prim_all,diff_all)
!    call create_PairReduced(System%OrbPairReduced,&
!         sumprim_pairs,prim_all,diff_all,0,isUsed)
! 
!    call check_consistency_PairReduced(System%NexpSQ,&
!         System%OrbPairReduced,System%PairReduced)
! 
!    call mem_dealloc(expSQ_restrict)
! 
!    call mem_dealloc(diff_all)
!    call mem_dealloc(diff)
!    call mem_dealloc(prim_all)
!    call mem_dealloc(prim)
! 
! endif
! 
! call mem_dealloc(exp_reorder)
! 
! call print_System(System,Control%LPRINT)
! 
! call print_EndSection
! 
end subroutine create_System

subroutine sort_prim(n,prim,diff,compare_prim)
use misc, only : swap
implicit none
integer :: n
integer :: prim(:,:),diff(:)
interface
function compare_prim(prim1,prim2) result(diff_type)
integer :: diff_type
integer,intent(in) :: prim1(:),prim2(:)
end function compare_prim
end interface
integer :: i,j

do j=2,n
   i = j
   do while(i>1)
      if(compare_prim(prim(:,i-1),prim(:,i))<0) then
         call swap(prim(:,i-1),prim(:,i))
         i = i - 1
      else
         exit
      endif
   enddo
enddo

diff = 1
do i=2,n
   diff(i) = compare_prim(prim(:,i-1),prim(:,i))
enddo
if(any(diff<0)) then
   write(LOUT,'(a)') &
        'ERROR!!! Primitives have not been sorted properly in sort_prim!'
   stop
endif

end subroutine sort_prim

subroutine reduce_prim(n,prim,diff)
implicit none
integer :: n
integer :: prim(:,:),diff(:)
integer :: diff_type
integer :: i,j

j = 2
do while(j<=n)
   diff_type = diff(j)
   if(diff_type==0.or.diff_type>10) then

      if(diff_type>10) then
         do i=diff_type/10,mod(diff_type,10)
            prim(i,j-1) = max(prim(i,j-1),prim(i,j))
         enddo
      endif

      do i=j,n-1
         prim(:,i) = prim(:,i+1)
         diff(i)   = diff(i+1)
      enddo

      n = n - 1

   else

      j = j + 1

   endif
enddo

if(any(diff(1:n)==0).or.any(diff(1:n)>10)) then
   write(LOUT,'(a)') &
        'ERROR!!! Primitives have not been reduced properly in reduce_prim!'
   stop
endif

end subroutine reduce_prim

subroutine create_OrbSystem(OrbSystem,n,prim)
implicit none
type(OrbSystemData) :: OrbSystem
integer,intent(in) :: n
integer,intent(in) :: prim(:,:)
integer :: offset
integer :: i

OrbSystem%n_prim = n

allocate(OrbSystem%OrbSpec(OrbSystem%n_prim))

offset = 0
do i=1,OrbSystem%n_prim
   associate(OrbSpec => OrbSystem%OrbSpec(i))

     OrbSpec%iexp   = prim(1,i)
     OrbSpec%irange = prim(2,i)
     OrbSpec%nbas   = OrbSpec%irange + 1
     OrbSpec%offset = offset

     offset = offset + OrbSpec%nbas

   end associate
enddo

OrbSystem%nbas = offset

end subroutine create_OrbSystem

subroutine create_OrbReduced(OrbReduced,n,prim,Nexp,isUsed)
implicit none
type(OrbReducedData) :: OrbReduced(:)
integer,intent(in) :: n
integer,intent(in) :: prim(:,:)
integer,intent(in) :: Nexp
logical,intent(in) :: isUsed(:)
integer :: i,iexp

do i=1,n
   iexp = prim(1,i)
   associate(OrbSpec => OrbReduced(iexp))

     OrbSpec%isUsed   = .true.
     OrbSpec%iexp     = iexp
     OrbSpec%maxrange = prim(2,i)

   end associate
enddo

do i=1,Nexp
   if(OrbReduced(i)%isUsed.neqv.isUsed(i)) then
      write(LOUT,'(a)') &
           'ERROR!!! Used orbital exponents have not been predicted properly!'
      stop
   endif
enddo

end subroutine create_OrbReduced

subroutine add_OrbReduced(n,prim,Nexp,OrbReduced)
implicit none
integer :: n
integer :: prim(:,:)
integer,intent(in) :: Nexp
type(OrbReducedData),intent(in) :: OrbReduced(:)
integer :: i,j

do j=1,Nexp
   do i=1,j
      associate(OrbReduced1 => OrbReduced(i), OrbReduced2 => OrbReduced(j))
        if(OrbReduced1%isUsed.and.OrbReduced2%isUsed) then
           n = n + 1
           prim(1,n) = OrbReduced1%iexp
           prim(2,n) = OrbReduced2%iexp
           prim(3,n) = 3
           prim(4,n) = OrbReduced1%maxrange
           prim(5,n) = OrbReduced2%maxrange
           prim(6,n) = 0
        endif
      end associate
   enddo
enddo

end subroutine add_OrbReduced

subroutine create_PairSystem(PairSystem,n,prim,diff,expSQ_restrict)
implicit none
type(PairSystemData) :: PairSystem
integer,intent(in) :: n
integer,intent(in) :: prim(:,:)
integer,intent(in) :: diff(:)
logical,intent(in) :: expSQ_restrict(:)
integer :: icnt,istart,iend
integer :: iexp1,iexp2,iexpSQ,itype
integer :: irange,maxrange
integer :: offset,nbas
integer :: i

PairSystem%n_prim = count(diff(1:n)==1)

allocate(PairSystem%PairSpec(PairSystem%n_prim))

offset = 0
icnt = 0
istart = 1
do while(istart<=n)
   icnt = icnt + 1
   iend = find_end(istart,n,diff,10)

   iexp1 = prim(1,istart)
   iexp2 = prim(2,istart)
   if(iexp1>iexp2) then
      write(LOUT,'(a)') 'ERROR!!! &
           &Incorrect exponents order in create_PairSystem!'
      stop
   endif
   iexpSQ = iexp1 + (iexp2 - 1)*iexp2/2

   itype = prim(3,istart)

   do i=istart+1,iend
      if(prim(1,i)/=iexp1.or.prim(2,i)/=iexp2.or.prim(3,i)/=itype) then
         write(LOUT,'(a)') 'ERROR!!! &
              &Incorrect recognition of pair primitive in create_PairSystem!'
         stop
      endif
   enddo

   associate(PairSpec => PairSystem%PairSpec(icnt))

     PairSpec%iexp1    = iexp1
     PairSpec%iexp2    = iexp2
     PairSpec%iexpSQ   = iexpSQ
     PairSpec%itype    = itype
     PairSpec%restrict = expSQ_restrict(PairSpec%iexpSQ)
     PairSpec%n_range  = iend - istart + 1
     PairSpec%irange1  = -1

     nbas = 0

     select case(PairSpec%itype)
     case(1)

        if(PairSpec%n_range>1) then
           write(LOUT,'(a)') 'ERROR!!! &
                &More than one range specification for itype = 1 !'
           stop
        endif

        PairSpec%irange1 = prim(4,istart)

        irange = PairSpec%irange1
        if(PairSpec%restrict) then
           nbas = ((irange + 2)*(irange + 4)*(2*irange + 3))/24
        else
           nbas = ((irange + 1)*(irange + 2)*(irange + 3))/6
        endif

     case(2)

        maxrange = 0
        do i=istart,iend
           maxrange = max(maxrange,prim(5,i))
        enddo
        if(PairSpec%n_range/=maxrange+1) then
           write(LOUT,'(a)') 'WARNING!!! &
                &Some r12 powers may be missing!'
           call print_pair_control(PairSpec%n_range,prim(:,istart:iend))
        endif

        call mem_alloc(PairSpec%irange23,2,PairSpec%n_range)
        do i=1,PairSpec%n_range
           PairSpec%irange23(:,i) = prim(4:5,istart+i-1)
        enddo

        do i=1,PairSpec%n_range
           irange = PairSpec%irange23(1,i)
           if(PairSpec%restrict) then
              nbas = nbas + (irange + 2)**2/4
           else
              nbas = nbas + ((irange + 1)*(irange + 2))/2
           endif
        enddo

     case(3)

        maxrange = 0
        do i=istart,iend
           maxrange = max(maxrange,prim(6,i))
        enddo
        if(PairSpec%n_range/=maxrange+1) then
           write(LOUT,'(a)') 'WARNING!!! &
                &Some r12 powers may be missing!'
           call print_pair_control(PairSpec%n_range,prim(:,istart:iend))
        endif
        if(PairSpec%restrict) then
           do i=istart,iend
              if(prim(4,i)/=prim(5,i)) then
                 write(LOUT,'(a)') 'ERROR!!! &
                      &In restricted range type 3, &
                      &r1 and r2 powers must be equal!'
                 call print_pair_control(PairSpec%n_range,prim(:,istart:iend))
                 stop
              endif
           enddo
        endif

        call mem_alloc(PairSpec%irange23,3,PairSpec%n_range)
        do i=1,PairSpec%n_range
           PairSpec%irange23(:,i) = prim(4:6,istart+i-1)
        enddo

        do i=1,PairSpec%n_range
           if(PairSpec%restrict) then
              irange = PairSpec%irange23(2,i)
              nbas = nbas + ((irange + 1)*(irange + 2))/2
           else
              nbas = nbas + &
                   (PairSpec%irange23(1,i) + 1)*(PairSpec%irange23(2,i) + 1)
           endif
        enddo

     case default

        write(LOUT,'(a)') &
             'ERROR!!! Incorrect range type specification in create_PairSystem!'
        stop

     end select

     PairSpec%nbas   = nbas
     PairSpec%offset = offset

     offset = offset + PairSpec%nbas

   end associate

   istart = iend + 1
enddo
if(icnt/=PairSystem%n_prim) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Number of prim pairs was incorrectly predicted in create_PairSystem!'
   stop
endif

PairSystem%nbas = offset

end subroutine create_PairSystem

subroutine create_PairReduced(PairReduced,n,prim,diff,Nexp,isUsed)
implicit none
type(PairReducedData) :: PairReduced(:)
integer,intent(in) :: n
integer,intent(in) :: prim(:,:),diff(:)
integer,intent(in) :: Nexp
logical,intent(in) :: isUsed(:)
integer :: icnt,istart,iend
integer :: jcnt,jstart,jend
integer :: iexp1,iexp2,iexpSQ,itype
integer :: n_range
integer :: i,j,k

icnt = 0
istart = 1
do while(istart<=n)
   icnt = icnt + 1
   iend = find_end(istart,n,diff,10)

   iexp1 = prim(1,istart)
   iexp2 = prim(2,istart)
   if(iexp1>iexp2) then
      write(LOUT,'(a)') 'ERROR!!! &
           &Incorrect exponents order in create_PairReduced!'
      stop
   endif
   iexpSQ = iexp1 + (iexp2 - 1)*iexp2/2

   do i=istart+1,iend
      if(prim(1,i)/=iexp1.or.prim(2,i)/=iexp2) then
         write(LOUT,'(a)') 'ERROR!!! &
              &Incorrect recognition of pair primitive in create_PairReduced!'
         stop
      endif
   enddo

   associate(PairSpec => PairReduced(iexpSQ))

     PairSpec%isUsed = .true.
     PairSpec%iexp1  = iexp1
     PairSpec%iexp2  = iexp2
     PairSpec%iexpSQ = iexpSQ

     jcnt = 0
     jstart = istart
     do while(jstart<=iend)
        jcnt = jcnt + 1
        jend = find_end(jstart,iend,diff,5)

        itype = prim(3,jstart)

        do j=jstart+1,jend
           if(prim(3,j)/=itype) then
              write(LOUT,'(a)') 'ERROR!!! &
                   &Incorrect recognition of type in create_PairReduced!'
              stop
           endif
        enddo

        n_range = jend - jstart + 1

        select case(itype)
        case(1)

           if(PairSpec%istype1) then
              write(LOUT,'(a)') 'ERROR!!! &
                   &Second occurence of itype = 1 in create_PairReduced!'
              stop
           endif
           if(n_range>1) then
              write(LOUT,'(a)') 'ERROR!!! &
                   &More than one range specification for itype = 1 !'
              stop
           endif

           PairSpec%istype1   = .true.
           PairSpec%maxrange1 = prim(4,jstart)

        case(2)

           if(PairSpec%istype2) then
              write(LOUT,'(a)') 'ERROR!!! &
                   &Second occurence of itype = 2 in create_PairReduced!'
              stop
           endif

           PairSpec%istype2  = .true.
           PairSpec%n_range2 = n_range

           call mem_alloc(PairSpec%maxrange2,2,PairSpec%n_range2)
           do j=1,PairSpec%n_range2
              PairSpec%maxrange2(:,j) = prim(4:5,jstart+j-1)
           enddo

        case(3)

           if(PairSpec%istype3) then
              write(LOUT,'(a)') 'ERROR!!! &
                   &Second occurence of itype = 3 in create_PairReduced!'
              stop
           endif

           PairSpec%istype3  = .true.
           PairSpec%n_range3 = n_range

           call mem_alloc(PairSpec%maxrange3,3,PairSpec%n_range3)
           do j=1,PairSpec%n_range3
              PairSpec%maxrange3(:,j) = prim(4:6,jstart+j-1)
           enddo

        case default

           write(LOUT,'(a)') 'ERROR!!! &
                &Incorrect range type specification in create_PairReduced!'
           stop

        end select

        jstart = jend + 1
     enddo
     if(jcnt/=count(mod(diff(istart:iend),5)==1)) then
        write(LOUT,'(a)') 'ERROR!!! &
             &Number of types was incorrectly predicted in create_PairReduced!'
        stop
     endif

   end associate

   istart = iend + 1
enddo
if(icnt/=count(diff(1:n)==1)) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Number of prim pairs was incorrectly predicted in create_PairReduced!'
   stop
endif

i=0
do k=1,Nexp
   do j=1,k
      i = i + 1
      if(PairReduced(i)%isUsed) then
         if(.not.(isUsed(j).and.isUsed(k))) then
            write(LOUT,'(a)') 'ERROR!!! &
                 &Used pair exponents have not been predicted properly!'
            stop
         endif
      endif
   enddo
enddo

end subroutine create_PairReduced

function find_end(istart,n,diff,divisor) result(iend)
implicit none
integer :: iend
integer,intent(in) :: istart,n,diff(:),divisor

iend = istart
do while(iend<n)
   if(mod(diff(iend+1),divisor)==1) exit
   iend = iend + 1
enddo

end function find_end

end module systemdef
