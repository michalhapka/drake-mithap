module SCFint
use file_OUT, only : LOUT
use precision, only : prec
use memory
use commontypes, only : OrbSystemData,OrbReducedData
use int12core
implicit none

private
public create_SCFint,free_SCFint
public SCFint_matH,SCFint_matS,SCFint_matJ,SCFint_matK

integer,parameter :: smallestOne = -1
integer,parameter :: offsetOne   = 1 - smallestOne

integer,parameter :: smallestTwo_uv = 0
integer,parameter :: offsetTwo_uv   = 1 - smallestTwo_uv
!integer,parameter :: onlyTwo_t      = -1
integer,parameter :: smallestTwo_t      = -1

type OneIntData
logical :: isUsed
integer :: iexp,jexp
real(prec) :: ijalphaPL,ijalphaMI
integer :: ijrange
real(prec),allocatable :: elms(:)
end type OneIntData

type TwoInt_t_Data
real(prec),allocatable :: elms(:,:)
integer :: shift
integer :: two_t
end type TwoInt_t_Data

type TwoIntData
logical :: isUsed
integer :: iexp,jexp,kexp,lexp
logical :: same_exp
real(prec) :: ijalpha,klalpha
integer :: ijrange,klrange
!real(prec),allocatable :: elms(:,:)
type(TwoInt_t_Data),allocatable :: t_val(:) 
!type(TwoInt_t_Data),allocatable :: lam_val(:) 
end type TwoIntData

type(OneIntData),allocatable :: OneInt(:)
type(TwoIntData),allocatable :: TwoInt(:)

contains

subroutine create_SCFint(Nexp,exponents,OrbReduced,LPRINT)
implicit none
integer,intent(in) :: Nexp
real(prec),intent(in) :: exponents(:)
type(OrbReducedData),intent(in) :: OrbReduced(:)
integer,intent(in) :: LPRINT

call init_OneInt(Nexp)
call create_OneInt(Nexp,exponents,OrbReduced)
call print_OneInt(LPRINT)

call init_TwoInt(Nexp)
call create_TwoInt(Nexp,exponents,OrbReduced)
call print_TwoInt(LPRINT)

end subroutine create_SCFint

subroutine free_SCFint
implicit none

call free_TwoInt

call free_OneInt

end subroutine free_SCFint

subroutine SCFint_matH(matH,iOrbSystem,jOrbSystem,nucZ)
implicit none
real(prec) :: matH(:,:)
type(OrbSystemData),intent(in) :: iOrbSystem,jOrbSystem
real(prec),intent(in) :: nucZ
integer :: iOrbl,jOrbl
integer :: i_prim,j_prim
integer :: iexp,jexp,ijmin,ijmax
logical :: ijorder
integer :: i,j,ipos,jpos,ijelms
real(prec) :: alphaPL,alphaMI,A0,A1_pre,A1
integer :: aPL,aMI,A2
real(prec) :: tmp

do j_prim=1,jOrbSystem%n_prim
   do i_prim=1,iOrbSystem%n_prim
      associate(&
           iOrbSpec => iOrbSystem%OrbSpec(i_prim), &
           jOrbSpec => jOrbSystem%OrbSpec(j_prim))

        iOrbl = iOrbSystem%lorb
        jOrbl = jOrbSystem%lorb

        iexp = iOrbSpec%iexp
        jexp = jOrbSpec%iexp

      ! write(*,*) iexp,jexp,i_prim,j_prim

        ijmin = min(iexp,jexp)
        ijmax = max(iexp,jexp)

        ijorder = (iexp<=jexp)

        associate(One => OneInt(ijmin + ijmax*(ijmax-1)/2))
          if(.not.One%isUsed) then
             write(LOUT,'(a)') 'ERROR!!! &
                  &Non-existing entry in SCFint_matH!'
             stop
          endif
          if(ijmin/=One%iexp.or.ijmax/=One%jexp) then
             write(LOUT,'(a)') 'ERROR!!! &
                  &Wrong entry chosen in SCFint_matH!'
             stop
          endif

          alphaPL = One%ijalphaPL
          alphaMI = One%ijalphaMI
          if(.not.ijorder) alphaMI = -alphaMI

          A0     = -(alphaPL**2 + alphaMI**2)/2._prec
          A1_pre = -4._prec*nucZ

          jpos = jOrbSpec%offset
          do j=0,jOrbSpec%irange
             jpos = jpos + 1
             ipos = iOrbSpec%offset
             do i=0,iOrbSpec%irange
                ipos = ipos + 1

                aPL = i + j + iOrbl + jOrbl
                aMI = i - j + iOrbl - jOrbl

                ijelms = offsetOne + aPL  

                A1  = A1_pre + alphaPL*(aPL+2) + alphaMI*aMI
                tmp = A0*One%elms(ijelms) + A1*One%elms(ijelms - 1)

                if(aPL>0) then
                   A2  = -(aPL**2 + aMI**2)/2 - aPL & 
                       & + iOrbl*(iOrbl + 1) + jOrbl*(jOrbl + 1)
                   tmp = tmp + A2*One%elms(ijelms - 2)
                endif

                matH(ipos,jpos) = tmp/4._prec

             enddo
          enddo

        end associate

      end associate
   enddo
enddo

end subroutine SCFint_matH

subroutine SCFint_matS(matS,iOrbSystem,jOrbSystem)
implicit none
real(prec) :: matS(:,:)
type(OrbSystemData),intent(in) :: iOrbSystem,jOrbSystem
integer :: iOrbl, jOrbl
integer :: i_prim,j_prim
integer :: iexp,jexp,ijmin,ijmax
integer :: i,j,ipos,jpos


do j_prim=1,jOrbSystem%n_prim
   do i_prim=1,iOrbSystem%n_prim
      associate(&
           iOrbSpec => iOrbSystem%OrbSpec(i_prim), &
           jOrbSpec => jOrbSystem%OrbSpec(j_prim))

        iOrbl = iOrbSystem%lorb
        jOrbl = jOrbSystem%lorb

        iexp = iOrbSpec%iexp
        jexp = jOrbSpec%iexp

        ijmin = min(iexp,jexp)
        ijmax = max(iexp,jexp)

        associate(One => OneInt(ijmin + ijmax*(ijmax-1)/2))
          if(.not.One%isUsed) then
             write(LOUT,'(a)') 'ERROR!!! &
                  &Non-existing entry in SCFint_matS!'
             stop
          endif
          if(ijmin/=One%iexp.or.ijmax/=One%jexp) then
             write(LOUT,'(a)') 'ERROR!!! &
                  &Wrong entry chosen in SCFint_matS!'
             stop
          endif

          jpos = jOrbSpec%offset
          do j=0,jOrbSpec%irange
             jpos = jpos + 1
             ipos = iOrbSpec%offset
             do i=0,iOrbSpec%irange
                ipos = ipos + 1

                matS(ipos,jpos) = One%elms(offsetOne + i + j + iOrbl + jOrbl)

             enddo
          enddo

        end associate

      end associate
   enddo
enddo

end subroutine SCFint_matS

subroutine SCFint_matJ(matJ,iOrbSystem,jOrbSystem,kOrbSystem,lOrbSystem,prefac)
use misc, only : swap
implicit none
real(prec) :: matJ(:,:)
real(prec),intent(in) :: prefac
type(OrbSystemData),intent(in) :: iOrbSystem,jOrbSystem,kOrbSystem,lOrbSystem
integer :: i_prim,j_prim,k_prim,l_prim
integer :: iexp,jexp,kexp,lexp,ijmin,ijmax,klmin,klmax,ijexp,klexp
logical :: ijklorder
integer :: i_nbas,k_nbas
integer :: i,j,k,l,ijpos,jpos,klpos,lpos,klelms
integer :: lorb1,lorb2,lorb3,lorb4

i_nbas = iOrbSystem%nbas
k_nbas = kOrbSystem%nbas

do l_prim=1,lOrbSystem%n_prim
   do k_prim=1,kOrbSystem%n_prim
      do j_prim=1,jOrbSystem%n_prim
         do i_prim=1,iOrbSystem%n_prim
            associate(&
                 iOrbSpec => iOrbSystem%OrbSpec(i_prim), &
                 jOrbSpec => jOrbSystem%OrbSpec(j_prim), &
                 kOrbSpec => kOrbSystem%OrbSpec(k_prim), &
                 lOrbSpec => lOrbSystem%OrbSpec(l_prim))

              lorb1 = iOrbSystem%lorb
              lorb2 = jOrbSystem%lorb
              lorb3 = kOrbSystem%lorb
              lorb4 = lOrbSystem%lorb

              iexp = iOrbSpec%iexp
              jexp = jOrbSpec%iexp
              kexp = kOrbSpec%iexp
              lexp = lOrbSpec%iexp

              ijmin = min(iexp,jexp)
              ijmax = max(iexp,jexp)
              klmin = min(kexp,lexp)
              klmax = max(kexp,lexp)

              ijexp = ijmin + ijmax*(ijmax-1)/2
              klexp = klmin + klmax*(klmax-1)/2

              ijklorder = (ijexp<=klexp)
              if(.not.ijklorder) then
                 call swap(ijmin,klmin)
                 call swap(ijmax,klmax)
                 call swap(ijexp,klexp)
              endif

              !associate(Two => TwoInt(ijexp + klexp*(klexp-1)/2))
              associate(Two => TwoInt(ijexp + klexp*(klexp-1)/2))
                if(.not.Two%isUsed) then
                   write(LOUT,'(a)') 'ERROR!!! &
                        &Non-existing entry in SCFint_matJ!'
                   stop
                endif
                if(ijmin/=Two%iexp.or.ijmax/=Two%jexp.or.&
                     klmin/=Two%kexp.or.klmax/=Two%lexp) then
                   write(LOUT,'(a)') 'ERROR!!! &
                        &Wrong entry chosen in SCFint_matJ!'
                   stop
                endif

                if(ijklorder) then

                   lpos = lOrbSpec%offset*k_nbas
                   do l=0,lOrbSpec%irange
                      klpos = lpos + kOrbSpec%offset
                      do k=0,kOrbSpec%irange
                         klpos = klpos + 1

                         klelms = offsetTwo_uv + k + l + lorb3 + lorb4
                         !klelms = offsetTwo_uv + k + l

                         jpos = jOrbSpec%offset*i_nbas
                         do j=0,jOrbSpec%irange
                            ijpos = jpos + iOrbSpec%offset
                            do i=0,iOrbSpec%irange
                               ijpos = ijpos + 1

                               matJ(ijpos,klpos) = prefac*Two%t_val(1)%elms( &
                                    offsetTwo_uv + i + j + lorb1 + lorb2, &
                                    klelms)

!                               write(*,*) matJ(ijpos,klpos), i+1, j+1, k+1, l+1
!                               matJ(ijpos,klpos) = Two%elms( &
!                                    offsetTwo_uv + i + j, &
!                                    klelms)

                            enddo
                            jpos = jpos + i_nbas
                         enddo

                      enddo
                      lpos = lpos + k_nbas
                   enddo

                else

                   lpos = lOrbSpec%offset*k_nbas
                   do l=0,lOrbSpec%irange
                      klpos = lpos + kOrbSpec%offset
                      do k=0,kOrbSpec%irange
                         klpos = klpos + 1

                         klelms = offsetTwo_uv + k + l + lorb3 + lorb4
                         !klelms = offsetTwo_uv + k + l

                         jpos = jOrbSpec%offset*i_nbas
                         do j=0,jOrbSpec%irange
                            ijpos = jpos + iOrbSpec%offset
                            do i=0,iOrbSpec%irange
                               ijpos = ijpos + 1

                               matJ(ijpos,klpos) = prefac*Two%t_val(1)%elms( &
                                    klelms, &
                                    offsetTwo_uv + i + j + lorb1 + lorb2)
!                               write(*,*) matJ(ijpos,klpos), i+1, j+1, k+1, l+1

!                               matJ(ijpos,klpos) = Two%elms( &
!                                    klelms, &
!                                    offsetTwo_uv + i + j)

                            enddo
                            jpos = jpos + i_nbas
                         enddo

                      enddo
                      lpos = lpos + k_nbas
                   enddo

                endif

              end associate

            end associate
         enddo
      enddo
   enddo
enddo

end subroutine SCFint_matJ

subroutine SCFint_matK(matK,iOrbSystem,jOrbSystem,kOrbSystem,lOrbSystem,prefac)
use misc, only : swap
implicit none
real(prec) :: matK(:,:)
real(prec) :: prefac
type(OrbSystemData),intent(in) :: iOrbSystem,jOrbSystem,kOrbSystem,lOrbSystem
integer :: i_prim,j_prim,k_prim,l_prim
integer :: iexp,jexp,kexp,lexp,ikmin,ikmax,jlmin,jlmax,ikexp,jlexp
logical :: ikjlorder
integer :: i_nbas,k_nbas
integer :: i,j,k,l,ijpos,jpos,klpos,lpos,jlelms
integer :: lorb1,lorb2,lorb3,lorb4

i_nbas = iOrbSystem%nbas
k_nbas = kOrbSystem%nbas

do l_prim=1,lOrbSystem%n_prim
   do k_prim=1,kOrbSystem%n_prim
      do j_prim=1,jOrbSystem%n_prim
         do i_prim=1,iOrbSystem%n_prim
            associate(&
                 iOrbSpec => iOrbSystem%OrbSpec(i_prim), &
                 jOrbSpec => jOrbSystem%OrbSpec(j_prim), &
                 kOrbSpec => kOrbSystem%OrbSpec(k_prim), &
                 lOrbSpec => lOrbSystem%OrbSpec(l_prim))

              lorb1 = iOrbSystem%lorb
              lorb2 = jOrbSystem%lorb
              lorb3 = kOrbSystem%lorb
              lorb4 = lOrbSystem%lorb

              iexp = iOrbSpec%iexp
              jexp = jOrbSpec%iexp
              kexp = kOrbSpec%iexp
              lexp = lOrbSpec%iexp

              ikmin = min(iexp,kexp)
              ikmax = max(iexp,kexp)
              jlmin = min(jexp,lexp)
              jlmax = max(jexp,lexp)

              ikexp = ikmin + ikmax*(ikmax-1)/2
              jlexp = jlmin + jlmax*(jlmax-1)/2

              ikjlorder = (ikexp<=jlexp)
              if(.not.ikjlorder) then
                 call swap(ikmin,jlmin)
                 call swap(ikmax,jlmax)
                 call swap(ikexp,jlexp)
              endif

              associate(Two => TwoInt(ikexp + jlexp*(jlexp-1)/2))
                if(.not.Two%isUsed) then
                   write(LOUT,'(a)') 'ERROR!!! &
                        &Non-existing entry in SCFint_matK!'
                   stop
                endif
                if(ikmin/=Two%iexp.or.ikmax/=Two%jexp.or.&
                     jlmin/=Two%kexp.or.jlmax/=Two%lexp) then
                   write(LOUT,'(a)') 'ERROR!!! &
                        &Wrong entry chosen in SCFint_matK!'
                   stop
                endif

                if(ikjlorder) then

                   lpos = lOrbSpec%offset*k_nbas
                   do l=0,lOrbSpec%irange
                      klpos = lpos + kOrbSpec%offset
                      do k=0,kOrbSpec%irange
                         klpos = klpos + 1

                         jpos = jOrbSpec%offset*i_nbas
                         do j=0,jOrbSpec%irange
                            ijpos = jpos + iOrbSpec%offset
                            jlelms = offsetTwo_uv + j + lorb2 + l + lorb4
                            !jlelms = offsetTwo_uv + j + l
                            do i=0,iOrbSpec%irange
                               ijpos = ijpos + 1

                               matK(ijpos,klpos) = prefac*Two%t_val(1)%elms( &
                                    offsetTwo_uv + i + k + lorb1 + lorb3, &
                                    jlelms)


                               write(*,*) matK(ijpos,klpos), i+1,j+1,k+1,l+1
!                               matK(ijpos,klpos) = Two%elms( &
!                                    offsetTwo_uv + i + k, &
!                                    jlelms)

                            enddo
                            jpos = jpos + i_nbas
                         enddo

                      enddo
                      lpos = lpos + k_nbas
                   enddo

                else

                   lpos = lOrbSpec%offset*k_nbas
                   do l=0,lOrbSpec%irange
                      klpos = lpos + kOrbSpec%offset
                      do k=0,kOrbSpec%irange
                         klpos = klpos + 1

                         jpos = jOrbSpec%offset*i_nbas
                         do j=0,jOrbSpec%irange
                            ijpos = jpos + iOrbSpec%offset
                            jlelms = offsetTwo_uv + j + lorb2 + l + lorb4
                            do i=0,iOrbSpec%irange
                               ijpos = ijpos + 1

                               matK(ijpos,klpos) = Two%t_val(1)%elms( &
                                    jlelms, &
                                    offsetTwo_uv + i + k + lorb1 + lorb3)
                                  
                               !write(*,*) matK(ijpos,klpos), i+1,j+1,k+1,l+1
!                               matK(ijpos,klpos) = Two%elms( &
!                                    jlelms, &
!                                    offsetTwo_uv + i + k)

                            enddo
                            jpos = jpos + i_nbas
                         enddo

                      enddo
                      lpos = lpos + k_nbas
                   enddo

                endif

              end associate

            end associate
         enddo
      enddo
   enddo
enddo

end subroutine SCFint_matK

subroutine init_OneInt(Nexp)
implicit none
integer,intent(in) :: Nexp
integer :: iOne

allocate(OneInt(Nexp*(Nexp+1)/2))

do iOne=1,size(OneInt)
   associate(One => OneInt(iOne))

     One%isUsed    = .false.
     One%iexp      = 0
     One%jexp      = 0
     One%ijalphaPL = 0._prec
     One%ijalphaMI = 0._prec
     One%ijrange   = smallestOne - 1

   end associate
enddo

end subroutine init_OneInt

subroutine free_OneInt
implicit none
integer :: iOne

do iOne=1,size(OneInt)
   associate(One => OneInt(iOne))

     if(One%isUsed) call mem_dealloc(One%elms)

   end associate
enddo

deallocate(OneInt)

end subroutine free_OneInt

subroutine init_TwoInt(Nexp)
implicit none
integer,intent(in) :: Nexp
integer :: iTwo

iTwo = Nexp*(Nexp+1)/2
allocate(TwoInt(iTwo*(iTwo+1)/2))

do iTwo=1,size(TwoInt)
   associate(Two => TwoInt(iTwo))

     Two%isUsed  = .false.
     Two%iexp    = 0
     Two%jexp    = 0
     Two%kexp    = 0
     Two%lexp    = 0
     Two%ijalpha = 0._prec
     Two%klalpha = 0._prec
     Two%ijrange = smallestTwo_uv - 1
     Two%klrange = smallestTwo_uv - 1

   end associate
enddo

end subroutine init_TwoInt

subroutine free_TwoInt
implicit none
integer :: iTwo, j

do iTwo=1,size(TwoInt)
   associate(Two => TwoInt(iTwo))

      if(Two%isUsed) then

         do j=1,size(Two%t_val)
            call mem_dealloc(Two%t_val(j)%elms)
         enddo

         deallocate(Two%t_val)
      endif

! hapka: old_drake
!     if(Two%isUsed) call mem_dealloc(Two%elms)

   end associate
enddo

deallocate(TwoInt)

end subroutine free_TwoInt

subroutine choose_maxrange(iOrbSpec,jOrbSpec,maxrange,doint)
implicit none
type(OrbReducedData) :: iOrbSpec,jOrbSpec
logical :: doint
integer :: maxrange, tmp
integer :: ival, jval

maxrange=0
doint=.false.

do jval=1,jOrbSpec%nlang
   do ival=1,iOrbSpec%nlang
   if(jOrbSpec%lang(jval).eq.iOrbSpec%lang(ival)) then
     tmp=maxrange
     maxrange = jOrbSpec%max_lrange(jval) + iOrbSpec%max_lrange(ival)
     if(tmp.gt.maxrange) maxrange=tmp
     doint=.true.
   endif

   enddo
enddo

end subroutine choose_maxrange

subroutine create_OneInt(Nexp,exponents,OrbReduced)
implicit none
integer,intent(in) :: Nexp
real(prec),intent(in) :: exponents(:)
type(OrbReducedData),intent(in) :: OrbReduced(:)
logical :: doint
integer :: maxrange
integer :: iOne
integer :: iexp,jexp
integer :: i

iOne = 0
do jexp=1,Nexp
   do iexp=1,jexp
      iOne = iOne + 1
      associate(&
           iOrbSpec => OrbReduced(iexp), &
           jOrbSpec => OrbReduced(jexp))
           call choose_maxrange(iOrbSpec,jOrbSpec,maxrange,doint)
        if(iOrbSpec%isUsed.and.jOrbSpec%isUsed.and.doint) then
           associate(One => OneInt(iOne))

             One%isUsed    = .true.
             One%iexp      = iOrbSpec%iexp
             One%jexp      = jOrbSpec%iexp
             One%ijalphaPL = exponents(One%iexp) + exponents(One%jexp)
             One%ijalphaMI = exponents(One%iexp) - exponents(One%jexp)
            !One%ijrange   = iOrbSpec%maxrange   + jOrbSpec%maxrange
             One%ijrange   = maxrange 

             call mem_alloc(One%elms,One%ijrange-smallestOne+1)

             do i=smallestOne,One%ijrange
                One%elms(offsetOne + i) = int1_slater(i,One%ijalphaPL)
             enddo

           end associate
        endif
      end associate
   enddo
enddo

if(iOne/=size(OneInt)) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Number of entries in SCF OneInt has not been predicted correctly!'
   stop
endif

! hapka: old_drake
!iOne = 0
!do jexp=1,Nexp
!   do iexp=1,jexp
!      iOne = iOne + 1
!      associate(&
!           iOrbSpec => OrbReduced(iexp), &
!           jOrbSpec => OrbReduced(jexp))
!        if(iOrbSpec%isUsed.and.jOrbSpec%isUsed) then
!           associate(One => OneInt(iOne))
!
!             One%isUsed    = .true.
!             One%iexp      = iOrbSpec%iexp
!             One%jexp      = jOrbSpec%iexp
!             One%ijalphaPL = exponents(One%iexp) + exponents(One%jexp)
!             One%ijalphaMI = exponents(One%iexp) - exponents(One%jexp)
!             One%ijrange   = iOrbSpec%maxrange   + jOrbSpec%maxrange
!
!             call mem_alloc(One%elms,One%ijrange-smallestOne+1)
!
!             do i=smallestOne,One%ijrange
!                One%elms(offsetOne + i) = int1_slater(i,One%ijalphaPL)
!             enddo
!
!           end associate
!        endif
!      end associate
!   enddo
!enddo


end subroutine create_OneInt

subroutine choose_maxrange_max4l(iOrbSpec,jOrbSpec,kOrbSpec,lOrbSpec,&
                               & Max4lOdd,Max4lEven,maxrange1,maxrange2,doint)

implicit none
type(OrbReducedData) :: iOrbSpec,jOrbSpec
type(OrbReducedData) :: kOrbSpec,lOrbSpec
logical :: doint
integer :: Max4lOdd, Max4lEven
integer :: TmpOdd, TmpEven
integer :: maxrange1, maxrange2
integer :: tmp0, tmp1, tmp2
integer :: ival, jval, kval, lval
integer :: ijM, klM, ijP, klP

 Max4lEven = -1
 Max4lOdd  = -1
 maxrange1=0
 maxrange2=0
 doint=.false.
 do lval=1,lOrbSpec%nlang
    do kval=1,kOrbSpec%nlang
       do jval=1,jOrbSpec%nlang
          do ival=1,iOrbSpec%nlang

             ijM = abs(iOrbSpec%lang(ival)-jOrbSpec%lang(jval))
             klM = abs(kOrbSpec%lang(kval)-lOrbSpec%lang(lval))
             ijP = iOrbSpec%lang(ival)+jOrbSpec%lang(jval)
             klP = kOrbSpec%lang(kval)+lOrbSpec%lang(lval)

             ! check parity condition
             if(mod(ijP+klP,2).eq.0) then

                ! check coupling condition
                if(max(ijM,klM).le.min(ijP,klP)) then
                   tmp0=ijP+klP
                   if(mod(tmp0/2,2).eq.0) then
                      TmpEven = tmp0
                      if(TmpEven.gt.Max4lEven) Max4lEven = TmpEven
                   else
                      TmpOdd = tmp0
                      if(TmpOdd.gt.Max4lOdd) Max4lOdd = TmpOdd
                   endif



                   tmp1=maxrange1
                   tmp2=maxrange2
                   maxrange1 = jOrbSpec%max_lrange(jval) + iOrbSpec%max_lrange(ival)
                   maxrange2 = kOrbSpec%max_lrange(kval) + lOrbSpec%max_lrange(lval)
   !               write(*,'(a,*(i7))') '(ij|kl):', iexp, jexp, kexp, lexp 
   !               print*, "max_ij, max_kl", maxrange1, maxrange2
                   if(tmp1.gt.maxrange1) maxrange1 = tmp1
                   if(tmp2.gt.maxrange2) maxrange2 = tmp2
                   doint=.true.
                endif

             endif

          enddo
       enddo
    enddo
 enddo

end subroutine choose_maxrange_max4l

subroutine get_range_shift(Max4lEven,Max4lOdd,shift)
implicit none
integer :: Max4lEven, Max4lOdd
integer :: Max4l
integer :: shift(:)
integer,allocatable :: ShiftEven(:), ShiftOdd(:)
integer :: i

 Max4l=max(Max4lEven,Max4lOdd)
 allocate(ShiftEven(Max4l/2+1),ShiftOdd(Max4l/2+1))

 ShiftEven = -huge(0)
 ShiftOdd  = -huge(0)

 if(Max4lEven.ge.0) then
    do i=0,Max4lEven/2

       ShiftEven(i+1) = -(i/2)*2

    enddo
 endif

 if(Max4lOdd.ge.0) then
    do i=0,Max4lOdd/2

       ShiftOdd(i+1) = 1 - ((1+i)/2)*2

    enddo
 endif

 do i=1,Max4l/2+1
    shift(i) = max(ShiftEven(i),ShiftOdd(i))
 enddo

 deallocate(ShiftEven,ShiftOdd)

end subroutine get_range_shift

subroutine create_TwoInt(Nexp,exponents,OrbReduced)
implicit none
integer,intent(in) :: Nexp
real(prec),intent(in) :: exponents(:)
type(OrbReducedData),intent(in) :: OrbReduced(:)
integer :: iTwo
integer :: iexp,jexp,kexp,lexp
logical :: doint
integer :: max4l,Max4lEven,Max4lOdd
integer :: maxrange_ij,maxrange_kl
integer, allocatable :: shift(:)
integer :: currentTwo_t
integer :: i,j,k

iTwo = 0
do lexp=1,Nexp
   do kexp=1,lexp
      do jexp=1,lexp
         do iexp=1,merge(kexp,jexp,jexp==lexp)
            iTwo = iTwo + 1
            associate(&
                 iOrbSpec => OrbReduced(iexp), &
                 jOrbSpec => OrbReduced(jexp), &
                 kOrbSpec => OrbReduced(kexp), &
                 lOrbSpec => OrbReduced(lexp))
                 call choose_maxrange_max4l(iOrbSpec,jOrbSpec,kOrbSpec,lOrbSpec, &
                      & Max4lOdd,Max4lEven,maxrange_ij,maxrange_kl,doint)
                ! if(doint) then
                !    write(*,'(a,*(i7))') '(ij|kl) ', iexp,jexp,kexp,lexp
                !    write(*,'(a,*(i7))') 'max4l ', max4l,maxrange_ij,maxrange_kl
                ! endif
                ! write(*,*) 'doint ', doint
              if(iOrbSpec%isUsed.and.jOrbSpec%isUsed.and.&
                   kOrbSpec%isUsed.and.lOrbSpec%isUsed.and.doint) then
                 associate(Two => TwoInt(iTwo))

                   Two%isUsed   = .true.
                   Two%iexp     = iOrbSpec%iexp
                   Two%jexp     = jOrbSpec%iexp
                   Two%kexp     = kOrbSpec%iexp
                   Two%lexp     = lOrbSpec%iexp
                   Two%same_exp = (Two%iexp==Two%kexp).and.(Two%jexp==Two%lexp)
                   Two%ijalpha  = exponents(Two%iexp) + exponents(Two%jexp)
                   Two%klalpha  = exponents(Two%kexp) + exponents(Two%lexp)
                   Two%ijrange  = maxrange_ij
                   Two%klrange  = maxrange_kl
                  !Two%ijrange  = iOrbSpec%maxrange   + jOrbSpec%maxrange
                  !Two%klrange  = kOrbSpec%maxrange   + lOrbSpec%maxrange

                   max4l = max(Max4lOdd,Max4lEven)
                   allocate(shift(max4l/2+1))

                   call get_range_shift(Max4lEven,Max4lOdd,shift)                   
                   allocate(Two%t_val(1+max4l/2))

                   do i=1,max4l/2+1  
                      associate( tval => Two%t_val(i) )
                      tval%shift = shift(i)
                      tval%two_t = 2*i-3
                      call mem_alloc(tval%elms, &
                           Two%ijrange+tval%shift-smallestTwo_uv+1,&
                           Two%klrange+tval%shift-smallestTwo_uv+1)
                      end associate
                   enddo

                   do k=1,max4l/2+1

                      currentTwo_t = 2*k-3

                      if(Two%same_exp) then
   
                         do j=smallestTwo_uv,Two%klrange+shift(k)
                            do i=smallestTwo_uv,min(j,Two%ijrange+shift(k))
                               Two%t_val(k)%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
                                    !int2_slater(i,j,onlyTwo_t,&
                                   int2_slater(i,j,currentTwo_t,&
                                   Two%ijalpha,Two%klalpha)
                            enddo
                         enddo
                         do j=smallestTwo_uv,min(Two%ijrange+shift(k),Two%klrange+shift(k))
                            do i=j+1,min(Two%ijrange+shift(k),Two%klrange+shift(k))
                               Two%t_val(k)%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
                                    Two%t_val(k)%elms(offsetTwo_uv + j,offsetTwo_uv + i)
                            enddo
                         enddo
                         if(Two%ijrange+shift(k)>Two%klrange+shift(k)) then
                            do j=smallestTwo_uv,Two%klrange+shift(k)
                               do i=Two%klrange+shift(k)+1,Two%ijrange+shift(k)
                                  Two%t_val(k)%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
                                       !int2_slater(i,j,onlyTwo_t,&
                                      int2_slater(i,j,currentTwo_t,&
                                      Two%ijalpha,Two%klalpha)
                               enddo
                            enddo
                         endif
   
                      else
   
                         do j=smallestTwo_uv,Two%klrange+shift(k)
                            do i=smallestTwo_uv,Two%ijrange+shift(k)
                               Two%t_val(k)%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
                                    !int2_slater(i,j,onlyTwo_t,&
                                    int2_slater(i,j,currentTwo_t,&
                                    Two%ijalpha,Two%klalpha)
                            enddo
                         enddo
   
                      endif
                 
                   enddo

                   deallocate(shift)

                 end associate
              endif
            end associate
         enddo
      enddo
   enddo
enddo
if(iTwo/=size(TwoInt)) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Number of entries in SCF TwoInt has not been predicted correctly!'
   stop
endif

! hapka: old_drake
! iTwo = 0
! do lexp=1,Nexp
!    do kexp=1,lexp
!       do jexp=1,lexp
!          do iexp=1,merge(kexp,jexp,jexp==lexp)
!             iTwo = iTwo + 1
!             associate(&
!                  iOrbSpec => OrbReduced(iexp), &
!                  jOrbSpec => OrbReduced(jexp), &
!                  kOrbSpec => OrbReduced(kexp), &
!                  lOrbSpec => OrbReduced(lexp))
!               if(iOrbSpec%isUsed.and.jOrbSpec%isUsed.and.&
!                    kOrbSpec%isUsed.and.lOrbSpec%isUsed) then
!                  associate(Two => TwoInt(iTwo))
! 
!                    Two%isUsed   = .true.
!                    Two%iexp     = iOrbSpec%iexp
!                    Two%jexp     = jOrbSpec%iexp
!                    Two%kexp     = kOrbSpec%iexp
!                    Two%lexp     = lOrbSpec%iexp
!                    Two%same_exp = (Two%iexp==Two%kexp).and.(Two%jexp==Two%lexp)
!                    Two%ijalpha  = exponents(Two%iexp) + exponents(Two%jexp)
!                    Two%klalpha  = exponents(Two%kexp) + exponents(Two%lexp)
!                    Two%ijrange  = iOrbSpec%maxrange   + jOrbSpec%maxrange
!                    Two%klrange  = kOrbSpec%maxrange   + lOrbSpec%maxrange
! 
!                    call mem_alloc(Two%elms,&
!                         Two%ijrange-smallestTwo_uv+1,&
!                         Two%klrange-smallestTwo_uv+1)
! 
!                    if(Two%same_exp) then
! 
!                       do j=smallestTwo_uv,Two%klrange
!                          do i=smallestTwo_uv,min(j,Two%ijrange)
!                             Two%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
!                                  !int2_slater(i,j,onlyTwo_t,&
!                                  int2_slater(i,j,smallestTwo_t,&
!                                  Two%ijalpha,Two%klalpha)
!                          enddo
!                       enddo
!                       do j=smallestTwo_uv,min(Two%ijrange,Two%klrange)
!                          do i=j+1,min(Two%ijrange,Two%klrange)
!                             Two%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
!                                  Two%elms(offsetTwo_uv + j,offsetTwo_uv + i)
!                          enddo
!                       enddo
!                       if(Two%ijrange>Two%klrange) then
!                          do j=smallestTwo_uv,Two%klrange
!                             do i=Two%klrange+1,Two%ijrange
!                                Two%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
!                                     !int2_slater(i,j,onlyTwo_t,&
!                                     int2_slater(i,j,smallestTwo_t,&
!                                     Two%ijalpha,Two%klalpha)
!                             enddo
!                          enddo
!                       endif
! 
!                    else
! 
!                       do j=smallestTwo_uv,Two%klrange
!                          do i=smallestTwo_uv,Two%ijrange
!                             Two%elms(offsetTwo_uv + i,offsetTwo_uv + j) = &
!                                  !int2_slater(i,j,onlyTwo_t,&
!                                  int2_slater(i,j,smallestTwo_t,&
!                                  Two%ijalpha,Two%klalpha)
!                          enddo
!                       enddo
! 
!                    endif
! 
!                  end associate
!               endif
!             end associate
!          enddo
!       enddo
!    enddo
! enddo
! if(iTwo/=size(TwoInt)) then
!    write(LOUT,'(a)') 'ERROR!!! &
!         &Number of entries in SCF TwoInt has not been predicted correctly!'
!    stop
! endif

end subroutine create_TwoInt

subroutine print_OneInt(LPRINT)
implicit none
integer,intent(in) :: LPRINT
integer :: iOne
integer :: i

if(LPRINT>=10) then

   write(LOUT,'()')
   write(LOUT,'(5x,a)') '--- SCF one-electron integrals ---'

   write(LOUT,'()')
   write(LOUT,'(6x,a,2x,2(2x,a),8x,a,9x,a,4x,a)') &
        'no.','i','j','ijalphaPL','ijalphaMI','ijrange'
   do iOne=1,size(OneInt)
      associate(One => OneInt(iOne))
        write(LOUT,'(5x,i3,a)',advance='no') iOne,' : '
        if(One%isUsed) then
           write(LOUT,'(2i3,2f18.8,3x,i5)') &
                One%iexp,One%jexp,One%ijalphaPL,One%ijalphaMI,One%ijrange
           if(LPRINT>=100) then
              do i=smallestOne,One%ijrange
                 write(LOUT,'(10x,i3,a,es45.33)') i,' : ',&
                      One%elms(offsetOne + i)
              enddo
           endif
        else
           write(LOUT,'(a)') 'NOT USED'
        endif
      end associate
   enddo

endif

end subroutine print_OneInt

subroutine print_TwoInt(LPRINT)
implicit none
integer,intent(in) :: LPRINT
integer :: iTwo
integer :: i,j,k

if(LPRINT>=10) then

   write(LOUT,'()')
   write(LOUT,'(5x,a)') '--- SCF two-electron integrals ---'

   write(LOUT,'()')
   write(LOUT,'(6x,a,2x,4(2x,a),9x,a,3x,a,6x,a,3x,a)') &
        'no.','i','j','k','l','ijalpha','ijrange','klaplha','klrange'
   do iTwo=1,size(TwoInt)
      associate(Two => TwoInt(iTwo))
        write(LOUT,'(5x,i3,a)',advance='no') iTwo,' : '
        if(Two%isUsed) then
           write(LOUT,'(4i3,2(f18.8,i5))') &
                Two%iexp,Two%jexp,Two%kexp,Two%lexp,&
                Two%ijalpha,Two%ijrange,Two%klalpha,Two%klrange
           if(LPRINT>=100) then
              do k=1,size(Two%t_val)
                 associate(tval => Two%t_val(k))
                 do j=smallestTwo_uv,Two%klrange+tval%shift
                    do i=smallestTwo_uv,Two%ijrange+tval%shift
                       write(LOUT,'(10x,3i3,a,es30.22)') i,j,tval%two_t,' : ',&
                            tval%elms(offsetTwo_uv + i,offsetTwo_uv + j)
                    enddo
                 enddo
                 end associate
              enddo

           endif
        else
           write(LOUT,'(a)') 'NOT USED'
        endif
      end associate
   enddo

endif

end subroutine print_TwoInt

end module SCFint
