module CCint2
use file_OUT, only : LOUT
use precision, only : prec
use memory
use commontypes
use int12core
implicit none

private
public create_CCint2,free_CCint2
public CCint2_matH,CCint2_matS,CCint2_vec,CCint2_val

integer,parameter :: smallest_uv = -1
integer,parameter :: offset_uv   = 1 - smallest_uv
integer,parameter :: smallest_uv_vec(0:2) = [-1,0,1]
integer,parameter :: offset_uv_vec(0:2) = [1,1,1] - smallest_uv_vec

integer,parameter :: smallest_t = -1
integer,parameter :: offset_t   = 1 - smallest_t

type r12Data
real(prec),allocatable :: elms(:,:)
end type r12Data

type ILambdaData
integer :: range_t
integer,allocatable :: range_uv(:,:)
type(r12Data),allocatable :: r12(:)
end type ILambdaData

type TwoInt_t_Data
real(prec),allocatable :: elms(:,:)
integer :: shift
integer,allocatable :: level_range_uv(:,:)
integer :: two_t
integer :: lambda
end type TwoInt_t_Data

type TwoIntData
logical :: isUsed
integer :: iexp,jexp,kexp,lexp
logical :: same_exp
real(prec) :: ijalphaPL,ijalphaMI,klalphaPL,klalphaMI
type(ILambdaData) :: ILambda(0:2)
type(TwoInt_t_Data),allocatable :: t_val(:) 
end type TwoIntData

type(TwoIntData),allocatable :: TwoInt(:)

type InfoGenPair
logical :: isUsed
integer :: add_to_1
integer :: add_to_2
real(prec) :: coeff
end type InfoGenPair

contains

subroutine create_CCint2(Nexp,exponents,PairReduced,LPRINT)
implicit none
integer,intent(in) :: Nexp
real(prec),intent(in) :: exponents(:)
type(PairReducedData),intent(in) :: PairReduced(:)
integer,intent(in) :: LPRINT

call init1_TwoInt(Nexp,exponents)
call init2_TwoInt(PairReduced)
call init3_TwoInt(PairReduced)

call create_TwoInt

call print_TwoInt(LPRINT)

end subroutine create_CCint2

subroutine free_CCint2
implicit none

call free_TwoInt

end subroutine free_CCint2

subroutine CCint2_matH(ABtype,matH,iPairSystem,jPairSystem,nucZ)
implicit none
character(1),intent(in) :: ABtype
real(prec) :: matH(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: nucZ
integer :: i_prim,j_prim
integer :: iexp,jexp,kexp,lexp,iTwo
logical :: ijorder,klorder,ijklorder
integer :: ui,vi,ti,ipos,iflag
integer :: uj,vj,tj,jpos,jflag
real(prec) :: alphaPL,alphaMI,betaPL,betaMI
real(prec) :: A0_pre1,A0_pre2,B0_pre1,B0_pre2,X1_pre
integer :: aPL,aMI,bPL,bMI,cPL,cMI

do j_prim=1,jPairSystem%n_prim
   do i_prim=1,iPairSystem%n_prim
      associate(&
           iPairSpec => iPairSystem%PairSpec(i_prim), &
           jPairSpec => jPairSystem%PairSpec(j_prim))

        select case(ABtype)
        case('A','a')

           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp1
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp2
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(.not.Two%isUsed) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Non-existing entry in CCint2_matH!'
                stop
             endif
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in CCint2_matH!'
                stop
             endif

             call pre_element(Two)

             jpos  = jPairSpec%offset
             jflag = -1
             do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                jpos = jpos + 1

                ipos  = iPairSpec%offset
                iflag = -1
                do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                   ipos = ipos + 1

                   aPL = ui + uj; aMI = ui - uj
                   bPL = vi + vj; bMI = vi - vj
                   cPL = ti + tj; cMI = ti - tj

                   matH(ipos,jpos) = get_element(Two)

                enddo

             enddo

           end associate

        case('B','b')

           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp2
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp1
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(.not.Two%isUsed) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Non-existing entry in CCint2_matH!'
                stop
             endif
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in CCint2_matH!'
                stop
             endif

             call pre_element(Two)

             jpos  = jPairSpec%offset
             jflag = -1
             do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                jpos = jpos + 1

                ipos  = iPairSpec%offset
                iflag = -1
                do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                   ipos = ipos + 1

                   aPL = ui + vj; aMI = ui - vj
                   bPL = vi + uj; bMI = vi - uj
                   cPL = ti + tj; cMI = ti - tj

                   matH(ipos,jpos) = get_element(Two)

                enddo

             enddo

           end associate

        case default

           write(LOUT,'(a)') 'ERROR!!! Incorrect ABtype in CCint2_matH!'
           stop

        end select

      end associate
   enddo
enddo

contains

  subroutine pre_element(Two)
  implicit none
  type(TwoIntData),intent(in) :: Two

  if(ijklorder) then

     alphaPL = Two%ijalphaPL
     alphaMI = Two%ijalphaMI
     if(.not.ijorder) alphaMI = -alphaMI

     betaPL = Two%klalphaPL
     betaMI = Two%klalphaMI
     if(.not.klorder) betaMI = -betaMI

  else

     alphaPL = Two%klalphaPL
     alphaMI = Two%klalphaMI
     if(.not.klorder) alphaMI = -alphaMI

     betaPL = Two%ijalphaPL
     betaMI = Two%ijalphaMI
     if(.not.ijorder) betaMI = -betaMI

  endif

  A0_pre1 = -(alphaPL**2 + alphaMI**2)/2._prec
  A0_pre2 = alphaPL*alphaMI

  B0_pre1 = -(betaPL**2 + betaMI**2)/2._prec
  B0_pre2 = betaPL*betaMI

  X1_pre = -4._prec*nucZ

  end subroutine pre_element

  function get_element(Two) result(val)
  implicit none
  real(prec) :: val
  type(TwoIntData),intent(in) :: Two
  logical :: isA2,isB2,isC2
  real(prec) :: A0,A1,A2
  real(prec) :: B0,B1,B2
  real(prec) :: cRATIO,C2
  integer :: pos_u,pos_v,pos_t

  isA2 = (aPL>0)
  isB2 = (bPL>0)
  isC2 = (cPL>0)

  A0 = A0_pre1
  A1 = X1_pre + alphaPL*(aPL+2) + alphaMI*aMI
  if(isA2) A2 = -(aPL**2 + aMI**2)/2 - aPL

  B0 = B0_pre1
  B1 = X1_pre + betaPL*(bPL+2) + betaMI*bMI
  if(isB2) B2 = -(bPL**2 + bMI**2)/2 - bPL

  if(isC2) then

     cRATIO = real(cMI,prec)/cPL

     A0 = A0 + A0_pre2*cRATIO
     A1 = A1 - (alphaPL*aMI + alphaMI*(aPL+2))*cRATIO
     if(isA2) A2 = A2 + (aMI*(aPL+1))*cRATIO

     B0 = B0 + B0_pre2*cRATIO
     B1 = B1 - (betaPL*bMI + betaMI*(bPL+2))*cRATIO
     if(isB2) B2 = B2 + (bMI*(bPL+1))*cRATIO

     C2 = cPL**2 - cMI**2

  endif

  pos_t = offset_t  + cPL

!  associate(elms => Two%r12(pos_t)%elms)
!    if(ijklorder) then
!
!       pos_u = offset_uv + aPL
!       pos_v = offset_uv + bPL
!
!       val = (A0 + B0)*elms(pos_u,pos_v) &
!            + A1*elms(pos_u-1,pos_v) &
!            + B1*elms(pos_u,pos_v-1)
!       if(isA2) val = val + A2*elms(pos_u-2,pos_v)
!       if(isB2) val = val + B2*elms(pos_u,pos_v-2)
!
!    else
!
!       pos_u = offset_uv + bPL
!       pos_v = offset_uv + aPL
!
!       val = (B0 + A0)*elms(pos_u,pos_v) &
!            + B1*elms(pos_u-1,pos_v) &
!            + A1*elms(pos_u,pos_v-1)
!       if(isB2) val = val + B2*elms(pos_u-2,pos_v)
!       if(isA2) val = val + A2*elms(pos_u,pos_v-2)
!
!    endif
!  end associate
!
!  if(isC2) val = val + C2*Two%r12(pos_t-2)%elms(pos_u,pos_v)

  val = val/4._prec

  end function get_element

end subroutine CCint2_matH

subroutine CCint2_matS(ABtype,add_t,n1,n2,matS,iGen,jGen,iPairSystem,jPairSystem)
implicit none
character(1),intent(in) :: ABtype
integer,intent(in) :: add_t
integer,intent(in) :: n1,n2
integer,intent(in) :: iGen,jGen
!real(prec) :: matS(:,:)
real(prec) :: matS(n1,n2)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
integer :: i_prim,j_prim
integer :: iexp,jexp,kexp,lexp,iTwo
logical :: ijorder,klorder,ijklorder
integer :: ui,vi,ti,ipos,iflag
integer :: uj,vj,tj,jpos,jflag
integer :: ilam
type(InfoGenPair) :: info(0:2)
integer ::i 

if(add_t<smallest_t) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Additional t parameter in CCint2_matS is incorrect!'
   stop
endif

call get_genpair_info(iGen,jGen,ABtype,info)

matS = 0
do j_prim=1,jPairSystem%n_prim
   do i_prim=1,iPairSystem%n_prim
      associate(&
           iPairSpec => iPairSystem%PairSpec(i_prim), &
           jPairSpec => jPairSystem%PairSpec(j_prim))

        select case(ABtype)
        case('A','a')

           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp1
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp2
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(.not.Two%isUsed) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Non-existing entry in CCint2_matS!'
                stop
             endif
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in CCint2_matS!'
                stop
             endif

             if(ijklorder) then

                do ilam=0,2
                   if(info(ilam)%isUsed) then
                      associate(&
                           offset_uv => offset_uv_vec(ilam),&
                           coeff    => info(ilam)%coeff,&
                           add_to_1 => info(ilam)%add_to_1,&
                           add_to_2 => info(ilam)%add_to_2,&
                           TwoL => Two%ILambda(ilam))

                        !write(*,*) 'A-even',ilam
                        !write(*,*) iGen,jGen,coeff
                        !write(*,*) add_to_1,add_to_2
                        jpos  = jPairSpec%offset
                        jflag = -1
                        do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                           jpos = jpos + 1

                           ipos  = iPairSpec%offset
                           iflag = -1
                           do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                              ipos = ipos + 1
                               
                              matS(ipos,jpos) = matS(ipos,jpos) + & 
                                   coeff * TwoL%r12(&
                                   offset_t  + ti + tj + add_t)%elms( &
                                   offset_uv + ui + uj + add_to_1, &
                                   offset_uv + vi + vj + add_to_2)

                           !   matS(ipos,jpos) = Two%r12(&
                           !        offset_t  + ti + tj + add_t)%elms( &
                           !        offset_uv + ui + uj, &
                           !        offset_uv + vi + vj)

                           enddo

                        enddo
                        end associate
                   endif
                enddo

             else

                do ilam=0,2
                   if(info(ilam)%isUsed) then
                      associate(&
                           offset_uv => offset_uv_vec(ilam),&
                           coeff    => info(ilam)%coeff,&
                           add_to_1 => info(ilam)%add_to_1,&
                           add_to_2 => info(ilam)%add_to_2,&
                           TwoL => Two%ILambda(ilam))
                        
                        !write(*,*) 'A-odd',ilam
                        !write(*,*) iGen,jGen,coeff
                        !write(*,*) add_to_1,add_to_2

                        jpos  = jPairSpec%offset
                        jflag = -1
                        do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                           jpos = jpos + 1
        
                           ipos  = iPairSpec%offset
                           iflag = -1
                           do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                              ipos = ipos + 1

                              matS(ipos,jpos) = matS(ipos,jpos) + &
                                   coeff * TwoL%r12(&
                                   offset_t  + ti + tj + add_t)%elms( &
                                   offset_uv + vi + vj + add_to_2, &
                                   offset_uv + ui + uj + add_to_1)

                           !   matS(ipos,jpos) = Two%r12(&
                           !        offset_t  + ti + tj + add_t)%elms( &
                           !        offset_uv + vi + vj, &
                           !        offset_uv + ui + uj)
        
                           enddo
        
                        enddo
                      end associate      
                   endif
                enddo

             endif

           end associate

        case('B','b')

           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp2
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp1

           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(.not.Two%isUsed) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Non-existing entry in CCint2_matS!'
                stop
             endif
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in CCint2_matS!'
                stop
             endif
           
             if(ijklorder) then
              
                do ilam=0,2
                   if(info(ilam)%isUsed) then
                      associate(&
                           offset_uv => offset_uv_vec(ilam),&
                           coeff     => info(ilam)%coeff,&
                           add_to_1  => info(ilam)%add_to_1,&
                           add_to_2  => info(ilam)%add_to_2,&
                           TwoL      => Two%ILambda(ilam))

                       ! write(*,*) 'B-even',ilam
                       ! write(*,*) iGen,jGen,coeff
                       ! write(*,*) add_to_1,add_to_2

                        jpos  = jPairSpec%offset
                        jflag = -1
                        do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                           jpos = jpos + 1
        
                           ipos  = iPairSpec%offset
                           iflag = -1
                           do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                              ipos = ipos + 1
                             
                              matS(ipos,jpos) = matS(ipos,jpos) + & 
                                   coeff * TwoL%r12(&
                                   offset_t  + ti + tj + add_t)%elms( &
                                   offset_uv + ui + vj + add_to_1, &
                                   offset_uv + vi + uj + add_to_2)

                            !  matS(ipos,jpos) = Two%r12(&
                            !       offset_t  + ti + tj + add_t)%elms( &
                            !       offset_uv + ui + vj, &
                            !       offset_uv + vi + uj)
        
                           enddo
        
                        enddo
                      end associate
                   endif
                enddo

             else

                do ilam=0,2
                   if(info(ilam)%isUsed) then
                      associate(&
                           offset_uv => offset_uv_vec(ilam),&
                           coeff     => info(ilam)%coeff,&
                           add_to_1  => info(ilam)%add_to_1,&
                           add_to_2  => info(ilam)%add_to_2,&
                           TwoL      => Two%ILambda(ilam))

                       ! write(*,*) 'B-odd',ilam
                       ! write(*,*) iGen,jGen,coeff
                       ! write(*,*) add_to_1,add_to_2
 
                        jpos  = jPairSpec%offset
                        jflag = -1
                        do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                           jpos = jpos + 1
        
                           ipos  = iPairSpec%offset
                           iflag = -1
                           do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                              ipos = ipos + 1
        
                              matS(ipos,jpos) = matS(ipos,jpos) + &
                                   coeff * TwoL%r12(&
                                   offset_t  + ti + tj + add_t)%elms( &
                                   offset_uv + vi + uj + add_to_2, &
                                   offset_uv + ui + vj + add_to_1)
        
                            !  matS(ipos,jpos) = Two%r12(&
                            !       offset_t  + ti + tj + add_t)%elms( &
                            !       offset_uv + vi + uj, &
                            !       offset_uv + ui + vj)
        
                           enddo
        
                        enddo
                      end associate
                   endif
                enddo

             endif

           end associate

        case default

           write(LOUT,'(a)') 'ERROR!!! Incorrect ABtype in CCint2_matS!'
           stop

        end select

      end associate
   enddo
enddo

end subroutine CCint2_matS

subroutine CCint2_vec(add_t,vec,&
     PairSystem,i_vector,i_OrbSystem,j_vector,j_OrbSystem)
implicit none
integer,intent(in) :: add_t
real(prec) :: vec(:)
type(PairSystemData),intent(in) :: PairSystem
real(prec),intent(in) :: i_vector(:),j_vector(:)
type(OrbSystemData),intent(in) :: i_OrbSystem,j_OrbSystem
integer :: prim,i_prim,j_prim
integer :: iexp,jexp,kexp,lexp,iTwo
logical :: ijorder,klorder,ijklorder
integer :: u,v,t,pos,flag
integer :: i,j,i_pos,j_pos
integer :: pos_v
real(prec) :: rtmp

if(add_t<smallest_t) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Additional t parameter in CCint2_vec is incorrect!'
   stop
endif

vec = 0._prec

do j_prim=1,j_OrbSystem%n_prim
   do i_prim=1,i_OrbSystem%n_prim
      do prim=1,PairSystem%n_prim
         associate(&
              PairSpec  => PairSystem%PairSpec(prim),&
              i_OrbSpec => i_OrbSystem%OrbSpec(i_prim),&
              j_OrbSpec => j_OrbSystem%OrbSpec(j_prim))

           iexp = i_OrbSpec%iexp
           jexp = PairSpec%iexp1
           kexp = j_OrbSpec%iexp
           lexp = PairSpec%iexp2
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(.not.Two%isUsed) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Non-existing entry in CCint2_vec!'
                stop
             endif
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in CCint2_vec!'
                stop
             endif

             if(ijklorder) then

                pos  = PairSpec%offset
                flag = -1
                do while(increment_PairSpec(u,v,t,flag,PairSpec))
                   pos  = pos + 1
                   rtmp = 0._prec
            !       associate(elms => Two%r12(offset_t + t + add_t)%elms)

            !         j_pos = j_OrbSpec%offset
            !         do j=0,j_OrbSpec%irange
            !            j_pos = j_pos + 1

            !            pos_v = offset_uv + j + v

            !            i_pos = i_OrbSpec%offset
            !            do i=0,i_OrbSpec%irange
            !               i_pos = i_pos + 1

            !               rtmp = rtmp &
            !                    + i_vector(i_pos)*j_vector(j_pos) &
            !                    * elms(offset_uv + i + u,pos_v)

            !            enddo

            !         enddo

            !       end associate
                   vec(pos) = vec(pos) + rtmp
                enddo

             else

                pos  = PairSpec%offset
                flag = -1
                do while(increment_PairSpec(u,v,t,flag,PairSpec))
                   pos  = pos + 1
                   rtmp = 0._prec
              !     associate(elms => Two%r12(offset_t + t + add_t)%elms)

              !       i_pos = i_OrbSpec%offset
              !       do i=0,i_OrbSpec%irange
              !          i_pos = i_pos + 1

              !          pos_v = offset_uv + i + u

              !          j_pos = j_OrbSpec%offset
              !          do j=0,j_OrbSpec%irange
              !             j_pos = j_pos + 1

              !             rtmp = rtmp &
              !                  + j_vector(j_pos)*i_vector(i_pos) &
              !                  * elms(offset_uv + j + v,pos_v)

              !          enddo

              !       enddo

              !     end associate
                   vec(pos) = vec(pos) + rtmp
                enddo

             endif

           end associate

         end associate
      enddo
   enddo
enddo

end subroutine CCint2_vec

function CCint2_val(only_t,&
     i1_vector,i1_OrbSystem,j1_vector,j1_OrbSystem,&
     i2_vector,i2_OrbSystem,j2_vector,j2_OrbSystem) result(val)
implicit none
real(prec) :: val
integer,intent(in) :: only_t
real(prec),intent(in) :: i1_vector(:),j1_vector(:)
type(OrbSystemData),intent(in) :: i1_OrbSystem,j1_OrbSystem
real(prec),intent(in) :: i2_vector(:),j2_vector(:)
type(OrbSystemData),intent(in) :: i2_OrbSystem,j2_OrbSystem
integer :: i1_prim,j1_prim,i2_prim,j2_prim
integer :: iexp,jexp,kexp,lexp,iTwo
logical :: ijorder,klorder,ijklorder
integer :: i1,j1,i2,j2,i1_pos,j1_pos,i2_pos,j2_pos
integer :: pos_v
real(prec) :: rtmp

if(only_t<smallest_t) then
   write(LOUT,'(a)') 'ERROR!!! &
        &The t parameter in CCint2_val is incorrect!'
   stop
endif

val = 0._prec

do j2_prim=1,j2_OrbSystem%n_prim
   do i2_prim=1,i2_OrbSystem%n_prim
      do j1_prim=1,j1_OrbSystem%n_prim
         do i1_prim=1,i1_OrbSystem%n_prim
            associate(&
                 i1_OrbSpec => i1_OrbSystem%OrbSpec(i1_prim),&
                 j1_OrbSpec => j1_OrbSystem%OrbSpec(j1_prim),&
                 i2_OrbSpec => i2_OrbSystem%OrbSpec(i2_prim),&
                 j2_OrbSpec => j2_OrbSystem%OrbSpec(j2_prim))

              iexp = i1_OrbSpec%iexp
              jexp = i2_OrbSpec%iexp
              kexp = j1_OrbSpec%iexp
              lexp = j2_OrbSpec%iexp
              call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

              associate(Two => TwoInt(iTwo))
                if(.not.Two%isUsed) then
                   write(LOUT,'(a)') 'ERROR!!! &
                        &Non-existing entry in CCint2_val!'
                   stop
                endif
                if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                     kexp/=Two%kexp.or.lexp/=Two%lexp) then
                   write(LOUT,'(a)') 'ERROR!!! &
                        &Wrong entry chosen in CCint2_val!'
                   stop
                endif

               ! associate(elms => Two%r12(offset_t + only_t)%elms)

               !   if(ijklorder) then

               !      j2_pos = j2_OrbSpec%offset
               !      do j2=0,j2_OrbSpec%irange
               !         j2_pos = j2_pos + 1
               !         j1_pos = j1_OrbSpec%offset
               !         do j1=0,j1_OrbSpec%irange
               !            j1_pos = j1_pos + 1

               !            pos_v = offset_uv + j1 + j2
               !            rtmp  = 0._prec

               !            i2_pos = i2_OrbSpec%offset
               !            do i2=0,i2_OrbSpec%irange
               !               i2_pos = i2_pos + 1
               !               i1_pos = i1_OrbSpec%offset
               !               do i1=0,i1_OrbSpec%irange
               !                  i1_pos = i1_pos + 1

               !                  rtmp = rtmp &
               !                       + i1_vector(i1_pos)*i2_vector(i2_pos) &
               !                       * elms(offset_uv + i1 + i2,pos_v)

               !               enddo
               !            enddo

               !            val = val + rtmp &
               !                 * j1_vector(j1_pos)*j2_vector(j2_pos)

               !         enddo
               !      enddo

               !   else

               !      i2_pos = i2_OrbSpec%offset
               !      do i2=0,i2_OrbSpec%irange
               !         i2_pos = i2_pos + 1
               !         i1_pos = i1_OrbSpec%offset
               !         do i1=0,i1_OrbSpec%irange
               !            i1_pos = i1_pos + 1

               !            pos_v = offset_uv + i1 + i2
               !            rtmp  = 0._prec

               !            j2_pos = j2_OrbSpec%offset
               !            do j2=0,j2_OrbSpec%irange
               !               j2_pos = j2_pos + 1
               !               j1_pos = j1_OrbSpec%offset
               !               do j1=0,j1_OrbSpec%irange
               !                  j1_pos = j1_pos + 1

               !                  rtmp = rtmp &
               !                       + j1_vector(j1_pos)*j2_vector(j2_pos) &
               !                       * elms(offset_uv + j1 + j2,pos_v)

               !               enddo
               !            enddo

               !            val = val + rtmp &
               !                 * i1_vector(i1_pos)*i2_vector(i2_pos)

               !         enddo
               !      enddo

               !   endif

               ! end associate

              end associate

            end associate
         enddo
      enddo
   enddo
enddo

end function CCint2_val

subroutine order_exp(iTwo,ijorder,klorder,ijklorder,ijmin,ijmax,klmin,klmax)
use misc, only : swap
implicit none
integer,intent(out) :: iTwo
logical,intent(out) :: ijorder,klorder,ijklorder
integer,intent(inout) :: ijmin,ijmax,klmin,klmax
integer :: ijexp,klexp

ijorder = (ijmin<=ijmax)
if(.not.ijorder) call swap(ijmin,ijmax)

klorder = (klmin<=klmax)
if(.not.klorder) call swap(klmin,klmax)

ijexp = ijmin + ijmax*(ijmax-1)/2
klexp = klmin + klmax*(klmax-1)/2

ijklorder = (ijexp<=klexp)
if(.not.ijklorder) then
   call swap(ijmin,klmin)
   call swap(ijmax,klmax)
   call swap(ijorder,klorder)
   call swap(ijexp,klexp)
endif

iTwo = ijexp + klexp*(klexp-1)/2

end subroutine order_exp

subroutine init1_TwoInt(Nexp,exponents)
implicit none
integer,intent(in) :: Nexp
real(prec),intent(in) :: exponents(:)
integer :: iexp,jexp,kexp,lexp,iTwo
integer :: ilam

iTwo = Nexp*(Nexp+1)/2
allocate(TwoInt(iTwo*(iTwo+1)/2))

iTwo = 0
do lexp=1,Nexp
   do kexp=1,lexp
      do jexp=1,lexp
         do iexp=1,merge(kexp,jexp,jexp==lexp)
            iTwo = iTwo + 1
            associate(Two => TwoInt(iTwo))

              Two%isUsed    = .false.
              Two%iexp      = iexp
              Two%jexp      = jexp
              Two%kexp      = kexp
              Two%lexp      = lexp
              Two%same_exp  = (Two%iexp==Two%kexp).and.(Two%jexp==Two%lexp)
              Two%ijalphaPL = exponents(Two%iexp) + exponents(Two%jexp)
              Two%ijalphaMI = exponents(Two%iexp) - exponents(Two%jexp)
              Two%klalphaPL = exponents(Two%kexp) + exponents(Two%lexp)
              Two%klalphaMI = exponents(Two%kexp) - exponents(Two%lexp)
              !Two%range_t   = smallest_t - 1
              ! HERE
              do ilam=0,2
                 associate(TwoL => Two%ILambda(ilam))
                   TwoL%range_t   = smallest_t - 1
                 end associate
              enddo

            end associate
         enddo
      enddo
   enddo
enddo
if(iTwo/=size(TwoInt)) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Number of entries in CC TwoInt has not been predicted correctly!'
   stop
endif

end subroutine init1_TwoInt

subroutine init2_TwoInt(PairReduced)
implicit none
type(PairReducedData),intent(in) :: PairReduced(:)
integer :: i_prim,j_prim
integer :: iexp,jexp,kexp,lexp,iTwo
logical :: ijorder,klorder,ijklorder
integer :: range_ti,range_tj,range_t
integer :: igen,jgen,jGenMax
integer :: ilam,numlam
integer,allocatable :: jGenNum(:),lam(:)

call mem_alloc(jGenNum,2)
jGenNum = 0

do j_prim=1,size(PairReduced)
   do i_prim=1,j_prim
      associate(&
           iPairSpec => PairReduced(i_prim), &
           jPairSpec => PairReduced(j_prim))
        if(iPairSpec%isUsed.and.jPairSpec%isUsed) then

          ! range_ti = maxt_PairReduced(iPairSpec,smallest_t)
          ! range_tj = maxt_PairReduced(jPairSpec,smallest_t)
          ! if((range_ti<smallest_t).or.(range_tj<smallest_t)) then
          !    write(LOUT,'(a)') 'ERROR!!! &
          !         &It is impossible to determine max t range in init2_TwoInt!'
          !    stop
          ! endif

          ! range_t = range_ti + range_tj

           ! A part
           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp1
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp2

           !write(*,*) 'u:',iexp,kexp,jexp,lexp
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in init2_TwoInt!'
                stop
             endif

            ! loop over generators
            do igen=1,iPairSpec%n_gen
               associate(iPairSpecG => iPairSpec%PairReduced_G(igen),&
                         iPairGen => iPairSpec%PairReduced_G(igen)%gen_type)
!                 write(*,*) 'o:',iexp,kexp,jexp,lexp
!                 write(*,*) ''
                 call check_gen_pairs(iPairSpecG,jPairSpec,jGenMax,jGenNum)
                 if(jGenMax.eq.0) cycle
                 do jgen=1,jGenMax
                    associate(jPairSpecG => jPairSpec%PairReduced_G(jGenNum(jgen)),&
                              jPairGen => jPairSpec%PairReduced_G(jGenNum(jgen))%gen_type)
                             
                    ! write(*,*) 'iGen:',iPairSpecG%gen_type,'jGen',jPairSpecG%gen_type
                      range_ti = maxt_PairReduced(iPairSpecG,smallest_t)
                      range_tj = maxt_PairReduced(jPairSpecG,smallest_t)
                      if((range_ti<smallest_t).or.(range_tj<smallest_t)) then
                         write(LOUT,'(a)') 'ERROR!!! &
                         &It is impossible to determine max t range in init2_TwoInt!'
                         stop
                      endif

                      range_t = range_ti + range_tj
                      Two%isUsed  = .true.

                      ! ATTENTION!!!  
                      ! ordering of generators
                      ! according to possible_generators 
                      if(iPairGen==jPairGen) then
                         select case(iPairGen)
                         case(3,4,6)

                             numlam = 2
                             call mem_alloc(lam,numlam)
                             lam(1) = 0
                             lam(2) = 2
                         !   Two%isUsed  = .true.
                         !   Two%ILambda(0)%range_t = &
                         !       max(Two%ILambda(0)%range_t,range_t)
                         !   Two%ILambda(2)%range_t = &
                         !       max(Two%ILambda(2)%range_t,range_t)

                        !write(*,*) 'iGen:',iPairGen,'jGen',jPairGen
                        !write(*,*) 'Lambda: (0, 2)'
                         case(1,2,5)

                             numlam = 1
                             call mem_alloc(lam,numlam)
                             lam(1) = 0
                         !   Two%isUsed  = .true.
                         !   Two%ILambda(0)%range_t = &
                         !       max(Two%ILambda(0)%range_t,range_t)
                         !write(*,*) 'iGen:',iPairGen,'jGen',jPairGen
                         !write(*,*) 'Lambda: (0)'
                         case default
                            write(LOUT,'(a)') 'Incorrect generator number &
                                             & in init2_TwoInt'
                            stop
                         end select
                      else

                         numlam = 1
                         call mem_alloc(lam,numlam)
                         lam(1) = 1
                         !Two%isUsed  = .true.
                         !Two%ILambda(1)%range_t = &
                         !    max(Two%ILambda(1)%range_t,range_t)
                        !write(*,*) 'iGen:',iPairGen,'jGen',jPairGen
                        !write(*,*) 'Lambda: (1)'
                      endif

! hapka: old_drake
!                range_t = range_ti + range_tj
!                Two%isUsed  = .true.
!                Two%range_t = max(Two%range_t,range_t)

                      do ilam=1,numlam
                         associate(TwoL => Two%ILambda(lam(ilam)))
                           !Two%isUsed  = .true.
                           TwoL%range_t = &
                               max(TwoL%range_t,range_t)
                         end associate     
                      enddo
                      !write(*,*) 'iGen:',iPairGen,'jGen',jPairGen
                      !write(*,*) 'Lambda ',lam

                      call mem_dealloc(lam)

                    end associate
                 enddo
               end associate 
            enddo
           end associate

         ! B part
           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp2
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp1

 !          write(*,*) 'u:',iexp,kexp,jexp,lexp
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

          associate(Two => TwoInt(iTwo))
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in init2_TwoInt!'
                stop
             endif

            ! loop over generators
            do igen=1,iPairSpec%n_gen
               associate(iPairSpecG => iPairSpec%PairReduced_G(igen),&
                         iPairGen => iPairSpec%PairReduced_G(igen)%gen_type)
!                 write(*,*) 'o:',iexp,kexp,jexp,lexp
!                 write(*,*) ''
                 call check_gen_pairs(iPairSpecG,jPairSpec,jGenMax,jGenNum)
                 if(jGenMax.eq.0) cycle
                 do jgen=1,jGenMax
                    associate(jPairSpecG => jPairSpec%PairReduced_G(jGenNum(jgen)),&
                              jPairGen => jPairSpec%PairReduced_G(jGenNum(jgen))%gen_type)
                             
                    ! write(*,*) 'iGen:',iPairSpecG%gen_type,'jGen',jPairSpecG%gen_type
                      range_ti = maxt_PairReduced(iPairSpecG,smallest_t)
                      range_tj = maxt_PairReduced(jPairSpecG,smallest_t)
                      if((range_ti<smallest_t).or.(range_tj<smallest_t)) then
                         write(LOUT,'(a)') 'ERROR!!! &
                         &It is impossible to determine max t range in init2_TwoInt!'
                         stop
                      endif

                      range_t = range_ti + range_tj

                      ! ATTENTION!!!  
                      ! ordering of generators
                      ! according to possible_generators 
                      if(iPairGen==jPairGen) then
                        select case(iPairGen)
                        case(3,4,6)
                           numlam = 2
                           call mem_alloc(lam,numlam)
                           lam(1) = 0
                           lam(2) = 2

                        ! write(*,*) 'iGen:',iPairGen,'jGen',jPairGen
                        ! write(*,*) 'Lambda: (0, 2)'
                         case(1)
                            numlam = 1
                            call mem_alloc(lam,numlam)
                            lam(1) = 0

                         case(2)
                            numlam = 1
                            call mem_alloc(lam,numlam)
                            lam(1) = 1

                         case(5)
                            numlam = 1
                            call mem_alloc(lam,numlam)
                            lam(1) = 2

                         end select
                      else 
                         numlam = 1
                         call mem_alloc(lam,numlam)
                         lam(1) = 1

                      endif
! hapka: old_drake
!                range_t = range_ti + range_tj
!                Two%isUsed  = .true.
!                Two%range_t = max(Two%range_t,range_t)

                      do ilam=1,numlam
                         associate(TwoL => Two%ILambda(lam(ilam)))
                           !Two%isUsed  = .true.
                           TwoL%range_t = &
                               max(TwoL%range_t,range_t)
                         end associate     
                      enddo
                      !write(*,*) 'iGen:',iPairGen,'jGen',jPairGen
                      !write(*,*) 'Lambda ',lam

                      call mem_dealloc(lam)

                    end associate
                 enddo
               end associate 
            enddo
          end associate

        endif
      end associate
   enddo
enddo

call mem_dealloc(jGenNum)

end subroutine init2_TwoInt

subroutine init3_TwoInt(PairReduced)
implicit none
type(PairReducedData),intent(in) :: PairReduced(:)
integer :: i_prim,j_prim
integer :: iexp,jexp,kexp,lexp,iTwo
logical :: ijorder,klorder,ijklorder
integer :: ti,tj,pos_ti,pos_tj,range_ti,range_tj
integer :: t,pos_t,max_t
integer :: add_to_uv(2)
integer,allocatable :: range_uvi(:,:),range_uvj(:,:)
logical,allocatable :: used_uvi(:),used_uvj(:)
integer :: i
integer :: ilam
integer :: numlam
integer :: igen,jgen,jGenMax
integer,allocatable :: jGenNum(:),lam(:)

call mem_alloc(jGenNum,2)
jGenNum = 0

max_t = smallest_t - 1

do iTwo=1,size(TwoInt)
   associate(Two => TwoInt(iTwo))
     if(Two%isUsed) then

        do ilam=0,2
           associate(TwoL => Two%ILambda(ilam))
             max_t = max(max_t,TwoL%range_t)

             call mem_alloc(TwoL%range_uv,2,TwoL%range_t-smallest_t+1)
             TwoL%range_uv = smallest_uv - 1
           end associate
        enddo
     endif
   end associate
enddo

! hapka: old_drake
!do iTwo=1,size(TwoInt)
!   associate(Two => TwoInt(iTwo))
!     if(Two%isUsed) then
!
!        max_t = max(max_t,Two%range_t)
!
!        call mem_alloc(Two%range_uv,2,Two%range_t-smallest_t+1)
!        Two%range_uv = smallest_uv - 1
!
!     endif
!   end associate
!enddo

max_t = max_t - smallest_t + 1

call mem_alloc(range_uvi,2,max_t)
call mem_alloc(range_uvj,2,max_t)
call mem_alloc(used_uvi,max_t)
call mem_alloc(used_uvj,max_t)

do j_prim=1,size(PairReduced)
   do i_prim=1,j_prim
      associate(&
           iPairSpec => PairReduced(i_prim), &
           jPairSpec => PairReduced(j_prim))
        if(iPairSpec%isUsed.and.jPairSpec%isUsed) then
! hapka: old_drake
!           range_ti = maxt_PairReduced(iPairSpec,smallest_t)
!           range_tj = maxt_PairReduced(jPairSpec,smallest_t)
!           if((range_ti<smallest_t).or.(range_tj<smallest_t)) then
!              write(LOUT,'(a)') 'ERROR!!! &
!                   &It is impossible to determine max t range in init3_TwoInt!'
!              stop
!           endif
!
!           range_uvi = smallest_uv - 1
!           used_uvi  = .false.
!           do ti=0,range_ti
!              pos_ti = offset_t + ti
!              range_uvi(:,pos_ti) = maxuv_PairReduced(ti,iPairSpec,smallest_uv)
!              used_uvi(pos_ti)    = all(range_uvi(:,pos_ti)>=smallest_uv)
!           enddo
!
!           range_uvj = smallest_uv - 1
!           used_uvj  = .false.
!           do tj=0,range_tj
!              pos_tj = offset_t + tj
!              range_uvj(:,pos_tj) = maxuv_PairReduced(tj,jPairSpec,smallest_uv)
!              used_uvj(pos_tj)    = all(range_uvj(:,pos_tj)>=smallest_uv)
!           enddo
!
           ! A part  
           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp1
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp2
!           write(*,*) 'u:',iexp,kexp,jexp,lexp
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(.not.Two%isUsed) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Non-existing entry in init3_TwoInt!'
                stop
             endif
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in init3_TwoInt!'
                stop
             endif

!            write(*,*) 'o:',iexp,kexp,jexp,lexp,ijklorder
!            write(*,*) ''
            ! loop over generators
             do igen=1,iPairSpec%n_gen
                associate(iPairSpecG => iPairSpec%PairReduced_G(igen),&
                          iPairGen => iPairSpec%PairReduced_G(igen)%gen_type)
                  call check_gen_pairs(iPairSpecG,jPairSpec,jGenMax,jGenNum)
                  if(jGenMax.eq.0) cycle
                  do jgen=1,jGenMax
                     associate(jPairSpecG => jPairSpec%PairReduced_G(jGenNum(jgen)),&
                               jPairGen => jPairSpec%PairReduced_G(jGenNum(jgen))%gen_type)
                             
                       range_ti = maxt_PairReduced(iPairSpecG,smallest_t)
                       range_tj = maxt_PairReduced(jPairSpecG,smallest_t)
                       if((range_ti<smallest_t).or.(range_tj<smallest_t)) then
                          write(LOUT,'(a)') 'ERROR!!! &
                          &It is impossible to determine max t range in init3_TwoInt!'
                          stop
                       endif

                       range_uvi = smallest_uv - 1
                       used_uvi  = .false.
                       do ti=0,range_ti
                          pos_ti = offset_t + ti
                          range_uvi(:,pos_ti) = maxuv_PairReduced(ti,iPairSpecG,smallest_uv)
                          used_uvi(pos_ti)    = all(range_uvi(:,pos_ti)>=smallest_uv)
                       enddo

                       range_uvj = smallest_uv - 1
                       used_uvj  = .false.
                       do tj=0,range_tj
                          pos_tj = offset_t + tj
                          range_uvj(:,pos_tj) = maxuv_PairReduced(tj,jPairSpecG,smallest_uv)
                          used_uvj(pos_tj)    = all(range_uvj(:,pos_tj)>=smallest_uv)
                       enddo

                     ! ATTENTION!!!  
                     ! ordering of generators
                     ! according to possible_generators 
                       add_to_uv = 0
                       if(iPairGen==jPairGen) then
                          ! even generators
                          select case(iPairGen)
                          case(3,4,6)
                            numlam = 2
                            call mem_alloc(lam,numlam)
                            lam(1) = 0
                            lam(2) = 2
                            add_to_uv(1) = 2
                            add_to_uv(2) = 2

                          case(1,2,5) 
                            numlam = 1
                            call mem_alloc(lam,numlam)
                            lam(1) = 0

                            select case(iPairGen)
                            case(1)
                               add_to_uv = 0
                            case(2)
                               add_to_uv(1) = 0
                               add_to_uv(2) = 2
                            case(5)
                               add_to_uv(1) = 0
                               add_to_uv(2) = 4
                            case default
                               write(LOUT,'()') 'Incorrect even pair &
                                   & generator number in init3_TwoInt'
                               stop
                            end select  

                          case default
                             write(LOUT,'(a)') 'Incorrect even pair & 
                                  &  generator number in init3_TwoInt &
                                  & (A part) '
                             stop
                          end select

                       else ! odd generators
                          numlam = 1
                          call mem_alloc(lam,numlam)
                          lam(1) = 1

                          select case(iPairGen)
                          case(1,6)
                             add_to_uv = 1
                          case(4,5)
                             add_to_uv(1) = 3 !1 - here is it ok?
                             add_to_uv(2) = 3
                          case default
                             write(LOUT,'(a)') 'Incorrect odd pair &
                                & generator number in init3_TwoInt '
                          end select
                       endif  

                      !  check A part
                      ! write(*,*) 'iGen:',iPairSpecG%gen_type,'jGen',jPairSpecG%gen_type
                      ! write(*,*) 'Lambda   ', lam
                      ! write(*,*) 'add_to_uv', add_to_uv

                       if(ijklorder) then

                       do tj=0,range_tj
                          pos_tj = offset_t + tj
                          do ti=0,range_ti
                             pos_ti = offset_t + ti
                             if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
                                pos_t = offset_t + ti + tj

                                do ilam=1,numlam
                                   associate(TwoL => Two%ILambda(lam(ilam)))
 
                                     TwoL%range_uv(1,pos_t) = max(TwoL%range_uv(1,pos_t),&
                                          range_uvi(1,pos_ti) + range_uvj(1,pos_tj) + &
                                          add_to_uv(1))
                                     TwoL%range_uv(2,pos_t) = max(TwoL%range_uv(2,pos_t),&
                                          range_uvi(2,pos_ti) + range_uvj(2,pos_tj) + &
                                          add_to_uv(2))

                                   end associate
                                enddo     
 
                             endif
                          enddo
                       enddo

                       else

                       do tj=0,range_tj
                          pos_tj = offset_t + tj
                          do ti=0,range_ti
                             pos_ti = offset_t + ti
                             if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
                                pos_t = offset_t + ti + tj

                                do ilam=1,numlam
                                   associate(TwoL => Two%ILambda(lam(ilam)))

                                   TwoL%range_uv(1,pos_t) = max(TwoL%range_uv(1,pos_t),&
                                     range_uvi(2,pos_ti) + range_uvj(2,pos_tj) + &
                                     add_to_uv(2))
                                   TwoL%range_uv(2,pos_t) = max(TwoL%range_uv(2,pos_t),&
                                     range_uvi(1,pos_ti) + range_uvj(1,pos_tj) + &
                                     add_to_uv(1))

                                   end associate  
                                enddo
 
                             endif
                          enddo
                       enddo

                   endif

                        call mem_dealloc(lam)
                     end associate
                  enddo
                end associate
             enddo
           end associate

! hapka: old_drake (A part)
!             if(ijklorder) then
!
!                do tj=0,range_tj
!                   pos_tj = offset_t + tj
!                   do ti=0,range_ti
!                      pos_ti = offset_t + ti
!                      if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
!                         pos_t = offset_t + ti + tj
!
!                         Two%range_uv(1,pos_t) = max(Two%range_uv(1,pos_t),&
!                              range_uvi(1,pos_ti) + range_uvj(1,pos_tj))
!                         Two%range_uv(2,pos_t) = max(Two%range_uv(2,pos_t),&
!                              range_uvi(2,pos_ti) + range_uvj(2,pos_tj))
!
!                      endif
!                   enddo
!                enddo
!
!             else
!
!                do tj=0,range_tj
!                   pos_tj = offset_t + tj
!                   do ti=0,range_ti
!                      pos_ti = offset_t + ti
!                      if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
!                         pos_t = offset_t + ti + tj
!
!                         Two%range_uv(1,pos_t) = max(Two%range_uv(1,pos_t),&
!                              range_uvi(2,pos_ti) + range_uvj(2,pos_tj))
!                         Two%range_uv(2,pos_t) = max(Two%range_uv(2,pos_t),&
!                              range_uvi(1,pos_ti) + range_uvj(1,pos_tj))
!
!                      endif
!                   enddo
!                enddo
!
!             endif
!          end associate
!
!          B part
           iexp = iPairSpec%iexp1
           jexp = jPairSpec%iexp2
           kexp = iPairSpec%iexp2
           lexp = jPairSpec%iexp1
          
           call order_exp(iTwo,ijorder,klorder,ijklorder,iexp,jexp,kexp,lexp)

           associate(Two => TwoInt(iTwo))
             if(.not.Two%isUsed) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Non-existing entry in init3_TwoInt!'
                stop
             endif
             if(iexp/=Two%iexp.or.jexp/=Two%jexp.or.&
                  kexp/=Two%kexp.or.lexp/=Two%lexp) then
                write(LOUT,'(a)') 'ERROR!!! &
                     &Wrong entry chosen in init3_TwoInt!'
                stop
             endif

           ! write(*,*) 'o:',iexp,kexp,jexp,lexp,ijklorder
           ! write(*,*) ''
            ! loop over generators
             do igen=1,iPairSpec%n_gen
                associate(iPairSpecG => iPairSpec%PairReduced_G(igen),&
                          iPairGen => iPairSpec%PairReduced_G(igen)%gen_type)
                  call check_gen_pairs(iPairSpecG,jPairSpec,jGenMax,jGenNum)
                  if(jGenMax.eq.0) cycle
                  do jgen=1,jGenMax
                     associate(jPairSpecG => jPairSpec%PairReduced_G(jGenNum(jgen)),&
                               jPairGen => jPairSpec%PairReduced_G(jGenNum(jgen))%gen_type)
                             
                       range_ti = maxt_PairReduced(iPairSpecG,smallest_t)
                       range_tj = maxt_PairReduced(jPairSpecG,smallest_t)
                       if((range_ti<smallest_t).or.(range_tj<smallest_t)) then
                          write(LOUT,'(a)') 'ERROR!!! &
                          &It is impossible to determine max t range in init3_TwoInt!'
                          stop
                       endif

                       range_uvi = smallest_uv - 1
                       used_uvi  = .false.
                       do ti=0,range_ti
                          pos_ti = offset_t + ti
                          range_uvi(:,pos_ti) = maxuv_PairReduced(ti,iPairSpecG,smallest_uv)
                          used_uvi(pos_ti)    = all(range_uvi(:,pos_ti)>=smallest_uv)
                       enddo

                       range_uvj = smallest_uv - 1
                       used_uvj  = .false.
                       do tj=0,range_tj
                          pos_tj = offset_t + tj
                          range_uvj(:,pos_tj) = maxuv_PairReduced(tj,jPairSpecG,smallest_uv)
                          used_uvj(pos_tj)    = all(range_uvj(:,pos_tj)>=smallest_uv)
                       enddo

                     ! ATTENTION!!!  
                     ! ordering of generators
                     ! according to possible_generators 
                       add_to_uv = 0
                       if(iPairGen==jPairGen) then
                          ! even generators
                          select case(iPairGen)
                          case(3,4,6)
                             numlam = 2
                             call mem_alloc(lam,numlam)
                             lam(1) = 0
                             lam(2) = 2
                             add_to_uv(1) = 2
                             add_to_uv(2) = 2
                          case(1,2,5)
                             numlam = 1
                             call mem_alloc(lam,numlam)

                             select case(iPairGen)
                             case(1)
                                lam(1) = 0
                                add_to_uv = 0
                             case(2)
                                lam(1) = 1
                                add_to_uv(1) = 1
                                add_to_uv(2) = 1
                             case(5)
                                lam(1) = 2
                                add_to_uv(1) = 2
                                add_to_uv(2) = 2
                             case default
                                write(LOUT,'()') 'Incorrect even pair &
                                    & genererator number in init3_TwoInt &
                                    & (B part) '
                                stop
                             end select
                           case default
                              write(LOUT,'()') 'Incorrect even pair &
                                  & genererator number in init3_TwoInt &
                                  & (B part) '
                              stop
                          end select
                       else ! odd generators
                          numlam = 1
                          call mem_alloc(lam,numlam)
                          lam(1) = 1
                          select case(iPairGen)
                          case(1,6)
                             add_to_uv = 1
                          case(4,5)
                             add_to_uv(1) = 3
                             add_to_uv(2) = 3 !1
                          case default
                             write(LOUT,'(a)') 'Incorrect odd pair &
                                & generator number in init3_TwoInt &
                                & (B part)'
                          end select

                       endif
                       ! check B part
                       !write(*,*) 'iGen:',iPairSpecG%gen_type,'jGen',jPairSpecG%gen_type
                       !write(*,*) 'Lambda   ', lam
                       !write(*,*) 'add_to_uv', add_to_uv

                       if(ijklorder) then
          
                       do tj=0,range_tj
                          pos_tj = offset_t + tj
                          do ti=0,range_ti
                             pos_ti = offset_t + ti
                             if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
                                pos_t = offset_t + ti + tj

                                do ilam=1,numlam
                                   associate(TwoL => Two%ILambda(lam(ilam)))

                                     TwoL%range_uv(1,pos_t) = max(TwoL%range_uv(1,pos_t),&
                                         range_uvi(1,pos_ti) + range_uvj(2,pos_tj) + &
                                         add_to_uv(1))
                                     TwoL%range_uv(2,pos_t) = max(TwoL%range_uv(2,pos_t),&
                                         range_uvi(2,pos_ti) + range_uvj(1,pos_tj) + &
                                         add_to_uv(2))
                                   end associate
                                enddo
          
                             endif
                          enddo
                       enddo
          
                       else
          
                       do tj=0,range_tj
                          pos_tj = offset_t + tj
                          do ti=0,range_ti
                             pos_ti = offset_t + ti
                             if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
                                pos_t = offset_t + ti + tj
                                
                                do ilam=1,numlam
                                   associate(TwoL => Two%ILambda(lam(ilam)))
          
                                     TwoL%range_uv(1,pos_t) = max(TwoL%range_uv(1,pos_t),&
                                          range_uvi(2,pos_ti) + range_uvj(1,pos_tj) + &
                                          add_to_uv(2))
                                     TwoL%range_uv(2,pos_t) = max(TwoL%range_uv(2,pos_t),&
                                          range_uvi(1,pos_ti) + range_uvj(2,pos_tj) + & 
                                          add_to_uv(1))
                                   end associate
                                enddo
          
                             endif
                          enddo
                       enddo
          
                       endif

! hapka: old_drake (B part)
!             if(ijklorder) then
!
!                do tj=0,range_tj
!                   pos_tj = offset_t + tj
!                   do ti=0,range_ti
!                      pos_ti = offset_t + ti
!                      if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
!                         pos_t = offset_t + ti + tj
!
!                         Two%range_uv(1,pos_t) = max(Two%range_uv(1,pos_t),&
!                              range_uvi(1,pos_ti) + range_uvj(2,pos_tj))
!                         Two%range_uv(2,pos_t) = max(Two%range_uv(2,pos_t),&
!                              range_uvi(2,pos_ti) + range_uvj(1,pos_tj))
!
!                      endif
!                   enddo
!                enddo
!
!             else
!
!                do tj=0,range_tj
!                   pos_tj = offset_t + tj
!                   do ti=0,range_ti
!                      pos_ti = offset_t + ti
!                      if(used_uvi(pos_ti).and.used_uvj(pos_tj)) then
!                         pos_t = offset_t + ti + tj
!
!                         Two%range_uv(1,pos_t) = max(Two%range_uv(1,pos_t),&
!                              range_uvi(2,pos_ti) + range_uvj(1,pos_tj))
!                         Two%range_uv(2,pos_t) = max(Two%range_uv(2,pos_t),&
!                              range_uvi(1,pos_ti) + range_uvj(2,pos_tj))
!
!                      endif
!                   enddo
!                enddo
!
!             endif
!
                        call mem_dealloc(lam)
                     end associate
                  enddo
                end associate
             enddo
           end associate

        endif
      end associate
   enddo
enddo

call mem_dealloc(used_uvj)
call mem_dealloc(used_uvi)
call mem_dealloc(range_uvj)
call mem_dealloc(range_uvi)

call mem_dealloc(jGenNum)

! adapt uv-range to Drake (11.33) 
do iTwo=1,size(TwoInt)
   associate(Two => TwoInt(iTwo))
     if(Two%isUsed) then

        do ilam=0,2
           associate(TwoL => Two%ILambda(ilam))

             do t=smallest_t,TwoL%range_t-1
                pos_t = offset_t + t
     
                do i=1,min(2,TwoL%range_t-t)
                   TwoL%range_uv(:,pos_t) = &
                        max(TwoL%range_uv(:,pos_t),TwoL%range_uv(:,pos_t+i))
                enddo
     
             enddo

           end associate
        enddo

     endif
   end associate
enddo

end subroutine init3_TwoInt

subroutine create_Twoint
implicit none
integer :: iTwo
integer :: t,pos_t,u,v
integer :: ilam,maxt
integer :: max_lam
integer :: ctr

 do iTwo=1,size(TwoInt)
    associate(Two => TwoInt(iTwo))
      if(Two%isUsed) then
         ! allocate r12%elms
         do ilam=0,2 
            associate(TwoL => Two%ILambda(ilam))
              allocate(TwoL%r12(TwoL%range_t-smallest_t+1))

              do t=smallest_t,TwoL%range_t
                 pos_t = offset_t + t
                 associate(&
                      r12      => TwoL%r12(pos_t), &
                      range_uv => TwoL%range_uv(:,pos_t))

                      call mem_alloc(r12%elms,&
                           range_uv(1)-smallest_uv_vec(ilam)+1,&
                           range_uv(2)-smallest_uv_vec(ilam)+1)

                 end associate
              enddo

            end associate
         enddo 

!        write(LOUT,'(4i3)') Two%iexp,Two%jexp,Two%kexp,Two%lexp
         maxt = maxval(Two%ILambda(:)%range_t)
         ! loop over t values 
         do t=smallest_t,maxt
            pos_t = offset_t + t 
            max_lam = 0
            do ilam=0,2
               if(pos_t.gt.(Two%ILambda(ilam)%range_t+offset_t)) cycle
               max_lam = max(ilam,max_lam)
            enddo

            call create_tval(Two,t,pos_t,max_lam)
            call create_Izero(Two,max_lam)

! check I0 manually
!            write(LOUT,'(5x,2f15.8)') Two%ijalphaPL,Two%klalphaPL
!            do ctr=1,size(Two%t_val)
!                   associate(&
!                     tval     => Two%t_val(ctr),&
!                     two_t    => Two%t_val(ctr)%two_t,&
!                     range_uv => Two%t_val(ctr)%level_range_uv(1,:))
!                     do v=smallest_uv,range_uv(2)
!                        do u=smallest_uv,range_uv(1)
!                           write(LOUT,'(10x,3i3,a,es30.22)') two_t,u,v,' : ',&
!                                tval%elms(offset_uv + u,offset_uv + v)
!                        enddo
!                     enddo
!
!                    end associate
!            enddo
            if(max_lam/=0) then
               call create_Ilambda(Two,max_lam)
            endif  

            call rewrite_Ilambda(Two,t,pos_t,max_lam)
            call free_tval(Two,max_lam)
         enddo

      endif
    end associate
 enddo

end subroutine create_Twoint

subroutine create_tval(Two,t,pos_t,max_lam)
implicit none
type(TwoIntData) :: Two
integer,intent(in) :: t,pos_t
integer,intent(in) :: max_lam
integer,allocatable :: range_tmp(:,:)
integer :: range_uv_alloc(2)
integer :: ilam,ival,ilev,ilev2


allocate(Two%t_val(max_lam+1))
allocate(range_tmp(2,3))

     do ival=1,3
        if(pos_t.gt.(Two%ILambda(ival-1)%range_t+offset_t)) then 
           range_tmp(:,ival) = smallest_uv - 1
        else
           range_tmp(:,ival) = Two%ILambda(ival-1)%range_uv(:,pos_t)
        endif

     enddo

!get level_range_uv
do ilam=max_lam,0,-1

   call mem_alloc(Two%t_val(ilam+1)%level_range_uv,max_lam+1,2)
   associate(tval => Two%t_val(ilam+1))

! not general - should be done 
! differently for higher angular momenta
     select case(ilam)
     case(0)
        do ilev=0,max_lam
           if(ilev==max_lam) then
              ! for t, level=last: uv(I0)
              tval%level_range_uv(max_lam+1,:) = range_tmp(:,1)
           else
              ! for t, other levels: max_uv(I0,I1+1,I2)
              tval%level_range_uv(ilev+1,:) = & 
                   max(range_tmp(:,1),range_tmp(:,2)+1,range_tmp(:,3))
           endif
        enddo
        !write(*,*) 'Lam=0, level=0',tval%level_range_uv(1,:)
        !write(*,*) 'Lam=0, level=1',tval%level_range_uv(2,:)
        !write(*,*) 'Lam=0, level=2',tval%level_range_uv(3,:)
        !write(*,*) ''
          
     case(1)
        do ilev=0,max_lam
           if(ilev==max_lam) then
              ! for t+2, level=last: uv(I1)
              tval%level_range_uv(max_lam+1,:) = range_tmp(:,2)
           elseif(ilev==0) then
              ! for t+2, other levels: max_uv(I1-1,I2)
              tval%level_range_uv(1,:) = &
                   max(range_tmp(:,2)-1,range_tmp(:,3))
           else
              ! for t+2, other levels: max_uv(I1-1,I2)
              tval%level_range_uv(2,:) = &
                   range_tmp(:,2)-1
           endif
        enddo
        !write(*,*) 'Lam=1, level=0',tval%level_range_uv(1,:)
        !write(*,*) 'Lam=1, level=1',tval%level_range_uv(2,:)
        !write(*,*) 'Lam=1, level=2',tval%level_range_uv(3,:)
        !write(*,*) ''

     case(2)
        do ilev=0,max_lam
           if(ilev==max_lam) then
              ! for t+4, level=last: max_uv(I2)
              tval%level_range_uv(3,:) = range_tmp(:,3)
           else
              ! for t+4, level=0: max_uv(I2-2)
              tval%level_range_uv(ilev+1,:) = &
                   range_tmp(:,3)-2+ilev
           endif

        enddo
        !write(*,*) 'Lam=2, level=0',tval%level_range_uv(1,:)
        !write(*,*) 'Lam=2, level=1',tval%level_range_uv(2,:)
        !write(*,*) 'Lam=2, level=2',tval%level_range_uv(3,:)
        !write(*,*) ''

     case default
        write(LOUT,'(a)') 'ERROR!!! Wrong max lambda in &
                       & create_TwoInt (CCint2)!'
        stop
     end select

   end associate
enddo

! allocate tval's
do ilam=0,max_lam
   associate(tval => Two%t_val(ilam+1))

     tval%shift = 0 ! not necessary?
     tval%two_t = t + 2*ilam
     tval%lambda = 0
     select case(ilam)
     case(0)
        ! allocate range max(I0,I1+1,I2)
        call mem_alloc(tval%elms, &
             tval%level_range_uv(1,1)-smallest_uv+1,&
             tval%level_range_uv(1,2)-smallest_uv+1)
     case(1)
        ! allocate range max(I1,I2)
        range_uv_alloc = max(range_tmp(:,2),range_tmp(:,3))
        call mem_alloc(tval%elms, &
             range_uv_alloc(1)-smallest_uv+1,&
             range_uv_alloc(2)-smallest_uv+1)
     case(2)
        ! allocate range (-1:I2)
        call mem_alloc(tval%elms, &
             tval%level_range_uv(3,1)-smallest_uv+1,&
             tval%level_range_uv(3,2)-smallest_uv+1)
     case default 
        write(LOUT,'(a)') 'ERROR!!! Incorrect &
                    & tval allocation in CCint2!'
        stop
     end select
   end associate
enddo

deallocate(range_tmp)

end subroutine create_tval

subroutine free_tval(Two,max_lam)
implicit none
type(TwoIntData) :: Two
integer :: max_lam, j

if(Two%isUsed) then

   do j=1,max_lam+1!size(Two%t_val)
      call mem_dealloc(Two%t_val(j)%level_range_uv)
      call mem_dealloc(Two%t_val(j)%elms)
   enddo

    deallocate(Two%t_val)
endif

end subroutine free_tval

subroutine create_Izero(Two,max_lam)
implicit none

type(TwoIntData) :: Two
integer :: ival,max_lam
integer :: i,j,u,v

do ival=1,max_lam+1
   associate(&
        tval     => Two%t_val(ival), &
        two_t    => Two%t_val(ival)%two_t,& 
        range_uv => Two%t_val(ival)%level_range_uv(1,:) )

        if(Two%same_exp) then
 
           do v=smallest_uv,range_uv(2)
              do u=smallest_uv,min(v,range_uv(1))
                 tval%elms(offset_uv + u,offset_uv + v) = &
                      int2_slater(u,v,two_t,Two%ijalphaPL,Two%klalphaPL)
              enddo
           enddo
           do v=smallest_uv,minval(range_uv)
              do u=v+1,minval(range_uv)
                 tval%elms(offset_uv + u,offset_uv + v) = &
                      tval%elms(offset_uv + v,offset_uv + u)
              enddo
           enddo
           if(range_uv(1)>range_uv(2)) then
              do v=smallest_uv,range_uv(2)
                 do u=range_uv(2)+1,range_uv(1)
                    tval%elms(offset_uv + u,offset_uv + v) = &
                         int2_slater(u,v,two_t,Two%ijalphaPL,Two%klalphaPL)
                 enddo
              enddo
           endif

        else

           do v=smallest_uv,range_uv(2)
              do u=smallest_uv,range_uv(1)
                 tval%elms(offset_uv + u,offset_uv + v) = &
                      int2_slater(u,v,two_t,Two%ijalphaPL,Two%klalphaPL)
              enddo
           enddo

        endif

   end associate
enddo
!
! hapka: old_drake
!        do t=smallest_t,Two%range_t
 !           pos_t = offset_t + t
 !           associate(&
 !                r12      => Two%r12(pos_t), &
 !                range_uv => Two%range_uv(:,pos_t))
 !
 !             call mem_alloc(r12%elms,&
 !                  range_uv(1)-smallest_uv+1,&
 !                  range_uv(2)-smallest_uv+1)
 !
 !             if(Two%same_exp) then
 !
 !                do v=smallest_uv,range_uv(2)
 !                   do u=smallest_uv,min(v,range_uv(1))
 !                      r12%elms(offset_uv + u,offset_uv + v) = &
 !                           int2_slater(u,v,t,Two%ijalphaPL,Two%klalphaPL)
 !                   enddo
 !                enddo
 !                do v=smallest_uv,minval(range_uv)
 !                   do u=v+1,minval(range_uv)
 !                      r12%elms(offset_uv + u,offset_uv + v) = &
 !                           r12%elms(offset_uv + v,offset_uv + u)
 !                   enddo
 !                enddo
 !                if(range_uv(1)>range_uv(2)) then
 !                   do v=smallest_uv,range_uv(2)
 !                      do u=range_uv(2)+1,range_uv(1)
 !                         r12%elms(offset_uv + u,offset_uv + v) = &
 !                              int2_slater(u,v,t,Two%ijalphaPL,Two%klalphaPL)
 !                      enddo
 !                   enddo
 !                endif
 !
 !             else
 !
 !                do v=smallest_uv,range_uv(2)
 !                   do u=smallest_uv,range_uv(1)
 !                      r12%elms(offset_uv + u,offset_uv + v) = &
 !                           int2_slater(u,v,t,Two%ijalphaPL,Two%klalphaPL)
 !                   enddo
 !                enddo
 !
 !             endif
 !
 !           end associate
 !        enddo

end subroutine create_Izero

subroutine create_Ilambda(Two,max_lam)
implicit none

type(TwoIntData) :: Two
integer :: max_lam
real(prec) :: lam_tmp,c_tmp,prefac
integer :: i,j
integer :: ival,jval
integer :: range_uv(2)
integer :: span,start_uv,pos

! loop over recursion
do ival=1,max_lam

   ! create I_1
   range_uv(1) = Two%t_val(max_lam+2-ival)%level_range_uv(ival+1,1)
   range_uv(2) = Two%t_val(max_lam+2-ival)%level_range_uv(ival+1,2)
   do j=range_uv(2),smallest_uv_vec(1),-1
      do i=range_uv(1),smallest_uv_vec(1),-1
         Two%t_val(max_lam+2-ival)%elms(offset_uv+i,offset_uv+j)= 0.5_prec* ( &
         Two%t_val(max_lam+1-ival)%elms(offset_uv+i+1,offset_uv+j-1) + &
         Two%t_val(max_lam+1-ival)%elms(offset_uv+i-1,offset_uv+j+1) - &
         Two%t_val(max_lam+2-ival)%elms(offset_uv+i-1,offset_uv+j-1) )
      enddo
   enddo

   Two%t_val(max_lam+2-ival)%two_t = Two%t_val(max_lam+2-ival)%two_t - 2
   Two%t_val(max_lam+2-ival)%lambda = 1

   if(ival.eq.1) cycle

   ! create I_n
   pos = 0
   do jval=1,ival-1

   lam_tmp = Two%t_val(max_lam+2-ival+jval)%lambda
   c_tmp = Two%t_val(max_lam+2-ival+jval)%two_t - 2._prec 
   prefac = (2._prec*lam_tmp + 1._prec) / (c_tmp + 2._prec)

   start_uv = jval
   range_uv(1) = Two%t_val(max_lam+2-ival+jval)%level_range_uv(ival+1,1) 
   range_uv(2) = Two%t_val(max_lam+2-ival+jval)%level_range_uv(ival+1,2) 
   do j=range_uv(2),start_uv,-1 
      do i=range_uv(1),start_uv,-1 

           Two%t_val(max_lam+2-ival+jval)%elms(offset_uv+i,offset_uv+j) = prefac* ( &
           Two%t_val(max_lam+2-ival+jval)%elms(offset_uv+i+pos-1,offset_uv+j+pos-1)) + &
           Two%t_val(max_lam-ival+jval)%elms(offset_uv+i+pos,offset_uv+j+pos) 

      enddo
   enddo

   pos = pos - 1
   Two%t_val(max_lam+2-ival+jval)%lambda = Two%t_val(max_lam+2-ival+jval)%lambda + 1
   Two%t_val(max_lam+2-ival+jval)%two_t = Two%t_val(max_lam+2-ival+jval)%two_t - 2

   enddo 

enddo

end subroutine create_Ilambda

subroutine rewrite_Ilambda(Two,t,pos_t,max_lam)
implicit none

type(TwoIntData) :: Two
integer :: t,pos_t,max_lam
integer :: u,v
integer :: ilam

do ilam=0,2
   associate(&
          TwoL     => Two%ILambda(ilam),&
          tval     => Two%t_val(ilam+1),&
          lam      => Two%t_val(ilam+1)%lambda,&
          two_t    => Two%t_val(ilam+1)%two_t,& 
         range_uv => Two%t_val(ilam+1)%level_range_uv(max_lam+1,:))

     if(pos_t.gt.(TwoL%range_t+offset_t)) cycle
     associate(&
          r12      => TwoL%r12(pos_t), &
          range_uv2 => TwoL%range_uv(:,pos_t) )

       do v=smallest_uv_vec(ilam),range_uv2(2)
          do u=smallest_uv_vec(ilam),range_uv2(1)
             !r12%elms(offset_uv-ilam+u,offset_uv-ilam+v) = & 
             r12%elms(offset_uv_vec(ilam)+u,offset_uv_vec(ilam)+v) = & 
                 tval%elms(offset_uv + u,offset_uv + v)
          enddo
       enddo
     end associate     
   end associate
enddo

end subroutine rewrite_Ilambda

! hapka: old_drake
!subroutine create_Twoint
!implicit none
!integer :: iTwo
!integer :: t,pos_t,u,v
!
! do iTwo=1,size(TwoInt)
!    associate(Two => TwoInt(iTwo))
!      if(Two%isUsed) then
 
 !        allocate(Two%r12(Two%range_t-smallest_t+1))
 !
 !        do t=smallest_t,Two%range_t
 !           pos_t = offset_t + t
 !           associate(&
 !                r12      => Two%r12(pos_t), &
 !                range_uv => Two%range_uv(:,pos_t))
 !
 !             call mem_alloc(r12%elms,&
 !                  range_uv(1)-smallest_uv+1,&
 !                  range_uv(2)-smallest_uv+1)
 !
 !             if(Two%same_exp) then
 !
 !                do v=smallest_uv,range_uv(2)
 !                   do u=smallest_uv,min(v,range_uv(1))
 !                      r12%elms(offset_uv + u,offset_uv + v) = &
 !                           int2_slater(u,v,t,Two%ijalphaPL,Two%klalphaPL)
 !                   enddo
 !                enddo
 !                do v=smallest_uv,minval(range_uv)
 !                   do u=v+1,minval(range_uv)
 !                      r12%elms(offset_uv + u,offset_uv + v) = &
 !                           r12%elms(offset_uv + v,offset_uv + u)
 !                   enddo
 !                enddo
 !                if(range_uv(1)>range_uv(2)) then
 !                   do v=smallest_uv,range_uv(2)
 !                      do u=range_uv(2)+1,range_uv(1)
 !                         r12%elms(offset_uv + u,offset_uv + v) = &
 !                              int2_slater(u,v,t,Two%ijalphaPL,Two%klalphaPL)
 !                      enddo
 !                   enddo
 !                endif
 !
 !             else
 !
 !                do v=smallest_uv,range_uv(2)
 !                   do u=smallest_uv,range_uv(1)
 !                      r12%elms(offset_uv + u,offset_uv + v) = &
 !                           int2_slater(u,v,t,Two%ijalphaPL,Two%klalphaPL)
 !                   enddo
 !                enddo
 !
 !             endif
 !
 !           end associate
 !        enddo
 !
!      endif
!    end associate
! enddo
!
!end subroutine create_Twoint


subroutine free_TwoInt
implicit none
integer :: iTwo
integer :: ilam
integer :: t

do iTwo=1,size(TwoInt)
   associate(Two => TwoInt(iTwo))
     if(Two%isUsed) then
        do ilam=0,2
           associate(TwoL => Two%ILambda(ilam))
             do t=smallest_t,TwoL%range_t
                call mem_dealloc(TwoL%r12(offset_t + t)%elms)
             enddo

             deallocate(TwoL%r12)
           end associate
        enddo  
     endif
   end associate
enddo

do iTwo=1,size(TwoInt)
   associate(Two => TwoInt(iTwo))

     if(Two%isUsed) then 
        do ilam=0,2
           associate(TwoL => Two%ILambda(ilam))
             call mem_dealloc(TwoL%range_uv)
           end associate
        enddo
     endif

   end associate
enddo

! hapka: old_drake
!do iTwo=1,size(TwoInt)
!   associate(Two => TwoInt(iTwo))
!     if(Two%isUsed) then
!
!        do t=smallest_t,Two%range_t
!           call mem_dealloc(Two%r12(offset_t + t)%elms)
!        enddo
!
!        deallocate(Two%r12)
!
!     endif
!   end associate
!enddo
!
!do iTwo=1,size(TwoInt)
!   associate(Two => TwoInt(iTwo))
!
!     if(Two%isUsed) call mem_dealloc(Two%range_uv)
!
!   end associate
!enddo

deallocate(TwoInt)

end subroutine free_TwoInt

subroutine print_TwoInt(LPRINT)
implicit none
integer,intent(in) :: LPRINT
integer :: iTwo
integer :: ilam
integer :: t,pos_t,u,v
integer :: ctr

 !if(LPRINT>=10) then
 if(LPRINT==6) then
 
    write(LOUT,'()')
    write(LOUT,'(5x,a)') '--- CC two-electron integrals ---'
 
    write(LOUT,'()')
    write(LOUT,'(1x,a,2x,4(2x,a),2x,4(5x,a,1x))') &
         'no.','i','j','k','l','ijalphaPL','ijalphaMI','klalphaPL','klalphaMI'
    do iTwo=1,size(TwoInt)
       associate(Two => TwoInt(iTwo))
         write(LOUT,'(1x,i3,a)',advance='no') iTwo,' :'
         if(Two%isUsed) then
            write(LOUT,'(4i3,1x,a,4f15.8)') &
                 Two%iexp,Two%jexp,Two%kexp,Two%lexp,&
                 merge('*',' ',Two%same_exp),&
                 Two%ijalphaPL,Two%ijalphaMI,Two%klalphaPL,Two%klalphaMI
            do ilam=0,2
               associate(TwoL => Two%ILambda(ilam))
                 do t=smallest_t,TwoL%range_t
                    pos_t = offset_t + t
                    write(LOUT,'(11x,a,i3,a,3x,a,i3,a,6x,a,2i4)') &
                          'lam = ',ilam,',','t = ',t,',','max_u max_v = ',TwoL%range_uv(:,pos_t)
! check ILambda
!              if(LPRINT>=100) then
               if(LPRINT>=66) then
                  associate(&
                       r12       => TwoL%r12(pos_t), &
                       range_uv  => TwoL%range_uv(:,pos_t),&
                       offset_uv => offset_uv_veC(ilam))
                    do v=smallest_uv_vec(ilam),range_uv(2)
                       do u=smallest_uv_vec(ilam),range_uv(1)
                          write(LOUT,'(10x,2i3,a,es40.32)') u,v,' : ',&
                              !r12%elms(offset_uv-ilam+u,offset_uv-ilam+v)
                               r12%elms(offset_uv+u,offset_uv+v)
                       enddo
                    enddo
                  end associate
               endif
                 enddo
               end associate
            enddo
!
! hapka: old_drake
!            do t=smallest_t,Two%range_t
!               pos_t = offset_t + t
!               write(LOUT,'(11x,a,i3,a,6x,a,2i4)') &
!                    't = ',t,',','max_u max_v = ',Two%range_uv(:,pos_t)
!!               if(LPRINT>=100) then
!!                  associate(&
!!                       r12      => Two%r12(pos_t), &
!!                       range_uv => Two%range_uv(:,pos_t))
!!                    do v=smallest_uv,range_uv(2)
!!                       do u=smallest_uv,range_uv(1)
!!                          write(LOUT,'(10x,2i3,a,es30.22)') u,v,' : ',&
!!                               r12%elms(offset_uv + u,offset_uv + v)
!!                       enddo
!!                    enddo
!!                  end associate
!!               endif
!            enddo
         else
            write(LOUT,'(a)') 'NOT USED'
         endif
       end associate
    enddo
 
 endif

end subroutine print_TwoInt

subroutine get_genpair_info(iPairGen,jPairGen,ABtype,info)
implicit none

integer,intent(in) :: iPairGen,jPairGen
character(1),intent(in) :: ABtype
type(InfoGenPair) :: info(0:2)

! init
info(:)%isUsed = .false.
info(:)%coeff = 0
info(:)%add_to_1 = 0
info(:)%add_to_2 = 0

select case(ABtype)
case('A','a')

   if(iPairGen==jPairGen) then
      ! even generators
      select case(iPairGen)
      case(1) ! S^e
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
      case(2) ! P^o
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
         info(0)%add_to_1 = 0
         info(0)%add_to_2 = 2
      case(3) ! P^e
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
         info(0)%add_to_1 = 2
         info(0)%add_to_2 = 2
         info(2)%isUsed = .true.
         info(2)%coeff = -0.5_prec
         info(2)%add_to_1 = 2
         info(2)%add_to_2 = 2
      case(4) ! D^e(1)
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
         info(0)%add_to_1 = 2
         info(0)%add_to_2 = 2
         info(2)%isUsed = .true.
         info(2)%coeff = 0.1_prec
         info(2)%add_to_1 = 2
         info(2)%add_to_2 = 2
      case(5) ! D^e(2)
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
         info(0)%add_to_1 = 0
         info(0)%add_to_2 = 4
      case default
         write(LOUT,'(a)') 'Incorrect generator &
          & passed to get_genpair_info (A-even))!!!'
         stop
      end select
   else
      ! odd generators
      select case(iPairGen)
      case(4,5)
         info(1)%isUsed = .true.
         info(1)%coeff = sqrt(0.3_prec)
         info(1)%add_to_1 = 1
         info(1)%add_to_2 = 3
      case default
         write(LOUT,'(a)') 'Incorrect generator &
          & passed to get_genpair_info (A-odd))!!!'
         stop
      end select
   endif

case('B','b')

   if(iPairGen==jPairGen) then
      ! even generators
      select case(iPairGen)
      case(1) ! S^e
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
      case(2) ! P^o
         info(1)%isUsed = .true.
         info(1)%coeff = 0.5_prec
         info(1)%add_to_1 = 1
         info(1)%add_to_2 = 1
      case(3) ! P^e
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
         info(0)%add_to_1 = 2
         info(0)%add_to_2 = 2
         info(2)%isUsed = .true.
         info(2)%coeff = -0.5_prec
         info(2)%add_to_1 = 2
         info(2)%add_to_2 = 2
      case(4) ! D^e(1)
         info(0)%isUsed = .true.
         info(0)%coeff = 0.5_prec
         info(0)%add_to_1 = 2
         info(0)%add_to_2 = 2
         info(2)%isUsed = .true.
         info(2)%coeff = 0.1_prec
         info(2)%add_to_1 = 2
         info(2)%add_to_2 = 2
      case(5) ! D^e(2)
         info(2)%isUsed = .true.
         info(2)%coeff = 0.5_prec
         info(2)%add_to_1 = 2
         info(2)%add_to_2 = 2
      case default
         write(LOUT,'(a)') 'Incorrect generator &
          & passed in get_genpair_info (B-even)!!!'
         stop
      end select
   else
      ! odd generators
      select case(iPairGen)
      case(4,5)
         info(1)%isUsed = .true.
         info(1)%coeff = sqrt(0.3_prec)
         info(1)%add_to_1 = 3
         info(1)%add_to_2 = 1
      case default
         write(LOUT,'(a)') 'Incorrect generator &
          & passed in get_genpair_info (B-odd)!!!'
         stop
      end select
   endif

case default

write(LOUT,'(a)') 'ERROR!!! Incorrect ABtype in get_genpair_info!'
stop

end select

end subroutine get_genpair_info

end module CCint2
