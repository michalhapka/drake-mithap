module CCint3
use file_OUT, only : LOUT
use precision, only : prec
use memory
use commontypes
implicit none

private
public create_CCint3,free_CCint3
public CCint3_matJ,CCint3_matK,CCint3_matM,CCint3_matP,CCint3_vecP

interface CCint3_matJ
module procedure CCint3_matJ1
module procedure CCint3_matJ2
end interface CCint3_matJ

interface CCint3_matK
module procedure CCint3_matK1
module procedure CCint3_matK2
end interface CCint3_matK

character :: int3_source

character(*),parameter :: script = 'script.sh'
character(*),parameter :: infile = 'ThreeElectron.in'
character(*),parameter :: core_J = 'intJ', core_K = 'intK'
integer,parameter :: name_length = 20

integer,parameter :: smallest_uv = 0
integer,parameter :: offset_uv = 1 - smallest_uv

real(prec),parameter :: &
     val_modifier = 5.039302255187420185065945732588076342666e-4_prec

type int3valData
real(prec),allocatable :: elms(:,:,:)
end type int3valData

type J_SpecData_Lam
logical :: isUsed
integer :: sum_ij,max_ij,max_k
end type J_SpecData_Lam

type J_SpecData
logical :: isUsed
character(2) :: symbol_ij
character(1) :: symbol_k
type(J_SpecData_Lam) :: J_SpecLam(0:2)
!integer :: sum_ij,max_ij,max_k
end type J_SpecData

type K_SpecData
logical :: isUsed
character(2) :: symbol_ij
integer :: sum_TOT,max_2,max_1
end type K_SpecData

type ThreeIntData
logical :: isUsed_J,isUsed_K
integer :: iexp,jexp,kexp,lexp,mexp,nexp
logical :: sameANY,sameALL,same12,same23
real(prec) :: ijalpha,klalpha,mnalpha
type(J_SpecData),allocatable :: J_Spec(:)
type(K_SpecData),allocatable :: K_Spec(:)
end type ThreeIntData

type(ThreeIntData),allocatable :: ThreeInt(:)

contains

subroutine create_CCint3(TOTexp,exponents,OrbReduced,PairReduced,INT3,LPRINT)
implicit none
integer,intent(in) :: TOTexp
real(prec),intent(in) :: exponents(:)
type(OrbReducedData),intent(in) :: OrbReduced(:)
type(PairReducedData),intent(in) :: PairReduced(:)
character,intent(in) :: INT3
integer,intent(in) :: LPRINT

int3_source = INT3

call init_ThreeInt(TOTexp,exponents)
call init_J1_ThreeInt(OrbReduced,PairReduced)
!call init_J2_ThreeInt(OrbReduced,PairReduced)
!call init_K_ThreeInt(OrbReduced,PairReduced)

call print_ThreeInt(LPRINT)

!call create_ThreeInt

end subroutine create_CCint3

subroutine free_CCint3
implicit none

call free_ThreeInt

int3_source = ''

end subroutine free_CCint3

subroutine CCint3_matJ1(ABtype,matJ,iPairSystem,jPairSystem,&
     vector,OrbSystem)
implicit none
character(1),intent(in) :: ABtype
real(prec) :: matJ(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: vector(:)
type(OrbSystemData),intent(in) :: OrbSystem

call classJ(ABtype,matJ,iPairSystem,jPairSystem,&
     vector,OrbSystem,vector,OrbSystem)

end subroutine CCint3_matJ1

subroutine CCint3_matJ2(ABtype,matJ,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)
implicit none
character(1),intent(in) :: ABtype
real(prec) :: matJ(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: i_vector(:),j_vector(:)
type(OrbSystemData),intent(in) :: i_OrbSystem,j_OrbSystem

call classJ(ABtype,matJ,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)

end subroutine CCint3_matJ2

subroutine CCint3_matK1(ABtype,matK,iPairSystem,jPairSystem,&
     vector,OrbSystem)
implicit none
character(1),intent(in) :: ABtype
real(prec) :: matK(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: vector(:)
type(OrbSystemData),intent(in) :: OrbSystem

call classK(1,ABtype,.true.,0,matK,iPairSystem,jPairSystem,&
     vector,OrbSystem,vector,OrbSystem)

end subroutine CCint3_matK1

subroutine CCint3_matK2(ABtype,matK,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)
implicit none
character(1),intent(in) :: ABtype
real(prec) :: matK(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: i_vector(:),j_vector(:)
type(OrbSystemData),intent(in) :: i_OrbSystem,j_OrbSystem

call classK(1,ABtype,.true.,0,matK,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)

end subroutine CCint3_matK2

subroutine CCint3_matM(ABtype,matM,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)
implicit none
character(1),intent(in) :: ABtype
real(prec) :: matM(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: i_vector(:),j_vector(:)
type(OrbSystemData),intent(in) :: i_OrbSystem,j_OrbSystem

call classK(1,ABtype,.false.,0,matM,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)

end subroutine CCint3_matM

subroutine CCint3_matP(ABtype,add_tj,matP,iPairSystem,jPairSystem,&
     vector,OrbSystem)
implicit none
character(1),intent(in) :: ABtype
integer,intent(in) :: add_tj
real(prec) :: matP(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: vector(:)
type(OrbSystemData),intent(in) :: OrbSystem

call classK(0,ABtype,.true.,add_tj,matP,iPairSystem,jPairSystem,&
     vector,OrbSystem,vector,OrbSystem)

end subroutine CCint3_matP

subroutine classJ(ABtype,matJ,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)
implicit none
character(1),intent(in) :: ABtype
real(prec) :: matJ(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: i_vector(:),j_vector(:)
type(OrbSystemData),intent(in) :: i_OrbSystem,j_OrbSystem
integer :: i_prim,j_prim,i_orb,j_orb
integer :: sum_12,max_12,max_3
integer :: iThree,perm1(3),perm2(3),order2(3),order_ij,r_xk
logical :: pos_uj_vj
integer :: pos_1,pos_2
integer :: ui,vi,ti,ipos,iflag
integer :: uj,vj,tj,jpos,jflag
integer :: i,i_pos
integer :: j,j_pos
real(prec) :: rtmp
type(int3valData),allocatable :: int3val(:)

matJ = 0._prec

do j_prim=1,jPairSystem%n_prim
   do i_prim=1,iPairSystem%n_prim
      associate(&
           iPairSpec => iPairSystem%PairSpec(i_prim), &
           jPairSpec => jPairSystem%PairSpec(j_prim))

        do j_orb=1,j_OrbSystem%n_prim
           do i_orb=1,i_OrbSystem%n_prim
              associate(&
                   i_OrbSpec => i_OrbSystem%OrbSpec(i_orb),&
                   j_OrbSpec => j_OrbSystem%OrbSpec(j_orb))

                sum_12 = &
                     maxOmega_PairSpec(iPairSpec) + &
                     maxOmega_PairSpec(jPairSpec)

                max_12 = &
                     maxt_PairSpec(iPairSpec) + &
                     maxt_PairSpec(jPairSpec)

                max_3 = i_OrbSpec%irange + j_OrbSpec%irange

                select case(ABtype)
                case('A','a')

                   order2 = [1,0,0]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        iPairSpec%iexp1,jPairSpec%iexp1,&
                        iPairSpec%iexp2,jPairSpec%iexp2,&
                        i_OrbSpec%iexp ,j_OrbSpec%iexp)
                   order_ij = maxloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   pos_uj_vj = .true.

                case('B','b')

                   order2 = [1,0,0]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        iPairSpec%iexp1,jPairSpec%iexp2,&
                        iPairSpec%iexp2,jPairSpec%iexp1,&
                        i_OrbSpec%iexp ,j_OrbSpec%iexp)
                   order_ij = maxloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   pos_uj_vj = .false.

                case default

                   write(LOUT,'(a)') 'ERROR!!! Incorrect ABtype in classJ!'
                   stop

                end select

                do r_xk=1,2

                   call read_ThreeInt_J(int3val,&
                        iThree,order_ij,r_xk,sum_12,max_12,max_3,perm1,perm2)

                   jpos  = jPairSpec%offset
                   jflag = -1
                   do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                      jpos = jpos + 1

                      ipos  = iPairSpec%offset
                      iflag = -1
                      do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                         ipos = ipos + 1

                         rtmp = 0._prec

                         if(pos_uj_vj) then
                            pos_1 = offset_uv + ui + uj
                            pos_2 = offset_uv + vi + vj
                         else
                            pos_1 = offset_uv + ui + vj
                            pos_2 = offset_uv + vi + uj
                         endif
                         associate(elms => int3val(ti + tj)%elms(:,pos_1,pos_2))

                           j_pos = j_OrbSpec%offset
                           do j=0,j_OrbSpec%irange
                              j_pos = j_pos + 1

                              i_pos = i_OrbSpec%offset
                              do i=0,i_OrbSpec%irange
                                 i_pos = i_pos + 1

                                 rtmp = rtmp &
                                      + elms(offset_uv + i + j) &
                                      * i_vector(i_pos)*j_vector(j_pos)

                              enddo

                           enddo

                         end associate

                         matJ(ipos,jpos) = matJ(ipos,jpos) + rtmp
                            
                      enddo

                   enddo

                   call free_ThreeInt_J(int3val)

                enddo

              end associate
           enddo
        enddo

      end associate
   enddo
enddo

end subroutine classJ

subroutine classK(extra,ABtype,add_2,add_tj,matK,iPairSystem,jPairSystem,&
     i_vector,i_OrbSystem,j_vector,j_OrbSystem)
implicit none
integer,intent(in) :: extra
character(1),intent(in) :: ABtype
logical,intent(in) :: add_2
integer,intent(in) :: add_tj
real(prec) :: matK(:,:)
type(PairSystemData),intent(in) :: iPairSystem,jPairSystem
real(prec),intent(in) :: i_vector(:),j_vector(:)
type(OrbSystemData),intent(in) :: i_OrbSystem,j_OrbSystem
integer :: i_prim,j_prim,i_orb,j_orb
integer :: imax,jmax,sum_TOT,max_2,max_1
integer :: iThree,perm1(3),perm2(3),order2(3),order_ij
integer :: pos_1,pos_2
integer :: ui,vi,ti,ipos,iflag
integer :: uj,vj,tj,jpos,jflag
integer :: i,i_pos
integer :: j,j_pos
real(prec) :: rtmp
type(int3valData),allocatable :: int3val(:,:)

if(add_tj<-1) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Additional t parameter in classK is incorrect!'
   stop
endif

matK = 0._prec

do j_prim=1,jPairSystem%n_prim
   do i_prim=1,iPairSystem%n_prim
      associate(&
           iPairSpec => iPairSystem%PairSpec(i_prim),&
           jPairSpec => jPairSystem%PairSpec(j_prim))

        do j_orb=1,j_OrbSystem%n_prim
           do i_orb=1,i_OrbSystem%n_prim
              associate(&
                   i_OrbSpec => i_OrbSystem%OrbSpec(i_orb),&
                   j_OrbSpec => j_OrbSystem%OrbSpec(j_orb))

                sum_TOT = &
                     maxOmega_PairSpec(iPairSpec) + &
                     maxOmega_PairSpec(jPairSpec) + &
                     i_OrbSpec%irange + &
                     j_OrbSpec%irange

                max_2 = max(&
                     maxt_PairSpec(iPairSpec),&
                     maxt_PairSpec(jPairSpec))

                imax = maxuv_PairSpec(iPairSpec)
                jmax = maxuv_PairSpec(jPairSpec)

                max_1 = max(&
                     imax + j_OrbSpec%irange,&
                     imax + jmax,&
                     i_OrbSpec%irange + jmax)

                select case(ABtype)
                case('A','a')

                   order2 = [0,-1,0]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        iPairSpec%iexp1,j_OrbSpec%iexp,&
                        iPairSpec%iexp2,jPairSpec%iexp2,&
                        i_OrbSpec%iexp ,jPairSpec%iexp1)
                   order_ij = minloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   call read_ThreeInt_K(int3val,&
                        iThree,order_ij,extra,sum_TOT,max_2,max_1,perm1,perm2)

                   jpos  = jPairSpec%offset
                   jflag = -1
                   do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                      jpos = jpos + 1

                      ipos  = iPairSpec%offset
                      iflag = -1
                      do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                         ipos = ipos + 1

                         rtmp = 0._prec

                         associate(elms => int3val(ti, tj + add_tj)%elms)
                           pos_2 = offset_uv + vi + vj

                           j_pos = j_OrbSpec%offset
                           do j=0,j_OrbSpec%irange
                              j_pos = j_pos + 1

                              i_pos = i_OrbSpec%offset
                              do i=0,i_OrbSpec%irange
                                 i_pos = i_pos + 1

                                 rtmp = rtmp + elms(&
                                      offset_uv + ui + j ,&
                                      pos_2,&
                                      offset_uv + i  + uj)&
                                      * i_vector(i_pos)*j_vector(j_pos)

                              enddo

                           enddo

                         end associate

                         matK(ipos,jpos) = matK(ipos,jpos) + rtmp

                      enddo

                   enddo

                   call free_ThreeInt_K(int3val)

                   order2 = [0,0,-1]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        iPairSpec%iexp1,jPairSpec%iexp1,&
                        iPairSpec%iexp2,j_OrbSpec%iexp,&
                        i_OrbSpec%iexp ,jPairSpec%iexp2)
                   order_ij = minloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   call read_ThreeInt_K(int3val,&
                        iThree,order_ij,extra,sum_TOT,max_2,max_1,perm1,perm2)

                   jpos  = jPairSpec%offset
                   jflag = -1
                   do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                      jpos = jpos + 1

                      ipos  = iPairSpec%offset
                      iflag = -1
                      do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                         ipos = ipos + 1

                         rtmp = 0._prec

                         associate(elms => int3val(ti, tj + add_tj)%elms)
                           pos_1 = offset_uv + ui + uj

                           j_pos = j_OrbSpec%offset
                           do j=0,j_OrbSpec%irange
                              j_pos = j_pos + 1

                              i_pos = i_OrbSpec%offset
                              do i=0,i_OrbSpec%irange
                                 i_pos = i_pos + 1

                                 rtmp = rtmp + elms(&
                                      pos_1,&
                                      offset_uv + vi + j ,&
                                      offset_uv + i  + vj)&
                                      * i_vector(i_pos)*j_vector(j_pos)

                              enddo

                           enddo

                         end associate

                         if(add_2) then
                            matK(ipos,jpos) = matK(ipos,jpos) + rtmp
                         else
                            matK(ipos,jpos) = matK(ipos,jpos) - rtmp
                         endif

                      enddo

                   enddo

                   call free_ThreeInt_K(int3val)

                case('B','b')

                   order2 = [0,-1,0]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        iPairSpec%iexp1,j_OrbSpec%iexp,&
                        iPairSpec%iexp2,jPairSpec%iexp1,&
                        i_OrbSpec%iexp ,jPairSpec%iexp2)
                   order_ij = minloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   call read_ThreeInt_K(int3val,&
                        iThree,order_ij,extra,sum_TOT,max_2,max_1,perm1,perm2)

                   jpos  = jPairSpec%offset
                   jflag = -1
                   do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                      jpos = jpos + 1

                      ipos  = iPairSpec%offset
                      iflag = -1
                      do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                         ipos = ipos + 1

                         rtmp = 0._prec

                         associate(elms => int3val(ti, tj + add_tj)%elms)
                           pos_2 = offset_uv + vi + uj

                           j_pos = j_OrbSpec%offset
                           do j=0,j_OrbSpec%irange
                              j_pos = j_pos + 1

                              i_pos = i_OrbSpec%offset
                              do i=0,i_OrbSpec%irange
                                 i_pos = i_pos + 1

                                 rtmp = rtmp + elms(&
                                      offset_uv + ui + j ,&
                                      pos_2,&
                                      offset_uv + i  + vj)&
                                      * i_vector(i_pos)*j_vector(j_pos)

                              enddo

                           enddo

                         end associate

                         matK(ipos,jpos) = matK(ipos,jpos) + rtmp

                      enddo

                   enddo

                   call free_ThreeInt_K(int3val)

                   order2 = [0,0,-1]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        iPairSpec%iexp1,jPairSpec%iexp2,&
                        iPairSpec%iexp2,j_OrbSpec%iexp,&
                        i_OrbSpec%iexp ,jPairSpec%iexp1)
                   order_ij = minloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   call read_ThreeInt_K(int3val,&
                        iThree,order_ij,extra,sum_TOT,max_2,max_1,perm1,perm2)

                   jpos  = jPairSpec%offset
                   jflag = -1
                   do while(increment_PairSpec(uj,vj,tj,jflag,jPairSpec))
                      jpos = jpos + 1

                      ipos  = iPairSpec%offset
                      iflag = -1
                      do while(increment_PairSpec(ui,vi,ti,iflag,iPairSpec))
                         ipos = ipos + 1

                         rtmp = 0._prec

                         associate(elms => int3val(ti, tj + add_tj)%elms)
                           pos_1 = offset_uv + ui + vj

                           j_pos = j_OrbSpec%offset
                           do j=0,j_OrbSpec%irange
                              j_pos = j_pos + 1

                              i_pos = i_OrbSpec%offset
                              do i=0,i_OrbSpec%irange
                                 i_pos = i_pos + 1

                                 rtmp = rtmp + elms(&
                                      pos_1,&
                                      offset_uv + vi + j ,&
                                      offset_uv + i  + uj)&
                                      * i_vector(i_pos)*j_vector(j_pos)

                              enddo

                           enddo

                         end associate

                         if(add_2) then
                            matK(ipos,jpos) = matK(ipos,jpos) + rtmp
                         else
                            matK(ipos,jpos) = matK(ipos,jpos) - rtmp
                         endif

                      enddo

                   enddo

                   call free_ThreeInt_K(int3val)

                case default

                   write(LOUT,'(a)') 'ERROR!!! &
                        &Incorrect ABtype in classK!'
                   stop

                end select

              end associate
           enddo
        enddo

      end associate
   enddo
enddo

end subroutine classK

subroutine CCint3_vecP(vecP,&
     PairSystem,A_vector,A_OrbSystem,B_vector,B_OrbSystem,&
     vector,OrbSystem)
implicit none
real(prec) :: vecP(:)
type(PairSystemData),intent(in) :: PairSystem
real(prec),intent(in) :: A_vector(:),B_vector(:)
type(OrbSystemData),intent(in) :: A_OrbSystem,B_OrbSystem
real(prec) :: vector(:)
type(OrbSystemData),intent(in) :: OrbSystem
integer :: prim,A_prim,B_prim,i_orb,j_orb
integer :: sum_12,max_12
integer :: iThree,perm1(3),perm2(3),order2(3),order_ij
integer :: u,v,t,pos,flag
integer :: A,A_pos
integer :: B,B_pos
integer :: i,i_pos
integer :: j,j_pos
real(prec) :: rtmp
type(int3valData),allocatable :: int3val(:)

vecP = 0._prec

do B_prim=1,B_OrbSystem%n_prim
   do A_prim=1,A_OrbSystem%n_prim
      do prim=1,PairSystem%n_prim
         associate(&
              PairSpec  => PairSystem%PairSpec(prim),&
              A_OrbSpec => A_OrbSystem%OrbSpec(A_prim),&
              B_OrbSpec => B_OrbSystem%OrbSpec(B_prim))

           do j_orb=1,OrbSystem%n_prim
              do i_orb=1,OrbSystem%n_prim
                 associate(&
                      i_OrbSpec => OrbSystem%OrbSpec(i_orb),&
                      j_OrbSpec => OrbSystem%OrbSpec(j_orb))

                   sum_12 = maxOmega_PairSpec(PairSpec)

                   max_12 = maxt_PairSpec(PairSpec)

                   order2 = [1,0,0]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        PairSpec%iexp1,A_OrbSpec%iexp,&
                        PairSpec%iexp2,j_OrbSpec%iexp,&
                        i_OrbSpec%iexp,B_OrbSpec%iexp)
                   order_ij = maxloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   call read_ThreeInt_J(int3val,&
                        iThree,order_ij,1,&
                        sum_12 + A_OrbSpec%irange + j_OrbSpec%irange,&
                        max_12,&
                        i_OrbSpec%irange + B_OrbSpec%irange,&
                        perm1,perm2)

                   pos  = PairSpec%offset
                   flag = -1
                   do while(increment_PairSpec(u,v,t,flag,PairSpec))
                      pos = pos + 1

                      j_pos = j_OrbSpec%offset
                      do j=0,j_OrbSpec%irange
                         j_pos = j_pos + 1
                         A_pos = A_OrbSpec%offset
                         do A=0,A_OrbSpec%irange
                            A_pos = A_pos + 1

                            rtmp = 0._prec

                            associate(elms => int3val(t)%elms(:,&
                                 offset_uv + u + A,&
                                 offset_uv + v + j))

                              B_pos = B_OrbSpec%offset
                              do B=0,B_OrbSpec%irange
                                 B_pos = B_pos + 1

                                 i_pos = i_OrbSpec%offset
                                 do i=0,i_OrbSpec%irange
                                    i_pos = i_pos + 1

                                    rtmp = rtmp &
                                         + elms(offset_uv + i + B) &
                                         * vector(i_pos)*B_vector(B_pos)

                                 enddo

                              enddo

                            end associate

                            vecP(pos) = vecP(pos) + rtmp &
                                 * A_vector(A_pos)*vector(j_pos)

                         enddo
                      enddo

                   enddo

                   call free_ThreeInt_J(int3val)

                   order2 = [1,0,0]
                   call perm_part1(iThree,perm1,perm2,order2,&
                        PairSpec%iexp1,j_OrbSpec%iexp,&
                        PairSpec%iexp2,B_OrbSpec%iexp,&
                        i_OrbSpec%iexp,A_OrbSpec%iexp)
                   order_ij = maxloc(order2,dim=1)
                   call perm_part2(perm1,perm2,order2)

                   call read_ThreeInt_J(int3val,&
                        iThree,order_ij,2,&
                        sum_12 + j_OrbSpec%irange + B_OrbSpec%irange,&
                        max_12,&
                        i_OrbSpec%irange + A_OrbSpec%irange,&
                        perm1,perm2)

                   pos  = PairSpec%offset
                   flag = -1
                   do while(increment_PairSpec(u,v,t,flag,PairSpec))
                      pos = pos + 1

                      B_pos = B_OrbSpec%offset
                      do B=0,B_OrbSpec%irange
                         B_pos = B_pos + 1
                         j_pos = j_OrbSpec%offset
                         do j=0,j_OrbSpec%irange
                            j_pos = j_pos + 1

                            rtmp = 0._prec

                            associate(elms => int3val(t)%elms(:,&
                                 offset_uv + u + j,&
                                 offset_uv + v + B))

                              A_pos = A_OrbSpec%offset
                              do A=0,A_OrbSpec%irange
                                 A_pos = A_pos + 1

                                 i_pos = i_OrbSpec%offset
                                 do i=0,i_OrbSpec%irange
                                    i_pos = i_pos + 1

                                    rtmp = rtmp &
                                         + elms(offset_uv + i + A) &
                                         * vector(i_pos)*A_vector(A_pos)

                                 enddo

                              enddo

                            end associate

                            vecP(pos) = vecP(pos) + rtmp &
                                 * vector(j_pos)*B_vector(B_pos)

                         enddo
                      enddo

                   enddo

                   call free_ThreeInt_J(int3val)

                 end associate
              enddo
           enddo

         end associate
      enddo
   enddo
enddo

end subroutine CCint3_vecP

subroutine init_ThreeInt(TOTexp,exponents)
implicit none
integer,intent(in) :: TOTexp
real(prec),intent(in) :: exponents(:)
integer :: iexp,jexp,kexp,lexp,mexp,nexp,iThree

iThree = TOTexp*(TOTexp+1)/2
allocate(ThreeInt(iThree*(iThree+1)*(iThree+2)/6))

iThree = 0
do nexp=1,TOTexp
   do mexp=1,nexp
      do lexp=1,nexp
         do kexp=1,merge(mexp,lexp,lexp==nexp)
            do jexp=1,lexp
               do iexp=1,merge(kexp,jexp,jexp==lexp)
                  iThree = iThree + 1
                  associate(Three => ThreeInt(iThree))

                    Three%isUsed_J = .false.
                    Three%isUsed_K = .false.
                    Three%iexp     = iexp
                    Three%jexp     = jexp
                    Three%kexp     = kexp
                    Three%lexp     = lexp
                    Three%mexp     = mexp
                    Three%nexp     = nexp
                    Three%same12   = &
                         (Three%iexp==Three%kexp) .and. &
                         (Three%jexp==Three%lexp)
                    Three%same23   = &
                         (Three%kexp==Three%mexp) .and. &
                         (Three%lexp==Three%nexp)
                    Three%sameANY  = Three%same12.or. Three%same23
                    Three%sameALL  = Three%same12.and.Three%same23
                    Three%ijalpha  = &
                         exponents(Three%iexp) + &
                         exponents(Three%jexp)
                    Three%klalpha  = &
                         exponents(Three%kexp) + &
                         exponents(Three%lexp)
                    Three%mnalpha  = &
                         exponents(Three%mexp) + &
                         exponents(Three%nexp)

                  end associate
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
if(iThree/=size(ThreeInt)) then
   write(LOUT,'(a)') 'ERROR!!! &
        &Number of entries in CC ThreeInt has not been predicted correctly!'
   stop
endif

end subroutine init_ThreeInt

subroutine init_J1_ThreeInt(OrbReduced,PairReduced)
implicit none
type(OrbReducedData),intent(in) :: OrbReduced(:)
type(PairReducedData),intent(in) :: PairReduced(:)
integer :: i_prim,j_prim,i_orb,j_orb
integer :: sum_12,max_12,max_3
integer :: iThree,perm1(3),perm2(3),order2(3)
integer :: igen,jgen,jGenMax
integer :: ilval, jlval, iadd, nadd
integer :: add_lval
integer,allocatable :: jGenNum(:)
integer :: add_A(2), add_B(2)

call mem_alloc(jGenNum,2)
jGenNum = 0

do j_prim=1,size(PairReduced)
   do i_prim=1,j_prim
      associate(&
           iPairSpec => PairReduced(i_prim),&
           jPairSpec => PairReduced(j_prim))
        if(iPairSpec%isUsed.and.jPairSpec%isUsed) then

           ! loop over generators
           do igen=1,iPairSpec%n_gen
               associate(iPairSpecG => iPairSpec%PairReduced_G(igen),&
                         iPairGen => iPairSpec%PairReduced_G(igen)%gen_type)

                 call check_gen_pairs(iPairSpecG,jPairSpec,jGenMax,jGenNum)
                 if(jGenMax.eq.0) cycle
                 do jgen=1,jGenMax
                    associate(jPairSpecG => jPairSpec%PairReduced_G(jGenNum(jgen)),&
                              jPairGen => jPairSpec%PairReduced_G(jGenNum(jgen))%gen_type)
                             
                    ! write(*,*) 'iGen:',iPairSpecG%gen_type,'jGen',jPairSpecG%gen_type
!
                   do j_orb=1,size(OrbReduced)
                      do i_orb=1,j_orb
                         associate(&
                              i_OrbSpec => OrbReduced(i_orb),&
                              j_OrbSpec => OrbReduced(j_orb))
                           if(i_OrbSpec%isUsed.and.j_OrbSpec%isUsed) then

                              do ilval=1,i_OrbSpec%nlang
                                 do jlval=1,j_OrbSpec%nlang

                                    if(j_OrbSpec%lang(jlval).eq.i_OrbSpec%lang(ilval)) then
                                      ! print*, i_OrbSpec%lang(ilval), j_OrbSpec%lang(jlval)

                                      add_A = 0
                                      add_B = 0 
                                      if(iPairGen==jPairGen) then
                                         select case(iPairGen)
                                         case(1)   !S^e
                                            nadd     = 1
                                            add_A(1) = 0 
                                            add_B(1) = 0 
                                            add_lval = 0
                                         case(2)   !P^o
                                            nadd     = 1
                                            add_A(1) = 0
                                            add_B(1) = 1
                                            add_lval = 2
                                         case(3,4) !P^e, D(1)
                                            nadd     = 2
                                            add_A(1) = 0
                                            add_A(2) = 2
                                            add_B(1) = 0
                                            add_B(2) = 2
                                            add_lval = 4
                                         case(5)   !D(2)
                                            nadd     = 1
                                            add_A(1) = 0
                                            add_B(1) = 2
                                            add_lval = 4
                                         case default
                                            write(LOUT,'(a)') 'ERROR! Wrong generator!'
                                         end select
                                      else
                                         select case(iPairGen)
                                         case(4,5)
                                            nadd     = 1
                                            add_A(1) = 1
                                            add_B(1) = 1
                                            add_lval = 4
                                         end select
                                      endif

                                      sum_12 = &
                                             maxOmega_PairReduced(iPairSpecG) + &
                                             maxOmega_PairReduced(jPairSpecG) + &
                                             add_lval 

                                       max_12 = &
                                              maxt_PairReduced(iPairSpecG) + &
                                              maxt_PairReduced(jPairSpecG)

                                     ! max_3 = i_OrbSpec%maxrange + j_OrbSpec%maxrange
                                       max_3 = & 
                                             i_OrbSpec%max_lrange(ilval) + &
                                             j_OrbSpec%max_lrange(jlval)

                                       ! A-part
                                       order2 = [1,0,0]
                                       call perm_part1(iThree,perm1,perm2,order2,&
                                            iPairSpec%iexp1,jPairSpec%iexp1,&
                                            iPairSpec%iexp2,jPairSpec%iexp2,&
                                            i_OrbSpec%iexp ,j_OrbSpec%iexp)
                                       do iadd=1,nadd     
                                          call add_J_Spec(iThree,add_A(iadd),maxloc(order2,dim=1),&
                                               sum_12-2*add_A(iadd),max_12,max_3)
                                       enddo 

                                       ! B-part
                                       order2 = [1,0,0]
                                       call perm_part1(iThree,perm1,perm2,order2,&
                                            iPairSpec%iexp1,jPairSpec%iexp2,&
                                            iPairSpec%iexp2,jPairSpec%iexp1,&
                                            i_OrbSpec%iexp ,j_OrbSpec%iexp)
                                       do iadd=1,nadd     
                                          call add_J_Spec(iThree,add_B(iadd),maxloc(order2,dim=1),&
                                               sum_12-2*add_B(iadd),max_12,max_3)
                                       enddo

                                    endif
                                 enddo
                              enddo

                               endif
                             end associate
                          enddo
                       enddo
                   end associate
                 enddo
              end associate
           enddo
        endif
      end associate
   enddo
enddo

call mem_dealloc(jGenNum)

end subroutine init_J1_ThreeInt

subroutine init_J2_ThreeInt(OrbReduced,PairReduced)
implicit none
type(OrbReducedData),intent(in) :: OrbReduced(:)
type(PairReducedData),intent(in) :: PairReduced(:)
integer :: prim,A_prim,B_prim,i_orb,j_orb
integer :: sum_12,max_12
integer :: iThree,perm1(3),perm2(3),order2(3)

!do B_prim=1,size(OrbReduced)
!   do A_prim=1,size(OrbReduced)
!      do prim=1,size(PairReduced)
!         associate(&
!              PairSpec  => PairReduced(prim),&
!              A_OrbSpec => OrbReduced(A_prim),&
!              B_OrbSpec => OrbReduced(B_prim))
!           if(PairSpec%isUsed.and.A_OrbSpec%isUsed.and.B_OrbSpec%isUsed) then
!
!              do j_orb=1,size(OrbReduced)
!                 do i_orb=1,size(OrbReduced)
!                    associate(&
!                         i_OrbSpec => OrbReduced(i_orb),&
!                         j_OrbSpec => OrbReduced(j_orb))
!                      if(i_OrbSpec%isUsed.and.j_OrbSpec%isUsed) then
!
!                         sum_12 = maxOmega_PairReduced(PairSpec)
!
!                         max_12 = maxt_PairReduced(PairSpec)
!
!                         order2 = [1,0,0]
!                         call perm_part1(iThree,perm1,perm2,order2,&
!                              PairSpec%iexp1,A_OrbSpec%iexp,&
!                              PairSpec%iexp2,j_OrbSpec%iexp,&
!                              i_OrbSpec%iexp,B_OrbSpec%iexp)
!                         call add_J_Spec(iThree,maxloc(order2,dim=1),&
!                              sum_12 + A_OrbSpec%maxrange + j_OrbSpec%maxrange,&
!                              max_12,&
!                              i_OrbSpec%maxrange + B_OrbSpec%maxrange)
!
!                         order2 = [1,0,0]
!                         call perm_part1(iThree,perm1,perm2,order2,&
!                              PairSpec%iexp1,j_OrbSpec%iexp,&
!                              PairSpec%iexp2,B_OrbSpec%iexp,&
!                              i_OrbSpec%iexp,A_OrbSpec%iexp)
!                         call add_J_Spec(iThree,maxloc(order2,dim=1),&
!                              sum_12 + j_OrbSpec%maxrange + B_OrbSpec%maxrange,&
!                              max_12,&
!                              i_OrbSpec%maxrange + A_OrbSpec%maxrange)
!
!                      endif
!                    end associate
!                 enddo
!              enddo
!
!           endif
!         end associate
!      enddo
!   enddo
!enddo

end subroutine init_J2_ThreeInt

subroutine init_K_ThreeInt(OrbReduced,PairReduced)
implicit none
type(OrbReducedData),intent(in) :: OrbReduced(:)
type(PairReducedData),intent(in) :: PairReduced(:)
integer :: i_prim,j_prim,i_orb,j_orb
integer :: imax,jmax,sum_TOT,max_2,max_1
integer :: iThree,perm1(3),perm2(3),order2(3)

!do j_prim=1,size(PairReduced)
!   do i_prim=1,j_prim
!      associate(&
!           iPairSpec => PairReduced(i_prim),&
!           jPairSpec => PairReduced(j_prim))
!        if(iPairSpec%isUsed.and.jPairSpec%isUsed) then
!
!           do j_orb=1,size(OrbReduced)
!              do i_orb=1,size(OrbReduced)
!                 associate(&
!                      i_OrbSpec => OrbReduced(i_orb),&
!                      j_OrbSpec => OrbReduced(j_orb))
!                   if(i_OrbSpec%isUsed.and.j_OrbSpec%isUsed) then
!
!                      sum_TOT = &
!                           maxOmega_PairReduced(iPairSpec) + &
!                           maxOmega_PairReduced(jPairSpec) + &
!                           i_OrbSpec%maxrange + &
!                           j_OrbSpec%maxrange
!
!                      max_2 = max(&
!                           maxt_PairReduced(iPairSpec),&
!                           maxt_PairReduced(jPairSpec))
!
!                      imax = maxuv_PairReduced(iPairSpec)
!                      jmax = maxuv_PairReduced(jPairSpec)
!
!                      max_1 = max(&
!                           imax + j_OrbSpec%maxrange,&
!                           imax + jmax,&
!                           i_OrbSpec%maxrange + jmax)
!
!                      order2 = [0,-1,0]
!                      call perm_part1(iThree,perm1,perm2,order2,&
!                           iPairSpec%iexp1,j_OrbSpec%iexp,&
!                           iPairSpec%iexp2,jPairSpec%iexp2,&
!                           i_OrbSpec%iexp ,jPairSpec%iexp1)
!                      call add_K_Spec(iThree,minloc(order2,dim=1),&
!                           sum_TOT,max_2,max_1)
!
!                      order2 = [0,-1,0]
!                      call perm_part1(iThree,perm1,perm2,order2,&
!                           iPairSpec%iexp1,j_OrbSpec%iexp,&
!                           iPairSpec%iexp2,jPairSpec%iexp1,&
!                           i_OrbSpec%iexp ,jPairSpec%iexp2)
!                      call add_K_Spec(iThree,minloc(order2,dim=1),&
!                           sum_TOT,max_2,max_1)
!
!                      order2 = [0,0,-1]
!                      call perm_part1(iThree,perm1,perm2,order2,&
!                           iPairSpec%iexp1,jPairSpec%iexp1,&
!                           iPairSpec%iexp2,j_OrbSpec%iexp,&
!                           i_OrbSpec%iexp ,jPairSpec%iexp2)
!                      call add_K_Spec(iThree,minloc(order2,dim=1),&
!                           sum_TOT,max_2,max_1)
!
!                      order2 = [0,0,-1]
!                      call perm_part1(iThree,perm1,perm2,order2,&
!                           iPairSpec%iexp1,jPairSpec%iexp2,&
!                           iPairSpec%iexp2,j_OrbSpec%iexp,&
!                           i_OrbSpec%iexp ,jPairSpec%iexp1)
!                      call add_K_Spec(iThree,minloc(order2,dim=1),&
!                           sum_TOT,max_2,max_1)
!
!                   endif
!                 end associate
!              enddo
!           enddo
!
!        endif
!      end associate
!   enddo
!enddo

end subroutine init_K_ThreeInt

subroutine add_J_Spec(iThree,lambda,ij,sum_ij,max_ij,max_k)
implicit none
integer,intent(in) :: iThree,ij
integer,intent(in) :: lambda
integer,intent(in) :: sum_ij,max_ij,max_k
integer :: ival,lval

associate(Three => ThreeInt(iThree))

  if(.not.Three%isUsed_J) then

     allocate(Three%J_Spec(3))

! hapka: old_drake
!     Three%J_Spec(1) = J_SpecData(.false.,'12','3',0,0,0)
!     Three%J_Spec(2) = J_SpecData(.false.,'13','2',0,0,0)
!     Three%J_Spec(3) = J_SpecData(.false.,'23','1',0,0,0)

     ! initialize
     do ival=1,size(Three%J_Spec)
        associate(Spec => Three%J_Spec(ival))
         
          ! init JSpec
          Spec%isUsed = .false.
          select case(ival)
          case(1)
             Spec%symbol_ij = '12'
             Spec%symbol_k  = '3'
          case(2)
             Spec%symbol_ij = '13'
             Spec%symbol_k  = '2'
          case(3)
             Spec%symbol_ij = '23'
             Spec%symbol_k  = '1' 
          end select

          ! init JSpecLam
          do lval=0,2
             associate(LSpec => Spec%J_SpecLam(lval))

               LSpec%isUsed = .false.
               LSpec%sum_ij = 0 
               LSpec%max_ij = 0 
               LSpec%max_k  = 0 

             end associate
          enddo

        end associate
     enddo

     Three%isUsed_J = .true.

  endif

  associate(Spec => Three%J_Spec(ij))
    associate(LSpec => Spec%J_SpecLam(lambda))

      LSpec%sum_ij = max(LSpec%sum_ij,sum_ij)
      LSpec%max_ij = max(LSpec%max_ij,max_ij)
      LSpec%max_k  = max(LSpec%max_k ,max_k)

      LSpec%isUsed = .true. 
      
    end associate

    Spec%isUsed = .true.
  end associate

end associate

end subroutine add_J_Spec

subroutine add_K_Spec(iThree,ij,sum_TOT,max_2,max_1)
implicit none
integer,intent(in) :: iThree,ij
integer,intent(in) :: sum_TOT,max_2,max_1

associate(Three => ThreeInt(iThree))

  if(.not.Three%isUsed_K) then

     allocate(Three%K_Spec(3))

     Three%K_Spec(1) = K_SpecData(.false.,'12',0,0,0)
     Three%K_Spec(2) = K_SpecData(.false.,'13',0,0,0)
     Three%K_Spec(3) = K_SpecData(.false.,'23',0,0,0)

     Three%isUsed_K = .true.

  endif

  associate(Spec => Three%K_Spec(ij))

    Spec%sum_TOT = max(Spec%sum_TOT,sum_TOT)
    Spec%max_2   = max(Spec%max_2  ,max_2)
    Spec%max_1   = max(Spec%max_1  ,max_1)

    Spec%isUsed = .true.

  end associate

end associate

end subroutine add_K_Spec

subroutine perm_part1(iThree,perm1,perm2,order2,iexp,jexp,kexp,lexp,mexp,nexp)
use misc, only : swap
implicit none
integer,intent(out) :: iThree,perm1(3),perm2(3)
integer,intent(inout) :: order2(3)
integer,intent(in) :: iexp,jexp,kexp,lexp,mexp,nexp
logical :: ijorder,klorder,mnorder
integer :: exp1,exp2,exp3
logical :: valid

ijorder = (iexp<=jexp)
klorder = (kexp<=lexp)
mnorder = (mexp<=nexp)

if(ijorder) then
   exp1 = iexp + jexp*(jexp-1)/2
else
   exp1 = jexp + iexp*(iexp-1)/2
endif

if(klorder) then
   exp2 = kexp + lexp*(lexp-1)/2
else
   exp2 = lexp + kexp*(kexp-1)/2
endif

if(mnorder) then
   exp3 = mexp + nexp*(nexp-1)/2
else
   exp3 = nexp + mexp*(mexp-1)/2
endif

perm1 = [1,2,3]
perm2 = [1,2,3]

if(exp1>exp2) then
   call swap(exp1,exp2)
   call swap(perm1(1),perm1(2))
   call swap(perm2(2),perm2(3))
   call swap(order2(2),order2(3))
endif

if(exp2>exp3) then
   call swap(exp2,exp3)
   call swap(perm1(2),perm1(3))
   call swap(perm2(1),perm2(2))
   call swap(order2(1),order2(2))
endif

if(exp1>exp2) then
   call swap(exp1,exp2)
   call swap(perm1(1),perm1(2))
   call swap(perm2(2),perm2(3))
   call swap(order2(2),order2(3))
endif

iThree = exp1 + exp2*(exp2-1)/2 + (exp3+1)*exp3*(exp3-1)/6

associate(Three => ThreeInt(iThree))

  valid = .true.
  call check_exp(perm1(1),Three%iexp,Three%jexp)
  call check_exp(perm1(2),Three%kexp,Three%lexp)
  call check_exp(perm1(3),Three%mexp,Three%nexp)
  if(.not.valid) then
     write(LOUT,'(a)') 'ERROR!!! Wrong entry chosen in perm_part1!'
     stop
  endif

  if(Three%sameANY) then

     if(Three%sameALL) then

        if(order2(1)<order2(2)) then
           call swap(perm1(2),perm1(3))
           call swap(perm2(1),perm2(2))
           call swap(order2(1),order2(2))
        endif

        if(order2(2)<order2(3)) then
           call swap(perm1(1),perm1(2))
           call swap(perm2(2),perm2(3))
           call swap(order2(2),order2(3))
        endif

        if(order2(1)<order2(2)) then
           call swap(perm1(2),perm1(3))
           call swap(perm2(1),perm2(2))
           call swap(order2(1),order2(2))
        endif

     elseif(Three%same12) then

        if(order2(2)<order2(3)) then
           call swap(perm1(1),perm1(2))
           call swap(perm2(2),perm2(3))
           call swap(order2(2),order2(3))
        endif

     elseif(Three%same23) then

        if(order2(1)<order2(2)) then
           call swap(perm1(2),perm1(3))
           call swap(perm2(1),perm2(2))
           call swap(order2(1),order2(2))
        endif

     else

        write(LOUT,'(a)') 'ERROR!!! &
             &Incorrect determination of identical exponents in perm_part1!'
        stop

     endif

  endif

end associate

contains

  subroutine check_exp(pos,check_1,check_2)
  implicit none
  integer,intent(in) :: pos
  integer,intent(in) :: check_1,check_2

  select case(pos)
  case(1)
     if(ijorder) then
        valid = valid.and.(iexp==check_1).and.(jexp==check_2)
     else
        valid = valid.and.(jexp==check_1).and.(iexp==check_2)
     endif
  case(2)
     if(klorder) then
        valid = valid.and.(kexp==check_1).and.(lexp==check_2)
     else
        valid = valid.and.(lexp==check_1).and.(kexp==check_2)
     endif
  case(3)
     if(mnorder) then
        valid = valid.and.(mexp==check_1).and.(nexp==check_2)
     else
        valid = valid.and.(nexp==check_1).and.(mexp==check_2)
     endif
  case default
     write(LOUT,'(a)') 'ERROR!!! Incorrect position in perm_part1!'
     stop
  end select

  end subroutine check_exp

end subroutine perm_part1

subroutine perm_part2(perm1,perm2,order2)
use misc, only : swap
implicit none
integer,intent(inout) :: perm1(3),perm2(3),order2(3)

if(order2(1)<order2(2)) then
   call swap(perm1(2),perm1(3))
   call swap(perm2(1),perm2(2))
   call swap(order2(1),order2(2))
endif

if(order2(2)<order2(3)) then
   call swap(perm1(1),perm1(2))
   call swap(perm2(2),perm2(3))
   call swap(order2(2),order2(3))
endif

if(order2(1)<order2(2)) then
   call swap(perm1(2),perm1(3))
   call swap(perm2(1),perm2(2))
   call swap(order2(1),order2(2))
endif

end subroutine perm_part2

function permuted_alpha(Three,perm1) result(alpha)
implicit none
real(prec) :: alpha(3)
type(ThreeIntData),intent(in) :: Three
integer,intent(in) :: perm1(3)
integer :: i

do i=1,3
   select case(perm1(i))
   case(1)
      alpha(i) = Three%ijalpha
   case(2)
      alpha(i) = Three%klalpha
   case(3)
      alpha(i) = Three%mnalpha
   case default
      write(LOUT,'(a)') 'ERROR!!! Incorrect perm1 in permuted_alpha!'
      stop
   end select
enddo

end function permuted_alpha

function permuted_same12(Three,perm2) result(same12)
implicit none
logical :: same12
type(ThreeIntData),intent(in) :: Three
integer,intent(in) :: perm2(3)

select case(perm2(1))
case(1)
   same12 = Three%same12
case(2)
   same12 = Three%sameALL
case(3)
   same12 = Three%same23
case default
   write(LOUT,'(a)') 'ERROR!!! Incorrect perm2 in permuted_same12!'
   stop
end select

end function permuted_same12

subroutine create_ThreeInt
use misc, only : swap
implicit none
integer :: iThree,ij
integer :: iunit
integer :: perm1(3),perm2(3),order2(3)
real(prec) :: alpha(3)
character(name_length) :: name_J1,name_J2,name_K0,name_K1

call execute_command_line('rm -f '//core_J//'* '//core_K//'*')

!select case(int3_source)
!case('L')
!
!   do iThree=1,size(ThreeInt)
!      associate(Three => ThreeInt(iThree))
!
!        if(Three%isUsed_J) then
!           do ij=1,3
!              associate(Spec => Three%J_Spec(ij))
!                if(Spec%isUsed) then
!
!                   name_J1 = prepare_name(core_J,iThree,ij,1)
!                   name_J2 = prepare_name(core_J,iThree,ij,2)
!
!                   perm1 = [1,2,3]
!                   perm2 = [1,2,3]
!                   order2 = 0
!                   order2(ij) = 1
!                   call perm_part2(perm1,perm2,order2)
!                   alpha = permuted_alpha(Three,perm1)
!
!                   open(newunit=iunit,file=infile)
!                   write(iunit,'(3es23.15,2x,a1,2x,3i5)') alpha,'J',&
!                        Spec%sum_ij,Spec%max_ij,Spec%max_k
!                   close(iunit)
!                   call execute_command_line(&
!                        'bash '//script//' L '//infile//' J '//trim(name_J1))
!
!                   if(permuted_same12(Three,perm2)) then
!
!                      call execute_command_line('cp '//&
!                           trim(name_J1)//' '//trim(name_J2))
!
!                   else
!
!                      call swap(alpha(1),alpha(2))
!
!                      open(newunit=iunit,file=infile)
!                      write(iunit,'(3es23.15,2x,a1,2x,3i5)') alpha,'J',&
!                           Spec%sum_ij,Spec%max_ij,Spec%max_k
!                      close(iunit)
!                      call execute_command_line(&
!                           'bash '//script//' L '//infile//' J '//trim(name_J2))
!
!                   endif
!
!                endif
!              end associate
!           enddo
!        endif
!
!        if(Three%isUsed_K) then
!           do ij=1,3
!              associate(Spec => Three%K_Spec(ij))
!                if(Spec%isUsed) then
!
!                   name_K0 = prepare_name(core_K,iThree,ij,0)
!                   name_K1 = prepare_name(core_K,iThree,ij,1)
!
!                   perm1 = [1,2,3]
!                   perm2 = [1,2,3]
!                   order2 = 0
!                   order2(ij) = -1
!                   call perm_part2(perm1,perm2,order2)
!                   alpha = permuted_alpha(Three,perm1)
!
!                   open(newunit=iunit,file=infile)
!                   write(iunit,'(3es23.15,2x,a1,2x,3i5)') alpha,'K',&
!                        Spec%sum_TOT,Spec%max_2,Spec%max_1
!                   close(iunit)
!                   call execute_command_line(&
!                        'bash '//script//' L '//infile//' K '//&
!                        trim(name_K0)//' '//trim(name_K1))
!
!                endif
!              end associate
!           enddo
!        endif
!
!      end associate
!   enddo
!
!case('P','Q')
!
!   do iThree=1,size(ThreeInt)
!      associate(Three => ThreeInt(iThree))
!
!        open(newunit=iunit,file=infile)
!        write(iunit,'(3es23.15)') Three%ijalpha,Three%klalpha,Three%mnalpha
!
!        if(Three%isUsed_J) then
!           write(iunit,'(i1)') count(Three%J_Spec(:)%isUsed)
!           do ij=1,3
!              associate(Spec => Three%J_Spec(ij))
!                if(Spec%isUsed) then
!
!                   name_J1 = prepare_name(core_J,iThree,ij,1)
!                   name_J2 = prepare_name(core_J,iThree,ij,2)
!
!                   write(iunit,'(a2,2x,3i5,2(2x,a))') Spec%symbol_ij,&
!                        Spec%sum_ij,Spec%max_ij,Spec%max_k,&
!                        trim(name_J1),trim(name_J2)
!
!                endif
!              end associate
!           enddo
!        else
!           write(iunit,'(i1)') 0
!        endif
!
!        if(Three%isUsed_K) then
!           write(iunit,'(i1)') count(Three%K_Spec(:)%isUsed)
!           do ij=1,3
!              associate(Spec => Three%K_Spec(ij))
!                if(Spec%isUsed) then
!
!                   name_K0 = prepare_name(core_K,iThree,ij,0)
!                   name_K1 = prepare_name(core_K,iThree,ij,1)
!
!                   write(iunit,'(a2,2x,3i5,2(2x,a))') Spec%symbol_ij,&
!                        Spec%sum_TOT,Spec%max_2,Spec%max_1,&
!                        trim(name_K0),trim(name_K1)
!
!                endif
!              end associate
!           enddo
!        else
!           write(iunit,'(i1)') 0
!        endif
!
!        close(iunit)
!        call execute_command_line(&
!             'bash '//script//' '//int3_source//' '//infile)
!
!      end associate
!   enddo
!
!case default
!
!   write(LOUT,'(a)') 'ERROR!!! Unrecognizable source of 3-el integrals!'
!   stop
!
!end select

end subroutine create_ThreeInt

subroutine read_ThreeInt_J(int3val,&
     iThree,ij,extra,sum_12,max_12,max_3,perm1,perm2)
use precision, only : shortint,dble
implicit none
type(int3valData),allocatable :: int3val(:)
integer,intent(in) :: iThree,ij,extra
integer,intent(in) :: sum_12,max_12,max_3
integer,intent(in) :: perm1(3),perm2(3)
integer :: t,range_uv
integer :: extra_name
character(name_length) :: name_J
integer :: iunit,stat
integer(shortint) :: idx(4)
real(dble) :: val2
real(prec) :: val4
integer :: i1,i2,i3,i12
procedure(),pointer :: permute_indices

!if(.not.ThreeInt(iThree)%isUsed_J) then
!   write(LOUT,'(a)') 'ERROR!!! &
!        &An entry incorrectly predicted as not necessary for J integrals!'
!   stop
!endif
!associate(Spec => ThreeInt(iThree)%J_Spec(ij))
!  if(.not.Spec%isUsed) then
!     write(LOUT,'(a)') 'ERROR!!! &
!          &An ordering incorrectly predicted as not necessary for J integrals!'
!     stop
!  endif
!  if(sum_12>Spec%sum_ij.or.max_12>Spec%max_ij.or.max_3>Spec%max_k) then
!     write(LOUT,'(a)') 'ERROR!!! &
!          &Ranges of 3-electron J integrals have been predicted incorrectly!'
!     stop
!  endif
!end associate
!if(perm1(3)/=3.or.perm2(1)/=1.or.&
!     (perm1(1)<perm1(2)).neqv.(perm2(2)<perm2(3))) then
!   write(LOUT,'(a)') 'ERROR!!! Unexpected permutation for J integrals!'
!   stop
!endif

select case(extra)
case(1)

   extra_name = merge(1,2,perm2(2)<perm2(3))
   permute_indices => permute12

case(2)

   extra_name = merge(2,1,perm2(2)<perm2(3))
   permute_indices => permute21

case default

   write(LOUT,'(a)') 'ERROR!!! Incorrect extra parameter for J integrals!'
   stop

end select

allocate(int3val(0:max_12))

do t=0,max_12
   range_uv = sum_12 - t
   call mem_alloc(int3val(t)%elms,&
        max_3    - smallest_uv + 1,&
        range_uv - smallest_uv + 1,&
        range_uv - smallest_uv + 1)
   int3val(t)%elms = 0._prec
enddo

name_J = prepare_name(core_J,iThree,ij,extra_name)

open(newunit=iunit,file=name_J,access='stream')
select case(int3_source)
case('L')

   do
      read(iunit,iostat=stat) idx,val2
      if(stat/=0) exit
      if(idx(4)>max_12) cycle
      if(idx(3)>max_3) cycle
      if(any(idx(1:2)>sum_12-idx(4))) cycle

      call permute_indices
      int3val(i12)%elms(&
           offset_uv + i3, &
           offset_uv + i1, &
           offset_uv + i2) = val2*val_modifier

   enddo

case('P','Q')

   do
      read(iunit,iostat=stat) idx,val4
      if(stat/=0) exit
      if(idx(4)>max_12) cycle
      if(idx(3)>max_3) cycle
      if(any(idx(1:2)>sum_12-idx(4))) cycle

      call permute_indices
      int3val(i12)%elms(&
           offset_uv + i3, &
           offset_uv + i1, &
           offset_uv + i2) = val4

   enddo

case default

   write(LOUT,'(a)') 'ERROR!!! Unrecognizable source of 3-el integrals!'
   stop

end select
close(iunit)

contains

  subroutine permute12
  implicit none
  i1  = idx(1)
  i2  = idx(2)
  i3  = idx(3)
  i12 = idx(4)
  end subroutine permute12

  subroutine permute21
  implicit none
  i2  = idx(1)
  i1  = idx(2)
  i3  = idx(3)
  i12 = idx(4)
  end subroutine permute21

end subroutine read_ThreeInt_J

subroutine free_ThreeInt_J(int3val)
implicit none
type(int3valData),allocatable :: int3val(:)
integer :: t

do t=lbound(int3val,dim=1),ubound(int3val,dim=1)
   call mem_dealloc(int3val(t)%elms)
enddo

deallocate(int3val)

end subroutine free_ThreeInt_J

subroutine read_ThreeInt_K(int3val,&
     iThree,ij,extra,sum_TOT,max_2,max_1,perm1,perm2)
use precision, only : shortint,dble
implicit none
type(int3valData),allocatable :: int3val(:,:)
integer,intent(in) :: iThree,ij,extra
integer,intent(in) :: sum_TOT,max_2,max_1
integer,intent(in) :: perm1(3),perm2(3)
integer :: start_2,range_1
character(name_length) :: name_K
integer :: iunit,stat
integer(shortint) :: idx(5)
real(dble) :: val2
real(prec) :: val4
integer :: i1,i2,i3,i12,ix3
character(3) :: permute_string
procedure(),pointer :: permute_indices_1,permute_indices_2

if(.not.ThreeInt(iThree)%isUsed_K) then
   write(LOUT,'(a)') 'ERROR!!! &
        &An entry incorrectly predicted as not necessary for K integrals!'
   stop
endif
associate(Spec => ThreeInt(iThree)%K_Spec(ij))
  if(.not.Spec%isUsed) then
     write(LOUT,'(a)') 'ERROR!!! &
          &An ordering incorrectly predicted as not necessary for K integrals!'
     stop
  endif
  if(sum_TOT>Spec%sum_TOT.or.max_2>Spec%max_2.or.max_1>Spec%max_1) then
     write(LOUT,'(a)') 'ERROR!!! &
          &Ranges of 3-electron K integrals have been predicted incorrectly!'
     stop
  endif
end associate
if(perm2(3)==1) then
   write(LOUT,'(a)') 'ERROR!!! Unexpected permutation for K integrals!'
   stop
endif

select case(extra)
case(0)

   start_2 = -1

case(1)

   start_2 = 0

case default

   write(LOUT,'(a)') 'ERROR!!! Incorrect extra parameter for K integrals!'
   stop

end select

if(perm2(1)==1) then
   permute_indices_2 => permute12
else
   permute_indices_2 => permute21
endif

write(permute_string,'(3i1)') perm1
select case(permute_string)
case('123')
   permute_indices_1 => permute123
case('213')
   permute_indices_1 => permute213
case('132')
   permute_indices_1 => permute132
case('312')
   permute_indices_1 => permute312
case('231')
   permute_indices_1 => permute231
case('321')
   permute_indices_1 => permute321
case default
   write(LOUT,'(a)') 'ERROR!!! Impossible permutation of K integrals!'
   stop
end select

allocate(int3val(start_2:max_2,start_2:max_2))

range_1 = max_1 - smallest_uv + 1
do ix3=start_2,max_2
   do i12=start_2,max_2
      call mem_alloc(int3val(i12,ix3)%elms,range_1,range_1,range_1)
      int3val(i12,ix3)%elms = 0._prec
   enddo
enddo

name_K = prepare_name(core_K,iThree,ij,extra)

open(newunit=iunit,file=name_K,access='stream')
select case(int3_source)
case('L')

   do
      read(iunit,iostat=stat) idx,val2
      if(stat/=0) exit
      if(any(idx(4:5)>max_2)) cycle
      if(any(idx(1:3)>max_1)) cycle

      call permute_indices_2
      call permute_indices_1
      int3val(i12,ix3)%elms(&
           offset_uv + i1, &
           offset_uv + i2, &
           offset_uv + i3) = val2*val_modifier

   enddo

case('P','Q')

   do
      read(iunit,iostat=stat) idx,val4
      if(stat/=0) exit
      if(any(idx(4:5)>max_2)) cycle
      if(any(idx(1:3)>max_1)) cycle

      call permute_indices_2
      call permute_indices_1
      int3val(i12,ix3)%elms(&
           offset_uv + i1, &
           offset_uv + i2, &
           offset_uv + i3) = val4

   enddo

case default

   write(LOUT,'(a)') 'ERROR!!! Unrecognizable source of 3-el integrals!'
   stop

end select
close(iunit)

contains

  subroutine permute12
  implicit none
  i12 = idx(4)
  ix3 = idx(5)
  end subroutine permute12

  subroutine permute21
  implicit none
  ix3 = idx(4)
  i12 = idx(5)
  end subroutine permute21

  subroutine permute123
  implicit none
  i1  = idx(1)
  i2  = idx(2)
  i3  = idx(3)
  end subroutine permute123

  subroutine permute213
  implicit none
  i2  = idx(1)
  i1  = idx(2)
  i3  = idx(3)
  end subroutine permute213

  subroutine permute132
  implicit none
  i1  = idx(1)
  i3  = idx(2)
  i2  = idx(3)
  end subroutine permute132

  subroutine permute312
  implicit none
  i3  = idx(1)
  i1  = idx(2)
  i2  = idx(3)
  end subroutine permute312

  subroutine permute231
  implicit none
  i2  = idx(1)
  i3  = idx(2)
  i1  = idx(3)
  end subroutine permute231

  subroutine permute321
  implicit none
  i3  = idx(1)
  i2  = idx(2)
  i1  = idx(3)
  end subroutine permute321

end subroutine read_ThreeInt_K

subroutine free_ThreeInt_K(int3val)
implicit none
type(int3valData),allocatable :: int3val(:,:)
integer :: i12,ix3

do ix3=lbound(int3val,dim=2),ubound(int3val,dim=2)
   do i12=lbound(int3val,dim=1),ubound(int3val,dim=1)
      call mem_dealloc(int3val(i12,ix3)%elms)
   enddo
enddo

deallocate(int3val)

end subroutine free_ThreeInt_K

function prepare_name(core,iThree,ij,extra) result(name)
implicit none
character(name_length) :: name
character(*),intent(in) :: core
integer,intent(in) :: iThree,ij,extra

write(name,'(a,"_",i6.6,"_",i1,"_",i1)') trim(core),iThree,ij,extra

end function prepare_name

subroutine free_ThreeInt
implicit none
integer :: iThree

do iThree=1,size(ThreeInt)
   associate(Three => ThreeInt(iThree))

!     if(Three%isUsed_K) deallocate(Three%K_Spec)
     if(Three%isUsed_J) deallocate(Three%J_Spec)

   end associate
enddo

deallocate(ThreeInt)

end subroutine free_ThreeInt

subroutine print_ThreeInt(LPRINT)
implicit none
integer,intent(in) :: LPRINT
integer :: iThree,ij
integer :: ilam

!if(LPRINT>=10) then

   write(LOUT,'()')
   write(LOUT,'(5x,a)') '--- CC three-electron integrals ---'

   write(LOUT,'()')
   write(LOUT,'(3x,a,2x,6(2x,a),3x,3(6x,a,2x),1x,a)') &
        'no.','i','j','k','l','m','n','ijalpha','klalpha','mnalpha','same'
   do iThree=1,size(ThreeInt)
      associate(Three => ThreeInt(iThree))
        write(LOUT,'(1x,i5,a)',advance='no') iThree,' :'
        if(Three%isUsed_J.or.Three%isUsed_K) then
           write(LOUT,'(6i3,3x,3f15.8,1x,4a1)') &
                Three%iexp,Three%jexp,&
                Three%kexp,Three%lexp,&
                Three%mexp,Three%nexp,&
                Three%ijalpha,Three%klalpha,Three%mnalpha,&
                merge('*',' ',Three%sameANY),&
                merge('*',' ',Three%sameALL),&
                merge('*',' ',Three%same12),&
                merge('*',' ',Three%same23)
           if(Three%isUsed_J) then
              write(LOUT,'(8x,a,i1)') 'J: ',count(Three%J_Spec(:)%isUsed)
              do ij=1,3
                 associate(Spec => Three%J_Spec(ij))
                   if(Spec%isUsed) then
!                      write(LOUT,'(12x,2a,2(2x,3a,i3),2x,3a,i3)') &
                      write(LOUT,'(12x,2a)') &
                           'symbol:',Spec%symbol_ij
                      do ilam=0,2
                         associate(LSpec => Spec%J_SpecLam(ilam))
                           if(LSpec%isUsed) then
                           write(LOUT,'(12x,a,i3,3(2x,3a,i3))') &
                           'ilam = ',ilam, & 
                           'sum_',Spec%symbol_ij,'  = ',LSpec%sum_ij,&
                           'max_',Spec%symbol_ij,' = ',LSpec%max_ij,&
                           'max_',Spec%symbol_k,' = ',LSpec%max_k

!                      write(LOUT,'(12x,2a,2(2x,3a,i3),2x,3a,i3)') &
!                           'symbol:',Spec%symbol_ij ,&
!                           'sum_',Spec%symbol_ij,'  = ',Spec%sum_ij,&
!                           'max_',Spec%symbol_ij,' = ',Spec%max_ij,&
!                           'max_',Spec%symbol_k,' = ',Spec%max_k
                           endif
                         end associate
                      enddo
                   endif
                 end associate
              enddo
           endif
!           if(Three%isUsed_K) then
!              write(LOUT,'(8x,a,i1)') 'K: ',count(Three%K_Spec(:)%isUsed)
!              do ij=1,3
!                 associate(Spec => Three%K_Spec(ij))
!                   if(Spec%isUsed) then
!                      write(LOUT,'(12x,2a,3(2x,a,i3))') &
!                           'symbol:',Spec%symbol_ij,&
!                           'sum_TOT = ',Spec%sum_TOT,&
!                           'max_ij = ',Spec%max_2,&
!                           'max_i = ',Spec%max_1
!                   endif
!                 end associate
!              enddo
!           endif
           write(LOUT,'()')
        else
           write(LOUT,'(a)') 'NOT USED'
        endif
      end associate
   enddo

!endif

end subroutine print_ThreeInt

end module CCint3
