module indices
implicit none

private
public init_indices,free_indices
public MAXF,MAXG
public idx6,idx3

integer,save :: MAXF,MAXG

integer,allocatable,save :: idx6_data(:,:)
integer,allocatable,save :: idx3_data(:,:)

contains

subroutine init_indices(sm)
implicit none
integer,intent(in) :: sm
integer(8) :: SHELL
integer(8) :: MAXF_long,MAXG_long

SHELL = sm

MAXF_long = ((1+SHELL)*(2+SHELL)*(3+SHELL)*(4+SHELL)*(5+SHELL)*(6+SHELL))/720
MAXG_long = ((7+SHELL)*(8+SHELL)*(9+SHELL))/6

MAXF = int(MAXF_long)
MAXG = int(MAXG_long)

call init_idx6(SHELL)
call init_idx3(SHELL)

end subroutine init_indices

subroutine free_indices
implicit none

call free_idx3
call free_idx6

end subroutine free_indices

integer function idx6(n1,n2,n3,n4,n5,n6)
implicit none
integer,intent(in) :: n1,n2,n3,n4,n5,n6
integer :: nx

!!$if((n1<0).or.(n2<0).or.(n3<0).or.(n4<0).or.(n5<0).or.(n6<0)) then
!!$   idx6 = 0
!!$else
   idx6 = 1
   nx = n1
   idx6 = idx6 + idx6_data(nx,1)
   nx = nx + n2
   idx6 = idx6 + idx6_data(nx,2)
   nx = nx + n3
   idx6 = idx6 + idx6_data(nx,3)
   nx = nx + n4
   idx6 = idx6 + idx6_data(nx,4)
   nx = nx + n5
   idx6 = idx6 + idx6_data(nx,5)
   nx = nx + n6
   idx6 = idx6 + idx6_data(nx,6)
!!$   if(idx6>MAXF) then
!!$      write(*,*) 'ZAKRES MAXF:',MAXF,idx6,n1,n2,n3,n4,n5,n6
!!$      stop
!!$   endif
!!$endif

end function idx6

integer function idx3(n1,n2,n3)
implicit none
integer :: n1,n2,n3
integer :: nx

!!$if((n1<-1).or.(n2<-1).or.(n3<-1))then
!!$   idx3 = 0
!!$else
   idx3 = 1
   nx = n1+1
   idx3 = idx3 + idx3_data(nx,1)
   nx = nx + n2+1
   idx3 = idx3 + idx3_data(nx,2)
   nx = nx + n3+1
   idx3 = idx3 + idx3_data(nx,3)
!!$   if(idx3>MAXG) then
!!$      write(*,*) 'ZAKRES MAXG:',MAXG,idx3,n1,n2,n3
!!$      stop      
!!$   endif
!!$endif

end function idx3

subroutine init_idx6(SHELL)
implicit none
integer(8),intent(in) :: SHELL
integer(8) :: nx,i,tmp

allocate(idx6_data(0:SHELL,6))

do nx=0,SHELL
   tmp = 1
   do i=1,6
      tmp = (tmp*(nx + i - 1))/i
      idx6_data(nx,i) = int(tmp)
   enddo
enddo

end subroutine init_idx6

subroutine free_idx6
implicit none

deallocate(idx6_data)

end subroutine free_idx6

subroutine init_idx3(SHELL)
implicit none
integer(8),intent(in) :: SHELL
integer(8) :: nx,i,tmp

allocate(idx3_data(0:SHELL+6,3))

do nx=0,SHELL+6
   tmp = 1
   do i=1,3
      tmp = (tmp*(nx + i - 1))/i
      idx3_data(nx,i) = int(tmp)
   enddo
enddo

end subroutine init_idx3

subroutine free_idx3
implicit none

deallocate(idx3_data)

end subroutine free_idx3

end module indices
