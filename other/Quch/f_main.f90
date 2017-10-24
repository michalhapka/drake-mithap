subroutine f_main(c_length,c_filename) bind(c,name="f_main")
use iso_c_binding, only : c_int,c_char
use qdmodule
use indices
use pachucki
implicit none
integer(c_int),intent(in),value :: c_length
character(kind=c_char),intent(in) :: c_filename(*)
integer,parameter :: shortint = 2
integer,parameter :: lengthMAX = 40, name_length = 20
character(lengthMAX) :: filename
integer :: iunit,status
type SpecData
character(2) :: symbol
integer :: range1,range2,range3
character(name_length) :: name_1,name_2
end type SpecData
real(16) :: a,b,c
type(qd_real) :: qd_a,qd_b,qd_c
integer :: nJ,nK,i
type(SpecData) :: SpecJ(3),SpecK(3)
integer :: sm,max_1,max_2
type(qd_real),allocatable :: f(:)
integer(shortint) :: i1,i2,i3,i12,i13
real(16) :: val1,val2,error
logical :: test
integer(4) :: old_cw

call f_fpu_fix_start(old_cw)

test = .false.

if(c_length>lengthMAX) then
   write(*,'(a)') 'PROBLEM WITH THE INPUT-FILE NAME! QUITTING!'
   stop
endif
filename = transfer(c_filename(1:c_length),filename)

open(newunit=iunit,file=trim(filename))
read(iunit,*) a,b,c
read(iunit,*) nJ
do i=1,nJ
   read(iunit,*) SpecJ(i)
enddo
read(iunit,*) nK
do i=1,nK
   read(iunit,*) SpecK(i)
enddo
close(iunit)

sm = 0
max_2 = 0
max_1 = 0
do i=1,nJ
   sm = max(sm,SpecJ(i)%range1+SpecJ(i)%range3)
   max_2 = max(max_2,SpecJ(i)%range2)
   max_1 = max(max_1,SpecJ(i)%range1,SpecJ(i)%range3)
enddo
do i=1,nK
   sm = max(sm,SpecK(i)%range1)
   max_2 = max(max_2,SpecK(i)%range2)
   max_1 = max(max_1,SpecK(i)%range3)
enddo
sm = sm + 6
max_2 = max_2 + 1
max_1 = max_1 + 1

call init_indices(sm)
allocate(f(0:MAXF))

qd_a = to_qd(a)
qd_b = to_qd(b)
qd_c = to_qd(c)

call p_f(f,qd_a,qd_b,qd_c,sm,max_2,max_1)

do i=1,nJ
   associate(Spec => SpecJ(i))
     select case(Spec%symbol)
     case('12')

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_1),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( 1, 0, i12+1, i1+1, i2+1, i3+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_1),error
        endif

        open(newunit=iunit,file=trim(Spec%name_1),access='stream')
        do i12=0,int(Spec%range2,shortint)
           do i2=0,int(Spec%range1,shortint)-i12
              do i1=0,int(Spec%range1,shortint)-i12-i2
                 do i3=0,int(Spec%range3,shortint)

                    write(iunit) i1,i2,i3,i12,&
                         to_16(f(idx6( 1, 0, i12+1, i1+1, i2+1, i3+1)))

                 enddo
              enddo
           enddo
        enddo
        close(iunit)

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_2),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( 0, 1, i12+1, i2+1, i1+1, i3+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_2),error
        endif

        open(newunit=iunit,file=trim(Spec%name_2),access='stream')
        do i12=0,int(Spec%range2,shortint)
           do i2=0,int(Spec%range1,shortint)-i12
              do i1=0,int(Spec%range1,shortint)-i12-i2
                 do i3=0,int(Spec%range3,shortint)

                    write(iunit) i1,i2,i3,i12,&
                         to_16(f(idx6( 0, 1, i12+1, i2+1, i1+1, i3+1)))

                 enddo
              enddo
           enddo
        enddo
        close(iunit)

     case('13')

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_1),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( 1, i12+1, 0, i1+1, i3+1, i2+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_1),error
        endif

        open(newunit=iunit,file=trim(Spec%name_1),access='stream')
        do i12=0,int(Spec%range2,shortint)
           do i2=0,int(Spec%range1,shortint)-i12
              do i1=0,int(Spec%range1,shortint)-i12-i2
                 do i3=0,int(Spec%range3,shortint)

                    write(iunit) i1,i2,i3,i12,&
                         to_16(f(idx6( 1, i12+1, 0, i1+1, i3+1, i2+1)))

                 enddo
              enddo
           enddo
        enddo
        close(iunit)

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_2),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( 0, i12+1, 1, i2+1, i3+1, i1+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_2),error
        endif

        open(newunit=iunit,file=trim(Spec%name_2),access='stream')
        do i12=0,int(Spec%range2,shortint)
           do i2=0,int(Spec%range1,shortint)-i12
              do i1=0,int(Spec%range1,shortint)-i12-i2
                 do i3=0,int(Spec%range3,shortint)

                    write(iunit) i1,i2,i3,i12,&
                         to_16(f(idx6( 0, i12+1, 1, i2+1, i3+1, i1+1)))

                 enddo
              enddo
           enddo
        enddo
        close(iunit)

     case('23')

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_1),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( i12+1, 1, 0, i3+1, i1+1, i2+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_1),error
        endif

        open(newunit=iunit,file=trim(Spec%name_1),access='stream')
        do i12=0,int(Spec%range2,shortint)
           do i2=0,int(Spec%range1,shortint)-i12
              do i1=0,int(Spec%range1,shortint)-i12-i2
                 do i3=0,int(Spec%range3,shortint)

                    write(iunit) i1,i2,i3,i12,&
                         to_16(f(idx6( i12+1, 1, 0, i3+1, i1+1, i2+1)))

                 enddo
              enddo
           enddo
        enddo
        close(iunit)

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_2),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( i12+1, 0, 1, i3+1, i2+1, i1+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_2),error
        endif

        open(newunit=iunit,file=trim(Spec%name_2),access='stream')
        do i12=0,int(Spec%range2,shortint)
           do i2=0,int(Spec%range1,shortint)-i12
              do i1=0,int(Spec%range1,shortint)-i12-i2
                 do i3=0,int(Spec%range3,shortint)

                    write(iunit) i1,i2,i3,i12,&
                         to_16(f(idx6( i12+1, 0, 1, i3+1, i2+1, i1+1)))

                 enddo
              enddo
           enddo
        enddo
        close(iunit)

     case default
        write(*,'(a)') 'UNKNOWN SYMBOL! QUITTING!'
        stop
     end select
   end associate
enddo

do i=1,nK
   associate(Spec => SpecK(i))
     select case(Spec%symbol)
     case('12')

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_1),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,i13,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( i13+1, i12+1, 1, i2+1, i3+1, i1+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,i13,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_1),error
        endif

        open(newunit=iunit,file=trim(Spec%name_1),access='stream')
        do i13=-1,int(Spec%range2,shortint)
           do i12=-1,int(Spec%range2,shortint)
              if(i12==-1.and.i13==-1) cycle
              do i3=0,int(Spec%range3,shortint)
                 do i2=0,int(Spec%range3,shortint)
                    do i1=0,int(Spec%range3,shortint)
                       if(i1+i2+i3+i12+i13>Spec%range1) cycle

                       write(iunit) i1,i2,i3,i12,i13,&
                            to_16(f(idx6( i13+1, i12+1, 1, i2+1, i3+1, i1+1)))

                    enddo
                 enddo
              enddo
           enddo
        enddo
        close(iunit)

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_2),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,i13,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( i13+1, i12+1, 0, i2+1, i3+1, i1+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,i13,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_2),error
        endif

        open(newunit=iunit,file=trim(Spec%name_2),access='stream')
        do i13=0,int(Spec%range2,shortint)
           do i12=0,int(Spec%range2,shortint)
              do i3=0,int(Spec%range3,shortint)
                 do i2=0,int(Spec%range3,shortint)
                    do i1=0,int(Spec%range3,shortint)
                       if(i1+i2+i3+i12+i13>Spec%range1) cycle

                       write(iunit) i1,i2,i3,i12,i13,&
                            to_16(f(idx6( i13+1, i12+1, 0, i2+1, i3+1, i1+1)))

                    enddo
                 enddo
              enddo
           enddo
        enddo
        close(iunit)

     case('13')

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_1),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,i13,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( i13+1, 1, i12+1, i2+1, i1+1, i3+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,i13,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_1),error
        endif

        open(newunit=iunit,file=trim(Spec%name_1),access='stream')
        do i13=-1,int(Spec%range2,shortint)
           do i12=-1,int(Spec%range2,shortint)
              if(i12==-1.and.i13==-1) cycle
              do i3=0,int(Spec%range3,shortint)
                 do i2=0,int(Spec%range3,shortint)
                    do i1=0,int(Spec%range3,shortint)
                       if(i1+i2+i3+i12+i13>Spec%range1) cycle

                       write(iunit) i1,i2,i3,i12,i13,&
                            to_16(f(idx6( i13+1, 1, i12+1, i2+1, i1+1, i3+1)))

                    enddo
                 enddo
              enddo
           enddo
        enddo
        close(iunit)

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_2),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,i13,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( i13+1, 0, i12+1, i2+1, i1+1, i3+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,i13,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_2),error
        endif

        open(newunit=iunit,file=trim(Spec%name_2),access='stream')
        do i13=0,int(Spec%range2,shortint)
           do i12=0,int(Spec%range2,shortint)
              do i3=0,int(Spec%range3,shortint)
                 do i2=0,int(Spec%range3,shortint)
                    do i1=0,int(Spec%range3,shortint)
                       if(i1+i2+i3+i12+i13>Spec%range1) cycle

                       write(iunit) i1,i2,i3,i12,i13,&
                            to_16(f(idx6( i13+1, 0, i12+1, i2+1, i1+1, i3+1)))

                    enddo
                 enddo
              enddo
           enddo
        enddo
        close(iunit)

     case('23')

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_1),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,i13,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( 1, i13+1, i12+1, i1+1, i2+1, i3+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,i13,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_1),error
        endif

        open(newunit=iunit,file=trim(Spec%name_1),access='stream')
        do i13=-1,int(Spec%range2,shortint)
           do i12=-1,int(Spec%range2,shortint)
              if(i12==-1.and.i13==-1) cycle
              do i3=0,int(Spec%range3,shortint)
                 do i2=0,int(Spec%range3,shortint)
                    do i1=0,int(Spec%range3,shortint)
                       if(i1+i2+i3+i12+i13>Spec%range1) cycle

                       write(iunit) i1,i2,i3,i12,i13,&
                            to_16(f(idx6( 1, i13+1, i12+1, i1+1, i2+1, i3+1)))

                    enddo
                 enddo
              enddo
           enddo
        enddo
        close(iunit)

        if(test) then
           error = 0
           open(newunit=iunit,file=trim(Spec%name_2),access='stream')
           do
              read(iunit,iostat=status) i1,i2,i3,i12,i13,val1
              if(status/=0) exit
              val2 = to_16(f(idx6( 0, i13+1, i12+1, i1+1, i2+1, i3+1)))
              val2 = abs((val1-val2)/(val1+val2))
              error = max(error,val2)
!              if(val2>1.e-30_16) write(*,*) i1,i2,i3,i12,i13,val2
           enddo
           close(iunit)
           write(*,'(a,es15.3)') trim(Spec%name_2),error
        endif

        open(newunit=iunit,file=trim(Spec%name_2),access='stream')
        do i13=0,int(Spec%range2,shortint)
           do i12=0,int(Spec%range2,shortint)
              do i3=0,int(Spec%range3,shortint)
                 do i2=0,int(Spec%range3,shortint)
                    do i1=0,int(Spec%range3,shortint)
                       if(i1+i2+i3+i12+i13>Spec%range1) cycle

                       write(iunit) i1,i2,i3,i12,i13,&
                            to_16(f(idx6( 0, i13+1, i12+1, i1+1, i2+1, i3+1)))

                    enddo
                 enddo
              enddo
           enddo
        enddo
        close(iunit)

     case default
        write(*,'(a)') 'UNKNOWN SYMBOL! QUITTING!'
        stop
     end select
   end associate
enddo

deallocate(f)
call free_indices

call f_fpu_fix_end(old_cw)

contains

type(qd_real) function to_qd(a) result(qd_a)
implicit none
real(16),intent(in) :: a
real(16) :: tmp
integer :: i

tmp = a
do i=1,3
   qd_a%re(i) = tmp
   tmp = tmp - qd_a%re(i)
enddo

qd_a%re(4) = 0.d0

end function to_qd

real(16) function to_16(qd_a) result(a)
implicit none
type(qd_real),intent(in) :: qd_a

a = sum(real(qd_a%re(1:3),16))

end function to_16

end subroutine f_main
