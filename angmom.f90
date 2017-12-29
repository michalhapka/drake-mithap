module angmom
use file_OUT, only : LOUT
implicit none

private
public get_angmom_name,get_angmom_l
public get_gen_n,get_gen_name

contains

function get_angmom_name(l) result(name)
implicit none
character :: name
integer,intent(in) :: l
select case(l)
case(0)
name='S'
case(1)
name='P'
case(2)
name='D'
case(3)
name='F'
case(4)
name='G'
case(5)
name='H'
case(6)
name='I'
case(7)
name='K'
case(8)
name='M'
case(9)
name='N'
case(10)
name='O'
case(11)
name='Q'
case(12)
name='R'
case(13)
name='T'
case(14)
name='U'
case(15)
name='V'
case(16)
name='W'
case(17)
name='X'
case(18)
name='Y'
case(19)
name='Z'
case default
write(LOUT,'(a,i3)') 'Incorrect (too large) angular momentum: ',l
stop
end select
end function get_angmom_name

function get_angmom_l(name) result(l)
implicit none
integer :: l
character(*),intent(in) :: name
select case(trim(name))
case('S')
l=0
case('P')
l=1
case('D')
l=2
case('F')
l=3
case('G')
l=4
case('H')
l=5
case('I')
l=6
case('K')
l=7
case('M')
l=8
case('N')
l=9
case('O')
l=10
case('Q')
l=11
case('R')
l=12
case('T')
l=13
case('U')
l=14
case('V')
l=15
case('W')
l=16
case('X')
l=17
case('Y')
l=18
case('Z')
l=19
case default
write(LOUT,'(2a)') 'Incorrect symbol of angular momentum: ',name
stop
end select
end function get_angmom_l

function get_gen_n(name) result(n)
implicit none
integer :: n
character(*),intent(in) :: name
select case(trim(name))
case('S^E')
n=1
case('P^E')
n=2
case('P^O')
n=3
case('D^E(1)')
n=4
case('D^E(2)')
n=5
case default
write(LOUT,'(2a)') 'Incorrect generator symbol : ',name
stop
end select
end function get_gen_n

function get_gen_name(n) result(name)
implicit none
character(6) :: name
integer,intent(in) :: n
select case(n)
case(1)
name='S^E'
case(2)
name='P^E'
case(3)
name='P^O'
case(4)
name='D^E(1)'
case(5)
name='D^E(2)'
case default
write(LOUT,'(a,i3)') 'Incorrect generator number: ',n
stop
end select
end function get_gen_name

end module angmom
