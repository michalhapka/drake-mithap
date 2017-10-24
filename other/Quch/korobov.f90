module korobov
use qdmodule
implicit none

private
public lenl
public prepare_korobov
public dilog,ygam

integer,parameter :: lenl = 200

type(qd_real),save :: zero,one,two,half,fourth,thresh,epsil,pi2_6,b(40)

contains

subroutine prepare_korobov
implicit none

zero   = qd_zero
one    = qd_one
two    = '2.0'
half   = '0.5'
fourth = '0.25'
thresh = '0.05'
epsil  = epsilon(one)
pi2_6  ='1.6449340668482264364724151666460251892189499012067984377355582294'

b( 1) =' 2.7777777777777777777777777777777777777777777777777777777777777778e-02'
b( 2) ='-2.7777777777777777777777777777777777777777777777777777777777777778e-04'
b( 3) =' 4.7241118669690098261526832955404383975812547241118669690098261527e-06'
b( 4) ='-9.1857730746619635508524397413286302175191064079952968841857730747e-08'
b( 5) =' 1.8978869988970999072009173019274029375039476049577059678069779080e-09'
b( 6) ='-4.0647616451442255268059093862919666745470571274397078222882048686e-11'
b( 7) =' 8.9216910204564525552179873167527488515142836130490451478105799093e-13'
b( 8) ='-1.9939295860721075687236443477937897056306947496538801470360150477e-14'
b( 9) =' 4.5189800296199181916504765528555932283968190144666184051991838850e-16'
b(10) ='-1.0356517612181247014483411542218656665960912381686505159641961307e-17'
b(11) =' 2.3952186210261867457402837430009803816789490019429742562511899344e-19'
b(12) ='-5.5817858743250093362830745056254199055670546676443980951358637371e-21'
b(13) =' 1.3091507554183212858123073991865923017498498387833038368544436042e-22'
b(14) ='-3.0874198024267402932422797648664624315955652561327456953256347812e-24'
b(15) =' 7.3159756527022034203579056092521485910334010636908750356925913770e-26'
b(16) ='-1.7408456572340007409890551477597025453408414217542712641710615283e-27'
b(17) =' 4.1576356446138997196178996207752266734882541595115638608246668632e-29'
b(18) ='-9.9621484882846221031940067024558388498548600173944887680619044799e-31'
b(19) =' 2.3940344248961653005211679878937495629342791569329157502206432977e-32'
b(20) ='-5.7683473553673900842917931618776542440723323179262751100622686433e-34'
b(21) =' 1.3931794796470079778278866039115483317324116256733995658061959100e-35'
b(22) ='-3.3721219654850894704684736352549309589797428916565393043858101374e-37'
b(23) =' 8.1782087775621026217647772148728342678761894624955032761980417558e-39'
b(24) ='-1.9870108311523859255648206692347865675418589958247432017904970119e-40'
b(25) =' 4.8357785180405508962870593731153782076944653694208278427320923312e-42'
b(26) ='-1.1786937248718384326695767537213903193540705623058985191453526287e-43'
b(27) =' 2.8770964081172571450019667396886617096883686589416665749345914616e-45'
b(28) ='-7.0320590981560280149649336675824257240019745278063228556910463362e-47'
b(29) =' 1.7208603145033146290899515161658919198419990913276265011572582782e-48'
b(30) ='-4.2160723905604454916800318192859916156937012298717379685552719382e-50'
b(31) =' 1.0340406405133039573902277553399774639907999441367470221742615354e-51'
b(32) ='-2.5386630625994653161632288930451123218275666294221390161755551563e-53'
b(33) =' 6.2385531769245908878361003503125478103452046261700167536143752658e-55'
b(34) ='-1.5344398069134650391696261221436482425147929342287479203962139366e-56'
b(35) =' 3.7772946355785502340013871251288878024916820684633969743942317948e-58'
b(36) ='-9.3058621248046865883933976705943761806774554153853988051934508988e-60'
b(37) =' 2.2943436822241873207151332543443075523754500239938637343783402548e-61'
b(38) ='-5.6606887394141478485716881165885209209448396398149038612790264885e-63'
b(39) =' 1.3975687219854008545365058198030725608701086163590964435421288496e-64'
b(40) ='-3.4526734733063388977836959915044772269931973257932824144076618996e-66'

end subroutine prepare_korobov

type(qd_real) function dilog(x)
implicit none
type(qd_real),intent(in) :: x
type(qd_real) :: y0,y1

if(x==zero) then

   dilog = pi2_6

elseif(x==one) then

   dilog = zero

elseif(x<half) then

   y0 =  log(x)
   y1 = -log(one-x)
   dilog = pi2_6 + y0*y1 - debay(y1)

elseif(x<=two) then

   y0 = -log(x)
   dilog = debay(y0)

else

   y0 = log(x)
   y1 = y0 - log(x-one)
   dilog = -(pi2_6 + half*y0**2 - y0*y1 - debay(y1))

endif

end function dilog

type(qd_real) function debay(z)
implicit none
type(qd_real),intent(in) :: z
type(qd_real) :: z2,yk,suma,prev
integer :: k

z2 = z*z

yk   = one
suma = zero
do k=1,40
   yk = yk*z2
   prev = suma
   suma = suma + b(k)*yk
   if(suma==prev) exit
enddo

debay = z*(one - fourth*z + suma)

end function debay

subroutine ygam(yg,alf,bet,gam,n)
implicit none
type(qd_real) :: yg(0:lenl,0:lenl)
type(qd_real),intent(in) :: alf,bet,gam
integer,intent(in) :: n
type(qd_real) :: gg(0:lenl)
type(qd_real) :: u,ab,ag,ainv,binv,temp
integer :: i,j
logical :: swap

swap = (bet<gam)

u  = bet-gam
ab = alf+bet
ag = alf+gam
if(swap) then
   u = -u
   temp = ab; ab = ag; ag = temp
endif
ainv = one/(bet+gam)
binv = one/ag
u    = u*binv

call dplog(u,n,gg)

yg(0,0) = gg(0)
temp = binv
do i=1,n
   yg(i,0) = gg(i)*temp
   temp = temp*binv
enddo
do j=1,n
   do i=0,n-j
      yg(i,j) = binv*((i+j-1)*yg(i,j-1) - ab*yg(i+1,j-1))
   enddo
enddo

yg(0,0) = binv*yg(0,0)
do i=1,n
   yg(i,0) = binv*yg(i,0)
end do
do j=1,n
   yg(0,j) = binv*(yg(0,j) + j*yg(0,j-1))
   do i=1,n-j
      yg(i,j) = binv*(yg(i,j) + j*yg(i,j-1))
   end do
end do

yg(0,0) = ainv*yg(0,0)
do i=1,n
   yg(i,0) = ainv*(yg(i,0) + i*yg(i-1,0))
enddo
do j=1,n
   yg(0,j) = ainv*(yg(0,j) + j*yg(0,j-1))
   do i=1,n-j
      yg(i,j) = ainv*(yg(i,j) + i*yg(i-1,j) + j*yg(i,j-1))
   enddo
enddo

if(swap) then
   do j=1,n
      do i=0,min(j-1,n-j)
         temp = yg(i,j)
         yg(i,j) = yg(j,i)
         yg(j,i) = temp
      enddo
   enddo
endif

end subroutine ygam

subroutine dplog(u,n,gg)
implicit none
type(qd_real),intent(in) :: u
integer,intent(in) :: n
type(qd_real) :: gg(0:lenl)
type(qd_real) :: au,uinv,u1inv
integer :: nadd,i

au = abs(u)

if(au>thresh) then

   uinv  = one/u
   u1inv = one/(one+u)

   gg(0) = log(one+u)
   gg(1) = -u1inv
   do i=2,n
      gg(i) = (i-1)*gg(i-1)*u1inv
   enddo

   gg(0) = gg(0)*uinv
   do i=1,n
      gg(i) = (gg(i) + i*gg(i-1))*uinv
   enddo

elseif(au==zero) then

   gg(0) = -one
   do i=1,n
      gg(i) = i*gg(i-1)
   enddo

   do i=n,0,-1
      gg(i) = -gg(i)/(i+1)
   enddo

else

   nadd = int(log(epsil)/log(au)) + 5
   if(n+nadd>lenl) then
      write (*,*) 'Error!!! Insufficient array length in DPLOG!'
      stop
   endif

   u1inv = one/(one+u)

   gg(0) = -u1inv
   do i=1,n+nadd-1
      gg(i) = i*gg(i-1)*u1inv
   enddo

   gg(n+nadd) = one
   do i=n+nadd-1,0,-1
      gg(i) = -(gg(i) - u*gg(i+1))/(i+1)
   enddo

endif

end subroutine dplog

end module korobov
