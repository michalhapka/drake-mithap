module pachucki
use qdmodule
use indices
use korobov
implicit none

private
public p_f

contains

subroutine p_f(f,w1,w2,w3,sm,max_2,max_1)
implicit none
type(qd_real) :: f(0:MAXF)
type(qd_real),intent(in) :: w1,w2,w3
integer,intent(in) :: sm,max_2,max_1
integer :: n1,n2,n3,n4,n5,n6,min_n,max_n
integer :: i,j,k,l,m,n,index
type(qd_real) :: nd1,nd2,nd3,nd4,nd5,nd6
type(qd_real) :: expr,expr2,comb,inv_w
type(qd_real) :: tg1(0:MAXG)
type(qd_real) :: tg2(0:MAXG)
type(qd_real) :: tg3(0:MAXG)
type(qd_real) :: tg4(0:MAXG)
type(qd_real) :: tg5(0:MAXG)
type(qd_real) :: tg6(0:MAXG)

call prepare_korobov

call p_g(tg1,w2+w3,w1,sm) 
call p_g(tg2,w1+w3,w2,sm)
call p_g(tg3,w1+w2,w3,sm)
call p_g(tg4,w2,w3,sm)
call p_g(tg5,w3,w1,sm)
call p_g(tg6,w1,w2,sm)

f(0) = 0

if(sm>=0) then 
   f(idx6(0,0,0,0,0,0))=-1/(2*w1*w2*w3)* &
        ( &
        log(w3/(w1+w2))*log(1+(w3/(w1+w2))) &
        +dilog(1+w3/(w1+w2))+dilog((w3/(w1+w2))) &
        +log(w2/(w3+w1))*log(1+(w2/(w3+w1))) &
        +dilog(1+w2/(w3+w1))+dilog((w2/(w3+w1))) &
        +log(w1/(w2+w3))*log(1+(w1/(w2+w3))) &
        +dilog(1+w1/(w2+w3))+dilog((w1/(w2+w3))) &
        )
endif

if(sm>=1) then
   f(idx6(1,0,0,0,0,0))=-1/(w2**2*w3**2)*log((w1*(w1+w2+w3))/((w1+w2)*(w1+w3)))
   f(idx6(0,1,0,0,0,0))=-1/(w1**2*w3**2)*log((w2*(w1+w2+w3))/((w2+w3)*(w2+w1))) 
   f(idx6(0,0,1,0,0,0))=-1/(w1**2*w2**2)*log((w3*(w1+w2+w3))/((w3+w1)*(w3+w2)))
endif

if(sm>=2) then
   f(idx6(1,1,0,0,0,0))=1/(w1*w2*(w1+w2)*w3**2)
   f(idx6(1,0,1,0,0,0))=1/(w1*w3*(w1+w3)*w2**2) 
   f(idx6(0,1,1,0,0,0))=1/(w2*w3*(w2+w3)*w1**2)
endif

if(sm>=3) f(idx6(1,1,1,0,0,0))=1/(w1**2*w2**2*w3**2)

do k=2,sm
   do j=0,k
      do i=0,j

         n1 = i
         n2 = j-i
         n3 = k-j

         min_n = min(n1,n2,n3)
         max_n = max(n1,n2,n3)

         if(max_n<2) cycle
         if(min_n>5) cycle
         if(max_n>max_2+2) cycle

         nd1 = n1
         nd2 = n2
         nd3 = n3

! expr=f(n1+2,n2,n3)
         if(n1==max_n) then

            n1 = n1-2
            nd1 = n1

            expr2=(nd1+2*nd2+nd3+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd1+1)*tg6(idx3(n2-1,n1+1,n3))+ &
                 1/(nd3+1)*tg4(idx3(n3+1,n2-1,n1))
            if(n2==0) then
               expr2=expr2- &
                    1/(nd1+1)*tg2(idx3(-1,0,n1+n3+1))- &
                    1/(nd3+1)*tg2(idx3(-1,0,n1+n3+1))
            elseif(n2>1) then
               expr2=expr2+ &
                    nd2*(-1+nd2)/(nd1+1)*f(idx6(n1+2,n2-2,n3,0,0,0))+ &
                    nd2*(-1+nd2)/(nd3+1)*f(idx6(n1,n2-2,n3+2,0,0,0))
            endif

            expr=1/(w3**2)*expr2

            expr2=(2*nd3+nd2+nd1+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd1+1)*tg5(idx3(n1+1,n3-1,n2))+ &
                 1/(nd2+1)*tg4(idx3(n3-1,n2+1,n1))
            if(n3==0) then
               expr2=expr2- &
                    1/(nd1+1)*tg3(idx3(-1,0,n1+n2+1))- &
                    1/(nd2+1)*tg3(idx3(-1,0,n1+n2+1))
            elseif(n3>1) then
               expr2=expr2+ &
                    nd3*(-1+nd3)/(nd1+1)*f(idx6(n1+2,n2,n3-2,0,0,0))+ &
                    nd3*(-1+nd3)/(nd2+1)*f(idx6(n1,n2+2,n3-2,0,0,0))
            endif

            expr=expr+1/(w2**2)*expr2  

            expr2=(nd3+nd2+2*nd1+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd2+1)*tg6(idx3(n2+1,n1-1,n3))+ &
                 1/(nd3+1)*tg5(idx3(n1-1,n3+1,n2))
            if(n1==0) then
               expr2=expr2- &
                    1/(nd2+1)*tg1(idx3(-1,0,n2+n3+1))- &
                    1/(nd3+1)*tg1(idx3(-1,0,n2+n3+1))
            elseif(n1>1) then
               expr2=expr2+ &
                    nd1*(-1+nd1)/(nd2+1)*f(idx6(n1-2,n2+2,n3,0,0,0))+ &
                    nd1*(-1+nd1)/(n3+1)*f(idx6(n1-2,n2,n3+2,0,0,0))
            endif

            expr=expr-w1**2/((w3*w2)**2)*expr2

            f(idx6(n1+2,n2,n3,0,0,0))=expr*(1+nd1)/2

            cycle
         endif

! expr=f(n1,n2+2,n3)
         if(n2==max_n) then

            n2=n2-2
            nd2=n2

            expr2=(nd1+2*nd3+nd2+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd2+1)*tg4(idx3(n3-1,n2+1,n1))+ &
                 1/(nd1+1)*tg5(idx3(n1+1,n3-1,n2)) 
            if(n3==0) then
               expr2=expr2- &
                    1/(nd2+1)*tg3(idx3(-1,0,n1+n2+1))- &
                    1/(nd1+1)*tg3(idx3(-1,0,n1+n2+1))
            elseif(n3>1) then
               expr2=expr2+ &
                    nd3*(-1+nd3)/(nd2+1)*f(idx6(n1,n2+2,n3-2,0,0,0))+ &
                    nd3*(-1+nd3)/(nd1+1)*f(idx6(n1+2,n2,n3-2,0,0,0))
            endif

            expr=1/(w1**2)*expr2

            expr2=(2*nd1+nd2+nd3+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd2+1)*tg6(idx3(n2+1,n1-1,n3))+ &
                 1/(nd3+1)*tg5(idx3(n1-1,n3+1,n2))
            if(n1==0) then
               expr2=expr2- &
                    1/(nd2+1)*tg1(idx3(-1,0,n3+n2+1))- &
                    1/(nd3+1)*tg1(idx3(-1,0,n3+n2+1))
            elseif(n1>1) then
               expr2=expr2+ &
                    nd1*(-1+nd1)/(nd2+1)*f(idx6(n1-2,n2+2,n3,0,0,0))+ &
                    nd1*(-1+nd1)/(nd3+1)*f(idx6(n1-2,n2,n3+2,0,0,0))
            endif

            expr=expr+1/(w3**2)*expr2

            expr2=(nd1+nd3+2*nd2+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd3+1)*tg4(idx3(n3+1,n2-1,n1))+ &
                 1/(nd1+1)*tg6(idx3(n2-1,n1+1,n3)) 
            if(n2==0) then
               expr2=expr2- &
                    1/(nd3+1)*tg2(idx3(-1,0,n1+n3+1))- &
                    1/(nd1+1)*tg2(idx3(-1,0,n1+n3+1))
            elseif(n2>1) then
               expr2=expr2+ &
                    nd2*(-1+nd2)/(nd3+1)*f(idx6(n1,n2-2,n3+2,0,0,0))+ &
                    nd2*(-1+nd2)/(nd1+1)*f(idx6(n1+2,n2-2,n3,0,0,0))
            endif

            expr=expr-w2**2/((w1*w3)**2)*expr2

            f(idx6(n1,n2+2,n3,0,0,0))=expr*(1+nd2)/2

            cycle
         endif

! expr=f(n1,n2,n3+2)
         if(n3==max_n) then

            n3=n3-2
            nd3=n3

            expr2=(nd1+2*nd2+nd3+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd3+1)*tg4(idx3(n3+1,n2-1,n1))+ &
                 1/(nd1+1)*tg6(idx3(n2-1,n1+1,n3)) 
            if(n2==0) then
               expr2=expr2- &
                    1/(nd3+1)*tg2(idx3(-1,0,n1+n3+1))- &
                    1/(nd1+1)*tg2(idx3(-1,0,n1+n3+1))
            elseif(n2>1) then
               expr2=expr2+ &
                    nd2*(-1+nd2)/(nd3+1)*f(idx6(n1,n2-2,n3+2,0,0,0))+ &
                    nd2*(-1+nd2)/(nd1+1)*f(idx6(n1+2,n2-2,n3,0,0,0))
            endif

            expr=1/(w1**2)*expr2

            expr2=(2*nd1+nd2+nd3+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd3+1)*tg5(idx3(n1-1,n3+1,n2))+ &
                 1/(nd2+1)*tg6(idx3(n2+1,n1-1,n3))
            if(n1==0) then
               expr2=expr2- &
                    1/(nd3+1)*tg1(idx3(-1,0,n3+n2+1))- &
                    1/(nd2+1)*tg1(idx3(-1,0,n3+n2+1))
            elseif(n1>1) then
               expr2=expr2+ &
                    nd1*(-1+nd1)/(nd3+1)*f(idx6(n1-2,n2,n3+2,0,0,0))+ &
                    nd1*(-1+nd1)/(nd2+1)*f(idx6(n1-2,n2+2,n3,0,0,0))
            endif

            expr=expr+1/(w2**2)*expr2

            expr2=(nd1+nd2+2*nd3+2)*f(idx6(n1,n2,n3,0,0,0))+ &
                 1/(nd2+1)*tg4(idx3(n3-1,n2+1,n1))+ &
                 1/(nd1+1)*tg5(idx3(n1+1,n3-1,n2))
            if(n3==0) then
               expr2=expr2- &
                    1/(nd2+1)*tg3(idx3(-1,0,n1+n2+1))- &
                    1/(nd1+1)*tg3(idx3(-1,0,n1+n2+1))
            elseif(n3>1) then
               expr2=expr2+ &
                    nd3*(-1+nd3)/(nd2+1)*f(idx6(n1,n2+2,n3-2,0,0,0))+ &
                    nd3*(-1+nd3)/(nd1+1)*f(idx6(n1+2,n2,n3-2,0,0,0))
            endif

            expr=expr-w3**2/((w1*w2)**2)*expr2

            f(idx6(n1,n2,n3+2,0,0,0))=expr*(1+nd3)/2

            cycle
         endif

      enddo
   enddo
enddo

inv_w = 1/(w1*w2*w3)

index=1

do n=1,sm
   do m=0,n
      do l=0,m
         do k=0,l
            do j=0,k
               do i=0,j

                  index = index+1

                  n1 = i
                  n2 = j-i
                  n3 = k-j

                  min_n = min(n1,n2,n3)
                  max_n = max(n1,n2,n3)

                  if(min_n>1) cycle
                  if(max_n>max_2) cycle

                  n4 = l-k
                  n5 = m-l
                  n6 = n-m

                  max_n = max(n4,n5,n6)

                  if(n4+n5+n6==0) cycle
                  if(max_n>max_1+1) cycle

                  nd1 = n1
                  nd2 = n2
                  nd3 = n3
                  nd4 = n4
                  nd5 = n5
                  nd6 = n6

! f(n1,n2,n3,n4+1,n5,n6)
                  if(n4==max_n) then

                     n4 = n4-1
                     nd4 = n4

                     comb=(-1+nd1-nd2-nd3-nd4)

                     expr=-comb*f(idx6(n1,n2,n3,n4,n5,n6))*w2*w3
                     if(n5>0) then
                        expr2=comb* &
                             f(idx6(n1,n2,n3,n4,-1+n5,n6))+ &
                             f(idx6(n1,n2,n3,1+n4,-1+n5,n6))*w1
                        expr=expr+nd5*w3*expr2
                     endif
                     if(n6>0) then
                        expr2=comb* &
                             f(idx6(n1,n2,n3,n4,n5,-1+n6))+ &
                             f(idx6(n1,n2,n3,1+n4,n5,-1+n6))*w1
                        expr=expr+nd6*w2*expr2
                     endif
                     if(n5*n6>0) then
                        expr2=comb* &
                             f(idx6(n1,n2,n3,n4,-1+n5,-1+n6))+ &
                             f(idx6(n1,n2,n3,1+n4,-1+n5,-1+n6))*w1
                        expr=expr-nd5*nd6*expr2
                     endif
                     if(n2==0) then
                        expr=expr+tg2(idx3(n4+n6,n5,n1+n3-1))*w2
                        if(n5>0) expr=expr- &
                             nd5*tg2(idx3(n4+n6,n5-1,n1+n3-1))
                     elseif(n2>1) then
                        expr2=-f(idx6(n1,-2+n2,n3,n4,n5,1+n6))*w2
                        if(n5>0) expr2=expr2+ &
                             nd5*f(idx6(n1,-2+n2,n3,n4,-1+n5,1+n6))
                        expr=expr+(-1+nd2)*nd2*expr2
                     endif
                     if(n3==0) then
                        expr=expr+tg3(idx3(n4+n5,n6,n1+n2-1))*w3
                        if(n6>0) expr=expr- &
                             nd6*tg3(idx3(n4+n5,n6-1,n1+n2-1))
                     elseif(n3>1) then
                        expr2=-f(idx6(n1,n2,-2+n3,n4,1+n5,n6))*w3
                        if(n6>0) expr2=expr2+ &
                             nd6*f(idx6(n1,n2,-2+n3,n4,1+n5,-1+n6))
                        expr=expr+(-1+nd3)*nd3*expr2
                     endif
                     if(n1==0) then
                        expr=expr-tg1(idx3(n5+n6,n4,n2+n3-1))*(w2+w3)
                        if(n5+n6>0) expr=expr+ &
                             (nd5+nd6)*tg1(idx3(n5+n6-1,n4,n2+n3-1))
                     elseif(n1>1) then
                        expr2= &
                             f(idx6(-2+n1,n2,n3,n4,n5,1+n6))*w2+ &
                             f(idx6(-2+n1,n2,n3,n4,1+n5,n6))*w3
                        if(n5>0) expr2=expr2- &
                             nd5*f(idx6(-2+n1,n2,n3,n4,-1+n5,1+n6))
                        if(n6>0) expr2=expr2- &
                             nd6*f(idx6(-2+n1,n2,n3,n4,1+n5,-1+n6))
                        expr=expr+(-1+nd1)*nd1*expr2
                     endif

                     f(index)=expr*inv_w

                     cycle
                  endif

! f(n1,n2,n3,n4,n5+1,n6)
                  if(n5==max_n) then         

                     n5 = n5-1
                     nd5 = n5

                     comb=(1+nd1-nd2+nd3+nd5)

                     expr=comb*f(idx6(n1,n2,n3,n4,n5,n6))*w1*w3
                     if(n4>0) then
                        expr2=-comb* &
                             f(idx6(n1,n2,n3,-1+n4,n5,n6))+ &
                             f(idx6(n1,n2,n3,-1+n4,1+n5,n6))*w2
                        expr=expr+nd4*w3*expr2
                     endif
                     if(n6>0) then
                        expr2=-comb* &
                             f(idx6(n1,n2,n3,n4,n5,-1+n6))+ &
                             f(idx6(n1,n2,n3,n4,1+n5,-1+n6))*w2
                        expr=expr+nd6*w1*expr2
                     endif
                     if(n4*n6>0) then
                        expr2=-comb* &
                             f(idx6(n1,n2,n3,-1+n4,n5,-1+n6))+ &
                             f(idx6(n1,n2,n3,-1+n4,1+n5,-1+n6))*w2
                        expr=expr-nd4*nd6*expr2
                     endif
                     if(n1==0) then
                        expr=expr+tg1(idx3(n5+n6,n4,n3+n2-1))*w1
                        if(n4>0) expr=expr- &
                             nd4*tg1(idx3(n5+n6,n4-1,n3+n2-1))
                     elseif(n1>1) then
                        expr2=-f(idx6(-2+n1,n2,n3,n4,n5,1+n6))*w1
                        if(n4>0) expr2=expr2+ &
                             nd4*f(idx6(-2+n1,n2,n3,-1+n4,n5,1+n6))
                        expr=expr+(-1+nd1)*nd1*expr2
                     endif
                     if(n3==0) then
                        expr=expr+tg3(idx3(n4+n5,n6,n1+n2-1))*w3
                        if(n6>0) expr=expr- &
                             nd6*tg3(idx3(n4+n5,n6-1,n1+n2-1))
                     elseif(n3>1) then
                        expr2=-f(idx6(n1,n2,-2+n3,1+n4,n5,n6))*w3
                        if(n6>0) expr2=expr2+ &
                             nd6*f(idx6(n1,n2,-2+n3,1+n4,n5,-1+n6))
                        expr=expr+(-1+nd3)*nd3*expr2
                     endif
                     if(n2==0) then
                        expr=expr-tg2(idx3(n4+n6,n5,n1+n3-1))*(w1+w3)
                        if(n4+n6>0) expr=expr+ &
                             (nd4+nd6)*tg2(idx3(n4+n6-1,n5,n1+n3-1))
                     elseif(n2>1) then
                        expr2= &
                             f(idx6(n1,-2+n2,n3,n4,n5,1+n6))*w1+ &
                             f(idx6(n1,-2+n2,n3,1+n4,n5,n6))*w3
                        if(n4>0) expr2=expr2- &
                             nd4*f(idx6(n1,-2+n2,n3,-1+n4,n5,1+n6))
                        if(n6>0) expr2=expr2- &
                             nd6*f(idx6(n1,-2+n2,n3,1+n4,n5,-1+n6))
                        expr=expr+(-1+nd2)*nd2*expr2
                     endif

                     f(index)=expr*inv_w

                     cycle
                  endif

! f(n1,n2,n3,n4,n5,n6+1)
                  if(n6==max_n) then

                     n6 = n6-1
                     nd6 = n6

                     comb=(1+nd1+nd2-nd3+nd6)

                     expr=comb*f(idx6(n1,n2,n3,n4,n5,n6))*w1*w2
                     if(n4>0) then
                        expr2=-comb* &
                             f(idx6(n1,n2,n3,-1+n4,n5,n6))+ &
                             f(idx6(n1,n2,n3,-1+n4,n5,1+n6))*w3
                        expr=expr+nd4*w2*expr2
                     endif
                     if(n5>0) then
                        expr2=-comb* &
                             f(idx6(n1,n2,n3,n4,-1+n5,n6))+ &
                             f(idx6(n1,n2,n3,n4,-1+n5,1+n6))*w3
                        expr=expr+nd5*w1*expr2
                     endif
                     if(n4*n5>0) then
                        expr2=-comb* &
                             f(idx6(n1,n2,n3,-1+n4,-1+n5,n6))+ &
                             f(idx6(n1,n2,n3,-1+n4,-1+n5,1+n6))*w3
                        expr=expr-nd4*nd5*expr2
                     endif
                     if(n1==0) then
                        expr=expr+tg1(idx3(n5+n6,n4,n3+n2-1))*w1
                        if(n4>0) expr=expr- &
                             nd4*tg1(idx3(n5+n6,n4-1,n3+n2-1))
                     elseif(n1>1) then
                        expr2=-f(idx6(-2+n1,n2,n3,n4,1+n5,n6))*w1
                        if(n4>0) expr2=expr2+ &
                             nd4*f(idx6(-2+n1,n2,n3,-1+n4,1+n5,n6))
                        expr=expr+(-1+nd1)*nd1*expr2
                     endif
                     if(n2==0) then
                        expr=expr+tg2(idx3(n4+n6,n5,n1+n3-1))*w2
                        if(n5>0) expr=expr- &
                             nd5*tg2(idx3(n4+n6,n5-1,n1+n3-1))
                     elseif(n2>1) then
                        expr2=-f(idx6(n1,-2+n2,n3,1+n4,n5,n6))*w2
                        if(n5>0) expr2=expr2+ &
                             nd5*f(idx6(n1,-2+n2,n3,1+n4,-1+n5,n6))
                        expr=expr+(-1+nd2)*nd2*expr2
                     endif
                     if(n3==0) then
                        expr=expr-tg3(idx3(n4+n5,n6,n1+n2-1))*(w1+w2)
                        if(n4+n5>0) expr=expr+ &
                             (nd4+nd5)*tg3(idx3(n4+n5-1,n6,n1+n2-1))
                     elseif(n3>1) then
                        expr2= &
                             f(idx6(n1,n2,-2+n3,n4,1+n5,n6))*w1+ &
                             f(idx6(n1,n2,-2+n3,1+n4,n5,n6))*w2
                        if(n4>0) expr2=expr2- &
                             nd4*f(idx6(n1,n2,-2+n3,-1+n4,1+n5,n6))
                        if(n5>0) expr2=expr2- &
                             nd5*f(idx6(n1,n2,-2+n3,1+n4,-1+n5,n6))
                        expr=expr+(-1+nd3)*nd3*expr2
                     endif

                     f(index)=expr*inv_w

                     cycle
                  endif

               enddo
            enddo
         enddo
      enddo
   enddo
enddo

end subroutine p_f

subroutine p_g(tg,a,b,sm)
implicit none
type(qd_real) :: tg(0:MAXG)
type(qd_real),intent(in) :: a,b
integer,intent(in) :: sm
type(qd_real) :: tb(0:MAXG)
integer :: sc
integer :: l,m,n
type(qd_real) :: inv_ab,zero
type(qd_real) :: yg(0:lenl,0:lenl)

call p_tb(tb,a,b,sm)

inv_ab = 1/(a+b)

if(sm>=0) tg(idx3(0,0,0)) = inv_ab/(a*b)

do sc=1,sm
   do l=0,sc
      do m=0,sc-l
         n=sc-l-m

         tg(idx3(l,m,n))=inv_ab*tb(idx3(l,m,n))
         if(l>0) tg(idx3(l,m,n))=tg(idx3(l,m,n))+inv_ab*l*tg(idx3(l-1,m,n))
         if(m>0) tg(idx3(l,m,n))=tg(idx3(l,m,n))+inv_ab*m*tg(idx3(l,m-1,n))

      enddo
   enddo
enddo

zero = 0

call ygam(yg,a,b,zero,sm+1)
do n=0,sm+1
   do m=0,sm+1-n
      tg(idx3(-1,m,n))=yg(m,n)
   enddo
enddo

call ygam(yg,b,a,zero,sm+1)
do n=0,sm+1
   do m=0,sm+1-n
      tg(idx3(m,-1,n))=yg(m,n)
   enddo
enddo

call ygam(yg,zero,b,a,sm+1)
do n=0,sm+1
   do m=0,sm+1-n
      tg(idx3(n,m,-1))=yg(m,n)
   enddo
enddo

end subroutine p_g

subroutine p_tb(tb,a,b,sm)
implicit none
type(qd_real) :: tb(0:MAXG)
type(qd_real),intent(in) :: a,b
integer,intent(in) :: sm
type(qd_real) :: ta(0:sm)
integer :: sc
integer :: l,m,n

if(sm>=0) then
   ta(0) = 1/b
   tb(idx3(0,0,0))=1/(a*b)
endif

do sc=1,sm

   ta(sc) = 1/b*sc*ta(sc-1)

   do l=0,sc
      do m=0,sc-l  
         n=sc-l-m

         if(n==0) then          
            tb(idx3(l,m,n))=0 
         else
            tb(idx3(l,m,n))=1/a*n*tb(idx3(l,m,n-1))
         endif

         if(l==0) then
            if(n==0) then
               tb(idx3(l,m,n))=tb(idx3(l,m,n))+1/b*m*tb(idx3(l,m-1,n))
            else          
               tb(idx3(l,m,n))=tb(idx3(l,m,n))+1/a*ta(m+n)
            endif
         else
            tb(idx3(l,m,n))=tb(idx3(l,m,n))+1/a*l*tb(idx3(l-1,m,n))
         endif

      enddo
   enddo

enddo

end subroutine p_tb

end module pachucki
