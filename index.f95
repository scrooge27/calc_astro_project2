module rhs
    implicit none
    contains
    subroutine dydx(y,f,neq)
        real*8,dimension(neq),intent(in)::y
        real*8,dimension(neq),intent(out)::f
        integer,intent(in)::neq
        real*8::x1,x2,dphi_x,dphi_y

        real*8,parameter:: m1=1.989d30,m2=1.898d27,mt=m1+m2, &
                    au=1.495978707d11,a=5.204d0*au, &
                    g=6.67d-11,pi=acos(-1.d0), &
                    p=2*pi*sqrt(a**3/(g*(m1+m2))),n=2*pi/p, &
                    p_red=sqrt(a*g*mt)
        !modifica le riga successive a seconda della funzione in esame

        x1=-a*m2/mt
    
        x2=a*m1/mt

        dphi_x= m1*(y(1)*a-x1)/sqrt((y(1)*a-x1)**2+(y(3)*a)**2)**(3)+ &
                m2*(y(1)*a-x2)/sqrt((y(1)*a-x2)**2+(y(3)*a)**2)**(3)

        dphi_y= m1/sqrt((y(1)*a-x1)**2+(y(3)*a)**2)**(3)+ &
                m2/sqrt((y(1)*a-x2)**2+(y(3)*a)**2)**(3)

        f(1)=y(2)                                                          !dx/dt=vx
        f(2)=(2*n*p_red*y(4)+a*y(1)*n**2)/(g*mt)-dphi_x/mt                 !dvx/dt=ax
        f(3)=y(4)                                                          !dy/dt=vy
        f(4)=(-2*n*p_red*y(2)+a*y(3)*n**2)/(g*mt)-a*y(3)*dphi_y/mt         !dvy/dt=ay
    
    end subroutine dydx
    
end module rhs

module ode_solver
    use rhs
    implicit none
    contains

    subroutine rk4(h,yold,ynew,neq)
        real*8,intent(in)::h
        integer,intent(in)::neq
        real*8,dimension(neq),intent(in)::yold
        real*8,dimension(neq),intent(out)::ynew
        real*8,dimension(neq)::k1,k2,k3,k4
        integer::i

        call dydx(yold,k1,neq)

        call dydx(yold+0.5d0*h*k1,k2,neq)
        
        call dydx(yold+0.5d0*h*k2,k3,neq)

        call dydx(yold+h*k3,k4,neq)

        do i=1,neq
            ynew(i)=yold(i)+h/6.d0*(k1(i)+2.d0*k2(i)+2.d0*k3(i)+k4(i))
        end do

    end subroutine rk4

    subroutine save_results(doc,x,y,n,neq)
        character (len=*),intent(in)::doc   
        integer,intent(in)::n,neq
        real*8,dimension(0:n),intent(in)::x
        real*8,dimension(neq,0:n),intent(in)::y
        integer::i,j

        !mi serve per dirgli su quante colonne formattare perchè sono stronzo e voglio formattarlo
        character(len=25):: fmt_string
        !scrivo nell'oggetto fmt_string (non esiste un cast migliore)
        write(fmt_string,'(a,i3,a)') '(' ,neq+1, '(1pe19.9))' !l'intero che individua le colonne occupa al più tre spazi-i3

        open(10,file=doc)

        do i=0,n 
            write(10,fmt_string) x(i),(y(j,i),j=1,neq)
        end do

        close(10)
    end subroutine save_results
    
end module ode_solver

program main
    use ode_solver
    implicit none

    character(len=10)::doc="3corpi.txt"

    integer,parameter::np=100,neq=4
    integer::i
    !real*8::x1,x2,y1,y2, &
    !        r(n),r1(n),r2(n),
    real*8::h
    real*8,dimension(0:np)::t
    real*8,dimension(neq,0:np)::y
    real*8,parameter:: m1=1.989d30,m2=1.898d27,mt=m1+m2, &
                    au=1.495978707d11,a=5.204d0*au, &
                    g=6.67d-11,pi=acos(-1.d0), &
                    p=2*pi*sqrt(a**3/(g*(m1+m2))),n=2*pi/p, &
                    p_red=sqrt(a*g*mt)
    h=60*pi*a/np                 !h tilde
    
    !x1=-a*m2/mt
    !y1=0.d0

    !x2=a*m1/mt
    !y2=0.d0

    !condizioni iniziali
    y(1,0)=-1.02745d0          !x(0)   tilde
    y(2,0)=0.d0                 !vx(0)  tilde
    y(3,0)=0.d0                 !y(0)   tilde
    y(4,0)=526.59484d0/p_red    !vy(0)  tilde
    
    do i=0,np
        t(i)=h*i                !t tilde
        
        !r1(i)=sqrt((y(1,i)-x1)**2+y(1,t)**2)
        !r2(i)=sqrt((y(1,i)-x2)**2+y(1,t)**2)
    
        !r(i)=sqrt(x(1,t)**2+y(1,t)**2)
    end do

    do i=1,np
        call rk4(h,y(:,i-1),y(:,i),neq)       !la subroutine agisce già sulle derivate di ogni ordine
    end do

    !tolgo la normalizzazione
    y(1,:)=a*y(1,:)             !x
    y(2,:)=p_red*y(2,:)         !vx

    y(3,:)=a*y(3,:)             !y
    y(4,:)=p_red*y(4,:)         !vy

    t=a*t/p_red

    call save_results(doc,t,y,np,neq)

end program main