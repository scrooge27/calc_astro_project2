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

        f(1)=y(2)                                                       !dx/dt=vx
        f(2)=4*(pi)**2*(y(4)/pi+y(1)-a**2/mt*dphi_x)           !dvx/dt=ax
        f(3)=y(4)                                                       !dy/dt=vy
        f(4)=4*(pi)**2*(-y(2)/pi+y(3)-a**3*y(3)/mt*dphi_y)          !dvy/dt=ay
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

        !mi serve per dirgli su quante colonne formattare
        character(len=25):: fmt_string
        !scrivo nell'oggetto fmt_string (non esiste un cast migliore)
        write(fmt_string,'(a,i3,a)') '(' ,neq+1, '(1pe19.9))' !l'intero che individua le colonne occupa al più tre spazi-i3

        open(10,file=doc)

        do i=0,n 
            write(10,fmt_string) x(i),(y(j,i),j=1,neq)
            if (doc=="nrg_ds1.txt") then
                print*,y(j,i)
                read(*,*)
            end if
        end do

        close(10)
    end subroutine save_results
    
end module ode_solver

program main
    use ode_solver
    implicit none

    character(len=11)::doc

    integer,parameter::np=1000,neq=4
    integer::i,j
    real*8::x1,x2,&
        m0,m1t,m2t
    real*8::h(2),v2,w2,u
    real*8,dimension(0:np)::t,rad
    real*8,dimension(neq,0:np)::y
    real*8,dimension(2,0:np)::r
    real*8,dimension(1,0:np)::cj,e
    real*8,parameter::  m1=1.989d30,m2=1.898d27,mt=m1+m2, &
                        au=1.495978707d11,a=5.204d0*au, &
                        g=6.67d-11,pi=acos(-1.d0), &
                        p=2*pi*sqrt(a**3/(g*mt)), &
                        n=2*pi/p

    doc(8:11)=".txt"
    
    
    m0=mt/(4*pi**2)
    m1t=m1/m0                       !m1 tilde
    m2t=m2/m0                       !m2 tilde

    h(1)=60.d0*pi/np                   !h tilde
    h(2)=h(1)*0.1
    do j=1,2

        write(doc(7:7),'(i0.0)') j
        !condizioni iniziali
        y(1,0)=-1.02745d0         !x(0)   tilde
        y(2,0)=0.d0               !vx(0)  tilde
        y(3,0)=0.d0               !y(0)   tilde
        y(4,0)=526.59484d0*p/a    !vy(0)  tilde
        
        do i=0,np
            t(i)=h(j)*i                !t tilde
        end do

        do i=1,np
            call rk4(h(j),y(:,i-1),y(:,i),neq)       !la subroutine agisce già sulle derivate di ogni ordine
        end do

        !step 3 e 4 perchè vogliono la adimensionalizzazione

        x1=-m2/mt                   !x1 tilde (ho già diviso per a)     la divisione tra masse è comunque adimensionale
        x2=m1/mt                    !x2 tilde (ho già diviso per a)     la divisione tra masse è comunque adimensionale

        do i=0,np  
            r(1,i)=sqrt((y(1,i)-x1)**2+y(3,i)**2)
            r(2,i)=sqrt((y(1,i)-x2)**2+y(3,i)**2)    
            rad(i)=sqrt(y(1,i)**2+y(3,i)**2)
    
            v2=y(2,i)**2+y(4,i)**2
            !uso m1t ed m2t,
            !pongo g=1,
            !pongo n=2pi
    
            u=0.5d0*(2*pi)**2*(y(1,i)**2+y(3,i)**2)+m1t/r(1,i)+m2t/r(2,i)
            cj(1,i)=2.d0*u-v2
            
            w2=v2+(2*pi*rad(i))**2+4*pi*(y(4,i)*y(1,i)-y(2,i)*y(3,i))
            e(1,j)=w2*0.5d0-m1t/r(1,i)-m2t/r(2,i)
        end do

        doc(1:6)="jacobi"
        call save_results(doc,t,cj,np,1)

        doc(1:6)="nrg_tm"
        call save_results(doc,t,e,np,1)
        doc(1:6)="nrg_ds"
        call save_results(doc,r(2,:),e,np,1)

        !tolgo la adimensionalizzazione e mi porto in AU
        y(1,:)=a*y(1,:)/au           !x
        y(2,:)=a/p*y(2,:)/au         !vx

        y(3,:)=a*y(3,:)/au           !y
        y(4,:)=a/p*y(4,:)/au         !vy

        !non tolgo la adimensionalizzazione perché t è già in unità di p
        !t=t

        doc(1:6)="pos_tm"
        call save_results(doc,t,y,np,neq)

        !dimensiono e porto in AU le distanze
       
        r(1,:)=r(1,:)*a/au
        r(2,:)=r(2,:)*a/au

        doc(1:6)="dis_tm"
        call save_results(doc,t,r,np,2)
    end do

end program main