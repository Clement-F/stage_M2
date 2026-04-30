
function flux(u_) result(f)
   implicit none
   real, intent(in)  :: u_
   real              :: f

   ! f = u_
   f = 0.5 * u_**2
   ! f = 4*u_**2 /(4*u_**2 + (1-u_)**2)
   ! print *,f,u_, 'f'
   return
end function flux

function flux_p(u_) result(f)
   implicit none
   real, intent(in)  :: u_
   real              :: f
   ! f = 1
   f = u_
   ! f = 8*u_* (4*u_**2 + (1-u_)**2 - u_*(4*u_-(1-u_)))/(4*u_**2 + (1-u_)**2)**2 

   return
end function flux_p


function flux_pp(u_) result(f)
   implicit none
   real, intent(in)  :: u_
   real              :: f
   f = 0
   return
end function flux_pp

function godunov(ug,ud, not_convex)
   implicit none
   real, intent(in)  :: ug,ud
   real              :: godunov
   real              :: flux

   logical, intent(in), optional :: not_convex
   real              :: step
   integer           :: i 
   real, dimension(0:100):: X

   if( (.not. present(not_convex)) .or. (.not. not_convex) ) then
      if(ug<ud) then
         godunov = min(flux(ug), flux(ud))
         if(ug*ud<=0) godunov =0
      else 
         godunov = max(flux(ug), flux(ud))
      end if
   else 
      if(ug<ud) then ! -> min_{u \in [ug,ud]} f(u)
         step = (ud-ug)/100
         X = (/  (ug+ (i)*step, i = 0,100)  /)
         godunov = flux(ug)

         do i=0,100
            ! print *,ug,ud,X(i), flux(X(i)), godunov, 'min'
            if(flux(X(i))<godunov) godunov = flux(X(i))
         end do

      else           !  -> max_{u \in [ug,ud]} f(u)
         step = -(ud-ug)/100
         X = (/  (ug+ (i)*step, i = 0,100)  /)
         godunov = flux(ug)

         do i=0,100
            ! if(godunov>0) print *,ug,ud,X(i), flux(X(i)), godunov, 'max'
            if(flux(X(i))>godunov) godunov = flux(X(i))
         end do

      end if
   end if
   

   return
end function godunov

function rusanov(ug,ud)
   implicit none
   real, intent(in)  :: ug,ud
   real              :: rusanov
   real              :: flux, flux_p

   rusanov = 0.5*(flux(ug)+flux(ud))-0.5* max(abs(flux_p(ug)), abs(flux_p(ud))) *(ud-ug)

   return
end function rusanov

function Roe(ug,ud)
   implicit none
   real, intent(in)  :: ug,ud
   real              :: Roe
   real              :: flux, flux_p

    Roe = 0.5*(flux(ug)+flux(ud))
    if(abs(ug-ud)<10e-10) then
    Roe = Roe -0.5* (abs(flux_p(ug))*(ud-ug))
    else
    Roe = Roe - 0.5*sign(abs(flux(ud)-flux(ug)),ud-ug)
    end if
   return
end function Roe

function Lax_Friedrichs(ug,ud)
   implicit none
   real, intent(in)  :: ug,ud
   real              :: Lax_Friedrichs
   real              :: flux, flux_p
   real              :: beta

   
   beta = max(abs(flux_p(ug)), abs(flux_p(ud)))

   Lax_Friedrichs = 0.5* (flux(ug)+flux(ud) - beta*(ud -ug ))
   return 

end function Lax_Friedrichs
   
function Q_init(x) result(Q)
    implicit none
    real, intent(in)    :: x
    real                :: Q
    real,parameter      :: pi = acos(-1.)
    

    Q = sin(2*pi* x)
   ! Q = 0
   ! if(x<0.5) Q = x
   ! if(x>=0.5)Q = 1 - x
   
   ! if(x<0.) Q=0.
   ! if(x>1.) Q=0.
   

   return 
end function Q_init

function Q_init_p(x) result(Q)
    implicit none
    real, intent(in)    :: x
    real                :: Q
    real,parameter      :: pi = acos(-1.)

    Q = 2*pi* cos(2*pi* x)
   ! Q = 0 
   ! if(x<0.5) Q = 1
   ! if(x>=0.5)Q = -1
    return 
end function Q_init_p

function Newton_search(x,t) result(q)
    implicit none
    real, intent(in)    :: x,t
    real                :: Q_init, Q_init_p
    real                :: flux_p, flux_pp
    real                :: xk, err, q, epsi = 1e-20
    real,parameter      :: pi = acos(-1.)
    integer             :: n=0

    n = 0
    err = abs(flux_p(Q_init(xk))*t+ xk-x)
    ! print *, err, epsi
    xk = x

    do while(err>epsi .and. n<50)
        ! print *, n, err
        xk = xk -   (flux_p(Q_init(xk))*t + xk-x)/(flux_pp(Q_init(xk))*Q_init_p(xk)*t +1)
        err =    abs(flux_p(Q_init(xk))*t + xk-x)
      !   if(t>1./(2*pi)) then
            ! if(x<0.5 .and. xk>0.5) then 
            !     xk =0.5 - 1e-6
            !     ! print *, "rectification gauche "
            ! end if
            ! if(x>0.5 .and. xk<0.5) then
            !     xk =0.5 + 1e-6
            !     ! print *, "rectification droite "
            ! end if
      !   end if
        n = n+1
    end do
    

    q = Q_init(xk)
    return

end function Newton_search

function dicho (fct,xd,xf)

   real, intent(in) :: xd,xf
   real :: dicho

   interface
   function fct(x);   real,intent(in) :: x; real fct; end function fct
   end interface

   real :: a,b,c
   real :: t
   integer :: n = 0

   if(abs(fct(xf))<1e-20) then 
      dicho = xf
      return
   else if(abs(fct(xd))<1e-20) then
      dicho = xd
      return
   end if


   a=xd; b=xf; 

   n =0
   do while(abs(fct(c))>1e-20 .and. n<100 )
      c=(a+b)/2.
      t = sign(1.,fct(c)*fct(a))
      if(t <0) then 
         b=c
      else 
         a=c
      end if
      n = n+1
   end do

   dicho = c
   return

end function dicho

subroutine pied_charact(x,t,sol)

   real, intent(in) :: x,t
   real, intent(out):: sol
   real  :: Q_init
   ! print *,'pied caract'

   sol = Q_init(dicho(g,0.,1.))

   ! if(t<1./(2*pi)-1e-2) then
   !    if(x<=0.5)  sol = Q_init(dicho(g,0.,0.5))
   !    if(x>0.5)   sol = Q_init(dicho(g,0.5,1.))
   ! else       
   !    if(x<=0.5)  sol = Q_init(dicho(g,0.,0.5-1e-3))
   !    if(x>0.5)   sol = Q_init(dicho(g,0.5+1e-3,1.))
   ! end if
   contains
   function g(x_)
      real,intent(in) :: x_
      real  :: flux_p, Q_init
      real  :: g
      g = flux_p(Q_init(x_))*t + x_ -x
      return 
   end function g

end subroutine pied_charact

function minmod(x,y,z)

   implicit none
   real, intent(in)    :: x,y,z
   real                :: minmod

   minmod = min(0.,max(x,y,z)) + max(0.,min(x,y,z))
   return

end function


subroutine Update(Q,X,dt,nx, arg_string)
   implicit none
   integer,                intent(in)  :: nx
   real, dimension(0:nx-1),intent(inout)  :: Q
   real, dimension(0:nx+1),intent(in)  :: X

   real,                   intent(in) ::  dt
   character(len=32),      intent(in),optional :: arg_string

   character(len=32) :: methode_update
   integer :: i 
   real, dimension(0:nx) :: F
   real, dimension(:), Pointer :: Q_int
   real, dimension(:), Pointer :: delta
   real                  :: dx
   real  :: minmod
   real  :: alpha 

    interface
      function godunov(ug,ud, not_convex)
      real, intent (in) :: ug,ud
      logical, intent(in), optional :: not_convex
      real              :: godunov
      end function godunov

      
      function Rusanov(ug,ud)
      real, intent (in) :: ug,ud
      real              :: Rusanov
      end function Rusanov

      
      function Roe(ug,ud)
      real, intent (in) :: ug,ud
      real              :: Roe
      end function Roe

      function Lax_Friedrichs(ug,ud)
         implicit none
         real, intent(in)  :: ug,ud
         real              :: Lax_Friedrichs
      end function Lax_Friedrichs

   end interface




   dx = X(2)-X(1)

   if(present(arg_string)) methode_update = arg_string

   if(methode_update == "godunov") then
   

      do i=0,nx
        if(i==0) then;        F(0)  = godunov(Q(0),Q(0),       not_convex=.false.)
        else if (i==nx) then; F(nx) = godunov(Q(nx-1),Q(nx-1), not_convex=.false.)
        else;                 F(i)  = godunov(Q(i-1),Q(i),     not_convex=.false.)
        end if
      end do

      do i=0,nx-1
         Q(i) = Q(i) - ((dt/dx)* (X(i+1)-X(i)) )* (F(i+1)-F(i))
      end do

   else if(methode_update == "Rusanov") then
   
      do i=0,nx
        if(i==0) then;        F(0)  = Rusanov(Q(0),Q(0))
        else if (i==nx) then; F(nx) = Rusanov(Q(nx-1),Q(nx-1))
        else;                 F(i)  = Rusanov(Q(i-1),Q(i))
        end if
      end do

      Q(:) = Q(:) - ((dt/dx)* (F(1:nx)-F(0:nx-1)) )

   else if(methode_update == "Roe") then 
   
      do i=0,nx
        if(i==0) then;        F(0)  = Roe(Q(0),Q(0))
        else if (i==nx) then; F(nx) = Roe(Q(nx-1),Q(nx-1))
        else;                 F(i)  = Roe(Q(i-1),Q(i))
        end if
      end do
      
      Q(:) = Q(:) - ((dt/dx)* (F(1:nx)-F(0:nx-1)) )

   else if(methode_update == "Lax") then
      
      do i=0,nx
        if(i==0) then;        F(0)  = Lax_Friedrichs(Q(0),Q(0))
        else if (i==nx) then; F(nx) = Lax_Friedrichs(Q(nx-1),Q(nx-1))
        else;                 F(i)  = Lax_Friedrichs(Q(i-1),Q(i))
        end if
      end do
      
      Q(:) = Q(:) - ((dt/dx)* (F(1:nx)-F(0:nx-1)) )

   else if(methode_update == "reconstruction") then

      allocate(Q_int(0:nx-1))
      allocate(delta(0:nx-1))

      ! print *,"methode : reconstruction"

      delta(0) = (Q(1)-Q(0))/(2*dx);  delta(nx-1) = (Q(nx-1)-Q(nx-2))/(2*dx); 
      do i=1,nx-2
         delta(i) = (Q(i+1)-Q(i-1))/(2*dx)
      end do

      ! print *,"delta construit"

      do i=0,nx
         ! print *,i
      !   if(i==0) then;        F(0)  = Lax_Friedrichs(Q(i)+0.5*dx*delta(i)        ,Q(i)-0.5*dx*delta(i))
      !   else if (i==nx) then; F(nx) = Lax_Friedrichs(Q(i-1)+0.5*dx*delta(i-1)    ,Q(i-1)-0.5*dx*delta(i-1))


        if(i==0) then;        F(0)  = Lax_Friedrichs(Q(nx-1)+0.5*dx*delta(nx-1)  ,Q(i)-0.5*dx*delta(i))
        else if (i==nx) then; F(nx) = Lax_Friedrichs(Q(i-1)+0.5*dx*delta(i-1)    ,Q(0)-0.5*dx*delta(0))



        else;                 F(i)  = Lax_Friedrichs(Q(i-1)+0.5*dx*delta(i-1)    ,Q(i)-0.5*dx*delta(i))
        end if
      end do

      ! print *,"Flux construit"

      Q_int(:) = Q(:) - ((dt/(2*dx))* (F(1:nx)-F(0:nx-1)))
      
      ! print *,"sol 1/2 update"

      delta(0) = (Q(1)-Q(0))/(2*dx);  delta(nx-1) = (Q(nx-1)-Q(nx-2))/(2*dx); 
      do i=1,nx-2
         delta(i) = (Q(i+1)-Q(i-1))/(2*dx)
      end do

      do i=0,nx
      !   if(i==0) then;        F(0)  = Lax_Friedrichs(Q_int(i)+0.5*dx*delta(i)        ,Q_int(i)-0.5*dx*delta(i))
      !   else if (i==nx) then; F(nx) = Lax_Friedrichs(Q_int(i-1)+0.5*dx*delta(i-1)    ,Q_int(i)-0.5*dx*delta(i))

        
        if(i==0) then;        F(0)  = Lax_Friedrichs(Q_int(nx-1)+0.5*dx*delta(nx-1)  ,Q_int(i)-0.5*dx*delta(i))
        else if (i==nx) then; F(nx) = Lax_Friedrichs(Q_int(i-1)+0.5*dx*delta(i-1)    ,Q_int(0)-0.5*dx*delta(0))

        else;                 F(i)  = Lax_Friedrichs(Q_int(i-1)+0.5*dx*delta(i-1)    ,Q_int(i)-0.5*dx*delta(i))
        end if
      end do

      Q(:) = Q(:) - ((dt/dx)* (F(1:nx)-F(0:nx-1)))

      ! print *,"full sol update"

   else if(methode_update == "limitation") then

      allocate(Q_int(0:nx-1))
      allocate(delta(0:nx-1))
      alpha = 0.5

      ! print *,"methode : reconstruction"

      delta(0)    = (Q(1)-Q(0))/(2*dx);  
      delta(nx-1) = (Q(nx-1)-Q(nx-2))/(2*dx); 
      do i=1,nx-2
         delta(i) = (Q(i+1)-Q(i-1))/(2*dx)
      end do

      delta(0)    = minmod(delta(0),   2*alpha* (Q(0)-Q(0))/dx,      2*alpha* (Q(1)-Q(0))/dx)
      delta(nx-1) = minmod(delta(nx-1),2*alpha* (Q(nx-1)-Q(nx-2))/dx,2*alpha* (Q(nx-1)-Q(nx-1))/dx)

      do i=1,nx-2; delta(i) = minmod(delta(i),2*alpha* (Q(i)-Q(i-1))/dx,2*alpha* (Q(i+1)-Q(i))/dx); end do

      ! print *,"delta construit"

      do i=0,nx
         ! print *,i
      !   if(i==0) then;        F(0)  = Lax_Friedrichs(Q(i)+0.5*dx*delta(i)        ,Q(i)-0.5*dx*delta(i))
      !   else if (i==nx) then; F(nx) = Lax_Friedrichs(Q(i-1)+0.5*dx*delta(i-1)    ,Q(i-1)-0.5*dx*delta(i-1))


        if(i==0) then;        F(0)  = Lax_Friedrichs(Q(nx-1)+0.5*dx*delta(nx-1)  ,Q(i)-0.5*dx*delta(i))
        else if (i==nx) then; F(nx) = Lax_Friedrichs(Q(i-1)+0.5*dx*delta(i-1)    ,Q(0)-0.5*dx*delta(0))



        else;                 F(i)  = Lax_Friedrichs(Q(i-1)+0.5*dx*delta(i-1)    ,Q(i)-0.5*dx*delta(i))
        end if
      end do

      ! print *,"Flux construit"

      Q_int(:) = Q(:) - ((dt/(2*dx))* (F(1:nx)-F(0:nx-1)))
      
      ! print *,"sol 1/2 update"

      delta(0) = (Q(1)-Q(0))/(2*dx);  delta(nx-1) = (Q(nx-1)-Q(nx-2))/(2*dx); 
      do i=1,nx-2
         delta(i) = (Q(i+1)-Q(i-1))/(2*dx)
      end do
      
      delta(0)    = minmod(delta(0),   2*alpha* (Q(0)-Q(0))/dx,      2*alpha* (Q(1)-Q(0))/dx)
      delta(nx-1) = minmod(delta(nx-1),2*alpha* (Q(nx-1)-Q(nx-2))/dx,2*alpha* (Q(nx-1)-Q(nx-1))/dx)

      do i=1,nx-2; delta(i) = minmod(delta(i),2*alpha* (Q(i)-Q(i-1))/dx,2*alpha* (Q(i+1)-Q(i))/dx); end do

      do i=0,nx
      !   if(i==0) then;        F(0)  = Lax_Friedrichs(Q_int(i)+0.5*dx*delta(i)        ,Q_int(i)-0.5*dx*delta(i))
      !   else if (i==nx) then; F(nx) = Lax_Friedrichs(Q_int(i-1)+0.5*dx*delta(i-1)    ,Q_int(i)-0.5*dx*delta(i))

        
        if(i==0) then;        F(0)  = Lax_Friedrichs(Q_int(nx-1)+0.5*dx*delta(nx-1)  ,Q_int(i)-0.5*dx*delta(i))
        else if (i==nx) then; F(nx) = Lax_Friedrichs(Q_int(i-1)+0.5*dx*delta(i-1)    ,Q_int(0)-0.5*dx*delta(0))

        else;                 F(i)  = Lax_Friedrichs(Q_int(i-1)+0.5*dx*delta(i-1)    ,Q_int(i)-0.5*dx*delta(i))
        end if
      end do

      Q(:) = Q(:) - ((dt/dx)* (F(1:nx)-F(0:nx-1)))

   
   
   end if

end subroutine