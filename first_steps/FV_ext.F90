
function flux(u_) result(f)
   implicit none
   real, intent(in)  :: u_
   real              :: f
   f = 0.5 * u_**2
   return
end function flux

function flux_p(u_) result(f)
   implicit none
   real, intent(in)  :: u_
   real              :: f
   f = u_
   return
end function flux_p


function flux_pp(u_) result(f)
   implicit none
   real, intent(in)  :: u_
   real              :: f
   f = 1
   return
end function flux_pp

function godunov(u_,v_)
   implicit none
   real, intent(in)  :: u_,v_
   real              :: godunov
   real              :: flux

   if(u_<v_) then
      godunov = min(flux(u_), flux(v_))
      if(u_*v_<=0) godunov =0
   else 
      godunov = max(flux(u_), flux(v_))
   end if

   return
end function godunov


function rusanov(u_,v_)
   implicit none
   real, intent(in)  :: u_,v_
   real              :: rusanov
   real              :: flux, flux_p

   rusanov = 0.5*(flux(u_)+flux(v_))-0.5* max(abs(flux_p(u_)), abs(flux_p(v_))) *(v_-u_)

   return
end function rusanov

function Roe(u_,v_)
   implicit none
   real, intent(in)  :: u_,v_
   real              :: Roe
   real              :: flux, flux_p

    Roe = 0.5*(flux(u_)+flux(v_))
    if(abs(u_-v_)<10e-10) then
    Roe = Roe -0.5* (abs(flux_p(u_))*(v_-u_))
    else
    Roe = Roe - 0.5*sign(abs(flux(v_)-flux(u_)),v_-u_)
    end if
   return
end function Roe

function U_init(x) result(U)
    implicit none
    real, intent(in)    :: x
    real                :: U

    U = sin(2*3.1415* x)
    return 
end function U_init


function U_init_p(x) result(U)
    implicit none
    real, intent(in)    :: x
    real                :: U

    U = 2*3.1415* cos(2*3.1415* x)
    return 
end function U_init_p

function Newton_search(x,t) result(u)
    implicit none
    real, intent(in)    :: x,t
    real                :: U_init, U_init_p
    real                :: flux_p, flux_pp
    real                :: xk, err, u, epsi = 1e-10
    real,parameter      :: pi = acos(-1.)
    integer             :: n=0

    n = 0
    err = abs(flux_p(U_init(xk))*t+ xk-x)
    ! print *, err, epsi
    xk = x

    do while(err>epsi .and. n<50)
        ! print *, n, err
        xk = xk -   (flux_p(U_init(xk))*t + xk-x)/(flux_pp(U_init(xk))*U_init_p(xk)*t +1)
        err =    abs(flux_p(U_init(xk))*t + xk-x)
        ! if(t>1./(2*pi)) then
            if(x<0.5 .and. xk>0.5) then 
                xk =0.5 - 1e-6
                ! print *, "rectification gauche "
            end if
            if(x>0.5 .and. xk<0.5) then
                xk =0.5 + 1e-6
                ! print *, "rectification droite "
            end if
        ! end if
        n = n+1
    end do
    
    ! print *, n, err

    u = U_init(xk)
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

   if(abs(fct(xf))<1e-6) then 
      dicho = xf
      return
   else if(abs(fct(xd))<1e-6) then
      dicho = xd
      return
   end if


   a=xd; b=xf; 

   n =0
   ! print*, 'new guess'
   do while(abs(fct(c))>1e-6 .and. n<100 )
      c=(a+b)/2.
      t = sign(1.,fct(c)*fct(a))
      ! print *,t
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
   real  :: U_init
   ! print *,'pied caract'
   if(t<1./(2*pi)-1e-2) then
      if(x<=0.5)  sol = U_init(dicho(g,0.,0.5))
      if(x>0.5)   sol = U_init(dicho(g,0.5,1.))
   else       
      if(x<=0.5)  sol = U_init(dicho(g,0.,0.5-1e-3))
      if(x>0.5)   sol = U_init(dicho(g,0.5+1e-3,1.))
   end if
   contains
   function g(x_)
      real,intent(in) :: x_
      real  :: flux_p, U_init
      real  :: g
      g = flux_p(U_init(x_))*t + x_ -x
      ! print *,'g=',g
      return 
   end function g

end subroutine pied_charact
