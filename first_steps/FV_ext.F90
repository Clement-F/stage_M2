
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

! 0.5*(f(u)+f(v))-0.5*(np.abs(f_p(u))*(v-u)*(u==v)+np.sign(v-u)*np.abs(f(v)-f(u))*(u!=v))

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
