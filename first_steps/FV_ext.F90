
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