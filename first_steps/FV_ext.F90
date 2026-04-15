
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
