
function flux(u_) result(f)
   implicit none
   real, intent(in)  :: u_
   real              :: f
   f = 0.5 * u_**2
   return
end function flux



function godunov(u_,v_)
   implicit none
   real, intent(in)  :: u_,v_
   real              :: godunov
   real              :: flux

   if(u_>v_) then
      godunov = min(flux(u_), flux(v_))
   else 
      godunov = max(flux(u_), flux(v_))
      if(u_*v_<=0) godunov =0
   end if
   return
end function godunov