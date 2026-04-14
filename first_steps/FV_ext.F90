
function flux(u_)
   implicit none
   real, intent(in)  :: u_
   real :: flux
   flux = u_**2
end function flux



function godunov(u_,v_)
   implicit none
   real, intent(in)  :: u_,v_
   if(u_>v_) then
      godunov = min(flux(u_), flux(v_))
   else 
      godunov = max(flux(u_), flux(v_))
   end if
end function godunov