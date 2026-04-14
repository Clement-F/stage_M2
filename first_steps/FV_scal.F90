program FiniteVolume

!  This program solves a FiniteVolume scalar problem 
   implicit none

   interface
      function flux(u_)
      real, intent (in) :: u_
      real              :: flux
      end function flux
   end interface
   
   interface
      function flux_p(u_)
      real, intent (in) :: u_
      real              :: flux_p
      end function flux_p
   end interface


   interface
      function rusanov(u_,v_)
      real, intent (in) :: u_,v_
      real              :: rusanov
      end function rusanov
   end interface

! loop int   
   integer:: i,j,k


!  declaration Riemann problem
   real :: ud = 1, ug = -2
! print *, 'declarer les valeurs gauche et droite du probleme de Riemann'
! read *, ud,ug 

!  domaine spatial
   integer, parameter   :: nx=100,  L=2
   real, parameter      :: dx = real(L)/nx
   real, dimension(nx)  :: X(0:nx-1) = (/ ((i+0.5)*dx, i = 0,nx-1)  /)

!  domaine temporelle
   integer, parameter   :: nt=1000,   T=1.5
   real                 :: t_=0.0,  dt=0.0

!  Flux Numeriques
   real, dimension(nx)  :: FG(0:nx-1),FD(0:nx-1)

!  Solution scalaire
   real, dimension(nx)  :: U(0:nx-1)

!  diverses valeurs numériques nécessaire
   real     :: vitesse=0., cfl=1.
   integer  :: n =0

! =======================================================================================
! =======================================================================================
! =======================================================================================

!  INIT
   where(X<real(L)/2)
      U = ug
   elsewhere 
      U = ud
   end where

   ! print *, "init"


!  boucle while sur le temps 
   do while (t_<real(T) .and. n<nt)
      n = n +1

      do i=1,nx
         if(   abs(flux_p((U(i))))>vitesse )   vitesse = abs(flux_p((U(i))))
      end do

      ! print *, "vitesse max found : ", vitesse, ";"

      if(vitesse >0) then
         dt = min(cfl * dx /(2* vitesse), T-t_)
      else
         dt = cfl*dx 
      end if
      
      ! print *, "dt : ", dt, ";"
      print *, "loop : ",n," time :",t_," ; ","dt : ",dt, ";"

      do i=0,nx-1
         
        if(i==0) then   
               FG(0) = rusanov(U(0),U(0))
        else   
               FG(i) = rusanov(U(i-1),U(i))
        end if
        
        if(i==nx-1) then   
               FD(nx-1) = rusanov(U(nx-1),U(nx-1))
        else   
               FD(i)  = rusanov(U(i),U(i+1)) 
        end if
      end do

      ! print *, "flux :"
      ! print '(f8.2)', FG
      ! print '(f8.2)', FD

      U = U - ((dt/dx)* (FD-FG) )
      t_ = t_ +dt

   end do

print '(f10.2)',U 

end program FiniteVolume