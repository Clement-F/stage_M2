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
      function godunov(u_,v_)
      real, intent (in) :: u_,v_
      real              :: godunov
      end function godunov
   end interface

! loop int   
   integer:: i,j,k


!  declaration Riemann problem
   real :: ud = 1, ug = -1
! print *, 'declarer les valeurs gauche et droite du probleme de Riemann'
! read *, ud,ug 

!  domaine spatial
   integer, parameter   :: nx=100,  L=2
   real, parameter      :: dx = real(L)/nx
   real, dimension(nx) :: X = (/ ((i+0.5)*dx, i = 1,nx)  /)

!  domaine temporelle
   integer, parameter   :: nt=10,   T=1
   real,                :: t=0., dt_ref = real(T)/nt, dt

!  Flux Numeriques
   real, dimension(nx)  :: FG,FD

!  Solution scalaire
   real, dimension(nx)  :: U

!  INIT
   where(X<real(L)/2)
      U = ud
   elsewhere 
      U = ug
   end where

!  diverses valeurs numériques nécessaire
   real                 :: vitesse


!  boucle while sur le temps 
   do while (t<real(T))
      vitesse = max(abs(flux(U)))
   end do


end program FiniteVolume