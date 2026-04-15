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
      function godunov(u_,v_)
      real, intent (in) :: u_,v_
      real              :: godunov
      end function godunov
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
   integer, parameter   :: nt=500,   T=1.5
   real                 :: t_=0.0,  dt=0.0

!  Flux Numeriques
   real, dimension(nx)  :: FG(0:nx-1),FD(0:nx-1)

!  Solution scalaire
   real, dimension(nx)  :: U(0:nx-1)

!  diverses valeurs numériques nécessaire
   real     :: vitesse=0., cfl=1.
   integer  :: n =0

!  file parameter
   integer, parameter   :: numfile_sol=1, numfile_data=2
   character(len=14)    :: nomfile_sol = 'file_sol.txt',   nomfile_data = 'file_data.txt'
   character(len=nx +6) :: save_format = '(100(f8.4, 4x))'
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
      if(mod(n,10) ==0)  print *, "loop : ",n," time :",t_," ; ","dt : ",dt, ";"

      open(unit=numfile_sol, file=nomfile_sol, form ='formatted', status ='old')

      do i=0,nx-1
         
        if(i==0) then   
               FG(0) = godunov(U(0),U(0))
        else   
               FG(i) = godunov(U(i-1),U(i))
        end if
        
        if(i==nx-1) then   
               FD(nx-1) = godunov(U(nx-1),U(nx-1))
        else   
               FD(i)  = godunov(U(i),U(i+1)) 
        end if
      end do

      ! print *, "flux :"
      ! print '(f8.2)', FG
      ! print '(f8.2)', FD

      U = U - ((dt/dx)* (FD-FG) )
      write(unit=numfile_sol, fmt=save_format) U
      
      t_ = t_ +dt

   end do

print save_format,U 

end program FiniteVolume