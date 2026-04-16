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

!  domaine spatial
   integer, parameter   :: nx=8000,  L=1
   real, parameter      :: dx = real(L)/nx
   real, dimension(nx)  :: X(0:nx-1) = (/ ((i+0.5)*dx, i = 0,nx-1)  /)

!  domaine temporelle
   integer, parameter   :: nt=500,   T=1
   real                 :: t_=0.0,  dt=0.0

!  Flux Numeriques
   real, dimension(nx)  :: FG(0:nx-1),FD(0:nx-1)
   real, dimension(nx+1):: F(0:nx)

!  Solution scalaire
   real, dimension(nx)  :: U(0:nx-1)

!  diverses valeurs numériques nécessaire
   real     :: vitesse=0., cfl=1.
   integer  :: n =0

!  file parameter
   integer, parameter   :: numfile_sol=1, numfile_data=2
   integer              :: n_imp=0
   real,dimension (2,nx):: sol
   real                 :: t_imp=real(T)/10
   character(len=32)    :: nomfile_sol = 'file_sol.txt',   nomfile_data = 'file_data.txt'
   character(len=32)    :: str,save_format
! =======================================================================================
! =======================================================================================
! =======================================================================================

!  INIT
   where(X<0.3) U=0
   where(X>0.7) U=0.5
   where(X>0.3 .and. X<0.7) U=-1

   ! U = sin(2*3.1415 *X)

   ! print *, "init"
   
! print *, 'declarer les valeurs gauche et droite du probleme de Riemann'
! read *, ud,ug 

   write(save_format, '( "(" i5 "(f10.6, 1x, f12.8 /) )" )') nx 

   print *,save_format

   open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')
   
   write(unit= numfile_data, fmt='("nt = "i5)') nt
   write(unit= numfile_data, fmt='("nx = "i5)') nx
   write(unit= numfile_data, fmt='("save_max =" i5)')  int(T/t_imp)

!  boucle while sur le temps 
   do while (t_<real(T))
      n = n +1

      vitesse =0
      do i=1,nx
         if(   abs(flux_p((U(i))))>vitesse )   vitesse = abs(flux_p((U(i))))
      end do

      ! print *, "vitesse max found : ", vitesse, ";"

      if(vitesse >10e-10) then
         dt = min(cfl * dx /(2* vitesse), T-t_)
      else
         dt = cfl*dx 
      end if
      

      open(unit=numfile_sol, file=nomfile_sol, form ='formatted', status ='old')

      do i=0,nx
         
        if(i==0) then   
               F(0)  = godunov(U(0),U(0))
        else if (i==nx) then   
               F(nx) = godunov(U(nx-1),U(nx-1))
        else   
               F(i)  = godunov(U(i-1),U(i))
        end if
        
      end do

      U(:) = U(:) - ((dt/dx)* (F(1:nx)-F(0:nx-1)) )

      if(t_>   n_imp*t_imp)  then
         print *, "loop : ",n,", n_imp",n_imp,", time :",t_," ; ","dt : ",dt, ";"
         n_imp = n_imp +1
         sol(1,:)=X
         sol(2,:)=U
         write(unit=numfile_sol,  fmt=save_format) sol
         write(unit=numfile_data, fmt='("time_save =" f10.6)')  t_
      end if
      
      t_ = t_ +dt


   end do

   print *, "program complete !"

end program FiniteVolume