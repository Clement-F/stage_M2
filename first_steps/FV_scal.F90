program FiniteVolume

!  This program solves a FiniteVolume scalar problem 
   implicit none

   interface
      function flux(u_)
      real, intent (in) :: u_
      real              :: flux
      end function flux

      function flux_p(u_)
      real, intent (in) :: u_
      real              :: flux_p
      end function flux_p 

      function Q_init(x_)
      real,dimension(:), intent (in) :: x_
      real              :: Q_init
      end function Q_init

      subroutine Update(Q,X,dt,nx, arg_string)
         implicit none
         integer,intent(in)  :: nx
         real, dimension(0:nx-1), intent(inout)  :: Q
         real, dimension(0:nx-1), intent(in)  :: X
         real, intent(in) ::  dt
         character(len=32), intent(in), optional :: arg_string
      end subroutine

      function Newton_search(x,t)
      real, intent (in) :: x,t
      real              :: Newton_search
      end function Newton_search
   end interface




! loop int   
   integer:: i,j,k
   character(len=20) :: c1

!  constants
   real, parameter      :: pi = acos(-1.)

!  domaine spatial
   integer              :: nx
   real                 :: dx, xd, xf
   real, dimension(:),  Pointer :: X
!  domaine temporelle
   integer, parameter   :: nt=500
   real       :: T
   real                 :: t_=0.0,  dt=0.0

!  Flux Numeriques
   real, dimension(:),  Pointer :: F

!  Solution scalaire
   real, dimension(:),  Pointer :: Q, Q_ex
   real,dimension (:,:),Pointer :: sol

!  diverses valeurs numériques nécessaire
   real     :: vitesse=0., cfl
   integer  :: n =0

!  error variables
   real     :: err_L1=0, err_L2=0, err_Li=0

!  probleme variable
   character(len=32)    :: methode_update = "godunov" 


!  file parameter
   integer, parameter   :: numfile_sol=1, numfile_data=2, numfile_param=3, numfile_err =4, numfile_conv =5
   integer              :: n_imp=0, n_imp_max
   real                 :: t_imp
   character(len=32)    :: nomfile_sol = 'file_sol.txt',   nomfile_data = 'file_data.txt', nomfile_param = 'param.txt', nomfile_err = 'error_file.txt', &
                           nomfile_conv = 'convergence_err.txt'
   character(len=32)    :: str,save_format
! =======================================================================================
! =======================================================================================
! =======================================================================================
   
   open(unit=numfile_param, file=nomfile_param, form ='formatted', status ='old')

   read(numfile_param,  *) nx;  
   read(numfile_param,  *) xd;  
   read(numfile_param,  *) xf;      
   read(numfile_param,  *) T;     
   read(numfile_param,  *) cfl;   
   read(numfile_param,  *) n_imp_max;     

   dx = real((xf-xd))/nx

   allocate(X(0:nx+1));   X(1:nx) = (/  (xd+ i*dx, i = 1,nx)  /); X(0) = xd; X(nx+1) = xf
   allocate(F(0:nx))
   allocate(Q(0:nx-1),  Q_ex(0:nx-1))
   allocate(sol (3,nx))

   t_imp=T/real(n_imp_max)

   ! Q =0
   ! where (X>-0.5 .and. X<0) Q =1
   Q = sin(2*pi*( X(1:nx)+X(0:nx-1) )/2 ); Q_ex = 0

   ! print *, "init"
   
   write(save_format, '("(" i5 "(f10.6, f12.8, f12.8 /))")') nx 

   open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')
   
   write(unit= numfile_data, fmt='("nt = "i5)') nt
   write(unit= numfile_data, fmt='("nx = "i5)') nx
   write(unit= numfile_data, fmt='("save_max =" i5)')  int(T/t_imp)
   
   open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')
   open(unit=numfile_err,  file=nomfile_err, form ='formatted', status ='old')

!  boucle while sur le temps 
   do while (t_<real(T))
      n = n +1

      vitesse =0
      do i=1,nx
         if(   abs(flux_p((Q(i))))>vitesse )   vitesse = abs(flux_p((Q(i))))
      end do

      ! print *, "vitesse max found : ", vitesse, ";"

      if(vitesse >10e-10) then
         dt = min(cfl * dx /(2* vitesse), T-t_)
      else
         dt = cfl*dx 
      end if
      
      call Update(Q=Q, X=X,dt=dt,nx =nx, arg_string = methode_update) 

      if(t_ >=  n_imp*t_imp)  then

         print *, "exact sol calcul"

         ! do i=0,nx-1
         !    ! Q_ex(i) = Newton_search(X(i),t_)
         !    call pied_charact(X(i),t_,Q_ex(i))
         ! end do

         print *, "exact sol calculated"

         print *, "loop : ",n,", n_imp",n_imp,", time :",t_," ; ","dt : ",dt, ";"
         n_imp = n_imp +1
         sol(1,:)=X(0:nx-1)
         sol(2,:)=Q(0:nx-1)
         sol(3,:)=Q_ex(0:nx-1)
         write(unit=numfile_sol,  fmt=save_format) sol

         if(sum( abs(Q-Q_ex))*dx > err_L1) err_L1 = sum( abs(Q-Q_ex))*dx 
         if(sum( (Q-Q_ex)**2)*dx > err_L2) err_L2 = sum( (Q-Q_ex)**2)*dx 


         write(unit=numfile_err, fmt='(" --------------- at time : "f10.6" ----------------- ")') t_ 
         write(unit=numfile_err, fmt='("err_L1 :" f16.10 )') sum( abs(Q-Q_ex))*dx 
         write(unit=numfile_err, fmt='("err_L2 :" f16.10 )') sum( (Q-Q_ex)**2)*dx   
         write(unit=numfile_data, fmt='("time_save =" f10.6)')  t_
      end if
      
      t_ = t_ +dt


   end do
   
   close(unit=numfile_data)
   close(unit=numfile_err)
   close(unit=numfile_param)
   close(unit=numfile_sol)

   open(unit=numfile_conv,  file=nomfile_conv, form ='formatted', status ='old', position='append')
   write(unit=numfile_conv, fmt='("=====================")') 
   write(unit=numfile_conv, fmt='("for nx = "i5" we have error :")' ) nx
   write(unit=numfile_conv, fmt='("err_L1 :" f16.10 )') sum( abs(Q-Q_ex))*dx 
   write(unit=numfile_conv, fmt='("err_L2 :" f16.10 )') sum( (Q-Q_ex)**2)*dx  
   write(unit=numfile_conv, fmt='("=====================")') 

   print *, "program complete !"

end program FiniteVolume