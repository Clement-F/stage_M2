
MODULE mod_declaration
  USE precis
  IMPLICIT NONE

  TYPE var_type
     REAL(prec),DIMENSION(:),POINTER :: base_poly
     REAL(prec),DIMENSION(2) :: inter
  END TYPE var_type

  TYPE(var_type), DIMENSION(:), POINTER :: sol, sol_step, flux_h

  REAL(prec), DIMENSION(:),   POINTER :: x_cell, x_middle, cell_size
  REAL(prec), DIMENSION(:),   POINTER :: x_quad, w_quad
  REAL(prec), DIMENSION(:,:), POINTER :: coeff_DG

  REAL(prec), DIMENSION(:,:), POINTER :: coeff_Taylor, coeff_legendre

  REAL(prec), DIMENSION(:,:), POINTER :: RK_alpha
  REAL(prec), DIMENSION(:),   POINTER :: RK_time,RK_beta
  REAL(prec), DIMENSION(:), POINTER :: L_step
  REAL(prec), DIMENSION(:),   POINTER :: Time_stemp

  REAL(prec), DIMENSION(:,:), POINTER :: Masse, Masse_inv, Rigid, Rigid_inv

  INTEGER :: nb_cell, order_x, order_t
  INTEGER :: n_imp, n_time, frame_rule, print_rule
  REAL(prec) :: time, tmax,t_ini,dt, t_imp
  REAL(prec) :: xL,xR,CFL,max_dflux,dx
  REAL(prec) :: eps0 =0._prec
  REAL(prec) :: vit_adv 

  CHARACTER(LEN=32) :: DG_meth, quad_meth, sol_ini_name, flux_name

  INTEGER   :: numfile_sol=1, numfile_param=2, numfile_data = 3
  CHARACTER(len=32)    :: nomfile_sol = 'file_sol.txt', nomfile_param = 'param.txt', nomfile_data= 'file_data.txt'

  INTEGER :: i,j


  CONTAINS

  SUBROUTINE ALLOCATE_all
    IMPLICIT NONE

    ALLOCATE( sol(nb_cell), sol_step(nb_cell), flux_h(nb_cell))
    ALLOCATE( x_cell(nb_cell+1), x_middle(nb_cell), cell_size(nb_cell)) 
    ALLOCATE( x_quad(order_x), w_quad(order_x) ) 
    ALLOCATE( coeff_DG(order_x, order_x) ) 

    ALLOCATE( coeff_Taylor(10,10), coeff_legendre(10,10) )
    ALLOCATE( RK_alpha(order_t,order_t), RK_beta(order_t), RK_time(order_t) )
    ALLOCATE( L_step(order_x))
    ALLOCATE( Time_stemp(print_rule+1))

    ALLOCATE( Masse(order_x,order_x), Masse_inv(order_x,order_x) ) 
    ALLOCATE( Rigid(order_x,order_x), Rigid_inv(order_x,order_x) ) 
    
    DO i =1,nb_cell
      ALLOCATE(sol(i)%base_poly(order_x))
      ALLOCATE(sol_step(i)%base_poly(order_x))
      ALLOCATE(flux_h(i)%base_poly(order_x))
    END DO
    

  END SUBROUTINE ALLOCATE_all

  SUBROUTINE DEALLOCATE_all
    IMPLICIT NONE

    DO i =1,nb_cell
        DEALLOCATE(sol(i)%base_poly)
    END DO
    DEALLOCATE( sol)
    DEALLOCATE( x_cell, x_middle, cell_size) 
    DEALLOCATE( x_quad, w_quad) 
    DEALLOCATE( coeff_DG) 

    DEALLOCATE( coeff_Taylor, coeff_legendre )
    ! ALLOCATE( RK_alpha, RK_beta )

    DEALLOCATE( Masse, Masse_inv ) 
    
  END SUBROUTINE DEALLOCATE_all
  
  SUBROUTINE Coeff_RK_init(nrk)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nrk
    INTEGER :: nrk2

    RK_time   = 0._prec
    RK_alpha  = 0._prec
    RK_beta   = 0._prec

    SELECT CASE (nrk)
    CASE (1)
       RK_time(1)=0._prec
       RK_alpha(1,1)=1._prec
       RK_beta(1)=1._prec
    CASE (2)
       RK_time(1)=0._prec
       RK_alpha(1,1)=1._prec
       RK_beta(1)=1._prec

       RK_time(2)=1._prec
       RK_alpha(2,1)=0.5_prec
       RK_alpha(2,2)=0.5_prec
       RK_beta(2)=0.5_prec
    CASE (3)
       RK_time(1)=0._prec
       RK_alpha(1,1)=1._prec
       RK_beta(1)=1._prec

       RK_time(2)=0.5_prec
       RK_alpha(2,1)=0.75_prec
       RK_alpha(2,2)=0.25_prec
       RK_beta(2)=0.25_prec

       RK_time(3)=1._prec
       RK_alpha(3,1)=1._prec/3._prec
       RK_alpha(3,2)=2._prec/3._prec
       RK_beta(3)=2._prec/3._prec
    CASE DEFAULT
       RK_time=0._prec; RK_alpha = 0._prec; RK_beta = 0._prec
      !  DO nrk2=1,nrk
      !     RK_time(nrk2)=1._prec
      !     RK_beta(nrk2)=1._prec/REAL(nrk+1-nrk2,prec)
      !  END DO
    END SELECT


  END SUBROUTINE Coeff_RK_init

END MODULE mod_declaration 