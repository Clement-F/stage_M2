
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

  REAL(prec), DIMENSION(:,:), POINTER :: RK_alpha, RK_beta
  REAL(prec), DIMENSION(:),   POINTER :: RK_time
  REAL(prec), DIMENSION(:,:), POINTER :: L_step

  REAL(prec), DIMENSION(:,:), POINTER :: Masse, Masse_inv, Rigid, Rigid_inv

  INTEGER :: nb_cell, order_x, order_t, n_time
  REAL(prec) :: time, tmax,t_ini,dt
  REAL(prec) :: xL,xR,CFL,max_dflux,dx
  REAL(prec) :: eps0 =0._prec

  CHARACTER(LEN=32) :: DG_meth, quad_meth

  integer, parameter   :: numfile_sol=1, numfile_param=2, numfile_data = 3
  character(len=32)    :: nomfile_sol = 'file_sol.txt', nomfile_param = 'param.txt', nomfile_data= 'file_data.txt'

  INTEGER :: i,j


  CONTAINS

  SUBROUTINE ALLOCATE_all
    IMPLICIT NONE

    ALLOCATE( sol(nb_cell), flux_h(nb_cell), sol_step(order_t+1))
    ALLOCATE( x_cell(nb_cell+1), x_middle(nb_cell), cell_size(nb_cell)) 
    ALLOCATE( x_quad(order_x), w_quad(order_x) ) 
    ALLOCATE( coeff_DG(order_x, order_x) ) 

    ALLOCATE( coeff_Taylor(10,10), coeff_legendre(10,10) )
    ALLOCATE( RK_alpha(order_t,order_t), RK_beta(order_t,order_t), RK_time(order_t+1) )
    ALLOCATE( L_step(order_t, order_x))

    ALLOCATE( Masse(order_x,order_x), Masse_inv(order_x,order_x) ) 
    ALLOCATE( Rigid(order_x,order_x), Rigid_inv(order_x,order_x) ) 
    
    DO i =1,nb_cell
      ALLOCATE(sol(i)%base_poly(order_x))
      ALLOCATE(flux_h(i)%base_poly(order_x))
    END DO
    DO i =1,order_t+1
      ALLOCATE(sol_step(i)%base_poly(order_x))
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

END MODULE mod_declaration 