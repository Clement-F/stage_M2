
MODULE mod_declaration
  USE precis
  IMPLICIT NONE

  TYPE var_type
     REAL(prec),DIMENSION(:),POINTER :: base_poly
     REAL(prec),DIMENSION(2) :: interface
  END TYPE var_type

  TYPE(var_type), DIMENSION(:), POINTER :: sol

  REAL(prec), DIMENSION(:), POINTER :: x_cell, x_middle, cell_size
  REAL(prec), DIMENSION(:), POINTER :: x_quad, w_quad
  REAL(prec), DIMENSION(:,:), POINTER :: coeff_DG
  REAL(prec), DIMENSION(:), POINTER :: fct_h

  REAL(prec), DIMENSION(:,:), POINTER :: coeff_Taylor, coeff_legendre
  REAL(prec), DIMENSION(:,:), POINTER :: RK_alpha, RK_beta

  REAL(prec), DIMENSION(:,:), POINTER :: Masse, Masse_inv

  INTEGER :: nb_cell, order_x, order_t
  REAL(prec) :: time, tmax,t_ini,dt
  REAL(prec) :: xL,xR,CFL,max_flux,dx
  REAL(prec) :: eps0

  CHARACTER(LEN=32) :: DG_meth, quad_meth

  integer, parameter   :: numfile_sol=1, numfile_param=2, numfile_data = 3
  character(len=32)    :: nomfile_sol = 'file_sol.txt', nomfile_param = 'param.txt', nomfile_data= 'file_data.txt'

  INTEGER :: i,j


  CONTAINS

  SUBROUTINE ALLOCATE_all
    IMPLICIT NONE
    ALLOCATE( sol(nb_cell))
    ALLOCATE( x_cell(nb_cell+1), x_middle(nb_cell), cell_size(nb_cell)) 
    ALLOCATE( x_quad(order_x), w_quad(order_x) ) 
    ALLOCATE( coeff_DG(order_x, order_x) ) 

    ALLOCATE( coeff_Taylor(order_x,order_x), coeff_legendre(order_x,order_x) )
    ! ALLOCATE( RK_alpha, RK_beta )

    ALLOCATE( Masse(order_x,order_x), Masse_inv(order_x,order_x) ) 
    
    DO i =1,nb_cell
        ALLOCATE(sol(i)%base_poly(order_x))
    END DO
    

  END SUBROUTINE ALLOCATE_all

END MODULE mod_declaration 