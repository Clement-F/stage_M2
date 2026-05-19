
MODULE mod_declaration
  USE precis
  IMPLICIT NONE

  TYPE var_type
     REAL(prec),DIMENSION(:),POINTER :: base_poly
     REAL(prec),DIMENSION(2) :: interface
  END TYPE var_type

  TYPE(var_type), DIMENSION(:), POINTER :: sol

  INTEGER :: nb_cell, order_x, order_t

  REAL(prec), DIMENSION(:), POINTER :: x_cell, x_middle, cell_size
  REAL(prec), DIMENSION(:), POINTER :: x_quad, w_quad
  REAL(prec), DIMENSION(:,:), POINTER :: coeff_DG

  REAL(prec), DIMENSION(:,:), POINTER :: coeff_Taylor, coeff_legendre
  REAL(prec), DIMENSION(:,:), POINTER :: RK_alpha, RK_beta

  REAL(prec), DIMENSION(:,:), POINTER :: Masse, Masse_inv

  REAL(prec) :: time, tmax,t_ini,dt
  REAL(prec) :: xL,xR,CFL,max_flux,dx
  REAL(prec) :: eps0

  CHARACTER(LEN=32) :: DG_meth, quad_meth

END MODULE mod_declaration 