PROGRAM tests

   use mod_RKDG
   use mod_init
   implicit none

    REAL(prec) :: pi= acos(-1._prec)
    REAL(prec) :: t1, t2, err
    INTEGER    :: ii 
    
    print *, "init"
    xL = 0._prec ; xR = pi
    order_x = 3; DG_meth = "Legendre"; quad_meth = "Legendre"
    nb_cell =300


    CALL ALLOCATE_all
    dx = REAL((xR-xL)/nb_cell, prec)

    DO i =1,nb_cell+1 
        x_cell(i) = xL + (i-1)*dx
    END DO
    
    DO i=1,nb_cell
        cell_size(i)= x_cell(i+1)-x_cell(i)  
        x_middle(i)= x_cell(i) + cell_size(i)/2._prec
    END DO
    
    CALL Coeff_quad_init
    CALL Coeff_DG_init
    CALL Matrice_Masse_init
    CALL Matrice_Rigid_init

    ! print *, coeff_DG

    DO ii=1,nb_cell
        CALL Projection_Pk(f,sol(ii)%base_poly,ii)
    END DO

    DO ii=1,nb_cell

        t1 = quadrature(ii,eval_sol,ii,unit,0)
        t2 = quadrature(ii,f,0,unit,0)
        err = err + abs(t1-t2)
    
    END DO

    print *, err

    CONTAINS

    FUNCTION f(x,ni)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: x
        INTEGER,    INTENT(IN) :: ni
        REAL(prec) :: f

        f = sin(x)
    END FUNCTION f

END PROGRAM tests