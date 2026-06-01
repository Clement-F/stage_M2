PROGRAM tests

   use mod_RKDG
   use mod_init
   implicit none

    REAL(prec) :: xi, t1, t2, err
    INTEGER    :: ii 
    
    print *, "init"
    xL = 0._prec ; xR = pi
    order_x = 3; DG_meth = "Lobatto"; quad_meth = "Lobatto"
    nb_cell =10


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
    
    open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')
    write(unit= numfile_data, fmt='("ordre x = "i5)') order_x
    write(unit= numfile_data, fmt='("ordre t = "i5)') order_t
    write(unit= numfile_data, fmt='("nx = "i5)') nb_cell
    write(unit= numfile_data, fmt='("nt = "i5)') 1
    close(unit= numfile_data)

    open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')
    DO i=1,nb_cell
        DO j=1,order_x
            xi = Ref_to_loc(i,x_quad(j))
            t1 = eval_sol(xi,i)
            t2 = sinus(xi - time*vit_adv,0)
            write(unit=numfile_sol,  fmt='(f10.6, f16.6, f16.6)') xi,t1, t2
        END DO
    END DO
    close(unit=numfile_sol)

    CONTAINS

    FUNCTION f(x,ni)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: x
        INTEGER,    INTENT(IN) :: ni
        REAL(prec) :: f

        f = sin(2*pi*x)
    END FUNCTION f

END PROGRAM tests