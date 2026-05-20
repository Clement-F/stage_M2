PROGRAM MAIN

   use mod_polynome
   implicit none

    REAL(prec) :: t,xi

    xL = 0._prec ; xR = 5._prec
    order_x =2
    nb_cell =10

    print *, "init"

    CALL ALLOCATE_all

    dx = REAL((xR-xL)/nb_cell, prec)

    DO i =1,nb_cell+1 
        x_cell(i) = xL + (i-1)*dx
    END DO
    
    DO i=1,nb_cell
        cell_size(i)= x_cell(i+1)-x_cell(i)  
        x_middle(i)= x_cell(i) + cell_size(i)/2._prec
    END DO

    print *, x_cell(nb_cell+1)

    print *,"init polynome and quadrature"

    DG_meth = "Legendre"
    quad_meth = "Legendre"

    print *,"quad"
    CALL Coeff_quad_init
    print *,"poly"
    CALL Coeff_DG_init

    print *,"init mass"

    CALL Matrice_Masse_init

    print *,"projection"

    DO i=1,nb_cell
        CALL Projection_Pk(test,sol(i)%base_poly,i)
        print *, sol(i)%base_poly
    END DO

    print *,"evaluation"


   open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')
   
   write(unit= numfile_data, fmt='("nx = "i5)') nb_cell

   open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')

    DO i=1,nb_cell
        DO j=1,10
            xi = x_cell(i) + REAL(j, prec)/10._prec*cell_size(i)
            ! print *, i,j,xi
            t = eval_sol(i,xi)
            ! print *,t
            write(unit=numfile_sol,  fmt='(f10.6, f12.6)') xi,t
        END DO
    END DO

   close(unit=numfile_sol)
   close(unit=numfile_data)

    print *, "end"

END PROGRAM