PROGRAM MAIN

   use mod_RKDG
   implicit none

    REAL(prec) :: t,xi

    xL = 0._prec ; xR = 2* 3.1415_prec
    order_x =2; order_t = 1
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

    time = 0._prec; tmax= 0.5_prec; CFL=1._prec
    n_time =1
    print *,"init polynome and quadrature"

    DG_meth = "Legendre"
    quad_meth = "Legendre"

    print *,"quad"
    CALL Coeff_quad_init
    print *,"poly"
    CALL Coeff_DG_init
    print *,"RK"
    CALL Coeff_RK_init(order_t)

    print *,"init mass"

    CALL Matrice_Masse_init
    CALL Matrice_Rigid_init

    print *,"projection"

    DO i=1,nb_cell
        CALL Projection_Pk(test,sol(i)%base_poly,i)
        ! print *, sol(i)%base_poly
    END DO

    print *,"evaluation"


   open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')
   open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')
   
   write(unit= numfile_data, fmt='("nx = "i5)') nb_cell


    DO i=1,nb_cell
        DO j=1,2
            xi = x_cell(i) + REAL(j, prec)/2._prec*cell_size(i)
            ! write (*,"(i3,i3)",advance ='no') i,j
            t = eval_sol(i,xi)
            write(unit=numfile_sol,  fmt='(f10.6, f12.6)') xi,t
        END DO
    END DO

    max_dflux = 2._prec

    print *,"time_step"
    write(unit=numfile_sol, fmt='("------------------------")' ) 

    CALL Time_step

    DO i=1,nb_cell
        DO j=1,2
            xi = x_cell(i) + REAL(j, prec)/2._prec*cell_size(i)
            ! write (*,"(i3,i3)",advance ='no') i,j
            t = eval_sol(i,xi)
            write(unit=numfile_sol,  fmt='(f10.6, f12.6)') xi,t
        END DO
    END DO
    print *,"time_step"
    write(unit=numfile_sol, fmt='("------------------------")' ) 


    ! DO WHILE (time .LT. tmax)
    !     CALL Time_step
        
    !     DO i=1,nb_cell
    !         DO j=1,2
    !             xi = x_cell(i) + REAL(j, prec)/2._prec*cell_size(i)
    !             ! write (*,"(i3,i3)",advance ='no') i,j
    !             t = eval_sol(i,xi)
    !             write(unit=numfile_sol,  fmt='(f10.6, f12.6)') xi,t
    !         END DO
    !     END DO
    !     write(unit=numfile_sol, fmt='("------------------------")' ) 
    ! END DO

   write(unit= numfile_data, fmt='("nt = "i5)') n_time





    ! print *,"_"
   close(unit=numfile_sol)
   close(unit=numfile_data)
    print *,"closed"

    CALL DEALLOCATE_all


    print *, "end"

END PROGRAM