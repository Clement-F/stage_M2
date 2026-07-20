PROGRAM MAIN

   use mod_RKDG
   use mod_init
   implicit none

    INTEGER :: i



    CALL SYSTEM_CLOCK(COUNT_RATE=nb_prd_sec,COUNT_MAX=nb_prd_max)
    CALL SYSTEM_CLOCK(COUNT=nb_prd_ini)
    
    CALL INIT_ALL

    open(unit=numfile_data,     file=nomfile_data,      form ='formatted', status ='unknown')
    open(unit=numfile_sol,      file=nomfile_sol,       form ='formatted', status ='unknown')
    open(unit=numfile_solex,    file=nomfile_solex,     form ='formatted', status ='unknown')
    open(unit=numfile_meshout,  file=nomfile_meshout,   form ='formatted', status ='unknown')
   

    write(unit= numfile_data, fmt='("nb var  = ",i5)') nb_var
    write(unit= numfile_data, fmt='("ordre x = ",i5)') order_x
    write(unit= numfile_data, fmt='("ordre t = ",i5)') order_t

    IF(subcell_use) THEN
        write(unit= numfile_data, fmt='("nb_cell = ",i5)') nb_cell
        write(unit= numfile_data, fmt='("nb_subcell = ",i5)') nb_subcell
    ELSE
        write(unit= numfile_data, fmt='("nb_cell = ",i5)') nb_cell
        write(unit= numfile_data, fmt='("blablab = 00000000 ")') 
    END IF

    close(unit=numfile_data)


    !! Print the time lapse of the computation !!
    CALL eval_time(nb_prd_ini,nb_prd_max,nb_prd_sec,0)
    WRITE(*,*) " "
    
    CALL writout

    DO WHILE (tmax-time  .GT. eps0 )
        CALL dt_calc

        IF(subcell_use) THEN;   CALL Time_step_subcell
        ELSE;                   CALL Time_step
        END IF

        time = time +dt
        n_time = n_time +1
        
        CALL writout 
    END DO
        
    CALL writout 

    open(unit=numfile_data,     file=nomfile_data,      form ='formatted', status ='old',  position='append')
    write(unit= numfile_data, fmt='("nt = ",i5)') n_imp
    
    DO i=1,n_imp
        write(unit= numfile_data, fmt='("time ",i5," = ",f16.6)') i, Time_stemp(i)
    END DO

    close(unit=numfile_sol);    close(unit=numfile_solex)
    close(unit=numfile_data)
    close(unit=numfile_meshout)

    IF(convergence) THEN
        open(unit=numfile_conv,  file=nomfile_conv, form ='formatted', status ='unknown', position='append')
        write(unit=numfile_conv, fmt='("=====================")') 
        write(unit=numfile_conv, fmt='("for elements P",i1," and RK SSP of order ",i1)' ) size_base-1, order_t
        write(unit=numfile_conv, fmt='("for nx = ",i5," we have error :")' ) nb_cell
        write(unit=numfile_conv, fmt='("err_L1 :", e20.12 )') err_L1
        write(unit=numfile_conv, fmt='("err_L2 :", e20.12 )') err_L2
        write(unit=numfile_conv, fmt='("err_Li :", e20.12 )') err_Li
        write(unit=numfile_conv, fmt='("=====================")') 
        close(unit=numfile_conv)
    END IF

    print *,"closed"
    print *, counter1, counter2
    

    !! Print the time lapse of the computation !!
    CALL eval_time(nb_prd_ini,nb_prd_max,nb_prd_sec,0)
    WRITE(*,*) " "

    CALL DEALLOCATE_all


    print *, "end"

END PROGRAM