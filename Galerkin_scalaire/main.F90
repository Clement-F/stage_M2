PROGRAM MAIN

   use mod_RKDG
   use mod_init
   implicit none

    REAL(prec) :: t,xi

    CALL INIT_ALL
    open(unit=numfile_data, file=nomfile_data, form ='formatted', status ='old')
    open(unit=numfile_sol,  file=nomfile_sol, form ='formatted', status ='old')
   
    write(unit= numfile_data, fmt='("nx = "i5)') nb_cell
        
    
            print *,n_time, time, dt
            DO i=1,nb_cell
                DO j=1,5
                    xi = x_cell(i) + REAL(j, prec)/5._prec*cell_size(i)
                    t = eval_sol(i,xi)
                    write(unit=numfile_sol,  fmt='(f10.6, f16.6)') xi,t
                END DO
            END DO
            n_imp = n_imp +1
            Time_stemp(n_imp) = time
            ! print *, sol(1)%base_poly
            ! print *, sol(nb_cell)%base_poly
            
            write(unit=numfile_sol, fmt='("------------------------")' ) 
            
    DO WHILE (time .LT. tmax)     
        ! print *,time   
        CALL dt_calc
        CALL Time_step
        IF(time >=  n_imp*t_imp)  THEN
        
            print *,"--------------",n_time, time, dt,"--------------"
            DO i=1,nb_cell
                DO j=1,5
                    xi = x_cell(i) + REAL(j, prec)/5._prec*cell_size(i)
                    t = eval_sol(i,xi)
                    write(unit=numfile_sol,  fmt='(f10.6, f16.6)') xi,t
                END DO
            END DO
            ! print *, sol(1)%base_poly
            ! print *, sol(nb_cell)%base_poly
            n_imp = n_imp +1
            Time_stemp(n_imp) = time
            
            write(unit=numfile_sol, fmt='("------------------------")' ) 
        END IF
        

    END DO

    write(unit= numfile_data, fmt='("nt = "i5)') n_imp
    DO i=1,n_imp
        write(unit= numfile_data, fmt='("time ",i5," = ",f16.6)') i, Time_stemp(i)
    END DO






    ! print *,"_"
   close(unit=numfile_sol)
   close(unit=numfile_data)
    print *,"closed"

    CALL DEALLOCATE_all


    print *, "end"

END PROGRAM