MODULE mod_init
    USE mod_polynome
    IMPLICIT NONE
    
CONTAINS

    SUBROUTINE INIT_ALL
        IMPLICIT NONE

        CALL READ_DATA
        
        print *, "init"

        IF(TRIM(DG_meth)=="Lobatto") size_base = order_x 
        IF(TRIM(DG_meth)=="Legendre") size_base = order_x 

        if(TRIM(quad_meth)=="Lobatto") size_quad_nodes = CEILING((size_base+2 )/2.)
        if(TRIM(quad_meth)=="Legendre") size_quad_nodes = size_base

        CALL ALLOCATE_all
        dx = REAL((xR-xL)/nb_cell, prec)

        DO i =1,nb_cell+1 
            x_cell(i) = xL + (i-1)*dx
        END DO
        
        DO i=1,nb_cell
            cell_size(i)= x_cell(i+1)-x_cell(i)  
            x_middle(i)= x_cell(i) + cell_size(i)/2._prec
        END DO
        
        time = 0._prec
        n_time =1; 
        n_imp = 0
        t_imp = tmax/Real(print_rule,prec)

        CALL Coeff_quad_init
        CALL Coeff_DG_init
        CALL Coeff_RK_init(order_t)
        CALL Matrice_Masse_init
        CALL Matrice_Rigid_init
        
        DO i=1,nb_cell
            CALL Projection_Pk(Q_init,sol(i)%base_poly,i)
        END DO

        IF(TRIM(flux_name) == "advection") THEN 
            max_dflux = abs(vit_adv)
        END IF

        err_L1 = 0._prec; err_L2 =0._prec;  err_Li=0._prec

    END SUBROUTINE INIT_ALL

    SUBROUTINE Skip_lines(unit,nlines)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: unit, nlines
        INTEGER :: nl
        CHARACTER(len=1) :: c1

        DO nl=1,nlines
            read(unit,*)    c1
        END DO
    END SUBROUTINE  Skip_lines
    
    SUBROUTINE READ_DATA
        IMPLICIT NONE
        
        open(unit=numfile_param, file=nomfile_param, form ='formatted', status ='old')
        
        CALL Skip_lines(numfile_param,3)

        read(numfile_param,  *) nb_cell;  
        read(numfile_param,  *) xL;  
        read(numfile_param,  *) xR;      
        read(numfile_param,  *) tmax;    
        CALL Skip_lines(numfile_param,1) 
        read(numfile_param,  *) cfl;   
        read(numfile_param,  *) frame_rule; 
        read(numfile_param,  *) print_rule;     

        CALL Skip_lines(numfile_param,3)
        
        read(numfile_param,  *) order_x;   
        read(numfile_param,  *) order_t; 
        read(numfile_param,  *) DG_meth; 
        read(numfile_param,  *) quad_meth;

        CALL Skip_lines(numfile_param,3)
        
        read(numfile_param,  *) sol_ini_name; 
        read(numfile_param,  *) flux_name;
        read(numfile_param,  *) vit_adv;

        
        close(unit=numfile_data)
        
    END SUBROUTINE READ_DATA



END MODULE mod_init