MODULE mod_init
    USE mod_polynome
    IMPLICIT NONE
    
CONTAINS

    SUBROUTINE INIT_ALL
        IMPLICIT NONE

        INTEGER :: ni,ii,jj

        CALL READ_DATA
        
        print *, "init"

        IF(TRIM(DG_meth)=="Lobatto")  size_base = order_x 
        IF(TRIM(DG_meth)=="Legendre") size_base = order_x 

        if(TRIM(quad_meth)=="Lobatto") size_quad_nodes = size_base+1        !CEILING((2*(size_base-1)+3 )/2.) 
        if(TRIM(quad_meth)=="Legendre")size_quad_nodes = size_base          !CEILING((2*(size_base-1)+1 )/2.) 

        CALL ALLOCATE_all
        dx = REAL((xR-xL)/nb_cell, prec)
        sub_dx = 2._prec/nb_subcell

        DO i =1,nb_cell+1 
            x_cell(i) = xL + (i-1)*dx
        END DO       

        DO i=1,nb_cell
            cell_size(i)= x_cell(i+1)-x_cell(i)  
            x_middle(i) = x_cell(i) + cell_size(i)/2._prec
        END DO
        
        
        time = 0._prec
        n_time =1; 
        n_imp = 0
        t_imp = tmax/Real(print_rule,prec)


        CALL Coeff_quad_init
        CALL Coeff_DG_init
        
        DO ii=1,size_base
          sig_1(ii) = DG_base(-1._prec,ii, LOC=LRef,ni=0); 
          sig_2(ii) = DG_base( 1._prec,ii, LOC=LRef,ni=0); 
          DO jj = 1,size_quad_nodes
            sig_quad(ii,jj) = DG_base(x_quad(jj),ii, LOC=LRef,ni=0)
          END DO
        END DO
        


        CALL Coeff_RK_init(order_t)
        CALL Matrice_Masse_init
        CALL Matrice_Rigid_init
        
        print *,"matrices"


        DO ni=1,nb_cell
            CALL Projection_Pk(Q_init,sol(ni)%base_poly, LOC=LLoc, ni= ni)

            DO ii=1,size_quad_nodes
                sol(ni)%val_nodes(ii)  = eval_sol(x_quad(ii),ni, ii, LOC= LRef)
            END DO
                        
            IF(TRIM(quad_meth)=="Lobatto") THEN
            sol(ni)%inter(1)      = sol(ni)%val_nodes(1)
            sol(ni)%inter(2)      = sol(ni)%val_nodes(size_quad_nodes)
            ELSE 
            sol(ni)%inter(1)      = eval_sol(x_cell(ni),ni,  LOC=LLoc)
            sol(ni)%inter(2)      = eval_sol(x_cell(ni+1),ni,LOC=LLoc)
            END IF

        END DO
        
        ! cALL sub_cells_init
        
        IF(TRIM(flux_name) == "advection") THEN 
            max_dflux = abs(vit_adv)
        END IF

        err_L1 = 0._prec; err_L2 =0._prec;  err_Li=0._prec
        print *,"end init"

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

        read(numfile_param,  *) subcell_use
        read(numfile_param,  *) nb_subcell
        ! read(numfile_param,  *) subcell_repartition
        CALL Skip_lines(numfile_param,1) 

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
        read(numfile_param,  *) bdry_cond;

        
        close(unit=numfile_data)
        
    END SUBROUTINE READ_DATA



END MODULE mod_init