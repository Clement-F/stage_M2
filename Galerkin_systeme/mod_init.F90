MODULE mod_init
    USE mod_polynome
    IMPLICIT NONE
    
CONTAINS

    SUBROUTINE INIT_ALL
        IMPLICIT NONE

        INTEGER :: ni,ii,jj, kk
        INTEGER :: i,j
        REAL(prec) :: YY

        print *,"read"

        CALL READ_DATA
        ! print *, (nomfile_sol)
        print *,"red"
    
        ! IF( TRIM(nomfile_sol) == " ") nomfile_sol = "file_sol"

        nomfile_sol = TRIM(nomfile_sol)//".txt"

        print *, "init"

        IF((DG_meth)=="Lobatto")  size_base = order_x 
        IF((DG_meth)=="Legendre") size_base = order_x 

        if((quad_meth)=="Lobatto") nb_nodes = size_base+1        !CEILING((2*(size_base-1)+3 )/2.) 
        if((quad_meth)=="Legendre")nb_nodes = size_base          !CEILING((2*(size_base-1)+1 )/2.) 


        IF((flux_name) == "advection") THEN 
            max_dflux = abs(vit_adv)
            nb_var=2
        ELSE IF((flux_name) == "burgers") THEN
            nb_var=2
        ELSE IF((flux_name) == "Shallow") THEN
            nb_var=2
        ELSE IF((flux_name) == "Euler") THEN
            nb_var=3;
            gamma_iso = 1.4_prec
            IF((sol_ini_name) == "Sod") gamma_iso = 1.4_prec
            IF((sol_ini_name) == "isentropique") gamma_iso = 3._prec
        ELSE
            print *," flux non reconnue "
            STOP
        END IF

        print *, "allocation"

        CALL ALLOCATE_all
        dx = REAL((xR-xL)/Real(nb_cell,prec), prec)
        sub_dx = 2._prec/Real(nb_subcell,prec)

        DO i =1,nb_cell+1 
            x_cell(i) = xL + Real(i-1,prec)*dx
        END DO       

        DO i=1,nb_cell
            cell_size(i)= x_cell(i+1)-x_cell(i)  
            x_middle(i) = x_cell(i) + cell_size(i)/2._prec
        END DO
        
        tmax = REAL(tmax,prec)
        time = 0._prec
        n_time =1; 
        n_imp = 0
        t_imp = tmax/Real(print_rule,prec)

        dt = 10._prec**(-8)



        CALL Coeff_quad_init
        CALL Coeff_DG_init
        
        DO ii=1,size_base
          sig_1(ii) = DG_base(-1._prec,ii, LOC=LRef,ni=0); 
          sig_2(ii) = DG_base( 1._prec,ii, LOC=LRef,ni=0); 
          DO jj = 1,nb_nodes
            sig_quad(ii,jj) = DG_base(x_quad(jj),ii, LOC=LRef,ni=0)
          END DO
        END DO

        print *, "matrice points quad "
        CALL print_mat(sig_quad, size_base, nb_nodes)
        
        CALL Coeff_RK_init(order_t)
        CALL Matrice_Masse_init

        print *, "masse matrix "
        CALL print_mat(Masse, size_base, size_base)

        CALL Matrice_Rigid_init
        
        print *,"matrices"

        IF(subcell_use) CALL sub_cells_init


        DO ni=1,nb_cell
            
            IF(subcell_use) THEN
                DO ii=1,nb_var

                DO j =1,nb_subcell; DO kk =1,nb_nodes
                    YY = Ref_to_loc(ni=ni, XX=Refsub_to_Ref(ZZ=x_quad(kk),n_sub =j))
                    sol(ni)%val_subcells(j,ii) = sol(ni)%val_subcells(j,ii) + Q_init(YY,ni=ni,nvar=ii)*w_quad(kk)/2._prec
                END DO; END DO
            
                DO j=1,size_base
                    sol(ni)%base_poly(j,ii) = DOT_PRODUCT(Projection_VF_inv(j,:), sol(ni)%val_subcells(:,ii))
                END DO

                END DO

            ELSE 
                DO ii=1,nb_var
                    CALL Projection_Pk(Q_init,sol(ni)%base_poly(:,ii), LOC=LLoc, ni= ni, nvar = ii)
                END DO
            END IF


            DO jj=1,nb_var; DO ii=1,nb_nodes
                
                sol(ni)%val_nodes(ii,jj)  = eval_sol(x_quad(ii),ni= ni,nvar=jj, kk=ii, LOC= LRef)
                        
                IF((quad_meth)=="Lobatto") THEN
                sol(ni)%inter(1,jj)      = sol(ni)%val_nodes(1,jj)
                sol(ni)%inter(2,jj)      = sol(ni)%val_nodes(nb_nodes,jj)
                ELSE 
                sol(ni)%inter(1,jj)      = eval_sol(x_cell(ni),ni= ni,nvar=jj  ,LOC=LLoc)
                sol(ni)%inter(2,jj)      = eval_sol(x_cell(ni+1),ni= ni,nvar=jj,LOC=LLoc)
                END IF
            
            END DO; END DO

        END DO
                

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

        CALL Skip_lines(numfile_param,3) 

        read(numfile_param,  *) subcell_use
        read(numfile_param,  *) monolithique
        read(numfile_param,  *) max_rule
        read(numfile_param,  *) coeff_smooth
        read(numfile_param,  *) smooth_extrema
        read(numfile_param,  *) nb_subcell
        read(numfile_param,  *) subcell_repartition

        CALL Skip_lines(numfile_param,3) 
        
        read(numfile_param,  *) cfl;   
        read(numfile_param,  *) error_calc; 
        read(numfile_param,  *) max_check; 
        read(numfile_param,  *) print_rule;     
        read(numfile_param,  *) nomfile_sol;     

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