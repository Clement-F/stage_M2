
MODULE mod_declaration
  USE precis
  IMPLICIT NONE

  TYPE var_type
    REAL(prec),DIMENSION(:,:),POINTER :: base_poly   ! decomposition dans la base DG
    REAL(prec),DIMENSION(:,:),POINTER :: deriv       ! decomposition dans la base DG
    REAL(prec),DIMENSION(:,:),POINTER :: val_quad   ! valeurs aux points de quadratures
    REAL(prec),DIMENSION(:,:),POINTER :: val_subcells   ! valeurs aux points de quadratures
    REAL(prec),DIMENSION(:,:),POINTER :: inter
  END TYPE var_type


  
  TYPE flux_
    REAL(prec) :: theta
    REAL(prec), DIMENSION(:,:), POINTER :: flux_VF
    REAL(prec), DIMENSION(:,:), POINTER :: flux_DG
    REAL(prec), DIMENSION(:,:), POINTER :: flux_subcells
    REAL(prec), DIMENSION(:,:), POINTER :: flux_tilde
  END TYPE flux_

  Type subcells
    INTEGER :: index_m, index_s
    INTEGER, DIMENSION(2) :: L, LL, R, RR
    LOGICAL :: extrema
  END TYPE subcells

  TYPE(var_type), DIMENSION(:), POINTER :: sol, sol_step, sol_exa
  TYPE(flux_), DIMENSION(:), POINTER :: flux_h
  TYPE(subcells), DIMENSION(:,:), POINTER :: subcells_

  REAL(prec), DIMENSION(:,:),   POINTER :: g

  REAL(prec), DIMENSION(:),   POINTER :: x_cell, x_middle, cell_size
  REAL(prec), DIMENSION(:),   POINTER :: x_subcell, x_submiddle, subcell_size
  REAL(prec), DIMENSION(:),   POINTER :: x_quad, w_quad
  REAL(prec), DIMENSION(:),   POINTER :: pts_DG, sig_1, sig_2
  REAL(prec), DIMENSION(:,:), POINTER :: coeff_DG, sig_quad

  REAL(prec), DIMENSION(:,:), POINTER :: coeff_Taylor, coeff_legendre

  REAL(prec), DIMENSION(:,:),   POINTER :: theta_
  REAL(prec), DIMENSION(:,:), POINTER :: RK_alpha
  REAL(prec), DIMENSION(:),   POINTER :: RK_time,RK_beta
  REAL(prec), DIMENSION(:),   POINTER :: L_step
  REAL(prec), DIMENSION(:),   POINTER :: Time_stemp

  REAL(prec), DIMENSION(:),   POINTER :: C_m, C_p
  REAL(prec), DIMENSION(:,:), POINTER :: Projection_VF, Projection_VF_inv
  REAL(prec), DIMENSION(:,:), POINTER :: Projection_VF_d,Projection_VF_dd, Projection_VF_inv_d, Projection_VF_inv_dd
  REAL(prec), DIMENSION(:,:), POINTER :: Masse, Masse_inv, Rigid, Rigid_inv

  REAL(prec) :: min_glob, max_glob

  LOGICAL*1 :: subcell_use, convex_flux, monolithique, error_calc, convergence
  INTEGER   :: max_rule, coeff_smooth, max_check, smooth_extrema, flux_num


  INTEGER :: nb_var
  INTEGER :: nb_cell, nb_subcell, order_x, order_t 
  INTEGER :: size_base, nb_nodes
  INTEGER :: n_imp, n_time, frame_rule, print_rule
  INTEGER :: counter1, counter2
  REAL(prec) :: time, tmax,t_ini,dt, dt_old, t_imp
  REAL(prec) :: xL,xR,CFL,max_dflux,dx, sub_dx
  REAL(prec), parameter :: pi = acos(-1._prec)
  REAL(prec), parameter :: eps0=0.1_prec**(2*prec-3)
  REAL(prec) :: vit_adv, gamma_iso 
  REAL(prec) :: err_L1, err_L2, err_Li

  CHARACTER(Len=8) :: LRef = "Ref", LLoc = "Loc", LSub = "SubRef" 

  CHARACTER(LEN=32):: DG_meth, quad_meth, sol_ini_name, flux_name, bdry_cond, subcell_repartition

  INTEGER   :: numfile_sol=1, numfile_param=2, numfile_data = 3, numfile_conv = 4, numfile_meshout = 5, numfile_solex = 7
  CHARACTER(len=32)    :: nomfile_sol = 'file_sol.txt', nomfile_param = 'param.txt', nomfile_data= 'file_data.txt', & 
                        & nomfile_conv= 'convergence_err.txt', nomfile_meshout = 'mesh_out.txt', nomfile_solex = 'file_solex.txt'



  CONTAINS

  SUBROUTINE ALLOCATE_all
    IMPLICIT NONE

    INTEGER :: i

    ALLOCATE( sol(nb_cell), sol_step(nb_cell), flux_h(nb_cell+1), sol_exa(nb_cell))
    ALLOCATE( g(nb_cell +1, nb_var))

    ALLOCATE (subcells_(nb_cell, nb_subcell))

    ALLOCATE( x_cell(nb_cell+1),        x_middle(nb_cell),        cell_size(nb_cell)) 
    ALLOCATE( x_subcell(nb_subcell+1),  x_submiddle(nb_subcell),  subcell_size(nb_subcell))
    ALLOCATE( theta_(nb_cell,nb_subcell+1))

    ALLOCATE( x_quad(nb_nodes), w_quad(nb_nodes) ) 
    ALLOCATE( sig_1(size_base), sig_2(size_base), sig_quad(size_base,nb_nodes))
    ALLOCATE( coeff_DG(size_base, size_base) ) 

    ALLOCATE( coeff_Taylor(10,10), coeff_legendre(10,10) )
    ALLOCATE( RK_alpha(order_t,2), RK_beta(order_t), RK_time(order_t) )
    ALLOCATE( Time_stemp(print_rule+2))
    
    ALLOCATE( pts_DG(size_base)) 
    ALLOCATE( C_m(nb_subcell+1), C_p(nb_subcell+1))
    ALLOCATE( Projection_VF(nb_subcell,size_base),    Projection_VF_inv(size_base,nb_subcell))
    ALLOCATE( Projection_VF_d(nb_subcell,size_base))
    ALLOCATE( Projection_VF_dd(nb_subcell,size_base))
    ALLOCATE( Masse(size_base,size_base), Masse_inv(size_base,size_base) ) 
    ALLOCATE( Rigid(size_base,size_base), Rigid_inv(size_base,size_base) ) 

    DO i =1,nb_cell
      ALLOCATE(sol(i)%base_poly(size_base,nb_var));       ALLOCATE(sol(i)%val_quad(nb_nodes,nb_var));      ALLOCATE(sol(i)%val_subcells(nb_subcell,nb_var));      ALLOCATE(sol(i)%inter(2,nb_var));
      ALLOCATE(sol_exa(i)%base_poly(size_base,nb_var));   ALLOCATE(sol_exa(i)%val_quad(nb_nodes,nb_var));  ALLOCATE(sol_exa(i)%val_subcells(nb_subcell,nb_var));  ALLOCATE(sol_exa(i)%inter(2,nb_var)); 
      ALLOCATE(sol_step(i)%base_poly(size_base,nb_var));  ALLOCATE(sol_step(i)%val_quad(nb_nodes,nb_var)); ALLOCATE(sol_step(i)%val_subcells(nb_subcell,nb_var)); ALLOCATE(sol_step(i)%inter(2,nb_var))
                                     

      ALLOCATE(flux_h(i)%flux_DG(size_base,nb_var))      
      ALLOCATE(flux_h(i)%flux_tilde(nb_subcell+1,nb_var));ALLOCATE(flux_h(i)%flux_VF(nb_subcell+1,nb_var)); ALLOCATE(flux_h(i)%flux_subcells(nb_subcell+1,nb_var));  
    
      ALLOCATE(sol_step(i)%deriv(size_base, nb_var))
    END DO
    

  END SUBROUTINE ALLOCATE_all

  SUBROUTINE DEALLOCATE_all
    IMPLICIT NONE

    INTEGER :: i
    
    DO i =1,nb_cell
      DEALLOCATE(sol(i)%base_poly);       DEALLOCATE(sol(i)%val_quad);       DEALLOCATE(sol(i)%val_subcells)
      DEALLOCATE(sol_exa(i)%base_poly);   DEALLOCATE(sol_exa(i)%val_quad);   DEALLOCATE(sol_exa(i)%val_subcells) 
      DEALLOCATE(sol_step(i)%base_poly);  DEALLOCATE(sol_step(i)%val_quad);  DEALLOCATE(sol_step(i)%val_subcells); DEALLOCATE(sol_step(i)%deriv)
      DEALLOCATE(flux_h(i)%flux_VF);      DEALLOCATE(flux_h(i)%flux_DG);      DEALLOCATE(flux_h(i)%flux_tilde); DEALLOCATE(flux_h(i)%flux_subcells)

      DEALLOCATE(sol(i)%inter);           DEALLOCATE(sol_exa(i)%inter);       DEALLOCATE(sol_step(i)%inter)
    END DO
    
    DEALLOCATE( sol, sol_step, flux_h, sol_exa, g)

    DEALLOCATE (subcells_)
    DEALLOCATE (theta_)

    DEALLOCATE( x_cell,     x_middle,     cell_size) 
    DEALLOCATE( x_subcell,  x_submiddle,  subcell_size)

    DEALLOCATE( x_quad, w_quad ) 
    DEALLOCATE( sig_1, sig_2, sig_quad)
    DEALLOCATE( coeff_DG ) 

    DEALLOCATE( coeff_Taylor, coeff_legendre )
    DEALLOCATE( RK_alpha, RK_beta, RK_time )
    DEALLOCATE( Time_stemp)
    
    DEALLOCATE( pts_DG) 
    DEALLOCATE( C_m, C_p)
    DEALLOCATE( Projection_VF,  Projection_VF_inv)
    DEALLOCATE( Masse, Masse_inv ) 
    DEALLOCATE( Rigid, Rigid_inv ) 
    
  END SUBROUTINE DEALLOCATE_all
  
  SUBROUTINE Coeff_RK_init(nrk)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nrk
    INTEGER :: nrk2

    RK_time   = 0._prec
    RK_alpha  = 0._prec
    RK_beta   = 0._prec

    SELECT CASE (nrk)
    CASE (1)
       RK_time(1)=0._prec
       RK_alpha(1,1)=1._prec
       RK_beta(1)=1._prec
    CASE (2)
       RK_time(1)=0._prec
       RK_alpha(1,1)=1._prec
       RK_beta(1)=1._prec

       RK_time(2)=1._prec
       RK_alpha(2,1)=0.5_prec
       RK_alpha(2,2)=0.5_prec
       RK_beta(2)=0.5_prec
    CASE (3)
       RK_time(1)=0._prec
       RK_alpha(1,2)=1._prec
       RK_beta(1)=1._prec

       RK_time(2)=0.5_prec
       RK_alpha(2,1)=0.75_prec
       RK_alpha(2,2)=0.25_prec
       RK_beta(2)=0.25_prec

       RK_time(3)=1._prec
       RK_alpha(3,1)=1._prec/3._prec
       RK_alpha(3,2)=2._prec/3._prec
       RK_beta(3)=2._prec/3._prec
    CASE DEFAULT
       RK_alpha=0._prec; RK_beta =0._prec
       DO nrk2=1,nrk
          RK_alpha(nrk2,1)=1._prec
          RK_beta(nrk2)=1._prec/REAL(nrk+1-nrk2,prec)
       END DO
    END SELECT


  END SUBROUTINE Coeff_RK_init

END MODULE mod_declaration 