MODULE mod_RKDG 
   use mod_polynome
  IMPLICIT NONE

CONTAINS

  FUNCTION flux(u)
    IMPLICIT NONE
    REAL(prec), INTENT(in) :: u
    REAL(prec) :: flux

    flux = max_dflux*u

  END FUNCTION flux

  FUNCTION flux_uh(x,ni)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ni
    REAL(prec), INTENT(IN) :: x
    REAL(prec) :: flux_uh

    flux_uh = flux(eval_sol(ni,x))

  END FUNCTION flux_uh

  SUBROUTINE flux_numerique
    IMPLICIT NONE
    INTEGER :: ii
    REAL(prec) :: ug,ud

    print *,"update flux_num"

    DO ii = 1,nb_cell
      sol(ii)%inter(1) = eval_sol(ii, x_cell(ii))
      sol(ii)%inter(2) = eval_sol(ii, x_cell(ii+1))
    END DO

    ! flux_h(1)       %inter(1) = 0._prec
    ! flux_h(nb_cell) %inter(2) = 0._prec

    DO ii = 1,nb_cell
      
      IF (ii ==1) THEN 
        ug = sol(nb_cell)  %inter(2)
        ud = sol(1)        %inter(1)
      ELSE 
        ug = sol(ii-1)%inter(2)
        ud = sol(ii)  %inter(1)
      END IF

      flux_h(ii)  %inter(1) = 0.5_prec * (flux(ug) + flux(ud) - max_dflux*(ud-ug))
      flux_h(ii-1)%inter(2) = flux_h(ii)%inter(1)      

      IF (ii ==1) flux_h(nb_cell)%inter(2) = flux_h(1)%inter(1)    
    END DO
    
    ! flux_h(1)       %inter(1) = 0._prec
    ! flux_h(nb_cell) %inter(2) = 0._prec

    print *,"flux_interface cree"

    DO ii=1,nb_cell
      CALL Projection_Pk(flux_uh,flux_h(ii)%base_poly,ii)
    END DO

    print *,"projection Pk du flux"
  END SUBROUTINE flux_numerique

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
       RK_beta(1,1)=1._prec
    CASE (2)
       RK_time(1)=0._prec
       RK_alpha(1,1)=1._prec
       RK_beta(1,1)=1._prec

       RK_time(2)=1._prec
       RK_alpha(2,1)=0.5_prec
       RK_alpha(2,2)=0.5_prec
       RK_beta(2,2)=0.5_prec
    CASE (3)
       RK_time(1)=0._prec
       RK_alpha(1,1)=1._prec
       RK_beta(1,1)=1._prec

       RK_time(2)=0.5_prec
       RK_alpha(2,1)=0.75_prec
       RK_alpha(2,2)=0.25_prec
       RK_beta(2,2)=0.25_prec

       RK_time(3)=1._prec
       RK_alpha(3,1)=1._prec/3._prec
       RK_alpha(3,3)=2._prec/3._prec
       RK_beta(3,3)=2._prec/3._prec
    CASE DEFAULT
       RK_time=0._prec; RK_alpha = 0._prec; RK_beta = 0._prec
       DO nrk2=1,nrk
          RK_time(nrk2)=1._prec
          RK_beta(nrk2,3)=1._prec/REAL(nrk+1-nrk2,prec)
       END DO
    END SELECT


  END SUBROUTINE Coeff_RK_init

  SUBROUTINE  Time_step
    IMPLICIT NONE
    
    INTEGER :: ni

    INTEGER :: ii,jj
    REAL(prec) :: ti
    REAL(prec), DIMENSION(order_x) :: sig_1, sig_2
    REAL(prec), DIMENSION(order_x) :: V_B, S_B, BB
    
    CALL flux_numerique

    dt = min(CFL*dx/(2*max_dflux),tmax-time)
    print *,"dt = ",dt

    DO ni =1,nb_cell

      print *,"-----------"

      DO ii=1,order_x
        sig_1(ii) = DG_base(Loc_to_Ref(ni,x_cell(ni)),ii); 
        sig_2(ii) = DG_base(Loc_to_Ref(ni,x_cell(ni+1)),ii); 
      END DO

      print *,sig_1,sig_2
      print *, "base poly bord"
      print *,flux_h(ni)%base_poly

      sol_step(1)%base_poly = sol(ni)%base_poly
      V_B = MATMUL(Rigid,flux_h(ni)%base_poly)/(cell_size(ni))
      S_B = -(flux_h(ni)%inter(2)*sig_2 - flux_h(ni)%inter(1)*sig_1)/(cell_size(ni))
      BB  = V_B + S_B 
      print *,"B", V_B, S_B, BB

      L_step(1,:) = MATMUL(Masse_inv, BB  ) 
      
      ! print *,"init du step fini"

      print *,sol_step(1)%base_poly
      DO ii = 1, order_t
        ti = time+ RK_time(ii)*dt
        sol_step(ii+1)%base_poly = 0._prec

        DO jj = 1,ii 
          print *,ii,jj
          ! print *,RK_alpha(ii,jj),RK_beta(ii,jj)
          ! print *,L_step(jj,:)
          sol_step(ii+1)%base_poly = sol_step(ii+1)%base_poly + RK_alpha(ii,jj)*sol_step(jj)%base_poly &
          & + RK_beta(ii,jj)*dt*L_step(jj,:)  
        END DO

        V_B = MATMUL(Rigid,flux_h(ni)%base_poly)/(cell_size(ni))
        S_B = -(flux_h(ni)%inter(2)*sig_2 - flux_h(ni)%inter(1)*sig_1)/(cell_size(ni))
        BB  = V_B + S_B         
        L_step(ii+1,:) = MATMUL(Masse_inv, BB )
        ! print *,"a"
        ! print *,sol_step(ii+1)%base_poly
      END DO
      print *,sol_step(order_t+1)%base_poly


      sol(ni)%base_poly  = sol_step(order_t+1)%base_poly
    END DO

    time = time +dt
    n_time = n_time +1
    ! print *,sol(2)%base_poly

  END SUBROUTINE Time_step

END MODULE