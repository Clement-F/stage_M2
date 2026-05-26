MODULE mod_RKDG 
   use mod_polynome
  IMPLICIT NONE

CONTAINS

  FUNCTION flux(u)
    IMPLICIT NONE
    REAL(prec), INTENT(in) :: u
    REAL(prec) :: flux

    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      flux = vit_adv*u
    CASE("burgers")
      flux = 0.5_prec * u**2  
    CASE DEFAULT
       WRITE(*,*) "flux non reconnu ",flux_name
       flux = 0._prec
       STOP
    END SELECT


  END FUNCTION flux

  FUNCTION flux_uh(x,ni)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ni
    REAL(prec), INTENT(IN) :: x
    REAL(prec) :: flux_uh

    flux_uh = flux(eval_step(ni,x))

  END FUNCTION flux_uh

  SUBROUTINE flux_numerique
    IMPLICIT NONE
    INTEGER :: ni
    REAL(prec) :: ug,ud

    ! print *, "flux"
    DO ni = 1,nb_cell
      sol_step(ni)%inter(1) = eval_step(ni, x_cell(ni))
      sol_step(ni)%inter(2) = eval_step(ni, x_cell(ni+1))
    END DO

    DO ni = 1,nb_cell
      IF (ni ==1) THEN 
        ug = sol_step(nb_cell)  %inter(2)
        ud = sol_step(1)        %inter(1)
      ELSE 
        ug = sol_step(ni-1)%inter(2)
        ud = sol_step(ni)  %inter(1)
      END IF

      flux_h(ni)  %inter(1) = 0.5_prec * (flux(ug) + flux(ud) - max_dflux*(ud-ug))  

      IF (ni ==1) THEN ;  flux_h(nb_cell)%inter(2) = flux_h(1) %inter(1)    
      ELSE ;              flux_h(ni-1)   %inter(2) = flux_h(ni)%inter(1)    
      END IF

      CALL Projection_Pk(flux_uh,flux_h(ni)%base_poly,ni)
    END DO

  END SUBROUTINE flux_numerique

  SUBROUTINE  Time_step
    IMPLICIT NONE
    
    INTEGER :: ni, tni

    INTEGER :: ii,jj
    REAL(prec) :: ti
    REAL(prec), DIMENSION(order_x) :: sig_1, sig_2
    REAL(prec), DIMENSION(order_x) :: V_B, S_B, BB

    ! print *,"calc time"


    DO tni =1,order_t
      
        ! print *,RK_alpha(tni,1),RK_alpha(tni,2), RK_beta(tni)

      DO ni=1,nb_cell
        sol_step(ni)%base_poly = sol(ni)%base_poly
      END DO 
      
      CALL flux_numerique

      DO ni =1,nb_cell

        DO ii=1,order_x
          sig_1(ii) = DG_base(Loc_to_Ref(ni,x_cell(ni)),ii); 
          sig_2(ii) = DG_base(Loc_to_Ref(ni,x_cell(ni+1)),ii); 
        END DO

        V_B = MATMUL(Rigid,flux_h(ni)%base_poly)
        S_B = -(flux_h(ni)%inter(2)*sig_2 - flux_h(ni)%inter(1)*sig_1)
        BB  = (V_B + S_B)
        L_step = MATMUL(Masse_inv, BB  )*(2._prec/(cell_size(ni))) 
        ! print *, V_B,S_b,BB
        ! print *,RK_alpha(tni,1),RK_alpha(tni,2), RK_beta(tni)

        sol_step(ni)%base_poly = RK_alpha(tni,1) * sol(ni)%base_poly + RK_alpha(tni,2) * sol_step(ni)%base_poly &
                             &+  RK_beta(tni) *dt * L_step
        
      END DO
    END DO

    ! print *,"---------------------"

    DO ni=1,nb_cell
        sol(ni)%base_poly  = sol_step(ni)%base_poly
    END DO

    time = time +dt
    n_time = n_time +1

  END SUBROUTINE Time_step

  SUBROUTINE dt_calc
    IMPLICIT NONE

    IF(TRIM(flux_name) == "advection") THEN
      max_dflux = vit_adv
    ELSE IF(TRIM(flux_name) == "burgers") THEN
      !boucle pour max de u ? 
      DO i=1,nb_cell
        IF(max_dflux .LT. sol(i)%base_poly(1)) THEN
          max_dflux = sol(i)%base_poly(1)
        END IF
      END DO
    ELSE 
      print *,"flux non reconnue"
      max_dflux = 1._prec
    END IF
    
    dt = min(CFL*dx/((2*(order_x-1)+1)*max_dflux),tmax-time)

  END SUBROUTINE dt_calc

END MODULE