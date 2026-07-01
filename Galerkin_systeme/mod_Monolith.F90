MODULE mod_Monolith
    use mod_flux
  IMPLICIT NONE

CONTAINS

  FUNCTION minmax_loc(mc,minmax,nvar)
    IMPLICIT NONE
    CHARACTER(len=3) :: minmax
    INTEGER, DIMENSION(2), INTENT(IN) :: mc
    INTEGER, INTENT(IN) :: nvar

    INTEGER, DIMENSION(2) :: voi_L, voi_R
    REAL(prec) :: sol_mc
    REAL(prec) :: minmax_loc

    ! print *, mc

    IF(mc(2) .ne. 0) THEN 
    sol_mc = sol_step(mc   (1))%val_subcells(mc   (2),nvar)
    voi_L = subcells_(mc(1),mc(2))%L
    voi_R = subcells_(mc(1),mc(2))%R
    ELSE 
    IF(TRIM(minmax) == "min") sol_mc = minval(sol_step(mc   (1))%val_subcells(:,nvar))
    IF(TRIM(minmax) == "max") sol_mc = maxval(sol_step(mc   (1))%val_subcells(:,nvar))
    voi_L = subcells_(mc(1),1)%L
    voi_R = subcells_(mc(1),nb_subcell)%R
    END IF 
    IF(TRIM(minmax) == "min") THEN
    IF(max_rule == 1)   minmax_loc= min( sol_mc, &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2),nvar), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2),nvar))

    IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = min_glob

    ELSE IF(TRIM(minmax) == "max") THEN
    IF(max_rule == 1)   minmax_loc= max( sol_mc, &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2),nvar), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2),nvar))

    IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = max_glob
    END IF


  END FUNCTION minmax_loc

  FUNCTION theta(mc,pv, DF)
    IMPLICIT NONE
    INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: DF
    REAL(prec), DIMENSION(nb_var) :: theta_
    REAL(prec), DIMENSION(nb_var) :: ug,ud, u_Riemann
    REAL(prec) :: gamma_mp, param
    REAL(prec) :: alpha,beta

    REAL(prec) :: A,B,M
    INTEGER :: ii

    REAL(prec) :: theta
    
    theta = 1._prec
    IF(max_rule == 0) THEN 
      theta = 1._prec
      return
    END IF
    
    IF(minval(abs(DF)) < eps0) THEN; theta = 0._prec; return; END IF

    ug = sol_step(mc(1))%val_subcells(mc(2),:)
    ud = sol_step(pv(1))%val_subcells(pv(2),:)
    gamma_mp = gamma_calc(ug,ud)
    u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

    IF(TRIM(flux_name)=="Euler") THEN 
      ! positivité de rho
      theta_(1) = min(1._prec, abs(gamma_mp/DF(1))*u_Riemann(1))

      ! positivité de E//P
      A = 1/(gamma_mp**2) *(0.5_prec*abs(DF(2))**2 - theta_(1)*DF(1)*DF(3))
      B = 1/(gamma_mp)    *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - theta_(1)*u_Riemann(3)*DF(1))
      M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2
      theta_(2) = min(1._prec, M/(abs(B)+ max(eps0,A)))

      ! positivité des deux 
      theta = theta_(1)* theta_(2)       
    END IF

    IF(( .not. subcells_(mc(1),mc(2))%extrema) .and. ( .not. subcells_(pv(1),pv(2))%extrema) ) THEN
      DO ii=1,nb_var
      
        IF(max_rule ==1) THEN
          IF((abs(DF(ii))) < eps0) THEN; theta = 0._prec; return; END IF
          IF(DF(ii) .LT. -eps0) THEN
            beta = minmax_loc(mc,"max",nvar=ii)
            alpha= minmax_loc(pv,"min",nvar=ii)

          ELSE IF(DF(ii) .GT. eps0) THEN
            beta = minmax_loc(pv,"max",nvar=ii)
            alpha= minmax_loc(mc,"min",nvar=ii)
          END IF

        ELSE IF(max_rule==2) THEN; alpha= min_glob; beta = max_glob
        END IF

        param = min(beta - u_Riemann(ii), u_Riemann(ii)- alpha); IF(param .LT. eps0) param = 0._prec
        theta_(ii) = max(min(1._prec, abs(gamma_mp/DF(ii)) * param),0._prec)           
        theta = min(theta,theta_(ii))

      END DO
    END IF

  END FUNCTION

  SUBROUTINE extrema_detect 
    IMPLICIT NONE

    REAL(prec), DIMENSION(nb_cell, nb_subcell, nb_var) :: du,ddu, vL, vR
    INTEGER, DIMENSION(2) :: voi_L, voi_R
    REAL(prec) :: v_min,v_max
    INTEGER :: ni, jj, ii
    LOGICAL :: x_L, x_R

    vL = 0._prec; vR = 0._prec

    ! print *, "-------------------",n_time,"-----------------------------"

    DO ni = 1,nb_cell; DO jj =1,nb_subcell; DO ii= 1,nb_var
        du( ni,jj,ii) = DOT_PRODUCT( sol_step(ni)%base_poly(:,ii), Projection_VF_d (jj,:))/cell_size(ni)
        ddu(ni,jj,ii) = DOT_PRODUCT( sol_step(ni)%base_poly(:,ii), Projection_VF_dd(jj,:))/cell_size(ni)

        vL(ni,jj,ii) = du(ni,jj,ii) - cell_size(ni)/2._prec * ddu(ni,jj,ii)
        vR(ni,jj,ii) = du(ni,jj,ii) + cell_size(ni)/2._prec * ddu(ni,jj,ii)


        ! write(*,fmt ="(e12.6,1x, e12.6,1x, e12.6)") sol_step(ni)%val_subcells(jj), du(ni,jj), ddu(ni,jj)
    END DO;  END DO; END DO


    subcells_(:,:)%extrema = .False.

    DO ni = 1,nb_cell; DO jj =1,nb_subcell; DO ii= 1,nb_var
      voi_L = Voisin_Face(ni,jj,'L'); voi_R = Voisin_Face(ni,jj,'R')

      v_min = min(du(ni,jj,ii), du(voi_L(1),voi_L(2),ii)) - eps0
      v_max = max(du(ni,jj,ii), du(voi_L(1),voi_L(2),ii)) + eps0

      IF(((vL(ni,jj,ii) .LT. v_min) .or. (vL(ni,jj,ii) .GT. v_max)) ) THEN  ! check Left
        v_min = min(du(ni,jj,ii), du(voi_R(1),voi_R(2),ii)) - eps0
        v_max = max(du(ni,jj,ii), du(voi_R(1),voi_R(2),ii)) + eps0
        IF(((vR(ni,jj,ii) .LT. v_min) .or. (vR(ni,jj,ii) .GT. v_max)) ) THEN ;
          subcells_(ni,jj)%extrema = .True.
        END IF
      END IF

      IF(subcells_(ni,jj)%extrema) EXIT
      
    END DO;   END DO;   END DO


  END SUBROUTINE extrema_detect

END MODULE mod_Monolith