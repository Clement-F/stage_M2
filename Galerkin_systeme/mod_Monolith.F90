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
    INTEGER :: ni, n_sub, n_var 
    INTEGER , DIMENSION(2) :: nL,nR
    REAL(prec) :: du_L,du_R, du,ddu, vL,vR

    LOGICAL :: face_L, face_R
    
    subcells_(:,:)%extrema = .FALSE. 

    IF(smooth_extrema == 1) THEN 
      DO ni=1,nb_cell;  DO n_var=1,nb_var
        sol_step(ni)%deriv(:,n_var)= MATMUL(Masse_inv, MATMUL(sol_step(ni)%base_poly(:,n_var), Rigid))
      END DO;           END DO

      DO ni=1,nb_cell; DO n_var=1,nb_var
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,1         , 'L')
        nR = Voisin_cell(ni,nb_subcell, 'R')

        du   = 1._prec/(cell_size(ni))      * (sol_step(ni   )%inter(2,n_var)-sol_step(ni   )%inter(1,n_var))
        du_L = 1._prec/(cell_size(nL(1)))   * (sol_step(nL(1))%inter(2,n_var)-sol_step(nL(1))%inter(1,n_var))
        du_R = 1._prec/(cell_size(nR(1)))   * (sol_step(nR(1))%inter(2,n_var)-sol_step(nR(1))%inter(1,n_var))

        ddu  = 2._prec/(cell_size(ni))**2  *&
        &(eval_poly(1._prec,ni,sol_step(ni)%deriv(:,n_var), LOC=LRef)  - eval_poly(-1._prec,ni,sol_step(ni)%deriv(:,n_var), LOC=LRef))

        vL = du - 0.5_prec*cell_size(ni)*ddu
        vR = du + 0.5_prec*cell_size(ni)*ddu

        IF((MIN(du,du_L)-eps0 .LT. vL) .AND. (MAX(du,du_L)+eps0 .GT. vL )) face_L = .TRUE.
        IF((MIN(du,du_R)-eps0 .LT. vR) .AND. (MAX(du,du_R)+eps0 .GT. vR )) face_R = .TRUE.

        ! print *,"----------------------------------------------------------"
        ! write(*,fmt='("cell : ",i3,", [",f10.6,",",f10.6,"], ")') ni, x_cell(ni), x_cell(ni+1)
        ! write(*,fmt='("ddu  : ",e12.6)') ddu
        ! write(*,fmt='("Left  : (du=",e12.6,", du_L=",e12.6,", vL =",e12.6,", check :",l1,")" )') du, du_L,vL, face_L
        ! write(*,fmt='("Right : (du=",e12.6,", du_R=",e12.6,", vR =",e12.6,", check :",l1,")" )') du, du_R,vR, face_R

        IF( face_L .AND. face_R) THEN
          subcells_(ni,:)%extrema = .TRUE.
        END IF
      END DO; END DO
    ELSE 
      DO ni=1,nb_cell;  DO n_var=1,nb_var
        sol_step(ni)%deriv(:,n_var)= MATMUL(Masse_inv, MATMUL(sol_step(ni)%base_poly(:,n_var), Rigid))
      END DO;           END DO

      DO ni=1,nb_cell;  DO n_sub =1,nb_subcell; DO n_var = 1,nb_var
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,n_sub, 'L') 
        nR = Voisin_cell(ni,n_sub, 'R') 

        du     = DOT_PRODUCT(sol_step(ni   )%deriv(:,n_var), Projection_VF(n_sub,:))/( cell_size(ni)*0.5_prec)
        du_L   = DOT_PRODUCT(sol_step(nL(1))%deriv(:,n_var), Projection_VF(nL(2),:))/( cell_size(ni)*0.5_prec)
        du_R   = DOT_PRODUCT(sol_step(nR(1))%deriv(:,n_var), Projection_VF(nR(2),:))/( cell_size(ni)*0.5_prec)

        ddu  = 2._prec/(cell_size(ni)**2 * subcell_size(n_sub))  *&
        &(eval_poly(1._prec,n_sub,sol_step(ni)%deriv(:,n_var), LOC=LSub)  - eval_poly(-1._prec,n_sub,sol_step(ni)%deriv(:,n_var), LOC=LRef))

        vL = du - 0.5_prec*cell_size(ni)*subcell_size(n_sub)*ddu
        vR = du + 0.5_prec*cell_size(ni)*subcell_size(n_sub)*ddu

        IF((MIN(du,du_L)-eps0 .LT. vL) .AND. (MAX(du,du_L)+eps0 .GT. vL )) face_L = .TRUE.
        IF((MIN(du,du_R)-eps0 .LT. vR) .AND. (MAX(du,du_R)+eps0 .GT. vR )) face_R = .TRUE.

        ! print *,"----------------------------------------------------------"
        ! write(*,fmt='("cell : (",i3,",",i3,"), [",f10.6,",",f10.6,"], ")') ni,n_sub, Ref_to_loc(ni,x_subcell(n_sub)), Ref_to_loc(ni,x_subcell(n_sub+1))
        ! write(*,fmt='("ddu  : ",e12.6)') ddu
        ! write(*,fmt='("Left  : (du=",e12.6,", du_L=",e12.6,", vL =",e12.6,", check :",l1,")" )') du, du_L,vL, face_L
        ! write(*,fmt='("Right : (du=",e12.6,", du_R=",e12.6,", vR =",e12.6,", check :",l1,")" )') du, du_R,vR, face_R

        IF( face_L .AND. face_R) THEN
          subcells_(ni,:)%extrema = .TRUE.
        END IF

      END DO;           END DO;                  END DO

    END IF

  END SUBROUTINE extrema_detect


END MODULE mod_Monolith