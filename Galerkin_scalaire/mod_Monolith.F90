MODULE mod_Monolith
   use mod_flux
   use mod_Divers
  IMPLICIT NONE

CONTAINS

  FUNCTION minmax_loc(mc,minmax)
    IMPLICIT NONE
    CHARACTER(len=3) :: minmax
    INTEGER, DIMENSION(2), INTENT(IN) :: mc

    INTEGER, DIMENSION(2) :: voi_L, voi_R
    REAL(prec) :: sol_mc
    REAL(prec) :: minmax_loc

    ! print *, mc

    IF(mc(2) .ne. 0) THEN 
    sol_mc = sol_step(mc   (1))%val_subcells(mc   (2))
    voi_L = subcells_(mc(1),mc(2))%L
    voi_R = subcells_(mc(1),mc(2))%R
    ELSE 
    IF(TRIM(minmax) == "min") sol_mc = minval(sol_step(mc   (1))%val_subcells(:))
    IF(TRIM(minmax) == "max") sol_mc = maxval(sol_step(mc   (1))%val_subcells(:))
    voi_L = subcells_(mc(1),1)%L
    voi_R = subcells_(mc(1),nb_subcell)%R
    END IF 
    IF(TRIM(minmax) == "min") THEN
    IF(max_rule == 1)   minmax_loc= min( sol_mc, &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2)), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2)))

    IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = min_glob

    ELSE IF(TRIM(minmax) == "max") THEN
    IF(max_rule == 1)   minmax_loc= max( sol_mc, &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2)), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2)))

    IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = max_glob
    END IF


  END FUNCTION minmax_loc

  FUNCTION theta(mc,pv,DF)
    IMPLICIT NONE
    INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
    REAL(prec), INTENT(IN) :: DF
    
    REAL(prec) :: theta_
    REAL(prec) :: ug,ud, u_Riemann
    REAL(prec) :: gamma_mp, param
    REAL(prec) :: alpha,beta

    REAL(prec) :: theta
    
    theta = 1._prec
    IF(max_rule == 0) THEN 
      theta = 1._prec
      return
    END IF
    
    IF(abs(DF) < eps0) THEN; theta = 0._prec; return; END IF

    ug = sol_step(mc(1))%val_subcells(mc(2))
    ud = sol_step(pv(1))%val_subcells(pv(2))
    gamma_mp = gamma_calc(ug,ud)
    u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

    IF(( .not. subcells_(mc(1),mc(2))%extrema) .and. ( .not. subcells_(pv(1),pv(2))%extrema) ) THEN

      IF(max_rule ==1) THEN
        IF((abs(DF)) < eps0) THEN; theta = 0._prec; return; END IF
        IF(DF .LT. -eps0) THEN
          beta = minmax_loc(mc,"max")
          alpha= minmax_loc(pv,"min")

        ELSE IF(DF .GT. eps0) THEN
          beta = minmax_loc(pv,"max")
          alpha= minmax_loc(mc,"min")
        END IF

      ELSE IF(max_rule==2) THEN; alpha= min_glob; beta = max_glob
      END IF

      param = min(beta - u_Riemann, u_Riemann- alpha); IF(param .LT. eps0) param = 0._prec
      theta_= max(min(1._prec, abs(gamma_mp/DF) * param),0._prec)           
      theta = min(theta,theta_)

    END IF

  END FUNCTION

  SUBROUTINE extrema_detect 
    IMPLICIT NONE
    INTEGER :: ni, n_sub 
    INTEGER , DIMENSION(2) :: nL,nR
    REAL(prec) :: du_L,du_R, du,ddu, vL,vR

    LOGICAL :: face_L, face_R
    
    subcells_(:,:)%extrema = .FALSE. 

    IF(smooth_extrema == 1) THEN 
      DO ni=1,nb_cell
        sol_step(ni)%deriv= MATMUL(Masse_inv, MATMUL(sol_step(ni)%base_poly, Rigid))
      END DO

      DO ni=1,nb_cell
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,1         , 'L')
        nR = Voisin_cell(ni,nb_subcell, 'R')

        du   = 1._prec/(cell_size(ni))      * (sol_step(ni   )%inter(2)-sol_step(ni   )%inter(1))
        du_L = 1._prec/(cell_size(nL(1)))   * (sol_step(nL(1))%inter(2)-sol_step(nL(1))%inter(1))
        du_R = 1._prec/(cell_size(nR(1)))   * (sol_step(nR(1))%inter(2)-sol_step(nR(1))%inter(1))

        ddu  = 2._prec/(cell_size(ni))**2  *&
        &(eval_poly(1._prec,ni,sol_step(ni)%deriv, LOC=LRef)  - eval_poly(-1._prec,ni,sol_step(ni)%deriv, LOC=LRef))

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
      END DO
    ELSE 
      DO ni=1,nb_cell
        sol_step(ni)%deriv= MATMUL(Masse_inv, MATMUL(sol_step(ni)%base_poly, Rigid))
      END DO    

      DO ni=1,nb_cell;  DO n_sub =1,nb_subcell
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,n_sub, 'L') 
        nR = Voisin_cell(ni,n_sub, 'R') 

        du     = DOT_PRODUCT(sol_step(ni   )%deriv, Projection_VF(n_sub,:))/( cell_size(ni)*0.5_prec)
        du_L   = DOT_PRODUCT(sol_step(nL(1))%deriv, Projection_VF(nL(2),:))/( cell_size(ni)*0.5_prec)
        du_R   = DOT_PRODUCT(sol_step(nR(1))%deriv, Projection_VF(nR(2),:))/( cell_size(ni)*0.5_prec)

        ddu  = 2._prec/(cell_size(ni)**2 * subcell_size(n_sub))  *&
        &(eval_poly(1._prec,n_sub,sol_step(ni)%deriv, LOC=LSub)  - eval_poly(-1._prec,n_sub,sol_step(ni)%deriv, LOC=LRef))

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

      END DO;           END DO 

    END IF

  END SUBROUTINE extrema_detect

END MODULE mod_Monolith