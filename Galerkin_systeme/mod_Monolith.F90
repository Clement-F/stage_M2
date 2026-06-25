MODULE mod_Monolith
    use mod_flux
  IMPLICIT NONE

CONTAINS

  ! FUNCTION minmax_loc(mc,minmax)
  !   IMPLICIT NONE
  !   CHARACTER(len=3) :: minmax
  !   INTEGER, DIMENSION(2), INTENT(IN) :: mc
  !   INTEGER, DIMENSION(2) :: voi_L, voi_R
  !   REAL(prec) :: minmax_loc

  !   voi_L = subcells_(mc(1),mc(2))%L
  !   voi_R = subcells_(mc(1),mc(2))%R
  !   IF(TRIM(minmax) == "min") THEN
  !   IF(max_rule == 1)   minmax_loc= min( sol_step(mc   (1))%val_subcells(mc   (2)), &
  !                                       &sol_step(voi_L(1))%val_subcells(voi_L(2)), &
  !                                       &sol_step(voi_R(1))%val_subcells(voi_R(2)))
  !   IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = min_glob

  !   ELSE IF(TRIM(minmax) == "max") THEN
  !   IF(max_rule == 1)   minmax_loc= max( sol_step(mc   (1))%val_subcells(mc   (2)), &
  !                                       &sol_step(voi_L(1))%val_subcells(voi_L(2)), &
  !                                       &sol_step(voi_R(1))%val_subcells(voi_R(2)))
  !   IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = max_glob
  !   END IF


  ! END FUNCTION minmax_loc


  ! FUNCTION theta(mc,pv, DF)
  !     IMPLICIT NONE
  !     INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
  !     REAL(prec), INTENT(IN) :: DF
  !     REAL(prec) :: theta

  !     REAL(prec) :: ug,ud, f_FV
  !     REAL(prec) :: gamma_mp, param
  !     REAL(prec) :: alpha,beta, u_Riemann
  !     INTEGER, DIMENSION(2) :: voi_L, voi_R
      
  !     IF(max_rule == 0) THEN 
  !       theta = 1._prec

  !     ELSE IF(max_rule == 2) THEN
          
  !       ug = sol_step(mc(1))%val_subcells(mc(2))
  !       ud = sol_step(pv(1))%val_subcells(pv(2))
  !       gamma_mp = gamma_calc(ug,ud)
   
  !       u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

  !       beta = max_glob; alpha = min_glob


  !     ELSE IF(max_rule == 1) THEN
      
  !       ug = sol_step(mc(1))%val_subcells(mc(2))
  !       ud = sol_step(pv(1))%val_subcells(pv(2))
  !       gamma_mp = gamma_calc(ug,ud)
  !       u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

  !       IF(DF .LT. -eps0) THEN
  !         beta = minmax_loc(mc,"max")
  !         alpha= minmax_loc(pv,"min")

  !       ELSE IF(DF .GT. eps0) THEN
  !         beta = minmax_loc(pv,"max")
  !         alpha= minmax_loc(mc,"min")

  !       ELSE 

  !         beta = minmax_loc(mc,"max")
  !         beta = max(beta,minmax_loc(pv,"max"))
  !         alpha= minmax_loc(pv,"min")
  !         alpha= min(alpha,minmax_loc(mc,"min"))

  !         theta = 0._prec
  !       END IF
        
      
    
  !     END IF

  !     param = min(beta - u_Riemann, u_Riemann- alpha); IF(param .LT. eps0) param = 0._prec

  !     ! IF(theta .LT. eps0) theta = 0._prec

  !     IF(abs(DF) .LT. eps0) THEN; theta = 1._prec
  !     ELSE;                       theta = max(min(1._prec, abs(gamma_mp/DF) * param),0._prec)
  !     END IF

  !     voi_L = subcells_(mc(1),mc(2))%L;
  !     voi_R = subcells_(pv(1),pv(2))%R; 

  !     IF(((theta .GT. 1._prec) .or. (theta .LT.  0._prec)) .or. (param .LT. 0._prec))THEN
  !       print *,"================" ,mc,"---------", pv ,"============================"
  !       print *,"voi_L : ", voi_L, "voi_R : ", voi_R
  !       write(*, fmt="( 'stencil = (', e12.6, 2x,e12.6, 2x ,e12.6, 2x,e12.6 ,')')") sol_step(voi_L(1))%val_subcells(voi_L(2)), ug, ud , sol_step(voi_R(1))%val_subcells(voi_R(2)) 
  !       write(*, fmt="( 'sol interface : ', e12.6,1x, e12.6 )") ug,ud
  !       write(*, fmt="( 'sol Riemann : ', e12.6 )") u_Riemann
  !       write(*, fmt="( 'alpha,beta = ',e12.6,2x, e12.6 )") alpha, beta
  !       write(*, fmt="( 'theta = ', f10.6, ' gamma = ', f10.6)") theta, gamma_mp
  !       write(* ,fmt="( 'param = ', e12.6, ' DF = ', e12.6)") param, DF
  !     END IF

  !     ! IF(coeff_smooth == 2) THEN
  !     !   subcells_(mc(1),mc(2))%theta_cm = subcells_(mc(1),mc(2))%theta_cm + theta/2._prec
  !     !   subcells_(pv(1),pv(2))%theta_cm = subcells_(pv(1),pv(2))%theta_cm + theta/2._prec
  !     ! ELSE IF(coeff_smooth==1) THEN
  !     !   subcells_(mc(1),mc(2))%theta_cm = min(subcells_(mc(1),mc(2))%theta_cm, theta)
  !     !   subcells_(pv(1),pv(2))%theta_cm = min(subcells_(pv(1),pv(2))%theta_cm, theta)
  !     ! END IF


  ! END FUNCTION

  ! SUBROUTINE extrema_detect 
  !   IMPLICIT NONE

  !   REAL(prec), DIMENSION(nb_cell,nb_subcell) :: du,ddu, vL, vR
  !   INTEGER, DIMENSION(2) :: voi_L, voi_R
  !   REAL(prec) :: v_min,v_max
  !   INTEGER :: ni, jj
  !   LOGICAL :: x_L, x_R

  !   vL = 0._prec; vR = 0._prec

  !   ! print *, "-------------------",n_time,"-----------------------------"

  !   DO ni = 1,nb_cell
  !     ! print *,"------------------"
  !     DO jj =1,nb_subcell
  !       du( ni,jj) = DOT_PRODUCT( sol_step(ni)%base_poly, Projection_VF_d (jj,:))/cell_size(ni)
  !       ddu(ni,jj) = DOT_PRODUCT( sol_step(ni)%base_poly, Projection_VF_dd(jj,:))/cell_size(ni)

  !       vL(ni,jj) = du(ni,jj) - cell_size(ni)/2._prec * ddu(ni,jj)
  !       vR(ni,jj) = du(ni,jj) + cell_size(ni)/2._prec * ddu(ni,jj)


  !       ! write(*,fmt ="(e12.6,1x, e12.6,1x, e12.6)") sol_step(ni)%val_subcells(jj), du(ni,jj), ddu(ni,jj)
  !     END DO
  !   END DO


  !   subcells_(:,:)%extrema = .False.

  !   DO ni = 1,nb_cell
  !     DO jj =1,nb_subcell
  !       voi_L = Voisin_Face(ni,jj,'L'); voi_R = Voisin_Face(ni,jj,'R')
  !       v_min = min(du(ni,jj), du(voi_L(1),voi_L(2))) - eps0
  !       v_max = max(du(ni,jj), du(voi_L(1),voi_L(2))) + eps0

  !       IF(((vL(ni,jj) .LT. v_min) .or. (vL(ni,jj) .GT. v_max)) ) THEN  ! check Left
  !         v_min = min(du(ni,jj), du(voi_R(1),voi_R(2))) - eps0
  !         v_max = max(du(ni,jj), du(voi_R(1),voi_R(2))) + eps0
  !         ! print *,"Left"
  !         ! print *, vR(ni,jj), v_min, v_max
  !         IF(((vR(ni,jj) .LT. v_min) .or. (vR(ni,jj) .GT. v_max)) ) THEN ;
  !           subcells_(ni,jj)%extrema = .True.
  !           ! print *,"Right"
  !         END IF
  !       END IF
        
  !     END DO
  !   END DO


  ! END SUBROUTINE extrema_detect

END MODULE mod_Monolith