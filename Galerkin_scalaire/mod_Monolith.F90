MODULE mod_Monolith
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
    CASE("Buckley")
      flux = 4._prec*u**2/((4._prec*u**2)+ (1._prec-u)**2 )
    CASE DEFAULT
      WRITE(*,*) "flux non reconnu ",flux_name
      flux = 0._prec
      STOP
    END SELECT
  END FUNCTION flux

  FUNCTION flux_d(u)
    IMPLICIT NONE
    REAL(prec), INTENT(in) :: u
    REAL(prec) :: flux_d

    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      flux_d = vit_adv
    CASE("burgers")
      flux_d = u  
    CASE("Buckley")
      flux_d = 8._prec*u* (4._prec*u**2 + (1._prec-u)**2 - u*(4._prec*u-(1._prec-u)))/(4._prec*u**2 + (1._prec-u)**2)**2
    CASE DEFAULT
      WRITE(*,*) "flux_d non reconnu ",flux_name
      flux_d = 0._prec
      STOP
    END SELECT


  END FUNCTION flux_d

  FUNCTION flux_uh(x,ni)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ni
    REAL(prec), INTENT(IN) :: x
    REAL(prec) :: flux_uh

    flux_uh = flux(eval_step(x,ni, LOC=LLoc))

  END FUNCTION flux_uh

  FUNCTION Flux_FV(u,v)
    IMPLICIT NONE
    REAL(prec), INTENT(IN) :: u,v
    REAL(prec) :: Flux_FV

    Flux_FV = (flux(u) + flux(v) - max_dflux*(v-u))  * 0.5_prec

  END FUNCTION Flux_FV

  FUNCTION minmax_loc(mc,minmax)
    IMPLICIT NONE
    CHARACTER(len=3) :: minmax
    INTEGER, DIMENSION(2), INTENT(IN) :: mc
    INTEGER, DIMENSION(2) :: voi_L, voi_R
    REAL(prec) :: minmax_loc

    voi_L = subcells_(mc(1),mc(2))%L
    voi_R = subcells_(mc(1),mc(2))%R
    IF(TRIM(minmax) == "min") THEN
    IF(max_rule == 1)   minmax_loc= min( sol_step(mc   (1))%val_subcells(mc   (2)), &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2)), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2)))
    IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = min_glob

    ELSE IF(TRIM(minmax) == "max") THEN
    IF(max_rule == 1)   minmax_loc= max( sol_step(mc   (1))%val_subcells(mc   (2)), &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2)), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2)))
    IF(max_rule == 2 .or. max_rule == 0)  minmax_loc = max_glob
    END IF


  END FUNCTION minmax_loc


  FUNCTION gamma_calc(u,v)
    IMPLICIT NONE
    REAL(prec), INTENT(IN) :: u,v
    REAL(prec) :: gamma_calc
    REAL(prec) :: gamma_temp

    REAL(prec) :: u_step
    INTEGER :: i 
    
    gamma_calc = 0._prec
    
    IF(TRIM(flux_name) == "advection") THEN
      gamma_calc = abs(vit_adv)

    ! ELSE IF(TRIM(flux_name) == "Buckley") THEN
    !   gamma_calc = 2.34_prec
 
    ELSE
      gamma_calc = max( abs(flux_d(u)), abs(flux_d(v)))

      IF(.NOT. convex_flux) THEN
        u_step = (max(u,v)-min(u,v))/10._prec

        DO i=1,10
          gamma_temp = abs(flux_d(min(u,v)+REAL(i,prec)*u_step))

          IF(gamma_calc .LT. gamma_temp ) THEN
             gamma_calc = gamma_temp
          END IF

        END DO
      END IF

    END IF
  END FUNCTION gamma_calc

  FUNCTION theta(mc,pv, DF)
      IMPLICIT NONE
      INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
      REAL(prec), INTENT(IN) :: DF
      REAL(prec) :: theta

      REAL(prec) :: ug,ud, f_FV
      REAL(prec) :: gamma_mp, param
      REAL(prec) :: alpha,beta, u_Riemann
      INTEGER, DIMENSION(2) :: voi_L, voi_R

      IF(max_rule == 0) THEN 
        theta = 1._prec
      ELSE IF(max_rule == 2) THEN
          
        ug = sol_step(mc(1))%val_subcells(mc(2))
        ud = sol_step(pv(1))%val_subcells(pv(2))
        f_FV = Flux_FV(ug,ud);   gamma_mp = gamma_calc(ug,ud)
   
        u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

        beta = max_glob; alpha = min_glob


      ELSE IF(max_rule == 1) THEN
      
        ug = sol_step(mc(1))%val_subcells(mc(2))
        ud = sol_step(pv(1))%val_subcells(pv(2))
        f_FV = Flux_FV(ug,ud);   gamma_mp = gamma_calc(ug,ud)
        u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

        IF(DF .LT. -eps0) THEN

          beta = minmax_loc(mc,"max")
          alpha= minmax_loc(pv,"min")

        ELSE IF(DF .GT. eps0) THEN
          
          beta = minmax_loc(pv,"max")
          alpha= minmax_loc(mc,"min")

        ELSE 
          
          theta = 0._prec
        END IF
        
      
    
      END IF

      param = min(beta - u_Riemann, u_Riemann- alpha); IF((param .LT. 10*eps0) .or.(DF .LT. 10*eps0)) param = 0._prec
      theta = min(1._prec, abs(gamma_mp/DF) * param)

      ! theta = max(theta, eps0)
      
      voi_L = subcells_(mc(1),mc(2))%L;
      voi_R = subcells_(pv(1),pv(2))%R; 

      IF(((theta .GT. 1._prec) .or. (theta .LT.  0._prec )) .or. (param .LT. 0._prec  ))THEN
        print *,"================" ,mc,"---------", pv ,"============================"
        print *,"voi_L : ", voi_L, "voi_R : ", voi_R
        write(*, fmt="( 'stencil = (', e12.6, 2x,e12.6, 2x ,e12.6, 2x,e12.6 ')')") sol_step(voi_L(1))%val_subcells(voi_L(2)), ug, ud , sol_step(voi_R(1))%val_subcells(voi_R(2)) 
        write(*, fmt="( 'sol Riemann : ', e12.6 )") u_Riemann
        write(*, fmt="( 'alpha,beta = ',e12.6,2x, e12.6 )") alpha, beta
        write(*, fmt="( 'theta = ', f10.6, ' gamma = ', f10.6)") theta, gamma_mp
        write(* ,fmt="( 'param = ', e12.6, ' DF = ', e12.6)") param, DF
      END IF

      IF(coeff_smooth == 2) THEN
        subcells_(mc(1),mc(2))%theta_cm = subcells_(mc(1),mc(2))%theta_cm + theta/2._prec
        subcells_(pv(1),pv(2))%theta_cm = subcells_(pv(1),pv(2))%theta_cm + theta/2._prec
      ELSE IF(coeff_smooth==1) THEN
        subcells_(mc(1),mc(2))%theta_cm = min(subcells_(mc(1),mc(2))%theta_cm, theta)
        subcells_(pv(1),pv(2))%theta_cm = min(subcells_(pv(1),pv(2))%theta_cm, theta)
      END IF


  END FUNCTION

END MODULE mod_Monolith