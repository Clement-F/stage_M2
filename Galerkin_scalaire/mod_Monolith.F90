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
      ! print *,abs(u),eps0, (abs(u) .GT. eps0)
      IF((abs(u) .GT. eps0) .and. (abs(u-1._prec) .GT. eps0) )THEN
        ! print *,"hey"
        flux_d = 8._prec*u* (4._prec*u**2 + (1._prec-u)**2 - u*(4._prec*u-(1._prec-u))) &
              &/(4._prec*u**2 + (1._prec-u)**2)**2
      ELSE 
        flux_d = eps0
      END IF

      ! print *,flux_d

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

    Flux_FV = (flux(u) + flux(v) - gamma_calc(u,v)*(v-u))  * 0.5_prec

  END FUNCTION Flux_FV


  FUNCTION minmax_cell(mc,minmax)
    IMPLICIT NONE
    CHARACTER(len=3) :: minmax
    INTEGER, DIMENSION(2), INTENT(IN) :: mc
    INTEGER, DIMENSION(2) :: voi_L, voi_R
    REAL(prec) :: minmax_cell

    voi_L = subcells_(mc(1),1)%L
    voi_R = subcells_(mc(1),nb_subcell)%R

    IF(TRIM(minmax) == "min") THEN
    IF(max_rule == 1)   minmax_cell= minval(min(sol_step(mc   (1))%val_subcells(:), &
                                               &sol_step(voi_L(1))%val_subcells(:), &
                                               &sol_step(voi_R(1))%val_subcells(:)))
    IF(max_rule == 2 .or. max_rule == 0)  minmax_cell = min_glob

    ELSE IF(TRIM(minmax) == "max") THEN
    IF(max_rule == 1)   minmax_cell= maxval(max( sol_step(mc   (1))%val_subcells(mc   (:)), &
                                                &sol_step(voi_L(1))%val_subcells(voi_L(:)), &
                                                &sol_step(voi_R(1))%val_subcells(voi_R(:))))
    IF(max_rule == 2 .or. max_rule == 0)  minmax_cell = max_glob
    END IF


  END FUNCTION minmax_cell


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

 
    ELSE
      gamma_calc = max(abs(flux_d(u)), abs(flux_d(v)))

      IF(.NOT. convex_flux) THEN
        u_step = (max(u,v)-min(u,v))/10._prec

        DO i=1,10
          gamma_temp = abs(flux_d(min(u,v)+REAL(i,prec)*u_step))
          gamma_calc = MAX(gamma_calc, gamma_temp)

        END DO
      END IF

    END IF
  END FUNCTION gamma_calc

  FUNCTION theta(mc,pv, DF, rule)
      IMPLICIT NONE
      INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
      INTEGER, INTENT(IN) :: rule
      REAL(prec), INTENT(IN) :: DF
      REAL(prec) :: theta

      REAL(prec) :: ug,ud
      REAL(prec) :: gamma_mp, param
      REAL(prec) :: alpha,beta, u_Riemann
      INTEGER, DIMENSION(2) :: voi_L, voi_R
      
      IF(max_rule == 0) THEN 
        theta = 1._prec
        return
      ELSE IF(rule == 2) THEN
          
        ug = sol_step(mc(1))%val_subcells(mc(2))
        ud = sol_step(pv(1))%val_subcells(pv(2))
        gamma_mp = gamma_calc(ug,ud)
   
        u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

        beta = max_glob; alpha = min_glob


      ELSE IF(rule == 1) THEN
        ug = sol_step(mc(1))%val_subcells(mc(2))
        ud = sol_step(pv(1))%val_subcells(pv(2))
        gamma_mp = gamma_calc(ug,ud)
        u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

        IF(DF .LT. -eps0) THEN
          beta = minmax_loc(mc,"max")
          alpha= minmax_loc(pv,"min")

        ELSE IF(DF .GT. eps0) THEN
          beta = minmax_loc(pv,"max")
          alpha= minmax_loc(mc,"min")

        ELSE 

          beta = minmax_loc(mc,"max")
          beta = max(beta,minmax_loc(pv,"max"))
          alpha= minmax_loc(pv,"min")
          alpha= min(alpha,minmax_loc(mc,"min"))

          theta = 0._prec
          return
        END IF
        
      
    
      END IF

      param = min(beta - u_Riemann, u_Riemann- alpha); IF(param .LT. eps0) param = 0._prec

      IF(abs(DF) .LT. eps0) THEN; theta = 1._prec
      ELSE;                       theta = max(min(1._prec, abs(gamma_mp/DF) * param),0._prec)
      END IF

      voi_L = subcells_(mc(1),mc(2))%L;
      voi_R = subcells_(pv(1),pv(2))%R; 

      IF(((theta .GT. 1._prec) .or. (theta .LT.  0._prec)) .or. (param .LT. 0._prec))THEN
        print *,"================" ,mc,"---------", pv ,"============================"
        print *,"voi_L : ", voi_L, "voi_R : ", voi_R
        write(*, fmt="( 'stencil = (', e12.6, 2x,e12.6, 2x ,e12.6, 2x,e12.6 ,')')") sol_step(voi_L(1))%val_subcells(voi_L(2)), ug, ud , sol_step(voi_R(1))%val_subcells(voi_R(2)) 
        write(*, fmt="( 'sol interface : ', e12.6,1x, e12.6 )") ug,ud
        write(*, fmt="( 'sol Riemann : ', e12.6 )") u_Riemann
        write(*, fmt="( 'alpha,beta = ',e12.6,2x, e12.6 )") alpha, beta
        write(*, fmt="( 'theta = ', f10.6, ' gamma = ', f10.6)") theta, gamma_mp
        write(* ,fmt="( 'param = ', e12.6, ' DF = ', e12.6)") param, DF
      END IF

      return 
      IF(coeff_smooth == 2) THEN
        subcells_(mc(1),mc(2))%theta_cm = subcells_(mc(1),mc(2))%theta_cm + theta/2._prec
        subcells_(pv(1),pv(2))%theta_cm = subcells_(pv(1),pv(2))%theta_cm + theta/2._prec
      ELSE IF(coeff_smooth==1) THEN
        subcells_(mc(1),mc(2))%theta_cm = min(subcells_(mc(1),mc(2))%theta_cm, theta)
        subcells_(pv(1),pv(2))%theta_cm = min(subcells_(pv(1),pv(2))%theta_cm, theta)
      END IF


  END FUNCTION

  SUBROUTINE extrema_detect 
    IMPLICIT NONE

    REAL(prec), DIMENSION(nb_cell,nb_subcell) :: du,ddu, vL, vR
    INTEGER, DIMENSION(2) :: voi_L, voi_R
    REAL(prec) :: v_min,v_max
    INTEGER :: ni, jj

    vL = 0._prec; vR = 0._prec

    ! print *, "-------------------",n_time,"-----------------------------"

    DO ni = 1,nb_cell
      ! print *,"------------------"
      DO jj =1,nb_subcell
        du( ni,jj) = DOT_PRODUCT( sol_step(ni)%base_poly, Projection_VF_d (jj,:))/cell_size(ni)
        ddu(ni,jj) = DOT_PRODUCT( sol_step(ni)%base_poly, Projection_VF_dd(jj,:))/(cell_size(ni)**2)

        vL(ni,jj) = du(ni,jj) - subcell_size(jj)*cell_size(ni)/2._prec * ddu(ni,jj)
        vR(ni,jj) = du(ni,jj) + subcell_size(jj)*cell_size(ni)/2._prec * ddu(ni,jj)


        ! write(*,fmt ="(e12.6,1x, e12.6,1x, e12.6, 2x,'x =', f10.6 ) ") sol_step(ni)%val_subcells(jj), du(ni,jj), ddu(ni,jj), Ref_to_loc(ni=ni,XX=x_submiddle(jj))
      END DO
    END DO


    subcells_(:,:)%extrema = .False.

    DO ni = 1,nb_cell
      DO jj =1,nb_subcell
        voi_L = Voisin_Face(ni,jj,'L'); voi_R = Voisin_Face(ni,jj,'R')
        v_min = min(du(ni,jj), du(voi_L(1),voi_L(2))) - eps0
        v_max = max(du(ni,jj), du(voi_L(1),voi_L(2))) + eps0

        IF(((vL(ni,jj) .LT. v_min) .or. (vL(ni,jj) .GT. v_max)) ) THEN  ! check Left
          v_min = min(du(ni,jj), du(voi_R(1),voi_R(2))) - eps0
          v_max = max(du(ni,jj), du(voi_R(1),voi_R(2))) + eps0
          ! print *,"Left"
          ! print *, vR(ni,jj), v_min, v_max
          IF(((vR(ni,jj) .LT. v_min) .or. (vR(ni,jj) .GT. v_max)) ) THEN ;
            subcells_(ni,jj)%extrema = .True.
            ! print *,"Right"
          END IF
        END IF
        
      END DO
    END DO


  END SUBROUTINE extrema_detect

END MODULE mod_Monolith