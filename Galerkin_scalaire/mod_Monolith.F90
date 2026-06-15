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


  FUNCTION gamma_calc(u,v)
    IMPLICIT NONE
    REAL(prec), INTENT(IN) :: u,v
    REAL(prec) :: gamma_calc

    REAL(prec) :: u_step
    INTEGER :: i 
    
    gamma_calc = 0._prec
    
    IF(TRIM(flux_name) == "advection") THEN
      gamma_calc = abs(vit_adv)

    ELSE
      gamma_calc = max( abs(flux_d(u)), abs(flux_d(v)))

      IF(.NOT. convex_flux) THEN
        u_step = (max(u,v)-min(u,v))/10._prec
        DO i=1,10
          IF(gamma_calc .LT. abs(flux_d(min(u,v)+REAL(i,prec)*u_step)) ) THEN
            gamma_calc = abs(flux_d(min(u,v)+REAL(i,prec)*u_step))
          END IF
        END DO
      END IF
    END IF
  END FUNCTION gamma_calc

  FUNCTION theta(mc,pv)
      IMPLICIT NONE
      INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
      REAL(prec) :: theta

      REAL(prec) :: ug,ud, f_FV
      REAL(prec) :: gamma_mp, DF, param
      REAL(prec) :: alpha,beta, u_Riemann
      INTEGER, DIMENSION(2) :: voi, voi_L, voi_R

      IF(max_rule == 0) THEN 
        theta = 1._prec
      ELSE IF(max_rule == 2) THEN
          
          ug = sol_step(mc(1))%val_subcells(mc(2))
          ud = sol_step(pv(1))%val_subcells(pv(2))
          f_FV = Flux_FV(ug,ud);   gamma_mp = gamma_calc(ug,ud)
          DF = flux_h(pv(1))%val_subcells(pv(2)) - f_FV    
          u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

          param = min(max_glob -u_Riemann, u_Riemann - min_glob)

          theta = min(1._prec, abs(gamma_mp/DF)*param)


      ELSE IF(max_rule == 1) THEN
        
          ! ug = sol_step(mc(1))%val_subcells(mc(2))
          ! ud = sol_step(pv(1))%val_subcells(pv(2))
          ! f_FV = Flux_FV(ug,ud);   gamma_mp = gamma_calc(ug,ud)
          ! u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)
          ! IF(DF < 0._prec) THEN

          !     voi = Voisin_Face(ni,jj,'L');
          !     voi_L = subcells_(voi(1),voi(2))%L;
          !     voi_R = subcells_(voi(1),voi(2))%R;  

          !     ! print *,'L',voi
          !     beta  = max(sol_step(voi  (1))%val_subcells(voi  (2)), &
          !                 &   sol_step(voi_L(1))%val_subcells(voi_L(2)), &
          !                 &   sol_step(voi_R(1))%val_subcells(voi_R(2)))

          !     voi = Voisin_Face(ni,jj,'R');           
          !     voi_L = subcells_(voi(1),voi(2))%L;
          !     voi_R = subcells_(voi(1),voi(2))%R;  

          !     ! print *,'R',voi
          !     alpha  = min(sol_step(voi  (1))%val_subcells(voi  (2)), &
          !                 &   sol_step(voi_L(1))%val_subcells(voi_L(2)), &
          !                 &   sol_step(voi_R(1))%val_subcells(voi_R(2)))

          !     param =     min(beta - u_Riemann, u_Riemann- alpha)
          ! ELSE

          !     voi = Voisin_Face(ni,jj,'L');           
          !     voi_L = subcells_(voi(1),voi(2))%L;
          !     voi_R = subcells_(voi(1),voi(2))%R;  
          !     ! print *,'L',voi
          !     alpha  = min(sol_step(voi  (1))%val_subcells(voi  (2)), &
          !                 &   sol_step(voi_L(1))%val_subcells(voi_L(2)), &
          !                 &   sol_step(voi_R(1))%val_subcells(voi_R(2)))

          !     voi = Voisin_Face(ni,jj,'R');           
          !     voi_L = subcells_(voi(1),voi(2))%L;
          !     voi_R = subcells_(voi(1),voi(2))%R;  
          !     ! print *,'R',voi
          !     beta  = max(sol_step(voi  (1))%val_subcells(voi  (2)), &
          !                 &   sol_step(voi_L(1))%val_subcells(voi_L(2)), &
          !                 &   sol_step(voi_R(1))%val_subcells(voi_R(2)))

          !     param =     max(min(beta - u_Riemann, u_Riemann- alpha),0._prec)
          ! END IF
          !     theta = min(1._prec, abs(gamma_mp/DF) * param)
      
      
      END IF


  END FUNCTION

END MODULE mod_Monolith