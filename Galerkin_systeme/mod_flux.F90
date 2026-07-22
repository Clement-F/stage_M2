MODULE mod_flux
    use mod_polynome
  IMPLICIT NONE

CONTAINS



  FUNCTION flux(u)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(in) :: u
    REAL(prec), DIMENSION(nb_var) :: flux

    flux = eps0
    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      flux = vit_adv*u
    CASE("burgers_SCL")
      flux = 0.5_prec * u**2 
    CASE("Buckley")
      flux = 4._prec*u**2/((4._prec*u**2)+ (1._prec-u)**2 )
    CASE("burgers")
      flux(1) = 0.5_prec * u(1)**2 
      flux(2) = u(1)*u(2)
    CASE("Shallow")
      flux(1) = u(2)
      flux(2) = u(1)*u(2) + (u(1)**2)/2._prec
    CASE("Euler")
      flux(1) = u(2)
      IF(u(1) .GT. eps0 ) flux(2) = u(2)**2/u(1)  + pression(u)
      IF(u(1) .GT. eps0 ) flux(3)  = (u(2)/u(1))*(u(3)+ pression(u)) 
    CASE DEFAULT
      WRITE(*,*) "flux non reconnu _flux",flux_name
      flux = 0._prec
      STOP
    END SELECT
  END FUNCTION flux

  FUNCTION pression(u) result(p)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(in) :: u
    REAL(prec) :: p

    p = (gamma_iso-1._prec)*(U(3) - 0.5_prec*(U(2)*U(2))/U(1))

  END FUNCTION pression

  FUNCTION flux_uh(x,ni, nvar)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ni,nvar
    REAL(prec), INTENT(IN) :: x
    REAL(prec) :: flux_uh
    REAL(prec), DIMENSION(nb_var) :: U, flux_temp
    INTEGER :: i

    DO i=1,nb_var
      U(i) = eval_step(x,nvar=i,ni=ni, LOC=LLoc)
    END DO

    flux_temp = flux(U)
    flux_uh = flux_temp(nvar)

  END FUNCTION flux_uh

  FUNCTION Flux_FV(u,v)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: u,v
    REAL(prec), DIMENSION(nb_var) :: Flux_FV

    IF(flux_num == 0) Flux_FV = (flux(u) + flux(v) - max_dflux*(v-u))  * 0.5_prec
    IF(flux_num == 1) Flux_FV = (flux(u) + flux(v) - gamma_calc(u,v)*(v-u))  * 0.5_prec

  END FUNCTION Flux_FV

  FUNCTION gamma_calc(u,v)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: u,v
    REAL(prec) :: gamma_calc
    REAL(prec) :: u_step, gamma_temp
    INTEGER :: i
    
    gamma_calc = 0._prec
    
    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      gamma_calc = abs(abs(vit_adv))
    CASE("burgers") 
      gamma_calc = max(abs(u(1)), abs(v(1)))
    CASE("burgers_SCL") 
      gamma_calc = max(abs(u(1)), abs(v(1)))
      
    CASE("Shallow") 
      gamma_calc = max(abs(u(2)/u(1)), abs(v(2)/v(1)))
      IF(ISNAN(gamma_calc)) gamma_calc = eps0

    CASE("Euler") 
      IF(u(1) .LT. eps0 .AND. v(1) .LT. eps0) THEN;
        gamma_calc = eps0
      ELSE 
        gamma_calc = max(abs(u(2)/u(1)) + sqrt(gamma_iso* pression(u)/u(1)) , abs(v(2)/v(1))+ sqrt(gamma_iso* pression(v)/v(1)))
      END IF
      IF(ISNAN(gamma_calc)) gamma_calc = eps0

    CASE("Buckley")
      u_step = (max(u(1),v(1))-min(u(1),v(1)))/10._prec

      DO i=1,10
        gamma_temp = abs(flux_d(min(u(1),v(1))+REAL(i,prec)*u_step))
        gamma_calc = MAX(gamma_calc, gamma_temp)
      END DO
    CASE DEFAULT
      WRITE(*,*) "flux non reconnu _gamma ",flux_name
      gamma_calc = 0._prec
      STOP
    END SELECT

  END FUNCTION gamma_calc
  
  FUNCTION flux_d(u)
    IMPLICIT NONE
    REAL(prec), INTENT(in) :: u
    REAL(prec) :: flux_d

    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      flux_d = vit_adv
    CASE("burgers_SCL")
      flux_d = u  
    CASE("Buckley")
      IF((abs(u) .GT. eps0) .and. (abs(u-1._prec) .GT. eps0) )THEN
        flux_d = 8._prec*u* (4._prec*u**2 + (1._prec-u)**2 - u*(4._prec*u-(1._prec-u))) &
              &/(4._prec*u**2 + (1._prec-u)**2)**2
      ELSE 
        flux_d = eps0
      END IF

    CASE DEFAULT
      WRITE(*,*) "flux_d non reconnu  _dflux",flux_name
      flux_d = 0._prec
      STOP
    END SELECT

  END FUNCTION flux_d

  FUNCTION entropie_numerique(u)
    IMPLICIT NONE
    REAL(prec) :: entropie_numerique
    REAL(prec), DIMENSION(nb_var) :: u
    entropie_numerique = 0.5_prec * DOT_PRODUCT(u,u)
  END FUNCTION entropie_numerique


  

END MODULE mod_flux