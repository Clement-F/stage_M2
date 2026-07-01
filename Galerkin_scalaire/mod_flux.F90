MODULE mod_Flux
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

    Flux_FV = (flux(u) + flux(v) - max_dflux*(v-u))  * 0.5_prec

  END FUNCTION Flux_FV


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
        u_step = (max(u,v)-min(u,v))/20._prec

        DO i=1,20
          gamma_temp = abs(flux_d(min(u,v)+REAL(i,prec)*u_step))
          gamma_calc = MAX(gamma_calc, gamma_temp)

        END DO
      END IF

    END IF
  END FUNCTION gamma_calc

END MODULE mod_Flux
