MODULE mod_flux
    use mod_polynome
  IMPLICIT NONE

CONTAINS



  FUNCTION flux(u)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(in) :: u
    REAL(prec), DIMENSION(nb_var) :: flux

    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      flux = vit_adv*u
    CASE("burgers")
      flux(1) = 0.5_prec * u(1)**2 
      flux(2) = u(1)*u(2)
    CASE("Shallow")
      flux(1) = u(2)
      flux(2) = u(1)*u(2) + (u(1)**2)/2._prec
    CASE("Euler")
      flux(1) = u(2)
      flux(2) = u(2)**2/u(1)  + pression(u)
      flux(3) = (u(2)/u(1))*(u(3)+ pression(u)) 
    CASE DEFAULT
      WRITE(*,*) "flux non reconnu ",flux_name
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


  ! FUNCTION flux_d(u)
    !   IMPLICIT NONE
    !   REAL(prec), INTENT(in) :: u
    !   REAL(prec) :: flux_d

    !   SELECT CASE (TRIM(flux_name))
    !   CASE("advection")
    !     flux_d = vit_adv
    !   CASE("burgers")
    !     flux_d = u  
    !   CASE("Buckley")
    !     ! print *,abs(u),eps0, (abs(u) .GT. eps0)
    !     IF((abs(u) .GT. eps0) .and. (abs(u-1._prec) .GT. eps0) )THEN
    !       ! print *,"hey"
    !       flux_d = 8._prec*u* (4._prec*u**2 + (1._prec-u)**2 - u*(4._prec*u-(1._prec-u))) &
    !             &/(4._prec*u**2 + (1._prec-u)**2)**2
    !     ELSE 
    !       flux_d = eps0
    !     END IF

    !     ! print *,flux_d

    !   CASE DEFAULT
    !     WRITE(*,*) "flux_d non reconnu ",flux_name
    !     flux_d = 0._prec
    !     STOP
    !   END SELECT


  ! END FUNCTION flux_d

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

    Flux_FV = (flux(u) + flux(v) - max_dflux*(v-u))  * 0.5_prec

  END FUNCTION Flux_FV

  FUNCTION gamma_calc(u,v)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: u,v
    REAL(prec) :: gamma_calc
    
    gamma_calc = 0._prec
    
    IF(TRIM(flux_name) == "advection") THEN
        gamma_calc = abs(abs(vit_adv))
    ELSE IF(TRIM(flux_name) == "burgers") THEN
        gamma_calc = max(abs(u(1)), abs(v(1)))
    ELSE IF(TRIM(flux_name) == "Shallow") THEN
        gamma_calc = max(abs(u(2)/u(1)), abs(v(2)/v(1)))
    ELSE IF(TRIM(flux_name) == "Euler") THEN
        gamma_calc = max(abs(u(2)/u(1) ) + sqrt(gamma_iso* pression(u)/u(1)) , abs(v(2)/v(1))+ sqrt(gamma_iso* pression(v)/v(1)))
        ! print *, pression(u),u(1)
    ELSE 
      print *,"flux non reconnue "
      Stop
    END IF

  END FUNCTION gamma_calc

  


  

END MODULE mod_flux