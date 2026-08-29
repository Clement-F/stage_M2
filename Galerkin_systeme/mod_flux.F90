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
    CASE("burgers")
      flux = 0.5_prec * u**2 
    CASE("Buckley")
      flux = 4._prec*u**2/((4._prec*u**2)+ (1._prec-u)**2 )
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
      STOP "flux"
    END SELECT
  END FUNCTION flux

  SUBROUTINE Projection_Flux
    IMPLICIT NONE
    INTEGER :: ni,jj,kk,ii
    REAL(prec), DIMENSION(nb_var) :: f_loc
    REAL(prec), DIMENSION(nb_var)            :: u_loc
    REAL(prec), DIMENSION(size_base, nb_var) :: fh_loc
    DO ni= 1,nb_cell

      fh_loc = 0._prec

      DO kk =1,nb_nodes
        u_loc = sol_step(ni)%val_quad(kk,:)
        f_loc = flux(u_loc)
      DO jj =size_base,1,-1
              fh_loc(jj,:) = fh_loc(jj,:) + f_loc*sig_quad(jj,kk)*w_quad(kk)
      END DO
      END DO

      flux_h(ni)%flux_DG = MATMUL(Masse_inv,  fh_loc)

    END DO

  END SUBROUTINE Projection_Flux

  FUNCTION pression(u) result(p)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(in) :: u
    REAL(prec) :: p
    ! print *,u
    p = (gamma_iso-1._prec)*(U(3) - 0.5_prec*(U(2)*U(2))/U(1))
    ! if(p .LT. eps0) STOP " negative pressure"

    ! IF(p .LT. -eps0) STOP "negative pressure"

  END FUNCTION pression

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
      gamma_calc = abs(vit_adv)
    CASE("burgers") 
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
      ! print*, pression(u)
      ! print *,"gamma",gamma_calc
      IF(ISNAN(gamma_calc)) gamma_calc = eps0
      gamma_calc = 1.5_prec* gamma_calc

    CASE("Buckley")
      u_step = (max(u(1),v(1))-min(u(1),v(1)))/10._prec

      DO i=1,10
        gamma_temp = abs(flux_d(min(u(1),v(1))+REAL(i,prec)*u_step))
        gamma_calc = MAX(gamma_calc, gamma_temp)
      END DO
    CASE DEFAULT
      WRITE(*,*) "flux non reconnu _gamma ",flux_name
      gamma_calc = 0._prec
      STOP "gamma calc"
    END SELECT

  END FUNCTION gamma_calc
  
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
      IF((abs(u) .GT. eps0) .and. (abs(u-1._prec) .GT. eps0) )THEN
        flux_d = 8._prec*u* (4._prec*u**2 + (1._prec-u)**2 - u*(4._prec*u-(1._prec-u))) &
              &/(4._prec*u**2 + (1._prec-u)**2)**2
      ELSE 
        flux_d = eps0
      END IF

    CASE DEFAULT
      WRITE(*,*) "flux_d non reconnu  _dflux",flux_name
      flux_d = 0._prec
      STOP "flux_d"
    END SELECT

  END FUNCTION flux_d

  !===============================================
  !===============================================

  FUNCTION entropie_numerique(u)
    IMPLICIT NONE
    REAL(prec) :: entropie_numerique
    REAL(prec), DIMENSION(nb_var) :: u
    REAL(prec) :: ke = 1.00001_prec
    REAL(prec) :: epsi = 0.25_prec

  entropie_numerique = 0.5_prec * u(1)**2

    ! IF(nb_var == 1) THEN
    !   IF(entropie_num == 0) entropie_numerique = 0.5_prec * u(1)**2
    !   IF(entropie_num == 1) entropie_numerique = (abs(u(1)-ke)**(1+epsi) )/(1+epsi)
    !   RETURN
    ! END IF

    ! print *,"entrop :",entropie_numerique

    IF(nb_var == 1) THEN
      IF(entropie_num == 0) entropie_numerique = 0.5_prec * u(1)**2
      IF(entropie_num == 1) entropie_numerique = (abs(u(1)-ke)**(1+epsi) )/(1+epsi)
      RETURN
    END IF

    ! print *,"entrop :",entropie_numerique

  END FUNCTION entropie_numerique

  FUNCTION Flux_entrop(u)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: u
    REAL(prec) :: Flux_entrop

    ! IF(nb_var == 1) THEN
      IF(flux_name == "advection" .AND. entropie_num == 0) Flux_entrop = vit_adv *0.5_prec * u(1)**2
      IF(flux_name == "burgers" .AND. entropie_num == 0) Flux_entrop = u(1)**3 /3._prec
    ! END IF

  END FUNCTION Flux_entrop

  FUNCTION Flux_entrop_VF(ul,ur)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: ul,ur
    REAL(prec), DIMENSION(nb_var) :: vl,vr,u_Riemann
    REAL(prec) :: Flux_entrop_VF

    ! u_Riemann = (ul+ur)/2._prec - (Flux(ur)-Flux(ul))/(2._prec*max_dflux)
    ! vl = Var_entrop(ul); vr = Var_entrop(ur)
    ! Flux_entrop_VF = DOT_PRODUCT((vl+vr)/2._prec , Flux_FV(ul,ur)) - (entrop_pot_flux(ul)+entrop_pot_flux(ur))/2._prec
    
    
    IF(flux_num == 0) Flux_entrop_VF = (Flux_entrop(ul)+Flux_entrop(ur) - max_dflux*        (entropie_numerique(ur)-entropie_numerique(ul)))/2._prec
    IF(flux_num == 1) Flux_entrop_VF = (Flux_entrop(ul)+Flux_entrop(ur) - gamma_calc(ul,ur)*(entropie_numerique(ur)-entropie_numerique(ul)))/2._prec

    ! Flux_entrop_VF = Flux_entrop(u_Riemann)
  END FUNCTION Flux_entrop_VF
  
  FUNCTION Var_entrop(u)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var) :: Var_entrop
    REAL(prec), DIMENSION(nb_var) :: u

    Var_entrop = u
  END FUNCTION Var_entrop

  ! FUNCTION entrop_pot_flux(u)
  !   IMPLICIT NONE
  !   REAL(prec) :: entrop_pot_flux
  !   REAL(prec), DIMENSION(nb_var), INTENT(IN) ::u
    
  !   entrop_pot_flux = (Var_entrop(u)*u(1)) - Flux_entrop(u)

  ! END FUNCTION entrop_pot_flux
END MODULE mod_flux