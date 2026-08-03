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
    REAL(prec), DIMENSION(nb_var) :: theta_temp
    REAL(prec), DIMENSION(nb_var) :: ug,ud, u_Riemann
    REAL(prec) :: gamma_mp, param
    REAL(prec) :: alpha,beta

    REAL(prec) :: A,B,M
    INTEGER :: ii
    INTEGER, DIMENSION(2) :: voi_L, voi_R
    LOGICAL :: extrema

    REAL(prec) :: theta
    
    theta = 1._prec
    IF(max_rule == 0) THEN 
      theta = 1._prec
      return
    END IF
    
    IF(minval(abs(DF)) < eps0) THEN; theta = 1._prec; return; END IF

    ug = sol_step(mc(1))%val_subcells(mc(2),:)
    ud = sol_step(pv(1))%val_subcells(pv(2),:)
    
    IF(flux_num == 0) gamma_mp = max_dflux
    IF(flux_num == 1) gamma_mp = gamma_calc(ug,ud)
    
    u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

    IF(TRIM(flux_name)=="Euler") THEN 
      ! positivité de rho
      theta_temp(1) = min(1._prec, abs(gamma_mp/DF(1))*u_Riemann(1))

      ! positivité de E//P
      A = 1/(gamma_mp**2) *(0.5_prec*abs(DF(2))**2 - theta_temp(1)*DF(1)*DF(3))
      B = 1/(gamma_mp)    *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - theta_temp(1)*u_Riemann(3)*DF(1))
      M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2
      theta_temp(2) = min(1._prec, M/(abs(B)+ max(eps0,A)))

      ! positivité des deux 
      theta = theta_temp(1)* theta_temp(2)       
    END IF

    extrema = subcells_(mc(1),mc(2))%extrema .AND.  subcells_(pv(1),pv(2))%extrema 

    IF(.not. extrema) THEN
      DO ii=1,nb_var
      
        IF(max_rule ==1) THEN
          IF((abs(DF(ii))) < eps0) THEN; theta = 1._prec; return; END IF
          IF(DF(ii) .LT. -eps0) THEN
            beta = minmax_loc(mc,"max",nvar=ii) 
            alpha= minmax_loc(pv,"min",nvar=ii)

          ELSE IF(DF(ii) .GT. eps0) THEN
            beta = minmax_loc(pv,"max",nvar=ii) 
            alpha= minmax_loc(mc,"min",nvar=ii)
          END IF

        ELSE IF(max_rule==2) THEN; alpha= (min_glob+eps0); beta = (max_glob-eps0)
        END IF

        param = min(beta - u_Riemann(ii), u_Riemann(ii)- alpha)

        !! ATTENTION !!
        ! 3 lignes du dessous sont peu justifier
        ! IF(abs(beta-alpha) .LT. 2*eps0) THEN; theta_temp(ii) = 1._prec;
        ! ELSE IF(param .LT. eps0) THEN ; theta_temp(ii) = 1._prec;
        ! ELSE;
        
        theta_temp(ii) =max(min(1._prec, abs(gamma_mp/DF(ii)) * param),0._prec);
        ! END IF      
        
        theta = min(theta,theta_temp(ii))

      END DO
    END IF

    theta = max(theta,0._prec)

    IF(ISNAN(theta)) print *,"theta nan"

    ! IF(theta .LT. eps0) THEN 
    !   voi_L = subcells_(mc(1),mc(2))%L;
    !   voi_R = subcells_(pv(1),pv(2))%R; 
    !   print *,"================" ,mc,"---------", pv ,"============================"
    !   write(*, fmt="( 'between : [',f10.6,',',f10.6,'], at t=',f10.6)") Ref_to_loc(mc(1),x_submiddle(mc(2))),Ref_to_loc(pv(1),x_submiddle(pv(2))), time
    !   write(*, fmt="( 'stencil = (', e12.6, 2x,e12.6, 2x ,e12.6, 2x,e12.6 ,')')") sol_step(voi_L(1))%val_subcells(voi_L(2),:), ug, ud , sol_step(voi_R(1))%val_subcells(voi_R(2),:) 
    !   write(*, fmt="( 'sol interface : ', e12.6,1x, e12.6 )") ug,ud
    !   write(*, fmt="( 'sol Riemann : ', e12.6 )") u_Riemann
    !   write(*, fmt="( 'alpha,beta = ',e12.6,2x, e12.6 )") alpha, beta
    !   print *,"theta_temp = ", theta_temp
    !   write(*, fmt="( 'theta = ', f10.6, ' gamma = ', f10.6)") theta, gamma_mp
    !   write(* ,fmt="( 'param = ', e12.6, ' DF = ', e12.6)") param, DF
    ! END IF
    ! theta =1._prec

  END FUNCTION

  SUBROUTINE extrema_detect 
    IMPLICIT NONE
    INTEGER :: ni, n_sub, n_var 
    INTEGER , DIMENSION(2) :: nL,nR
    REAL(prec) :: du_L,du_R, du,ddu, vL,vR

    LOGICAL :: face_L, face_R
    
    subcells_(:,:)%extrema = .FALSE. 

    DO ni=1,nb_cell;  DO n_var=1,nb_var
      sol_step(ni)%deriv(:,n_var)= MATMUL(Masse_inv, MATMUL(sol_step(ni)%base_poly(:,n_var), Rigid))
    END DO;           END DO

    IF(smooth_extrema == 1) THEN 

      DO ni=1,nb_cell; DO n_var=1,nb_var
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,1         , 'L')
        nR = Voisin_cell(ni,nb_subcell, 'R')

        du   = 1._prec/(cell_size(ni))    * (sol_step(ni   )%inter(2,n_var)-sol_step(ni   )%inter(1,n_var))
        du_L = 1._prec/(cell_size(nL(1))) * (sol_step(nL(1))%inter(2,n_var)-sol_step(nL(1))%inter(1,n_var))
        du_R = 1._prec/(cell_size(nR(1))) * (sol_step(nR(1))%inter(2,n_var)-sol_step(nR(1))%inter(1,n_var))

        ddu  = 2._prec/(cell_size(ni)**2)  *&
        &(eval_poly(-1._prec,ni,sol_step(ni)%deriv(:,n_var), LOC=LRef)  - eval_poly(1._prec,ni,sol_step(ni)%deriv(:,n_var), LOC=LRef))
        IF(abs(ddu) .LT. eps0) ddu = 0._prec

        vL = du - 0.5_prec*cell_size(ni)*ddu
        vR = du + 0.5_prec*cell_size(ni)*ddu

        IF((MIN(du,du_L)-eps0 .LT. vL) .AND. (MAX(du,du_L)+eps0 .GT. vL )) face_L = .TRUE.
        IF((MIN(du,du_R)-eps0 .LT. vR) .AND. (MAX(du,du_R)+eps0 .GT. vR )) face_R = .TRUE.

        IF( face_L .AND. face_R) THEN
          subcells_(ni,:)%extrema = .TRUE.
        END IF

        du_L = du; du = du_R
        ! IF(ni .GT. 1) THEN
        ! IF(subcells_(ni,1)%extrema .AND. subcells_(ni-1,1)%extrema) THEN
        !   print *,"================" ,ni,'----------',n_sub,"============================"
        !   print *,nL, nR
        !   write(*, fmt="( 'between : [',f10.6,',',f10.6,'], at t=',f10.6)") x_cell(ni),x_cell(ni+1), time
        !   write(*,fmt='("u = ",e12.6 )') eval_sol(0._prec,nvar =n_var,ni=ni,Loc=LRef)
        !   write(*,fmt='("du_L = ",e12.6,", du = ",e12.6,", du_R = ",e12.6)') du_L, du, du_R
        !   write(*,fmt='("vL   = ",e12.6,",ddu = ",e12.6,", vR   = ",e12.6)') vL, ddu,vR
        ! END IF; END IF


      END DO; END DO
    ELSE 

      DO ni=1,nb_cell;  DO n_sub =1,nb_subcell; DO n_var = 1,nb_var
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,n_sub, 'L') 
        nR = Voisin_cell(ni,n_sub, 'R') 

        du     = DOT_PRODUCT(sol_step(ni   )%deriv(:,n_var), Projection_VF(n_sub,:))/( cell_size(ni)*0.5_prec)
        du_L   = DOT_PRODUCT(sol_step(nL(1))%deriv(:,n_var), Projection_VF(nL(2),:))/( cell_size(ni)*0.5_prec)
        du_R   = DOT_PRODUCT(sol_step(nR(1))%deriv(:,n_var), Projection_VF(nR(2),:))/( cell_size(ni)*0.5_prec)

        ddu  = (4._prec/(subcell_size(n_sub) * cell_size(ni)**2))*&
        &(eval_poly(1._prec,n_sub,sol_step(ni)%deriv(:,n_var), LOC=LSub)  - eval_poly(-1._prec,n_sub,sol_step(ni)%deriv(:,n_var), LOC=LSub))

        vL = du - 0.25_prec*cell_size(ni)*subcell_size(n_sub)*ddu
        vR = du + 0.25_prec*cell_size(ni)*subcell_size(n_sub)*ddu

        IF((MIN(du,du_L)-eps0 .LT. vL) .AND. (MAX(du,du_L)+eps0 .GT. vL )) face_L = .TRUE.
        IF((MIN(du,du_R)-eps0 .LT. vR) .AND. (MAX(du,du_R)+eps0 .GT. vR )) face_R = .TRUE.

        IF( face_L .AND. face_R) THEN
          subcells_(ni,n_sub)%extrema = .TRUE.
        END IF

        ! IF(subcells_(ni,n_sub)%extrema) THEN
        !   print *,"================" ,ni,'----------',n_sub,"============================"
        !   print *,nL, nR
        !   write(*, fmt="( 'between : [',f10.6,',',f10.6,'], at t=',f10.6)") x_cell(ni),x_cell(ni+1), time
        !   write(*,fmt='("u = ",e12.6 )') eval_sol(0._prec,nvar =n_var,ni=ni,Loc=LRef)
        !   write(*,fmt='("du_L = ",e12.6,", du = ",e12.6,", du_R = ",e12.6)') du_L, du, du_R
        !   write(*,fmt='("vL   = ",e12.6,",ddu = ",e12.6,", vR   = ",e12.6)') vL, ddu,vR
        ! END IF

      END DO;           END DO;                  END DO

    END IF

  END SUBROUTINE extrema_detect


END MODULE mod_Monolith