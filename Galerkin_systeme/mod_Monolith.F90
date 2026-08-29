MODULE mod_Monolith
    use mod_flux
    use mod_Abgrall
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
    IF(max_rule == 2)   minmax_loc= min( sol_mc, &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2),nvar), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2),nvar))

    IF(max_rule == 1 .or. max_rule == 0)  minmax_loc = min_glob

    ELSE IF(TRIM(minmax) == "max") THEN
    IF(max_rule == 2)   minmax_loc= max( sol_mc, &
                                        &sol_step(voi_L(1))%val_subcells(voi_L(2),nvar), &
                                        &sol_step(voi_R(1))%val_subcells(voi_R(2),nvar))

    IF(max_rule == 1 .or. max_rule == 0)  minmax_loc = max_glob
    END IF


  END FUNCTION minmax_loc

  FUNCTION theta(mc,pv, DF)
    IMPLICIT NONE
    INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: DF
    REAL(prec), DIMENSION(nb_var) :: ug,ud, u_Riemann
    REAL(prec) :: gamma_mp

    REAL(prec), DIMENSION(nb_var) :: theta
    
    theta = 1._prec
    
    IF(minval(abs(DF)) < eps0) THEN; theta = 1._prec; return; END IF

    ug = sol_step(mc(1))%val_subcells(mc(2),:)
    ud = sol_step(pv(1))%val_subcells(pv(2),:)
    
    IF(flux_num == 0) gamma_mp = max_dflux
    IF(flux_num == 1) gamma_mp = gamma_calc(ug,ud)
    
    u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

    theta = min(1._prec,THETA_max(mc,pv, u_Riemann,gamma_mp,DF), THETA_ent(ug,ud,DF,gamma_mp))
    theta = min(1._prec,THETA_pos(u_Riemann,gamma_mp,DF, theta))
    theta = max(theta, 0._prec)
    IF(ISNAN(theta(1)))  STOP "theta nan"

  END FUNCTION

  FUNCTION THETA_pos(u_Riemann,gamma_mp,DF, theta_inc)
    IMPLICIT NONE
    REAL(prec), DIMENSION(nb_var) :: THETA_pos
    REAL(prec), INTENT(IN) :: gamma_mp
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: DF
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: u_Riemann, theta_inc
    REAL(prec), DIMENSION(nb_var) :: theta_temp
    REAL(prec) :: A,B,C,M, E,U,rho, B_e,B_u,B_r,A_r

    INTEGER :: pm

    THETA_pos = 1._prec
    theta_temp = theta_inc

    IF(TRIM(flux_name)=="Euler" .AND. positivity .GT. 0) THEN 
      ! positivité de rho
      IF(abs(DF(1)) .LT. eps0) THEN; 
        theta_temp(1) = 1._prec 
      ELSE ;
        theta_temp(1) = min(1._prec, abs(gamma_mp/DF(1))*(u_Riemann(1)-eps0))
      END IF

      ! positivité de E//P
      IF(positivity == 1) THEN 
        A = 1._prec/(gamma_mp**2)  *(0.5_prec*abs(DF(2))**2 - DF(1)*DF(3)*theta_temp(1))
        B = 1._prec/(gamma_mp)   *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - u_Riemann(3)*DF(1)*theta_temp(1))
        M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2

        ! IF(M .LT. -eps0) print *,"bug"

        theta_temp(2) = min(1._prec, max(M,eps0)/(max(abs(B),eps0)+ max(eps0,A)))

        THETA_pos = theta_temp(1)* theta_temp(2) 

      ELSE IF(positivity == 2) THEN  
        A = (1._prec/(gamma_mp))**2  *(0.5_prec*abs(DF(2))**2 - DF(1)*DF(3))
        B = 1._prec/(gamma_mp)        *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - u_Riemann(3)*DF(1))
        M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2

        ! IF(M .LT. -eps0) print *,"bug"

        theta_temp(2) = min(1._prec, max(M,eps0)/(max(abs(B),eps0)+ max(eps0,A)))

        THETA_pos = min(theta_temp(1), theta_temp(2)) 

      ELSE IF(positivity == 3) THEN 
        IF(DF(1) .LT. 0) theta_temp(1) = 0._prec
        THETA_pos(1) = theta_temp(1)
        ! print *, theta_temp(1)
        A = DF(1)*u_Riemann(3)/gamma_mp
        B = (DF(2)**2 /(2._prec* gamma_mp**2 ))- DF(2)*u_Riemann(2)/gamma_mp
        C = DF(3)*u_Riemann(1)/gamma_mp - theta_temp(1)*sign(1._prec,DF(1))*u_Riemann(1)*DF(3)/gamma_mp 
        
        M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2 !- max(A,eps0)
        ! print *,theta_temp(1)
        ! print *,"M", M
        A = abs(A); B= abs(B); C= abs(C)
        IF(M .LT. eps0) STOP "M"

        CALL Knapsack_greedy(THETA_pos,(/A,B,C /),M,(/THETA_pos(1),1._prec,1._prec/),3)

      ELSE IF(positivity == 4) THEN 
        M   = u_Riemann(1)* u_Riemann(3) - u_Riemann(2)**2 /(2._prec)
        B_r = u_Riemann(1)*DF(3)/gamma_mp
        IF(theta_temp(1) .GT. M/max(B_r,eps0))THEN 
          Theta_pos(2:3)= flat_positivity_criteria(theta_temp(1),u_Riemann,gamma_mp,DF)
        ELSE 
          Theta_pos(2:3)= flat_positivity_criteria(0._prec,u_Riemann,gamma_mp,DF)
          
          B_u = abs(DF(2)*u_Riemann(2))/(gamma_mp)
          B_e = abs(DF(3)*u_Riemann(1))/(gamma_mp)
          A   = DF(2)**2/(2._prec* gamma_mp**2)
          A_r = abs(DF(1)*DF(3)/gamma_mp**2)

          THETA_pos(1)  = min((M-B_u*Theta_pos(2)- B_e*THETA_pos(3)-A*Theta_pos(2)**2)/max(B_r+A_r*Theta_pos(3),eps0) ,1._prec)
        END IF       
      
      ELSE IF(positivity == 5) THEN
        ! heuristique

        Theta_pos(2:3)= flat_positivity_criteria(theta_temp(1),u_Riemann,gamma_mp,DF)

        rho = u_Riemann(1) - theta_temp(1)*abs(DF(1))/gamma_mp; 
        M   = rho* u_Riemann(3) - u_Riemann(2)**2 /(2._prec)
        if(M .LT. eps0) Theta_pos(1) = 0._prec
        
        return
      END IF
    END IF
  END FUNCTION THETA_pos

  FUNCTION flat_positivity_criteria(Theta_rho,u_Riemann,gamma_mp,DF) result(Theta)
    IMPLICIT NONE
    Real(prec), DIMENSION(2) :: Theta
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: u_Riemann,DF
    REAL(prec), INTENT(IN) :: Theta_rho, gamma_mp
    REAL(prec) :: rho,M,B_u,B_e,A
    INTEGER :: pm
    
        
    DO pm = 0,1

      if(pm == 0)  rho = u_Riemann(1) + Theta_rho*DF(1)/gamma_mp; 
      if(pm == 1)  rho = u_Riemann(1) - Theta_rho*DF(1)/gamma_mp; 

      IF(rho .LT. 0) print *,"rho<0",rho,DF(1),Theta_rho, u_Riemann(1)

      M   = rho* u_Riemann(3) - u_Riemann(2)**2 /(2._prec)

      IF(M .LT. eps0) THEN; Theta =0._prec;  return; END IF;

      if(pm == 0) B_u = DF(2)*u_Riemann(2)/(gamma_mp)
      if(pm == 1) B_u = -DF(2)*u_Riemann(2)/(gamma_mp)

      if(pm == 0) B_e = -(DF(3)/gamma_mp) * rho
      if(pm == 1) B_e =  (DF(3)/gamma_mp) * rho
        
      A   = DF(2)**2/(2._prec* gamma_mp**2)

      IF(M .LT. 0 .AND. pm ==1) return

      
      IF(A .LT. 0) print *,"A<0", A, gamma_mp, DF(2)
      
      IF(B_e .LT. eps0 .AND. B_u .LT. -A) THEN 
        Theta(1) = min(1._prec,Theta(1))
        Theta(2) = min(1._prec,Theta(2))
        return
      END IF

      IF(B_e .LT. eps0) THEN
        Theta(1) = min(1._prec, (-B_u + sqrt(B_u**2 +4._prec*A*(M-B_e)))/(2._prec*A) )
        Theta(2) = min(1._prec,Theta(2))
      END IF

      IF(B_u .LT. -A) THEN
        Theta(1) = min(1._prec,Theta(1))
        Theta(2) = min(1._prec,(M-B_u-A)/(max(B_e,eps0)), Theta(2))
      END IF

      IF(B_e .GT. eps0 .AND. B_u .GT. -A) THEN 
        ! IF(M .LT. 0) print *,"M",M
        IF(B_e .LT. 0) print *,"Be",B_e
        IF(B_u .LT. -A) print *,"Bu",B_u
        IF(A .LT. 0) print *,"A",A
        
        Theta(2) = min(1._prec,M*    B_e/max(max(B_e,eps0)**2 + max(B_u+A,eps0)**2,eps0),Theta(2))
        Theta(1) = min(1._prec, (-B_u + sqrt(B_u**2 +4._prec*A*(M-B_e*Theta(2))))/(2._prec*A) )
      END IF
    END DO



  END FUNCTION flat_positivity_criteria

  FUNCTION THETA_max(mc,pv, u_Riemann,gamma_mp,DF)
    IMPLICIT NONE
    REAL(prec) :: THETA_max
    INTEGER :: maxi
    INTEGER, DIMENSION(2), INTENT(IN) :: mc,pv
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: DF
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: u_Riemann
    REAL(prec), INTENT(IN) :: gamma_mp
    REAL(prec) :: alpha,beta, param
    LOGICAL :: extrema

    INTEGER :: ii

    THETA_max = 1._prec

    if(max_rule .LT. 3) maxi= 1
    if(max_rule .GE. 3) maxi= 3

    extrema = subcells_(mc(1),mc(2))%extrema .AND.  subcells_(pv(1),pv(2))%extrema 

    IF(.not. extrema .AND. max_rule .GT. 0) THEN
      DO ii = 1,maxi
        IF(max_rule==1) THEN; 
          alpha= 1.01_prec*(min_glob+eps0); 
          beta = 0.99_prec*(max_glob-eps0);

        ELSEIF(max_rule ==2) THEN
          IF((abs(DF(ii))) < eps0) THEN; THETA_max = 1._prec; return; END IF

          IF(DF(ii) .LT. -eps0) THEN
            beta = minmax_loc(mc,"max",nvar=ii)
            alpha= minmax_loc(pv,"min",nvar=ii)

          ELSE IF(DF(ii) .GT. eps0) THEN
            beta = minmax_loc(pv,"max",nvar=ii)
            alpha= minmax_loc(mc,"min",nvar=ii)
          END IF

        END IF

        param = min(beta - u_Riemann(ii), u_Riemann(ii)- alpha)
        
        THETA_max =max(min(1._prec, abs(gamma_mp/DF(1)) * param),0._prec); 
        IF(max_rule == 3) THEN 
          IF((abs(DF(3))) < eps0) return

          IF(DF(3) .LT. -eps0) THEN
            beta = minmax_loc(mc,"max",nvar=3)
            alpha= minmax_loc(pv,"min",nvar=3)

          ELSE IF(DF(3) .GT. eps0) THEN
            beta = minmax_loc(pv,"max",nvar=3)
            alpha= minmax_loc(mc,"min",nvar=3)
          END IF
          param = min(beta - u_Riemann(3), u_Riemann(3)- alpha)
          
          THETA_max =min(THETA_max, max(min(1._prec, abs(gamma_mp/DF(3)) * param),0._prec)); 
          ! print *,"aaa"

        END IF
      END DO
    END IF

  END FUNCTION THETA_max

  FUNCTION THETA_ent(ug,ud,DF,gamma_mp)
    IMPLICIT NONE
    REAL(prec) :: THETA_ent
    REAL(prec), DIMENSION(nb_var), INTENT(IN) :: ug,ud,DF
    REAL(prec), INTENT(IN) :: gamma_mp
    REAL(prec), DIMENSION(nb_var) :: vg, vd

    LOGICAL :: extrema

    THETA_ent = 1._prec

    IF(.not. extrema .AND. entropie_rule .GT. 0) THEN 

      IF(entropie_rule == 1) THEN 

        IF(DOT_PRODUCT(ud-ug, DF) .GT. eps0) THETA_ent =  &
          & min(1._prec, (max_dflux-gamma_mp)* minval((ud-ug)/(2._prec *DF)) )
        
      ELSEIF(entropie_rule==2) THEN 

        vg = Var_entrop(ug); vd = Var_entrop(ud)
        ! IF((DOT_PRODUCT(vd-vg, DF)) .GT. eps0) THETA_ent = (DOT_PRODUCT(vg,(flux(ug)-flux(ud)+gamma_mp*(ud-ug)))/2._prec  -Flux_entrop(ug) - &
        !                                                 & DOT_PRODUCT(vd,(flux(ug)-flux(ud)+gamma_mp*(ud-ug)))/2._prec  -Flux_entrop(ud)) / &
        !                                                 & (DOT_PRODUCT((vd-vg), DF))

        IF(DOT_PRODUCT(vd-vg,DF) .GT. eps0) THETA_ent = ((entrop_pot_flux(ud)-entrop_pot_flux(ug)) - DOT_PRODUCT(Flux_FV(ug,ud),vd-vg))/DOT_PRODUCT(vd-vg,DF)

          ! print *,"------------------------"
        ! write(*,fmt='(f10.6, f10.6, f10.6, f10.6)') vg, vd, vd-vg, DF
        ! ! write(*,fmt='(f10.6, f10.6, f10.6)')(entrop_pot_flux(ud)-entrop_pot_flux(ug)), -DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg)),(DOT_PRODUCT((vd-vg), DF))
        ! write(*,fmt='(f10.6)')(DOT_PRODUCT(vg,(flux(ug)-flux(ud)+gamma_mp/2._prec *(ud-ug))) -Flux_entrop(ug) - DOT_PRODUCT(vd,(flux(ug)-flux(ud)+gamma_mp/2._prec *(ud-ug))) -Flux_entrop(ud)) 
        ! ! write(*,fmt='(f10.6)') ( (entrop_pot_flux(ud)-entrop_pot_flux(ug)) -DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg)))
        ! ! write(*,fmt='(f10.6)')( (entrop_pot_flux(ud)-entrop_pot_flux(ug)) -DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg))) / (DOT_PRODUCT((vd-vg), DF))
        ! write(*,fmt='(f10.6)') THETA_ent

        ! print *, THETA_ent                                                
        THETA_ent = max(min(1._prec,THETA_ent),0._prec)
        ! print *,((entrop_pot_flux(ud)-entrop_pot_flux(ug))-DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg)))/(DOT_PRODUCT((vd-vg), DF))       
      END IF
    END IF
  END FUNCTION THETA_ent

  ! SUBROUTINE ENTROPI_CELL(ni)
  !   INTEGER, INTENT(IN) :: ni
  !   INTEGER, DIMENSION(2) :: voi_L,voi_R
  !   REAL(prec),DIMENSION(nb_var) :: vd,vg,ud,ug
  !   REAL(prec),DIMENSION(nb_subcell,nb_var) :: u_moy,v_moy
  !   REAL(prec),DIMENSION(size_base, nb_var):: poly_ent
  !   REAL(prec) :: D_c
  !   REAL(prec), DIMENSION(nb_subcell+1) :: Cell_c
  !   INTEGER :: ii,jj


  !   ! mettre vh_c le polynome entropie sur c
  !   u_moy = sol_step(ni)%val_subcells(:,:)
  !   DO ii =1,nb_subcell; v_moy(ii,:) = Var_entrop(u_moy(ii,:)); END DO
  !   poly_ent = MATMUL(Projection_VF_inv(:,:), v_moy)


  !   DO ii=1,nb_subcell
  !     voi_L = Voisin_Face(ni,jj,'L'); ug = sol_step(voi_L(1))%val_subcells(voi_L(2),:)
  !     voi_R = Voisin_Face(ni,jj,'R'); ud = sol_step(voi_R(1))%val_subcells(voi_R(2),:)

  !     vd = Var_entrop(ud); vg = Var_entrop(ug)
  !     D_c = entrop_pot_flux(vd)      

  !   END DO

  ! END SUBROUTINE ENTROPI_CELL

  SUBROUTINE extrema_detect 
    IMPLICIT NONE
    INTEGER :: ni, n_sub, n_var 
    INTEGER , DIMENSION(2) :: nL,nR
    REAL(prec) :: du_L,du_R, du,ddu, vL,vR
    LOGICAL :: face_L, face_R
    
    n_var =1
    subcells_(:,:)%extrema = .FALSE. 

    ! à optimiser!
    DO ni=1,nb_cell
      sol_step(ni)%deriv(:,n_var)= MATMUL(Masse_inv, MATMUL(sol_step(ni)%base_poly(:,n_var), Rigid))
    END DO;      

    IF(smooth_extrema == 1) THEN 

      DO ni=1,nb_cell; 
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,1         , 'L')
        nR = Voisin_cell(ni,nb_subcell, 'R')

        IF(ni == 1)  THEN 
        du_L = 1._prec/(cell_size(nL(1))) * (sol_step(nL(1))%inter(2,n_var)-sol_step(nL(1))%inter(1,n_var))
        du   = 1._prec/(cell_size(ni))    * (sol_step(ni   )%inter(2,n_var)-sol_step(ni   )%inter(1,n_var))
        END IF
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

      END DO;
    ELSE 

      DO ni=1,nb_cell;  DO n_sub =1,nb_subcell; 
        face_L = .FALSE.; face_R = .FALSE.
        nL = Voisin_cell(ni,n_sub, 'L') 
        nR = Voisin_cell(ni,n_sub, 'R') 

        IF(n_sub ==1) THEN
        du_L   = DOT_PRODUCT(sol_step(nL(1))%deriv(:,n_var), Projection_VF(nL(2),:))/( cell_size(ni)*0.5_prec)
        du     = DOT_PRODUCT(sol_step(ni   )%deriv(:,n_var), Projection_VF(n_sub,:))/( cell_size(ni)*0.5_prec)
        END IF
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

        du_L = du; du = du_R

      END DO;           END DO;                

    END IF

  END SUBROUTINE extrema_detect

  SUBROUTINE Construct_thetaMesh
    IMPLICIT NONE
    INTEGER, DIMENSION(2) :: voi_L,voi_R
    REAL(prec), DIMENSION(nb_var) :: theta_temp
    REAL(prec), DIMENSION(nb_var) :: ug,ud,u_Riemann
    REAL(prec), DIMENSION(nb_var) :: DF,pos, LMP, ent,out1
    REAL(prec) :: gamma_mp
    REAL(prec) :: xi
    Character(len=128) theta_outstring

    LOGICAL :: extrema
    INTEGER :: ni, jj

    extrema = .FALSE. 
    theta_(:,:,:) = 1._prec;  
    
    
    IF((mesh_out .AND. (time +dt .GE.  Time_stemp(n_imp+1)-eps0)).AND. outed_mesh ==0)  THEN
      write(unit=numfile_meshout, fmt='("---------",f10.6,"---------------")' ) time
      theta_outstring  = "(f7.4,1x"//Repeat(",f5.3,1x",nb_var)//",1x"//Repeat(",f5.3,1x",nb_var)//",1x"//Repeat(",f5.3,1x",nb_var)//",1x,l1 )"
    END IF

    DO ni=1,nb_cell; DO jj=1,nb_subcell
    IF(coeff_smooth ==0 .OR. coeff_smooth == 1)    subcells_(ni,jj)%theta = 1._prec 
    IF(coeff_smooth ==2)                           subcells_(ni,jj)%theta = 0._prec 
    ! subcells_(ni,jj)%theta(:) =1._prec
    END DO; END DO;

    DO ni=1,nb_cell; DO jj=1,nb_subcell+1
      ! print *,ni,jj
      pos = 1._prec; LMP = 1._prec; ent =1._prec

      theta_temp = 1._prec
      voi_L = Voisin_Face(ni,jj,'L'); ug = sol_step(voi_L(1))%val_subcells(voi_L(2),:)
      voi_R = Voisin_Face(ni,jj,'R'); ud = sol_step(voi_R(1))%val_subcells(voi_R(2),:)

      IF(TRIM(bdry_cond)=="Solid" .AND. TRIM(flux_name)=="Euler") THEN
        IF(ni == 1      .AND. jj ==1           )  ug(2) = -ug(2)
        IF(ni ==nb_cell .AND. jj ==nb_subcell+1)  ud(2) = -ud(2)
      END IF

      ! write(*,fmt="(3(f10.6),1x, 3(f10.6))") ug,ud
      IF(flux_num == 0) gamma_mp = max_dflux
      IF(flux_num == 1) gamma_mp = gamma_calc(ug,ud)

      flux_h(ni)%flux_vf(jj,:) = Flux_FV(ug,ud)
      
      DF = ( flux_h(ni)%flux_subcells(jj,:)- flux_h(ni)%flux_vf(jj,:))

      IF(ISNAN(DF(1)))   STOP "NAN DF"

      ! print *,Voi_L,Voi_R


      IF(.not. mesh_out) THEN
      theta_(ni,jj,:) = theta(voi_L,voi_R, DF)
      ELSE 
        
        IF(smooth_extrema .GT. 0) extrema = subcells_(voi_L(1),voi_L(2))%extrema .AND.  subcells_(voi_R(1),voi_R(2))%extrema 

        IF(maxval(abs(DF)) < eps0) THEN; theta_(ni,jj,:) = 1._prec;  ! check le besoin du calcul qui suit
        ELSE 
          ! print *,"--------------",ni,jj,"----------"
          ! print *, gamma_mp
          ! print *, ug,ud
          u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)


          ! positivité pour le pb  d'Euler
          ! check si l'interface se trouve entre deux sous-cellules qui captent un extrema
          LMP = THETA_max(voi_L,voi_R, u_Riemann,gamma_mp,DF)
          theta_(ni,jj,:) = min(theta_(ni,jj,:), LMP)

          pos = THETA_pos(u_Riemann,gamma_mp,DF,theta_(ni,jj,:))
          theta_(ni,jj,:) = min(theta_(ni,jj,:), pos)
          

          ! Entropie 
          ent = THETA_ent(ug,ud,DF,gamma_mp)
          IF(.not. extrema )theta_(ni,jj,:) = min(theta_(ni,jj,:), ent)
            

          theta_(ni,jj,:) = min(max(theta_(ni,jj,:),0._prec),1._prec)
          
          IF(ISNAN(theta_(ni,jj,1))) STOP 'theta nan'

          END IF
      END IF

      ! write(*,fmt="(f6.3)") theta_(ni,jj)

      IF(coeff_smooth == 1) THEN 
      subcells_(voi_L(1),voi_L(2))%theta = min(theta_(ni,jj,:), subcells_(voi_L(1),voi_L(2))%theta)
      subcells_(voi_R(1),voi_R(2))%theta = min(theta_(ni,jj,:), subcells_(voi_R(1),voi_R(2))%theta)
      ELSE IF(coeff_smooth == 2) THEN 
      subcells_(voi_L(1),voi_L(2))%theta = theta_(ni,jj,:)/2._prec + subcells_(voi_L(1),voi_L(2))%theta
      subcells_(voi_R(1),voi_R(2))%theta = theta_(ni,jj,:)/2._prec + subcells_(voi_R(1),voi_R(2))%theta
      ELSE IF(coeff_smooth == 0) THEN 
      ! print *,voi_L
      subcells_(voi_L(1),voi_L(2))%theta = theta_(ni,jj,:)
      END IF


      IF((mesh_out.AND. (time +dt .GE.  Time_stemp(n_imp+1)-eps0)).AND. outed_mesh ==0) THEN 
        xi = Ref_to_loc(voi_L(1),x_subcell(voi_L(2)))
        out1 = theta_(ni,jj,:)
        write(unit=numfile_meshout,  fmt=theta_outstring) xi,out1, pos,LMP,extrema
      END IF
      
    END DO; END DO

    IF(outed_mesh == 0)  outed_mesh = 1

  END SUBROUTINE Construct_thetaMesh

END MODULE mod_Monolith