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
    REAL(prec), DIMENSION(nb_var) :: theta_temp
    REAL(prec), DIMENSION(nb_var) :: ug,ud, u_Riemann
    REAL(prec) :: gamma_mp, param
    REAL(prec) :: alpha,beta

    REAL(prec) :: A,B,M
    INTEGER :: ii
    LOGICAL :: extrema

    REAL(prec) :: theta
    
    theta = 1._prec
    
    IF(minval(abs(DF)) < eps0) THEN; theta = 1._prec; return; END IF

    ug = sol_step(mc(1))%val_subcells(mc(2),:)
    ud = sol_step(pv(1))%val_subcells(pv(2),:)
    
    IF(flux_num == 0) gamma_mp = max_dflux
    IF(flux_num == 1) gamma_mp = gamma_calc(ug,ud)
    
    u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

    IF(TRIM(flux_name)=="Euler" .AND. positivity .GT. 0) THEN 
      ! positivité de rho
      theta_temp(1) = min(1._prec, abs(gamma_mp/DF(1))*(u_Riemann(1)-eps0))

      ! positivité de E//P
      IF(positivity == 1) THEN 
        A = 1/(gamma_mp**2)  *(0.5_prec*abs(DF(2))**2 - DF(1)*DF(3)*theta_temp(1))
        B = 1/(gamma_mp)   *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - u_Riemann(3)*DF(1)*theta_temp(1))
        M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2

        theta_temp(2) = min(1._prec, max(M,eps0)/(max(abs(B),eps0)+ max(eps0,A)))

        ! positivité des deux 
        theta = theta_temp(1)* theta_temp(2) 

      ELSE IF(positivity == 2) THEN  
        A = 1/(gamma_mp**2) *theta_temp(1)**2 *(0.5_prec*abs(DF(2))**2 - DF(1)*DF(3))
        B = 1/(gamma_mp) *theta_temp(1)    *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - u_Riemann(3)*DF(1))
        M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2

        theta_temp(2) = min(1._prec, max(M,eps0)/(max(abs(B),eps0)+ max(eps0,A)))

        ! positivité des deux 
        theta = min(theta_temp(1), theta_temp(2)) 
      END IF
    END IF

    theta_temp = 0._prec
    
    extrema = subcells_(mc(1),mc(2))%extrema .AND.  subcells_(pv(1),pv(2))%extrema 


    IF(.not. extrema .AND. max_rule .GT. 0) THEN
      DO ii=1,1

        IF(max_rule==1) THEN; 
          alpha= 1.01_prec*(min_glob+eps0); 
          beta = 0.99_prec*(max_glob-eps0);

        ELSEIF(max_rule ==2) THEN
          IF((abs(DF(ii))) < eps0) THEN; theta = 1._prec; return; END IF

          IF(DF(ii) .LT. -eps0) THEN
            beta = minmax_loc(mc,"max",nvar=ii)
            alpha= minmax_loc(pv,"min",nvar=ii)

          ELSE IF(DF(ii) .GT. eps0) THEN
            beta = minmax_loc(pv,"max",nvar=ii)
            alpha= minmax_loc(mc,"min",nvar=ii)
          END IF

        END IF

        param = min(beta - u_Riemann(ii), u_Riemann(ii)- alpha)
        
        theta_temp(ii) =max(min(1._prec, abs(gamma_mp/DF(ii)) * param),0._prec);
        ! END IF      
        
        theta = min(theta,theta_temp(ii))

      END DO
    END IF

    theta = max(theta,0._prec)
    
    IF(ISNAN(theta)) THEN; print *,"theta nan"; STOP; END IF

  END FUNCTION

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
    REAL(prec), DIMENSION(nb_var) :: ug,ud, vg,vd, u_Riemann
    REAL(prec), DIMENSION(nb_var) :: DF
    REAL(prec) :: gamma_mp, param
    REAL(prec) :: alpha,beta 
    REAL(prec) :: xi,out1, pos, LMP, ent

    REAL(prec) :: A,B,M
    LOGICAL :: extrema
    INTEGER :: ni, jj,ii

  extrema = .FALSE. 
  subcells_(:,:)%theta = 0._prec 
  theta_(:,:) = 1._prec;  
   
  
  IF((mesh_out .AND. (time +dt .GE.  Time_stemp(n_imp+1)-eps0)).AND. outed_mesh ==0)  THEN
    write(unit=numfile_meshout, fmt='("---------",f10.6,"---------------")' ) time
  END IF

  DO ni=1,nb_cell; DO jj=1,nb_subcell+1
    pos = 1._prec; LMP = 1._prec; ent =1._prec

    theta_temp = 1._prec
    ! print *,"------------------------------"
    ! print *, ni,jj
    IF(jj==1) THEN
    voi_L = Voisin_Face(ni,jj,'L'); ug = sol_step(voi_L(1))%val_subcells(voi_L(2),:)
    voi_R = Voisin_Face(ni,jj,'R'); ud = sol_step(voi_R(1))%val_subcells(voi_R(2),:)
    ELSE 
    voi_L = voi_R; ug = ud;
    voi_R = Voisin_Face(ni,jj,'R'); ud = sol_step(voi_R(1))%val_subcells(voi_R(2),:)
    END IF

    IF(TRIM(bdry_cond)=="Solid" .AND. TRIM(flux_name)=="Euler") THEN
      IF(ni == 1      .AND. jj ==1           )  ug(2) = -ug(2)
      IF(ni ==nb_cell .AND. jj ==nb_subcell+1)  ud(2) = -ud(2)
    END IF

    IF(flux_num == 0) gamma_mp = max_dflux
    IF(flux_num == 1) gamma_mp = gamma_calc(ug,ud)

    flux_h(ni)%flux_vf(jj,:) = (flux(ug) + flux(ud) - gamma_mp*(ud-ug))  * 0.5_prec
    
    DF = ( flux_h(ni)%flux_subcells(jj,:)- flux_h(ni)%flux_vf(jj,:))

    IF(ISNAN(DF(1)))  DF = 0._prec

    IF(.not. mesh_out) THEN
    theta_(ni,jj) = theta(voi_L,voi_R, DF)
    ELSE 
      IF(minval(abs(DF)) < eps0) THEN; theta_(ni,jj) = 1._prec;  ! check le besoin du calcul qui suit
      ELSE 
              
        u_Riemann = (ug+ud)/2._prec - (Flux(ud)-Flux(ug))/(2._prec*gamma_mp)

        ! positivité pour le pb  d'Euler
        IF(TRIM(flux_name)=="Euler" .AND. positivity .GT. 0) THEN 
          ! positivité de rho
          theta_temp(1) = min(1._prec, abs(gamma_mp/DF(1))*(u_Riemann(1)-eps0))

          ! positivité de E//P
          IF(positivity == 1) THEN 
            A = 1/(gamma_mp**2)  *(0.5_prec*abs(DF(2))**2 - DF(1)*DF(3)*theta_temp(1))
            B = 1/(gamma_mp)   *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - u_Riemann(3)*DF(1)*theta_temp(1))
            M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2

            theta_temp(2) = min(1._prec, max(M,eps0)/(max(abs(B),eps0)+ max(eps0,A)))

            ! positivité des deux 
            theta_(ni,jj) = theta_temp(1)* theta_temp(2) 

          ELSE IF(positivity == 2) THEN  
            A = 1/(gamma_mp**2) *theta_temp(1)**2 *(0.5_prec*abs(DF(2))**2 - DF(1)*DF(3))
            B = 1/(gamma_mp) *theta_temp(1)    *(u_Riemann(2)*DF(2) - u_Riemann(1)*DF(3) - u_Riemann(3)*DF(1))
            M = u_Riemann(1)*u_Riemann(3) - 0.5_prec * abs(u_Riemann(2))**2

            theta_temp(2) = min(1._prec, max(M,eps0)/(max(abs(B),eps0)+ max(eps0,A)))

            ! positivité des deux 
            theta_(ni,jj) = min(theta_temp(1), theta_temp(2)) 
          END IF
          pos = min(theta_temp(1), theta_temp(2))
        END IF

        theta_temp = 1._prec
        
        ! check si l'interface se trouve entre deux sous-cellules qui captent un extrema
        IF(smooth_extrema .GT. 0) extrema = subcells_(voi_L(1),voi_L(2))%extrema .AND.  subcells_(voi_R(1),voi_R(2))%extrema 

        ! Maximum Principle 
        IF(.not. extrema .AND. max_rule .GT. 0) THEN
          IF(max_rule==1) THEN; 
            alpha= 1.01_prec*(min_glob+eps0); 
            beta = 0.99_prec*(max_glob-eps0);

          ELSEIF(max_rule ==2) THEN
            IF((abs(DF(1))) < eps0) THEN; theta_(ni,jj) = 1._prec; return; END IF

            IF(DF(1) .LT. -eps0) THEN
              beta = minmax_loc(voi_L,"max",nvar=1)
              alpha= minmax_loc(voi_R,"min",nvar=1)

            ELSE IF(DF(1) .GT. eps0) THEN
              beta = minmax_loc(voi_R,"max",nvar=1)
              alpha= minmax_loc(voi_L,"min",nvar=1)
            END IF

          END IF

          param = min(beta - u_Riemann(1), u_Riemann(1)- alpha)
          
          theta_temp(1) =max(min(1._prec, abs(gamma_mp/DF(1)) * param),0._prec);
          ! END IF      
          
          theta_(ni,jj) = min(theta_(ni,jj),theta_temp(1))

          LMP = theta_temp(1)
        END IF

        theta_temp = 1._prec

        ! Entropie 
        IF(.not. extrema .AND. entropie_rule .GT. 0) THEN 

          IF(entropie_rule == 1) THEN 

            IF(DOT_PRODUCT(ud-ug, DF) .GT. eps0) theta_temp = min(1._prec, (max_dflux-gamma_mp)* minval((ud-ug)/(2._prec *DF)) )
            
            ent = theta_temp(1)

          ELSEIF(entropie_rule==2) THEN 

            vg = Var_entrop(ug); vd = Var_entrop(ud)
            ! IF((DOT_PRODUCT(vd-vg, DF)) .GT. eps0) theta_temp = (DOT_PRODUCT(vg,(flux(ug)-flux(ud)+gamma_mp/2._prec *(ud-ug))) -Flux_entrop(ug) - &
            !                                                 & DOT_PRODUCT(vd,(flux(ug)-flux(ud)+gamma_mp/2._prec *(ud-ug))) -Flux_entrop(ud)) / &
            !                                                 & (DOT_PRODUCT((vd-vg), DF))

            IF(DOT_PRODUCT(vd-vg, DF) .GT. eps0) theta_temp =& 
                & ( (entrop_pot_flux(ud)-entrop_pot_flux(ug))     &
                & -(flux_h(ni)%flux_vf(jj,:)*(vd-vg))) &
                & /((vd-vg) * DF)

              ! print *,"------------------------"
            ! write(*,fmt='(f10.6, f10.6, f10.6, f10.6)') vg, vd, vd-vg, DF
            ! write(*,fmt='(f10.6, f10.6, f10.6)')(entrop_pot_flux(ud)-entrop_pot_flux(ug)), -DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg)),(DOT_PRODUCT((vd-vg), DF))
            ! write(*,fmt='(f10.6)')(DOT_PRODUCT(vg,(flux(ug)-flux(ud)+gamma_mp/2._prec *(ud-ug))) -Flux_entrop(ug) - DOT_PRODUCT(vd,(flux(ug)-flux(ud)+gamma_mp/2._prec *(ud-ug))) -Flux_entrop(ud)) 
            ! ! write(*,fmt='(f10.6)') ( (entrop_pot_flux(ud)-entrop_pot_flux(ug)) -DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg)))
            ! ! write(*,fmt='(f10.6)')( (entrop_pot_flux(ud)-entrop_pot_flux(ug)) -DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg))) / (DOT_PRODUCT((vd-vg), DF))
            ! write(*,fmt='(f10.6)') theta_temp

            theta_temp = max(min(1._prec,theta_temp),0._prec)
            ent  =min(1._prec,theta_temp(1)) 
            ! print *,((entrop_pot_flux(ud)-entrop_pot_flux(ug))-DOT_PRODUCT(flux_h(ni)%flux_vf(jj,:),(vd-vg)))/(DOT_PRODUCT((vd-vg), DF))

          END IF
          theta_(ni,jj) = min(theta_(ni,jj),minval(theta_temp))
        END IF
          

        theta_(ni,jj) = max(theta_(ni,jj),0._prec)
        
        IF(ISNAN(theta_(ni,jj))) THEN; print *,"theta nan"; STOP; END IF

        END IF
    END IF

    IF(coeff_smooth == 1) THEN 
    subcells_(voi_L(1),voi_L(2))%theta = min(theta_(ni,jj), subcells_(voi_L(1),voi_L(2))%theta)
    subcells_(voi_R(1),voi_R(2))%theta = min(theta_(ni,jj), subcells_(voi_R(1),voi_R(2))%theta)
    ELSE IF(coeff_smooth == 2) THEN 
    subcells_(voi_L(1),voi_L(2))%theta = theta_(ni,jj)/2._prec + subcells_(voi_L(1),voi_L(2))%theta
    subcells_(voi_R(1),voi_R(2))%theta = theta_(ni,jj)/2._prec + subcells_(voi_R(1),voi_R(2))%theta
    ELSE IF(coeff_smooth == 0) THEN 
    subcells_(voi_L(1),voi_L(2))%theta = theta_(ni,jj)
    END IF


    IF((mesh_out.AND. (time +dt .GE.  Time_stemp(n_imp+1)-eps0)).AND. outed_mesh ==0) THEN 
      xi = Ref_to_loc(voi_L(1),x_subcell(voi_L(2)))
      out1 = theta_(ni,jj)
      write(unit=numfile_meshout,  fmt='(f10.6, f9.6, 1x, f9.6,  f9.6,1x,f9.6,1x, l1)') xi,out1, pos,LMP,ent,extrema
    END IF
    
  END DO; END DO

  IF(outed_mesh == 0)  outed_mesh = 1

  END SUBROUTINE Construct_thetaMesh

END MODULE mod_Monolith