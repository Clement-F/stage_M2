MODULE mod_Abgrall
    use mod_flux
  IMPLICIT NONE

CONTAINS



SUBROUTINE update_RF_entropie
    ! u_h^c est scalaire
    IMPLICIT NONE
    INTEGER :: ni,n_sub,ii,jj,kk
    REAL(prec), DIMENSION(nb_var) :: ug,ud,v_temp
    REAL(prec), DIMENSION(nb_subcell) :: vh, beta_entr, phi
    REAL(prec), DIMENSION(size_base)  :: v_tilde, sig
    REAL(prec), DIMENSION(nb_cell) :: alpha_entr
    REAL(prec), DIMENSION(nb_cell+1):: K
    REAL(prec), DIMENSION(nb_cell,nb_subcell) :: d
    REAL(prec) :: A,B,YY, vc,fh_L,fh_R
    REAL(prec) :: alpha_max

    beta_entr = 1._prec
    alpha_max = 0._prec
    d = 0._prec
    
    ! print *,"aaaaaaaaaaaaaaaaaaaaaaa"


    DO ni = 1,nb_cell+1
      IF (ni ==1) THEN 

        IF (TRIM(bdry_cond) == "period") THEN 
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(1)        %inter(1,:)
        ELSE IF(TRIM(bdry_cond) == "Sym") THEN
          ug = sol_step(1)        %inter(1,:)
          ud = sol_step(1)        %inter(1,:)
        ELSE IF(TRIM(bdry_cond) == "Solid" .AND. TRIM(flux_name)=="Euler") THEN
          ug = sol_step(1)        %inter(1,:); ug(2) = -ug(2)
          ud = sol_step(1)        %inter(1,:)

        ELSE 
          print *, "boundary condition non reconnue"
        END IF

      ELSE IF (ni ==nb_cell +1) THEN
        IF (TRIM(bdry_cond) == "period") THEN 
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(1)        %inter(1,:)
        ELSE IF(TRIM(bdry_cond) == "Sym") THEN
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(nb_cell)  %inter(2,:)
        ELSE IF(TRIM(bdry_cond) == "Solid" .AND. TRIM(flux_name)=="Euler") THEN
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(nb_cell)  %inter(2,:);  ud(2) = -ud(2)
        ELSE 
          print *, "boundary condition non reconnue"
        END IF

      ELSE 
        ug = sol_step(ni-1)%inter(2,:)
        ud = sol_step(ni)  %inter(1,:)
      END IF

      K(ni) = Flux_entrop_VF(ug,ud)
    END DO

    DO ni=1,nb_cell

      v_tilde = 0._prec; vh= 0._prec; vc = 0._prec
      DO ii=1,size_base; DO kk =1,nb_nodes
        YY =  x_quad(kk)
        v_temp = Var_entrop(eval_step(YY,ni,LOC=LRef))
        v_tilde(ii) = v_tilde(ii) + v_temp(1)*sig_quad(ii,kk)*w_quad(kk)
      END DO; END DO
      
      vh(:) = MATMUL(v_tilde,Projection_VF_inv)/(subcell_size)
      vc = sum(vh*beta_entr)/sum(beta_entr)

      A= 0._prec; B=0._prec
      d(ni,:) = (vh(:)-vc)*beta_entr(:)
      A = sum(vh(:)*(flux_h(ni)%flux_subcells(2:nb_subcell+1,1)-flux_h(ni)%flux_subcells(1:nb_subcell,1)))
      B = sum(vh*d(ni,:))

      alpha_entr(ni) = (A-(K(ni+1)-K(ni)))/(B+eps0)
    END DO


    DO ni = 1,nb_cell; DO ii = 1,nb_var
      fh_L = eval_poly(-1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)
      fh_R = eval_poly( 1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)

      DO jj =1,nb_subcell+1
        YY = x_subcell(jj)
        flux_h(ni)%flux_subcells(jj,ii) = eval_poly(YY,ni=ni, base_poly=flux_h(ni)%flux_DG(:,ii),LOC= LRef) &
                                        & - C_m(jj)*(fh_L-g(ni  ,ii)) &
                                        & - C_p(jj)*(fh_R-g(ni+1,ii)) &
                                        ! & - alpha_entr(ni)* DOT_PRODUCT(Adjacency(jj-1,:),d(ni,:))
                                        & - alpha_entr(ni)*sum(d(ni,1:jj-1))

      END DO
    END DO; END DO

END SUBROUTINE update_RF_entropie

SUBROUTINE calc_entrop(it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: it

    INTEGER :: ni,nk,i_max
    REAL(prec), DIMENSION(nb_var) :: ug,ud
    REAL(prec),allocatable, DIMENSION(:,:), save :: entr
    REAL(prec),allocatable, DIMENSION(:,:), save :: flux_k,diff_k
    REAL(prec) :: flux_t,max_diff

    IF (.not.allocated(entr)) THEN
        ALLOCATE(entr(nb_cell,2))
    END IF
    IF (.not.allocated(flux_k)) THEN
        ALLOCATE(flux_k(nb_cell+1,order_t))
        ALLOCATE(diff_k(nb_cell,order_t))
    END IF

    ! print *,it
    max_diff = 0._prec


    DO ni = 1,nb_cell+1
      IF (ni ==1) THEN 

        IF (TRIM(bdry_cond) == "period") THEN 
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(1)        %inter(1,:)
        ELSE IF(TRIM(bdry_cond) == "Sym") THEN
          ug = sol_step(1)        %inter(1,:)
          ud = sol_step(1)        %inter(1,:)
        ELSE IF(TRIM(bdry_cond) == "Solid" .AND. TRIM(flux_name)=="Euler") THEN
          ug = sol_step(1)        %inter(1,:); ug(2) = -ug(2)
          ud = sol_step(1)        %inter(1,:)

        ELSE 
          print *, "boundary condition non reconnue"
        END IF

      ELSE IF (ni ==nb_cell +1) THEN
        IF (TRIM(bdry_cond) == "period") THEN 
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(1)        %inter(1,:)
        ELSE IF(TRIM(bdry_cond) == "Sym") THEN
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(nb_cell)  %inter(2,:)
        ELSE IF(TRIM(bdry_cond) == "Solid" .AND. TRIM(flux_name)=="Euler") THEN
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(nb_cell)  %inter(2,:);  ud(2) = -ud(2)
        ELSE 
          print *, "boundary condition non reconnue"
        END IF

      ELSE 
        ug = sol_step(ni-1)%inter(2,:)
        ud = sol_step(ni)  %inter(1,:)
      END IF

      flux_k(ni,it) = Flux_entrop_VF(ug,ud)
    END DO

    DO ni=1,nb_cell
        diff_k(ni,it) = flux_k(ni+1,it)-flux_k(ni,it)
        if(abs(diff_k(ni,it)) .LT. eps0) diff_k(ni,it) = 0._prec
    END DO


    IF(it == 1) THEN
        entr = 0._prec
        DO ni =1,nb_cell;       DO nk = 1,nb_nodes
                entr(ni,1) = entr(ni,1)+ entropie_numerique(sol(ni)%val_quad(nk,:))*w_quad(nk)*cell_size(ni)/2._prec
            
        END DO;                 END DO        
    END IF

    IF(it == order_t) THEN
        DO ni =1,nb_cell;   DO nk = 1,nb_nodes
                entr(ni,2) = entr(ni,2)+ entropie_numerique(sol_step(ni)%val_quad(nk,:))*w_quad(nk)*cell_size(ni)/2._prec
        END DO;              END DO

        DO ni= 1,nb_cell
                
        IF(order_t==1) flux_t = -diff_k(ni,1)
        IF(order_t==2) flux_t = -(diff_k(ni,1)+diff_k(ni,2))/2._prec
        IF(order_t==3) flux_t = 2._prec/3._prec*(diff_k(ni,3)+0.5_prec*(diff_k(ni,1)+diff_k(ni,2)))

        ! print *,entr(ni,1),entr(ni,2)
        ! write(*,fmt='("entropie diff :",f10.6,"  entropie flux :",f10.6)') (entr(ni,2)-entr(ni,1))/dt , flux_t
        if(abs((entr(ni,2)-entr(ni,1))/dt - flux_t) .GT. max_diff) i_max = ni
        max_diff = max(max_diff,abs((entr(ni,2)-entr(ni,1))/dt - flux_t))

        print *,((entr(ni,2)-entr(ni,1))/dt - flux_t)
        END DO
        
    write(*,fmt='("max i    :",i5)')i_max
    write(*,fmt='("max diff :",f10.7)') max_diff

    END IF




END SUBROUTINE calc_entrop

END MODULE mod_Abgrall