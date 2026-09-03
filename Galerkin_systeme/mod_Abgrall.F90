MODULE mod_Abgrall
    use mod_flux
  IMPLICIT NONE

CONTAINS



SUBROUTINE update_RF_entropie
    ! u_h^c est scalaire
    IMPLICIT NONE
    INTEGER :: ni,n_sub,kk,ii
    INTEGER,DIMENSION(2) :: cL,cR
    REAL(prec), DIMENSION(nb_var) :: v_temp,ulr,ul,ur,url,diff
    REAL(prec), DIMENSION(nb_subcell) :: vh,d, beta_entr
    REAL(prec) :: A,B,YY, vc, alpha_entr,x_sub,fh_L,fh_R
    REAL(prec) :: alpha_max,chg_max

    beta_entr = 1._prec
    alpha_max = 0._prec


    DO ni=1,nb_cell

        n_sub=0;kk=0;ii=0
        v_temp=0._prec;ulr=0._prec;ul=0._prec;ur=0._prec;url=0._prec;diff=0._prec;
        vh=0._prec;d=0._prec;
        A=0._prec;B=0._prec;YY=0._prec; vc=0._prec; alpha_entr=0._prec;x_sub=0._prec;fh_L=0._prec;fh_R=0._prec;
        alpha_max=0._prec;chg_max=0._prec;

        ! print *,"a"

        ! vh(:) = 0._prec
        ! DO n_sub =1,nb_subcell;    DO kk =1,nb_nodes
        !     YY =  Refsub_to_Ref(ZZ=x_quad(kk),n_sub =n_sub)
        !     v_temp = Var_entrop(eval_step(YY,ni,LOC=LRef))
        !     vh(n_sub) = vh(n_sub) + v_temp(1)*w_quad(kk)*cell_size(ni)
        ! END DO;                     END DO

        ! DO n_sub=1,nb_subcell; 
        ! print *,subcell_size(n_sub)
        ! vh(n_sub) = sol_step(ni)%val_subcells(n_sub,1)*(subcell_size(n_sub)*cell_size(ni))/2._prec
        ! print *,vh(n_sub),sol_step(ni)%val_subcells(n_sub,1)*subcell_size(n_sub)*cell_size(ni)/2._prec

        ! END DO

        ! vh(:) = sol_step(ni)%val_subcells(:,1)*(subcell_size(:)/2._prec)

        ! print *,cell_size(ni)*subcell_size(:)/2._prec
        vh(:) = sol_step(ni)%val_subcells(:,1)!*(cell_size(ni)*subcell_size(:)/2._prec)

        ! write(*,fmt="(10(f10.6))")vh(:)
        vc =0._prec

        ! DO n_sub=1,nb_subcell
        !     vc = vc + vh(n_sub)*beta_entr(n_sub)
        ! END DO
        vc = sum(vh*beta_entr)/sum(beta_entr)

        d= 0._prec
        ! DO n_sub=1,nb_subcell
        !     d(n_sub)=  vh(n_sub)*(vh(n_sub)-vc)*beta_entr(n_sub)
        ! END DO
        d = sum(vh*(vh-vc)*beta_entr)

        A = 0._prec; B=0._prec

        DO n_sub = 1,nb_subcell;
            A = A + vh(n_sub)*(flux_h(ni)%flux_subcells(n_sub+1,1)- flux_h(ni)%flux_subcells(n_sub,1))
            B = B + vh(n_sub)*d(n_sub)
        END DO

        cL = Voisin_cell(ni,1,'L'); cR = Voisin_cell(ni,nb_subcell,'R')
        ulr =sol_step(cL(1))%inter(2,:);ul =sol_step(ni)   %inter(1,:)
        ur  =sol_step(ni)   %inter(2,:);url=sol_step(cR(1))%inter(1,:)
        
        ! diff= Flux_entrop_VF(ur,url) - Flux_entrop_VF(ulr,ul)
        alpha_entr = A - (Flux_entrop_VF(ur,url) - Flux_entrop_VF(ulr,ul))
        alpha_entr = alpha_entr /(B+eps0)
    END DO

    DO ni=1,nb_cell
    DO ii =1,nb_var
        ! print *,ii
    fh_L = eval_poly(-1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)
    fh_R = eval_poly( 1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)

    DO n_sub =1,nb_subcell+1
        x_sub = x_subcell(n_sub)
        flux_h(ni)%flux_subcells(n_sub,ii) = eval_poly(x_sub,ni=ni, base_poly=flux_h(ni)%flux_DG(:,ii),LOC= LRef) &
                                        & - C_m(n_sub)*(fh_L-g(ni  ,ii)) &
                                        & - C_p(n_sub)*(fh_R-g(ni+1,ii)) &
                                        & - alpha_entr*sum(d(1:n_sub-1))
    
        alpha_max = max(alpha_max,abs(alpha_entr))
        chg_max   = max(chg_max  ,abs(alpha_entr*sum(d(1:n_sub-1))) )
    END DO
    END DO
    END DO

    print *,"alpha max :",alpha_max
    print *,"max change:",chg_max

END SUBROUTINE update_RF_entropie

SUBROUTINE calc_entrop(it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: it

    INTEGER :: ni,nk
    INTEGER, DIMENSION(2) :: cl,cr
    REAL(prec), DIMENSION(nb_var) :: ulr,ul,ur,url
    REAL(prec),allocatable, DIMENSION(:,:), save :: entr
    REAL(prec),allocatable, DIMENSION(:,:), save :: flux_k
    REAL(prec) :: flux_t,max_diff

    IF (.not.allocated(entr)) THEN
        ALLOCATE(entr(nb_cell,2))
    END IF
    IF (.not.allocated(flux_k)) THEN
        ALLOCATE(flux_k(nb_cell,0:order_t))
    END IF

    max_diff = 0._prec

    DO ni=1,nb_cell
        cL = Voisin_cell(ni,1,'L'); cR = Voisin_cell(ni,nb_subcell,'R')
        ulr =sol_step(cL(1))%inter(2,:);ul =sol_step(ni)   %inter(1,:)
        ur  =sol_step(ni)   %inter(2,:);url=sol_step(cR(1))%inter(1,:)
        
        flux_k(ni,it-1) = Flux_entrop_VF(ur,url) - Flux_entrop_VF(ulr,ul)
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
            
        IF(order_t==1) flux_t = -flux_k(ni,0)
        IF(order_t==2) flux_t = -flux_k(ni,0)-flux_k(ni,1)
        IF(order_t==3) flux_t = -2._prec/3._prec * flux_k(ni,2) - 1._prec/6._prec *(flux_k(ni,1)+flux_k(ni,0))

        ! print *,entr(ni,1),entr(ni,2)
        ! write(*,fmt='("entropie diff :",f10.6,"  entropie flux :",f10.6)') (entr(ni,2)-entr(ni,1))/dt , flux_t
        max_diff = max(max_diff,abs((entr(ni,2)-entr(ni,1))/dt - flux_t))
        END DO
        
        write(*,fmt='("max diff :",f10.7)') max_diff
    END IF




END SUBROUTINE calc_entrop

END MODULE mod_Abgrall