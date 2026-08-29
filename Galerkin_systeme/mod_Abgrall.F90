MODULE mod_Abgrall
    use mod_flux
  IMPLICIT NONE

CONTAINS



SUBROUTINE update_RF_entropie
    ! u_h^c est scalaire
    IMPLICIT NONE
    INTEGER :: ni,n_sub,kk
    INTEGER,DIMENSION(2) :: cL,cR
    REAL(prec), DIMENSION(nb_var) :: v_temp,ulr,ul,ur,url,diff
    REAL(prec), DIMENSION(nb_subcell) :: vh,d, beta_entr
    REAL(prec) :: A,B,YY, vc, alpha_entr

    vc = 0._prec
    beta_entr = 1._prec

    DO ni=1,nb_cell
        vh(:) = 0._prec
        DO n_sub =1,nb_subcell;    DO kk =nb_nodes,1,-1
            YY = Ref_to_loc(ni=ni, XX=Refsub_to_Ref(ZZ=x_quad(kk),n_sub =n_sub))
            v_temp = Var_entrop(eval_step(YY,ni,LOC=LLoc))
            vh(n_sub) = vh(n_sub) + v_temp(1)*w_quad(kk)/2._prec
        END DO;                     END DO

        vh(:) = sol_step(ni)%val_subcells(:,1) 

        vc =0._prec

        vc = sum(vh*beta_entr)/sum(beta_entr)
        d= (vh(:)-vc)*beta_entr(:)

        A = 0._prec; B=0._prec
        A = sum(vh*(flux_h(ni)%flux_subcells(2:nb_subcell+1,1)-flux_h(ni)%flux_subcells(1:nb_subcell,1)))
        B = sum(vh*d)

        cL = Voisin_cell(ni,1,'L'); cR = Voisin_cell(ni,nb_subcell,'R')
        ulr =sol_step(cL(1))%inter(2,:);ul =sol_step(ni)   %inter(1,:)
        ur  =sol_step(ni)   %inter(2,:);url=sol_step(cR(1))%inter(1,:)
        
        diff= Flux_entrop_VF(ur,url) - Flux_entrop_VF(ulr,ul)
        alpha_entr = A - diff(1)
        alpha_entr = alpha_entr /(B+eps0)

        flux_h(ni)%flux_Abgrall(1,1) = g(ni,1)
        flux_h(ni)%flux_Abgrall(nb_subcell+1,1) = g(ni+1,1)


        flux_h(ni)%flux_Abgrall(2:nb_subcell,1) = flux_h(ni)%flux_subcells(2:nb_subcell,1)- alpha_entr*MATMUL(Adjacency,d)

    END DO

    DO ni =1,nb_cell
        flux_h(ni)%flux_subcells(:,1)= flux_h(ni)%flux_Abgrall(:,1) 
    END DO


END SUBROUTINE update_RF_entropie

END MODULE mod_Abgrall