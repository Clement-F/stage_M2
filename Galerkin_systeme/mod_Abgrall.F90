MODULE mod_Abgrall
    use mod_flux
  IMPLICIT NONE

CONTAINS



SUBROUTINE update_RF_entropie
    ! on suppose ici que v_m^c(x) = v(u_h^c(x))
    ! u_h^c est scalaire
    IMPLICIT NONE
    INTEGER :: ni,n_sub,kk
    INTEGER,DIMENSION(2) :: cL,cR
    REAL(prec), DIMENSION(nb_var) :: v_temp 
    REAL(prec), DIMENSION(nb_subcell) :: vh
    REAL(prec) :: A,B,YY

    vc_star = 0._prec

    DO ni=1,nb_cell

        DO n_sub =1,nb_subcell; DO kk =1,nb_nodes
            YY = Ref_to_loc(ni=ni, XX=Refsub_to_Ref(ZZ=x_quad(kk),n_sub =n_sub))
            v_temp = Var_entrop(eval_step(YY,ni,LOC=LLoc))
            vh(n_sub) = vh(n_sub) + v_temp(1)*w_quad(kk)/2._prec
        END DO; END DO

        DO n_sub =1,nb_subcell
            vc_star(ni) = vc_star(ni) + beta_entr(n_sub)*vh(n_sub)
        END DO
        vc_star(ni) = vc_star(ni)/sum(beta_entr)

        A = 0._prec; B=0._prec
        DO n_sub=1,nb_subcell
            A = A + vh(n_sub)*(flux_h(ni)%flux_subcells(n_sub+1,1)-flux_h(ni)%flux_subcells(n_sub,1))
            B = B + (vh(n_sub)-vc_star(ni))**2 *beta_entr(n_sub)
        END DO

        cL = Voisin_cell(ni,1,'L'); cR = Voisin_cell(ni,nb_subcell,'R')

        alpha_entr(ni) = A - entr_num_VF(sol_step(ni)%inter(2,:),sol_step(cR(1))%inter(1,:)) &
                    &  + entr_num_VF(sol_step(cL(1))%inter(2,:),sol_step(ni)%inter(1,:))
        alpha_entr(ni) = alpha_entr(ni) /(B+eps0)

        DO n_sub=1,nb_subcell
            flux_h(ni)%flux_subcells(n_sub,1) = flux_h(ni)%flux_subcells(n_sub,1) - alpha_entr(ni)/2._prec *(vh(n_sub)-vc_star(ni))*beta_entr(n_sub) 
        END DO
    END DO


END SUBROUTINE update_RF_entropie

END MODULE mod_Abgrall