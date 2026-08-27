MODULE mod_Abgrall
    use mod_flux
  IMPLICIT NONE

CONTAINS



SUBROUTINE update_RF_entropie
    ! u_h^c est scalaire
    IMPLICIT NONE
    INTEGER :: ni,n_sub,kk
    INTEGER,DIMENSION(2) :: cL,cR
    REAL(prec), DIMENSION(nb_var) :: v_temp,ulr,ul,ur,url,u_dum
    REAL(prec), DIMENSION(nb_subcell) :: vh,d
    REAL(prec) :: A,B,YY,diff

    vc_star = 0._prec
    alpha_entr = 0._prec

    print *,"============"

    DO ni=1,nb_cell

        print *,"------",ni,"--------"

        vh(:) = 0._prec; u_dum=0._prec
        DO n_sub =1,nb_subcell;             DO kk =1,nb_nodes
            YY = Ref_to_loc(ni=ni, XX=Refsub_to_Ref(ZZ=x_quad(kk),n_sub =n_sub))
            u_dum =eval_step(YY,ni,LOC=LLoc)
            v_temp = Var_entrop(u_dum(1))
            vh(n_sub) = vh(n_sub) + v_temp(1)*w_quad(kk)/2._prec
        END DO;                            
        ! print*,vh(n_sub),sol_step(ni)%val_subcells(n_sub,1) 
        END DO

        vc_star(ni) =0._prec
        DO n_sub =1,nb_subcell
            vc_star(ni) = vc_star(ni) + beta_entr(n_sub)*vh(n_sub)
        END DO
        vc_star(ni) = vc_star(ni)/sum(beta_entr)
        d= (vh(:)-vc_star(ni))*beta_entr(:)

        write(*,fmt="('u = ',10(f6.3))")sol_step(ni)%val_subcells(:,1)
        write(*,fmt="('vh =',10(f6.3))")vh
        write(*,fmt="('v*= ',10(f6.3))")vc_star(ni)
        write(*,fmt="('d = ',10(f6.3))")d

        A = 0._prec; B=0._prec

        A = sum(vh*(flux_h(ni)%flux_subcells(2:nb_subcell+1,1)-flux_h(ni)%flux_subcells(1:nb_subcell,1)))
        B = sum(vh*d)

        cL = Voisin_cell(ni,1,'L');         cR = Voisin_cell(ni,nb_subcell,'R')
        ur  =sol_step(ni)   %inter(2,:);    url=sol_step(cR(1))%inter(1,:)
        ulr =sol_step(cL(1))%inter(2,:);    ul =sol_step(ni)   %inter(1,:)
        
        ! print *,"-------------------"
        ! print *,ulr,ul
        ! print *,Flux_entrop_VF(ulr,ul) 
        ! print *,Flux_entrop(ulr),Flux_entrop(ul)
        ! print *,entropie_numerique(ulr),entropie_numerique(ul)
        ! print *,"---------"
        ! print *,ur,url
        ! print *,Flux_entrop_VF(ur,url) 
        ! print *,Flux_entrop(ur),Flux_entrop(url)
        ! print *,entropie_numerique(ur),entropie_numerique(url)
        ! print *,"-------------------"

        diff= Flux_entrop_VF(ur(1),url(1))- Flux_entrop_VF(ulr(1),ul(1))
        alpha_entr(ni) = (A - diff)/(B+eps0)

        print *,"l",Flux_entrop_VF(ulr(1),ul(1)), vit_adv*ulr**2/2._prec
        print *,"r",Flux_entrop_VF(ur(1),url(1)), vit_adv*ur**2/2._prec

        print *,"A",A
        print *,"B",B
        print *,"diff",diff
        print *,"alpha",alpha_entr(ni)

        ! flux_h(ni)%flux_Abgrall(1,1) = g(ni,1)
        ! flux_h(ni)%flux_Abgrall(nb_subcell+1,1) = g(ni+1,1)

        DO n_sub =2,nb_subcell
            ! print *,"modif",- alpha_entr(ni)*DOT_PRODUCT(Adjacency(n_sub-1,:),d)
            print *,"modif",  -alpha_entr(ni)*sum(d(1:n_sub-1))
            ! flux_h(ni)%flux_subcells(n_sub,1) = flux_h(ni)%flux_subcells(n_sub,1) - alpha_entr(ni)*DOT_PRODUCT(Adjacency(n_sub-1,:),d)
            flux_h(ni)%flux_subcells(n_sub,1) = flux_h(ni)%flux_subcells(n_sub,1) - alpha_entr(ni)*sum(d(1:n_sub-1))
        END DO


        ! print *,"-----------"
        ! DO n_sub=1,nb_subcell+1
        !     ! print *,sol_step(ni)%val_subcells(n_sub,1)
        !     print *,flux_h(ni)%flux_subcells(n_sub,1), flux_h(ni)%flux_Abgrall(n_sub,1) 
        ! END DO

        ! flux_h(ni)%flux_subcells(:,1)= flux_h(ni)%flux_Abgrall(:,1) 



    END DO


END SUBROUTINE update_RF_entropie

END MODULE mod_Abgrall