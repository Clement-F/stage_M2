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
    REAL(prec), DIMENSION(nb_var) :: v_temp,ulr,ul,ur,url,diff
    REAL(prec), DIMENSION(nb_subcell) :: vh,d
    REAL(prec) :: A,B,YY

    vc_star = 0._prec

    ! print *,"============"

    DO ni=1,nb_cell

        ! print *,"------",ni,"--------"

        vh(:) = 0._prec
        DO n_sub =1,nb_subcell;    DO kk =nb_nodes,1,-1
            YY = Ref_to_loc(ni=ni, XX=Refsub_to_Ref(ZZ=x_quad(kk),n_sub =n_sub))
            v_temp = Var_entrop(eval_step(YY,ni,LOC=LLoc))
            vh(n_sub) = vh(n_sub) + v_temp(1)*w_quad(kk)/2._prec
        END DO;                            
        ! print*,vh(n_sub),sol_step(ni)%val_subcells(n_sub,1) 
        END DO

        vh(:) = sol_step(ni)%val_subcells(:,1) 

        vc_star(ni) =0._prec
        ! DO n_sub =1,nb_subcell
           
        ! END DO

        vc_star(ni) = sum(vh*beta_entr)/sum(beta_entr)
        d= (vh(:)-vc_star(ni))*beta_entr(:)

        ! write(*,fmt="('u = ',10(f6.3))")sol_step(ni)%val_subcells(:,1)
        ! write(*,fmt="('vh =',10(f6.3))")vh
        ! write(*,fmt="('v*= ',10(f6.3))")vc_star(ni)
        ! write(*,fmt="('d = ',10(f6.3))")d

        A = 0._prec; B=0._prec
        ! vh(:)=1._prec
        DO n_sub=1,nb_subcell
            A = A + vh(n_sub)*(flux_h(ni)%flux_subcells(n_sub+1,1)-flux_h(ni)%flux_subcells(n_sub,1))
            B = B + (vh(n_sub)-vc_star(ni))**2 *beta_entr(n_sub)
        END DO

        A = sum(vh*(flux_h(ni)%flux_subcells(2:nb_subcell+1,1)-flux_h(ni)%flux_subcells(1:nb_subcell,1)))
        ! B = sum(vh*d)
        ! print *,"A comp",

        cL = Voisin_cell(ni,1,'L'); cR = Voisin_cell(ni,nb_subcell,'R')
        ulr =sol_step(cL(1))%inter(2,:);ul =sol_step(ni)   %inter(1,:)
        ur  =sol_step(ni)   %inter(2,:);url=sol_step(cR(1))%inter(1,:)
        
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

        ! print *,"flux",flux_h(ni)%flux_subcells(:,1)
        ! print *,(flux_h(ni)%flux_subcells(2:nb_subcell+1,1)-flux_h(ni)%flux_subcells(1:nb_subcell,1))

        diff= Flux_entrop_VF(ur,url) - Flux_entrop_VF(ulr,ul)
        alpha_entr(ni) = A - diff(1)
        alpha_entr(ni) = alpha_entr(ni) /(B+eps0)

        ! print *,"l",Flux_entrop_VF(ulr,ul)
        ! print *,"r",Flux_entrop_VF(ur,url)
        ! print *,"diff",diff
        ! print *,"A",A
        ! print *,"B",B
        ! write(*,fmt="(f10.6)")alpha_entr(ni)

        flux_h(ni)%flux_Abgrall(1,1) = g(ni,1)
        flux_h(ni)%flux_Abgrall(nb_subcell+1,1) = g(ni+1,1)

        ! DO n_sub =2,nb_subcell
            ! print *,"modif",- alpha_entr(ni)*DOT_PRODUCT(Adjacency(n_sub-1,:),d)
            ! print *,"modif",  -alpha_entr(ni)*sum(d(1:n_sub-1))
            ! flux_h(ni)%flux_Abgrall(n_sub,1) = flux_h(ni)%flux_subcells(n_sub,1) - alpha_entr(ni)*DOT_PRODUCT(Adjacency(n_sub-1,:),d)
            ! flux_h(ni)%flux_Abgrall(n_sub,1) = flux_h(ni)%flux_subcells(n_sub,1) - alpha_entr(ni)*sum(d(1:n_sub-1))
        ! END DO


        flux_h(ni)%flux_Abgrall(2:nb_subcell,1) = flux_h(ni)%flux_subcells(2:nb_subcell,1)- alpha_entr(ni)*MATMUL(Adjacency,d)
        ! write(*,fmt="(10(f10.6))") alpha_entr(ni)*MATMUL(Adjacency,d)
        ! write(*,fmt="(10(f10.6))") flux_h(ni)%flux_Abgrall(:,1) - flux_h(ni)%flux_subcells(:,1) 
        ! write(*,fmt="(2(f10.6))") flux_h(ni)%flux_subcells(1,1)-g(ni,1), flux_h(ni)%flux_subcells(nb_subcell+1,1)-g(ni+1,1)
        ! flux_h(ni)%flux_Abgrall(nb_subcell,1) = flux_h(ni)%flux_subcells(nb_subcell,1)+ alpha_entr(ni)*d(nb_subcell)


        ! print *,"-----------"
        ! DO n_sub=1,nb_subcell+1
        !     ! print *,sol_step(ni)%val_subcells(n_sub,1)
        !     print *,flux_h(ni)%flux_subcells(n_sub,1), flux_h(ni)%flux_Abgrall(n_sub,1) 
        ! END DO

        flux_h(ni)%flux_subcells(:,1)= flux_h(ni)%flux_Abgrall(:,1) 



    END DO


END SUBROUTINE update_RF_entropie

END MODULE mod_Abgrall