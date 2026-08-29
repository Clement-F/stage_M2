MODULE mod_RKDG 
   use mod_Monolith
   use mod_SolIni
   use mod_Divers
   use mod_Abgrall
  IMPLICIT NONE

CONTAINS

  
  SUBROUTINE flux_numerique
    IMPLICIT NONE
    INTEGER :: ni,ii,jj
    REAL(prec), DIMENSION(nb_var) :: ug,ud, DF,theta_mp

    REAL(prec) :: x_s
    REAL(prec) :: fh_L, fh_R, gamma_mp
    INTEGER, DIMENSION(2) :: voi_L,voi_R

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

      g(ni,:) = (flux(ug) + flux(ud) - max_dflux*(ud-ug))  * 0.5_prec


    END DO
    ! print *,"proj"

    CALL Projection_Flux
    ! print *,"projected"
    
    IF(subcell_use) THEN
      
    DO ni = 1,nb_cell; DO ii = 1,nb_var
      flux_h(ni)%flux_subcells(:,:)= 0._prec
      ! F(uh(x)) .ne. Fh(x)
      fh_L = eval_poly(-1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)
      fh_R = eval_poly( 1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)

      DO jj =1,nb_subcell+1
        x_s = x_subcell(jj)
        flux_h(ni)%flux_subcells(jj,ii) = eval_poly(x_s,ni=ni, base_poly=flux_h(ni)%flux_DG(:,ii),LOC= LRef) &
                                        & - C_m(jj)*(fh_L-g(ni  ,ii)) &
                                        & - C_p(jj)*(fh_R-g(ni+1,ii))
      END DO
    END DO; END DO

    IF(entropie_rule == 4) CALL update_RF_entropie
    
    IF(monolithique) THEN 
    IF(smooth_extrema .GT. 0) CALL extrema_detect


    CALL Construct_thetaMesh

    DO ni=1,nb_cell; DO jj=1,nb_subcell+1                    

      voi_L = Voisin_Face(ni,jj,'L'); ug = sol_step(voi_L(1))%val_subcells(voi_L(2),:)
      voi_R = Voisin_Face(ni,jj,'R'); ud = sol_step(voi_R(1))%val_subcells(voi_R(2),:)

      IF(flux_num == 0) gamma_mp = max_dflux
      IF(flux_num == 1) gamma_mp = gamma_calc(ug,ud)

      flux_h(ni)%flux_vf(jj,:) = (flux(ug) + flux(ud) - gamma_mp*(ud-ug))  * 0.5_prec

      DF = ( flux_h(ni)%flux_subcells(jj,:)- flux_h(ni)%flux_vf(jj,:))

      IF(coeff_smooth == 0)theta_mp =     theta_(ni,jj,:)
      IF(coeff_smooth == 1)theta_mp = Min(theta_(ni,jj,:), 0.5_prec*(subcells_(voi_L(1),voi_L(2))%theta + subcells_(voi_R(1),voi_R(2))%theta)  )
      IF(coeff_smooth == 2)theta_mp = Min(theta_(ni,jj,:), Min(      subcells_(voi_L(1),voi_L(2))%theta,  subcells_(voi_R(1),voi_R(2))%theta)  )

      flux_h(ni)%flux_subcells(jj,:) = flux_h(ni)%flux_vf(jj,:) +theta_mp*DF
    END DO; END DO

    END IF
    END IF

  END SUBROUTINE flux_numerique

  SUBROUTINE  Time_step
    IMPLICIT NONE
    
    INTEGER :: ni, tni

    INTEGER :: ii,jj
    REAL(prec), DIMENSION(size_base,nb_var) :: V_B, S_B, BB,L_step


    DO ni=1,nb_cell
      sol_step(ni)%base_poly  = sol(ni)%base_poly
      sol_step(ni)%val_quad  = sol(ni)%val_quad
      sol_step(ni)%inter      = sol(ni)%inter
    END DO
      
    ! print *,"-----------------"

    DO tni =1,order_t
      
      CALL flux_numerique
      ! print *,"========================="
      ! TO CHANGE
      DO ni =1,nb_cell;     

        V_B = MATMUL(Rigid,flux_h(ni)%flux_DG(:,:))

        DO jj=1,nb_var;   S_B(:,jj) = -(g(ni+1,jj)*sig_2 - g(ni,jj)*sig_1);   END DO
        BB  = (V_B + S_B)
                
        L_step = MATMUL(Masse_inv, BB)*(2._prec/(cell_size(ni))) 

        sol_step(ni)%base_poly(:,:) = RK_alpha(tni,1) * sol(ni)%base_poly(:,:)& 
                                    &+ RK_alpha(tni,2) * sol_step(ni)%base_poly(:,:) &
                                    &+ RK_beta(tni)    * dt * L_step
                          
        DO ii=1,nb_nodes          
          sol_step(ni)%val_quad(ii,:)  = eval_step(x_quad(ii),ni=ni,kk= ii,LOC= LRef )
        END DO

        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol_step(ni)%inter(1,:)      = sol_step(ni)%val_quad(1,jj)
          sol_step(ni)%inter(2,:)      = sol_step(ni)%val_quad(nb_nodes,jj)
        ELSE 
          sol_step(ni)%inter(1,:)      = eval_step(x_cell(ni),  ni=ni, LOC= LLoc)
          sol_step(ni)%inter(2,:)      = eval_step(x_cell(ni+1),ni=ni, LOC= LLoc)
        END IF

      END DO
    END DO

    ! print *,"---------------------"

    DO ni=1,nb_cell
        sol(ni)%base_poly  = sol_step(ni)%base_poly

        DO ii=1,nb_nodes          
          sol(ni)%val_quad(ii,:)  = eval_sol(x_quad(ii),ni=ni,kquad= ii,LOC= LRef )
        END DO

        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol(ni)%inter(1,:)      = sol(ni)%val_quad(1,jj)
          sol(ni)%inter(2,:)      = sol(ni)%val_quad(nb_nodes,jj)
        ELSE 
          sol(ni)%inter(1,:)      = eval_sol(x_cell(ni),  ni=ni,LOC=LLoc)
          sol(ni)%inter(2,:)      = eval_sol(x_cell(ni+1),ni=ni,LOC=LLoc)
        END IF
    END DO


  END SUBROUTINE Time_step

  SUBROUTINE Time_step_subcell
    IMPLICIT NONE
    INTEGER :: ni

    INTEGER :: ii,jj
    REAL(prec), DIMENSION(nb_var) :: L

    outed_mesh = 0

    DO ni=1,nb_cell
      sol_step(ni)%base_poly  = sol(ni)%base_poly
      sol_step(ni)%val_quad  = sol(ni)%val_quad
      sol_step(ni)%val_subcells=sol(ni)%val_subcells
      sol_step(ni)%inter      = sol(ni)%inter
    END DO

    ! print *,"==========="
    DO ii=1,order_t
      ! print *,"------------------"
      CALL flux_numerique
      
      DO ni = 1,nb_cell;DO jj =1,nb_subcell
          
          L = (flux_h(ni)%flux_subcells(jj+1,:)- flux_h(ni)%flux_subcells(jj,:))

          sol_step(ni)%val_subcells(jj,:)=  RK_alpha(ii,1) * sol(ni)     %val_subcells(jj,:) &
                                        &  + RK_alpha(ii,2) * sol_step(ni)%val_subcells(jj,:) &
                                        &  - L*RK_beta(ii)  *(2._prec *dt/(cell_size(ni)* subcell_size(jj)))
                                                                      
      END DO; END DO

      DO ni = 1,nb_cell;

        sol_step(ni)%base_poly(:,:) = MATMUL(Projection_VF_inv(:,:), sol_step(ni)%val_subcells(:,:))

        DO jj=1,nb_nodes
          sol_step(ni)%val_quad(jj,:)  = eval_step(x_quad(jj),ni=ni, kk=jj, LOC= LRef)
        END DO
        
        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol_step(ni)%inter(1,:)      = sol_step(ni)%val_quad(1       ,:)
          sol_step(ni)%inter(2,:)      = sol_step(ni)%val_quad(nb_nodes,:)
        ELSE 
          sol_step(ni)%inter(1,:)      = eval_step(x_cell(ni),   ni=ni,LOC= LLoc)
          sol_step(ni)%inter(2,:)      = eval_step(x_cell(ni+1), ni=ni,LOC= LLoc)
        END IF

      END DO

      ! IF(max_check .GT. 0) CALL Error_check

    END DO

    DO ni=1,nb_cell
      sol(ni)%base_poly    = sol_step(ni)%base_poly
      sol(ni)%val_quad     = sol_step(ni)%val_quad
      sol(ni)%val_subcells = sol_step(ni)%val_subcells
      sol(ni)%inter        = sol_step(ni)%inter
    END DO



  END SUBROUTINE Time_step_subcell

  SUBROUTINE dt_calc
    IMPLICIT NONE
    INTEGER :: i,j
    INTEGER, DIMENSION(2) :: nxt_
    REAL(prec) :: gamma_temp, dt_loc, gamma_bf
    REAL(prec), DIMENSION(nb_var) :: u_,v_

    ! write (*,*) "dt calc"

    IF(TRIM(flux_name) == "advection" .AND. (.NOT. monolithique)) THEN
      max_dflux = abs(vit_adv)
    ELSE 
      max_dflux = 0._prec

      IF(flux_num ==1 ) THEN;

        IF(monolithique) THEN
          gamma_bf =eps0; dt_loc = 2._prec

          DO i=1,nb_cell;    DO j=1,nb_subcell
            ! print *,i,j
            nxt_ = Voisin_Face(i,j,'R')

            u_=sol(i)%val_subcells(j,:); v_ = sol(nxt_(1))%val_subcells(nxt_(2),:)

            IF(flux_name == "Buckley") THEN; max_dflux =2.4_prec
              gamma_temp = 2.4_prec
              dt_loc = min(CFL* cell_size(i)*subcell_size(j)/(4._prec*max_dflux), dt_loc)
            ELSE; 
            gamma_temp = max(gamma_calc(u_, v_),eps0)
            dt_loc = min(CFL* cell_size(i)*subcell_size(j)/(2._prec*(gamma_bf + gamma_temp)), dt_loc)
            END IF



            max_dflux = max(max_dflux, gamma_temp)
            gamma_bf = gamma_temp
          END DO;END DO 

        ELSE 

          DO i=1,nb_cell;
            nxt_ = Voisin_quad(i,nb_nodes,'R')

            u_=sol(i)%val_quad(nb_nodes,:); v_ = sol(nxt_(1))%val_quad(nxt_(2),:)
            gamma_temp = gamma_calc(u_, v_)
                      
            max_dflux = max(max_dflux, gamma_temp)
          END DO
        END IF

      ELSE IF(flux_num == 0) THEN; 
        ! print *, max_dflux
        IF(flux_name == "Buckley") THEN; max_dflux =2.4_prec
        ELSE
        DO i=1,nb_cell;
          nxt_ = Voisin_quad(i,nb_nodes,'R')
          ! IF(i == 1) THEN
          ! u_=sol(i)%val_quad(nb_nodes,:); v_ = sol(nxt_(1))%val_quad(nxt_(2),:)  
          ! u_(2) = -u_(2)        
          ! ELSE
          u_=sol(i)%val_quad(nb_nodes,:); v_ = sol(nxt_(1))%val_quad(nxt_(2),:)
          ! END IF
          gamma_temp = gamma_calc(u_, v_)
                    
          max_dflux = max(max_dflux, gamma_temp)
        END DO
        dt_loc = CFL* minval(cell_size(:))*minval(subcell_size(:))/(4._prec*max(max_dflux,eps0))  
      END IF     
      END IF

    END IF

    if(.not. monolithique) dt_loc =  2._prec
    
    dt = min(CFL*dx/max(REAL(2*order_x-1,prec)*max_dflux,eps0), tmax-time, dt_loc)
    IF(exact_time) dt = min(dt, Time_stemp(n_imp+1)-time)

    IF((order_x .GT. order_t).AND. convergence ) THEN
      IF(      monolithique) dt = min(tmax-time  +10._prec**(-10),CFL*(dx/(REAL(2*order_x-1,prec)*max_dflux))**(Real(order_x,prec)/Real(order_t,prec)), dt_loc**(Real(order_x,prec)/Real(order_t,prec)))
      IF(.not. monolithique) dt = min(tmax-time  +10._prec**(-10),CFL*(dx/(REAL(2*order_x-1,prec)*max_dflux))**(Real(order_x,prec)/Real(order_t,prec)))

      IF(exact_time) dt = min(dt, Time_stemp(n_imp+1)-time  +10._prec**(-10))
    END IF 

    IF(dt .LT. 10._prec**(-20) .AND. (time+dt .LT.Time_stemp(n_imp+1))) THEN
      write(*, fmt ='("dt trop petit : dt =",e16.6, 1x,",max dflux =",e16.6 )') dt, max_dflux
      CALL Emergency_stop
    END IF
    
    IF(ISNAN(sol(1)%val_quad(1,1))) THEN
      write(*, fmt ='("NAN found emergency stop")')
      CALL Emergency_stop
    END IF

    ! write(*,*) "dt calculated"

  END SUBROUTINE dt_calc

  SUBROUTINE writout(switch)
    IMPLICIT NONE

    INTEGER :: i,j
    REAL(prec), DIMENSION(nb_var) :: out, out_ex
    REAL(prec) :: xi
    REAL(prec) :: err1 , err2, errLi, entropy, pression_,p_max,p_min

    LOGICAL, optional :: switch
    LOGICAL :: force
    
    Character(len=63) :: save_format

    IF(present(switch)) THEN; force = switch
    ELSE; force = .FALSE.
    END IF
    err1 = 0._prec; err2 =0._prec; errLi = 0._prec;entropy=0._prec
    IF(modulo(n_time,500) == 0)  THEN
      write(*,fmt='("---------------",i7,1x,f10.6,1x,e12.6,2x,f6.2, "% --------------")') n_time, time, dt, (time*100._prec)/tmax 
    END IF


    save_format = "(f10.6"//Repeat(",f16.6",nb_var)//")"

    IF((time .GE.  Time_stemp(n_imp+1)-eps0) .OR. force )  THEN
      ! IF(mesh_out)CALL Out_The_Mesh(time)
      write(*,fmt='("---------------",i7,1x,f10.6,1x,e12.6,2x,f6.2, "% --------------")') n_time, time, dt, (time*100._prec)/tmax 
      DO i=1,nb_cell
        IF(subcell_use .AND.(.not. convergence) ) THEN

          save_format = "(f10.6"//Repeat(",f16.6",nb_var)//")"

          DO j=1,nb_subcell; 
            xi = Ref_to_loc(i,x_submiddle(j))
            out = sol(i)%val_subcells(j,:)

            save_format = "(f10.6"//Repeat(",f16.6",nb_var)//", f10.6)"
            write(unit=numfile_sol,    fmt=save_format, advance="no") xi,out !, Ref_to_loc(i,x_subcell(j))
            IF(TRIM(flux_name)== "Euler" ) THEN
            ! u_ = sol(i)%val_quad(j,:); pression_ = pression(u_)
            ! p_max = max(pression_, p_max); p_min = min(pression_,p_min)
            write(unit=numfile_sol, fmt= '(1x,e12.6)') pression(out)
            ELSE;
               write(unit=numfile_sol,   fmt= '(1x)')
               write(unit=numfile_solex, fmt= '(1x)')
            END IF

          END DO
          IF(error_calc) THEN;          
          DO j=1,nb_nodes
            xi = Ref_to_loc(i,x_quad(j))

            IF(TRIM(flux_name) == "advection") THEN       
              out_ex =Q_init(xi - time*vit_adv,0,nb_var)

            ELSE IF(nb_var ==1) THEN
              call pied_charact(xi,time,out_ex)

            ELSE IF(sol_ini_name == "isentropique") THEN
              call pied_charact(xi,time,out_ex)

            ELSE 
              out_ex = 0._prec

            END IF

            errLi = max(errLi , (abs(sol(i)%val_quad(j,1)-out_ex(1))))
            err1 = err1 + (abs(sol(i)%val_quad(j,1)-out_ex(1)))*w_quad(j)    *cell_size(i)/2
            err2 = err2 + (((sol(i)%val_quad(j,1)-out_ex(1))*w_quad(j))**2)  *cell_size(i)/2

            write(unit=numfile_solex,  fmt=save_format) xi,out_ex

            IF(TRIM(flux_name)== "Euler" ) THEN
            ! u_ = sol(i)%val_quad(j,:); 
            pression_ = pression(sol(i)%val_subcells(j,:) )
            p_max = max(pression_, p_max); p_min = min(pression_,p_min)
            write(unit=numfile_sol, fmt= '(2x,f10.6)') pression_
            ELSE;
               write(unit=numfile_sol,   fmt= '(1x)')
               write(unit=numfile_solex, fmt= '(1x)')
            END IF

          END DO 
          END IF

          IF(entropie_rule .GT. 0) THEN
          DO j=1,size_base
            out = sol(i)%val_quad(j,:)
            entropy = entropy + entropie_numerique(out)*cell_size(i) *w_quad(j)/2._prec
          END DO  
          END IF
        ELSE
          
          DO j=1,nb_nodes
            xi = Ref_to_loc(i,x_quad(j))

            out = sol(i)%val_quad(j,:)

            IF(TRIM(flux_name) == "advection") THEN     
              out_ex =Q_init(xi - time*vit_adv,0,nb_var)

            ELSE IF(nb_var ==1) THEN
              call pied_charact(xi,time,out_ex)

            ELSE IF(sol_ini_name == "isentropique") THEN
              call pied_charact(xi,time,out_ex)

            ELSE 
              out_ex = 0._prec

            END IF

            write(unit=numfile_sol,    fmt=save_format,  advance="no") xi,out
            IF(error_calc) write(unit=numfile_solex,  fmt=save_format,  advance="no") xi,out_ex


            errLi = max(errLi , (abs(out(1)-out_ex(1))))
            err1 = err1 + (abs(out(1)-out_ex(1)))*w_quad(j)    *cell_size(i)/2
            err2 = err2 + (((out(1)-out_ex(1))*w_quad(j))**2)  *cell_size(i)/2

            IF(entropie_rule .GT. 0) entropy = entropy + entropie_numerique(out(1))*cell_size(i)/2 *w_quad(j) 

            IF(TRIM(flux_name)== "Euler" ) THEN
            ! u_ = sol(i)%val_quad(j,:); 
            pression_ = pression(sol(i)%val_quad(j,:) )
            p_max = max(pression_, p_max); p_min = min(pression_,p_min)
            write(unit=numfile_sol, fmt= '(f10.6)') pression_
            ELSE;
               write(unit=numfile_sol,   fmt= '(1x)')
               write(unit=numfile_solex, fmt= '(1x)')
            END IF

          END DO

        END IF
      END DO



      err2 = sqrt(err2)
      IF(error_calc) THEN 
      write(*, fmt ='("err L1 = ", e20.12)')  err1
      write(*, fmt ='("err L2 = ", e20.12)')  err2
      write(*, fmt ='("err Li = ", e20.12)')  errLi
      END IF

      IF(monolithique) THEN
      write(*, fmt ='("avg theta = ")',advance = "no")  
      DO i=1,nb_var; write(*, fmt ='(f10.6)',advance ="no") Sum(theta_(:,:,i))/REAL((nb_cell)*(nb_subcell),prec);        END DO
      write(*, fmt ='(1x)')

      write(*, fmt ='("max theta = ")',advance = "no")  
      DO i=1,nb_var; write(*, fmt ='(f10.6)',advance ="no") maxval(theta_(:,:,i));        END DO
      write(*, fmt ='(1x)')

      write(*, fmt ='("max theta = ")',advance = "no")  
      DO i=1,nb_var; write(*, fmt ='(f10.6)',advance ="no") minval(theta_(:,:,i));        END DO
      write(*, fmt ='(1x)')

      IF(entropie_rule .GT. 0) THEN 
      write(*, fmt ='("entropy = ", f12.6)')  entropy
      END IF

      END IF
      err_L1 = max(err1, err_L1); err_L2 = max(err2, err_L2); err_Li = max(errLi, err_Li)


      IF((time .GE.  Time_stemp(n_imp+1)-eps0)) THEN
      n_imp = n_imp +1
      Time_stemp(n_imp) = time    
      END IF
      
      write(unit=numfile_sol  , fmt='("----------",f10.6,"--------------")' ) time
      write(unit=numfile_solex, fmt='("----------",f10.6,"--------------")' ) time
    END IF
  END SUBROUTINE writout
  
  SUBROUTINE pied_charact(x,t,sol)

    REAL(prec), INTENT(IN) :: x,t
    REAL(prec), DIMENSION(nb_var) , intent(out):: sol
    REAL(prec), DIMENSION(nb_var) :: Up, Um
    REAL(prec) :: w_p,w_m
    REAL(prec) :: xd,xf, x_m,x_p, c,u


    IF(TRIM(flux_name)=="advection") THEN
      xd = xL-abs(vit_adv)*t; xf = xR + abs(vit_adv)*t
    ELSE IF(bdry_cond =="period") THEN
      xd = xL-max_dflux*t; xf = xR + max_dflux*t

    ELSE IF(bdry_cond =="Sym" .OR. bdry_cond =="Solid" ) THEN
      xd = xL-eps0; xf = xR+eps0
    END IF

    IF(nb_var ==1)    sol = Q_init(dicho(h,xd,xf),0,nb_var)

    IF(sol_ini_name == "isentropique") THEN
      x_m = dicho(Burg_Wm,xd,xf);   x_p = dicho(Burg_Wp,xd,xf)
      Um  = Q_init(x_m,0,nb_var);   Up  =Q_init(x_p,0,nb_var)

      w_m = Um(2)/Um(1) - 2/(gamma_iso-1)*sqrt(gamma_iso) *Um(1) 
      w_p = Up(2)/Up(1) + 2/(gamma_iso-1)*sqrt(gamma_iso) *Up(1) 

      c = (gamma_iso-1._prec)/4._prec * (w_p-w_m) 
      u = (w_m+w_p)/2._prec

      sol(1) = c/sqrt(gamma_iso); sol(2) = sol(1)*u
      sol(3) = (sol(1)**3) /(gamma_iso-1._prec) + (sol(2)**2)/(2._prec*sol(1))

    END IF 


    contains
    FUNCTION h(x_)
        use precis
        REAL(prec),INTENT(IN) :: x_
        REAL(prec)  :: h
        REAL(prec),DIMENSION(nb_var) :: Um
        Um  =Q_init(x_,0,nb_var)
        h = flux_d(Um(1))*t + x_ -x
        return 
    END FUNCTION h

    FUNCTION Burg_Wm(x_)
        use precis
        REAL(prec),INTENT(IN) :: x_
        REAL(prec),DIMENSION(nb_var) :: Um
        REAL(prec)  :: Burg_Wm, W_m
        Um  =Q_init(x_,0,nb_var)
        W_m = Um(2)/Um(1)- 2/(gamma_iso-1)*sqrt(gamma_iso) * Um(1)
        Burg_Wm = W_m*t + x_ -x

        return 
    END FUNCTION Burg_Wm

    FUNCTION Burg_Wp(x_)
        use precis
        REAL(prec),INTENT(IN) :: x_
        REAL(prec),DIMENSION(nb_var) :: Up
        REAL(prec)  :: Burg_Wp, W_p
        Up  =Q_init(x_,0,nb_var)
        W_p = Up(2)/Up(1)+ 2/(gamma_iso-1)*sqrt(gamma_iso) * Up(1)
        Burg_Wp = W_p*t + x_ -x

        return 
    END FUNCTION Burg_Wp

  END SUBROUTINE pied_charact


  ! SUBROUTINE Out_The_Mesh(ti)
  !   IMPLICIT NONE
  !   REAL(prec), INTENT(IN) :: ti
  !   INTEGER :: ni,n_sub
  !   REAL(prec) :: xi, out1

  !   write(unit=numfile_meshout, fmt='("---------",f10.6,"---------------")' ) ti
  !   DO ni =1,nb_cell;    DO n_sub=1,nb_subcell
  !     xi = Ref_to_loc(ni,x_subcell(n_sub))
  !     out1 = subcells_(ni,n_sub)%theta
  !     write(unit=numfile_meshout,  fmt='(f10.6, f16.6, f16.6, 2x, l1)') xi,out1
  !   END DO; END DO

  ! END SUBROUTINE Out_The_Mesh

  SUBROUTINE Emergency_stop
    INTEGER :: i

    CALL writout(.TRUE.)

    open(unit=numfile_data,     file=nomfile_data,      form ='formatted', status ='old',  position='append')
    write(unit= numfile_data, fmt='("nt = ",i5)') n_imp
    DO i=1,n_imp+1
        write(unit= numfile_data, fmt='("time ",i5," = ",f16.6)') i, Time_stemp(i)
    END DO

    close(unit=numfile_sol);    close(unit=numfile_solex)

    
    close(unit=numfile_data)

    open(unit=numfile_conv,  file=nomfile_conv, form ='formatted', status ='old', position='append')
    write(unit=numfile_conv, fmt='("=====================")') 
    write(unit=numfile_conv, fmt='("for elements P",i1," and RK SSP of order ",i1)' ) size_base-1, order_t
    write(unit=numfile_conv, fmt='("for nx = ",i5," we have error :")' ) nb_cell
    write(unit=numfile_conv, fmt='("err_L1 :", e20.12 )') err_L1
    write(unit=numfile_conv, fmt='("err_L2 :", e20.12 )') err_L2
    write(unit=numfile_conv, fmt='("err_Li :", e20.12 )') err_Li
    write(unit=numfile_conv, fmt='("=====================")') 
    close(unit=numfile_conv)

    print *,"closed"
    print *, counter1, counter2
    

    CALL DEALLOCATE_all

    STOP "EMERGENCY STOP"
  END SUBROUTINE Emergency_stop


  SUBROUTINE Error_check
    IMPLICIT NONE
    INTEGER ::i,j,k
    REAL(prec) :: min_,max_
    
    ! print *,"error check"
    ! minmax rule
    DO i=1,nb_cell;       DO k=1,nb_var
      min_ = minmax_loc((/i,0/),'min',k)
      max_ = minmax_loc((/i,0/),'max',k)

      IF(maxval(sol_step(i)%val_subcells(:,k)) .GT. max_) THEN
        write(*,fmt='("erreur (max), var=",i1," en cellule ",i4," : max ",e12.6," , val",e12.6)') k,i,max_,maxval(sol_step(i)%val_subcells(:,k))
      END IF 

      IF(minval(sol_step(i)%val_subcells(:,k)) .LT. min_) THEN
        write(*,fmt='("erreur (min), var=",i1," en cellule ",i4," : max ",e12.6," , val",e12.6)') k,i,min_,maxval(sol_step(i)%val_subcells(:,k))
      END IF 

      DO j=1,nb_subcell
        ! minmax rule on subcell 
        min_ = minmax_loc((/i,j/),'min',k)
        max_ = minmax_loc((/i,j/),'max',k)

        IF(sol_step(i)%val_subcells(j,k) .GT. max_) THEN
          write(*,fmt='("erreur (max), var=",i1," en cellule ",i4," sous-cellule ",i4," : max ",e12.6," , val",e12.6)') k,i,j,max_,sol_step(i)%val_subcells(j,k)
        END IF 

        IF(sol_step(i)%val_subcells(j,k) .LT. min_) THEN
          write(*,fmt='("erreur (min), var=",i1," en cellule ",i4," sous-cellule ",i4," : min ",e12.6," , val",e12.6)') k,i,j,min_,sol_step(i)%val_subcells(j,k)
        END IF 

      IF(trim(flux_name)=="Euler") THEN 
        IF(minval(sol_step(i)%val_subcells(:,1)) .LT. 0._prec) THEN
          write(*,fmt='("erreur (densité), en cellule ",i4," :  val",e12.6)') i,sol_step(i)%val_subcells(j,1)
        END IF 

        IF(pression(sol_step(i)%val_subcells(j,:)) .LT. 0._prec) THEN
          write(*,fmt='("erreur (pression), en cellule ",i4," :  val",e12.6)') i,pression(sol_step(i)%val_subcells(j,:))
        END IF 

      END IF



      END DO

    END DO;  END DO

  END SUBROUTINE Error_check

  
END MODULE