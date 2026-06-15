MODULE mod_RKDG 
   use mod_Monolith
   use mod_Divers
  IMPLICIT NONE

CONTAINS

  
  SUBROUTINE flux_numerique
    IMPLICIT NONE
    INTEGER :: ni,ii,jj
    REAL(prec), DIMENSION(size_quad_nodes) :: flux_uh_val

    REAL(prec) :: ug,ud
    REAL(prec) :: x_s
    REAL(prec) :: fh_L, fh_R, f_FV
    REAL(prec) :: gamma_loc, theta_loc, DF
    REAL(prec) :: u_Riemann, param
    REAL(prec) :: alpha_loc, beta_loc
    INTEGER, DIMENSION(2) :: voi, voi_L, voi_R

    DO ni = 1,nb_cell+1
      IF (ni ==1) THEN 

        IF (TRIM(bdry_cond) == "period") THEN 
          ug = sol_step(nb_cell)  %inter(2)
          ud = sol_step(1)        %inter(1)
        ELSE IF(TRIM(bdry_cond) == "Sym") THEN
          ug = sol_step(1)        %inter(1)
          ud = sol_step(1)        %inter(1)
        ELSE 
          print *, "boundary condition non reconnue"
        END IF

      ELSE IF (ni ==nb_cell +1) THEN
        IF (TRIM(bdry_cond) == "period") THEN 
          ug = sol_step(nb_cell)  %inter(2)
          ud = sol_step(1)        %inter(1)
        ELSE IF(TRIM(bdry_cond) == "Sym") THEN
          ug = sol_step(nb_cell)  %inter(2)
          ud = sol_step(nb_cell)  %inter(2)
        ELSE 
          print *, "boundary condition non reconnue"
        END IF

      ELSE 
        ug = sol_step(ni-1)%inter(2)
        ud = sol_step(ni)  %inter(1)
      END IF

      g(ni) = (flux(ug) + flux(ud) - max_dflux*(ud-ug))  * 0.5_prec


    END DO

    DO ni = 1,nb_cell
      DO ii=1,size_quad_nodes
        flux_uh_val(ii) = flux(sol_step(ni)%val_nodes(ii))
      END DO

      CALL Projection_Pk(flux_uh,flux_h(ni)%base_poly,LOC= LLoc,ni =ni,fct_val=flux_uh_val)
    END DO
    
    IF(subcell_use) THEN
    DO ni = 1,nb_cell
      fh_L = eval_poly(-1._prec,ni,flux_h(ni)%base_poly, LOC= LRef)
      fh_R = eval_poly( 1._prec,ni,flux_h(ni)%base_poly, LOC= LRef)

      DO jj =1,nb_subcell+1
        x_s = x_subcell(jj)

        flux_h(ni)%val_subcells(jj) = eval_poly(x_s,ni, flux_h(ni)%base_poly,LOC= LRef) &
                                  & - C_m(jj)*(fh_L-g(ni  )) &
                                  & - C_p(jj)*(fh_R-g(ni+1))
        IF(monolithique) THEN
          voi_L = Voisin_Face(ni,jj,'L'); ug = sol_step(voi_L(1))%val_subcells(voi_L(2))
          voi_R = Voisin_Face(ni,jj,'R'); ud = sol_step(voi_R(1))%val_subcells(voi_R(2))
          f_FV = Flux_FV(ug,ud);      

          DF = flux_h(ni)%val_subcells(jj) - f_FV          
          theta_loc = theta(voi_L,voi_R)

          flux_h(ni)%val_subcells(jj) = f_FV + theta_loc * DF 
        END IF                  
      END DO

      
    END DO
    END IF

  END SUBROUTINE flux_numerique

  SUBROUTINE  Time_step
    IMPLICIT NONE
    
    INTEGER :: ni, tni

    INTEGER :: ii,jj
    REAL(prec) :: ti
    REAL(prec), DIMENSION(size_base) :: V_B, S_B, BB


    DO ni=1,nb_cell
      sol_step(ni)%base_poly  = sol(ni)%base_poly
      sol_step(ni)%val_nodes  = sol(ni)%val_nodes
      sol_step(ni)%inter      = sol(ni)%inter
    END DO
      

    DO tni =1,order_t
      
      CALL flux_numerique

      DO ni =1,nb_cell

        V_B = MATMUL(Rigid,flux_h(ni)%base_poly)
        S_B = -(g(ni+1)*sig_2 - g(ni)*sig_1)
        BB  = (V_B + S_B)
        
        L_step = MATMUL(Masse_inv, BB  )*(2._prec/(cell_size(ni))) 

        sol_step(ni)%base_poly = RK_alpha(tni,1) * sol(ni)%base_poly + RK_alpha(tni,2) * sol_step(ni)%base_poly &
                             &+  RK_beta(tni) *dt * L_step
        
        sol_step(ni)%val_nodes = 0._prec 
                  
        DO ii=1,size_quad_nodes          
          sol_step(ni)%val_nodes(ii)  = eval_step(x_quad(ii),ni, ii,LOC= LRef )
        END DO

        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol_step(ni)%inter(1)      = sol_step(ni)%val_nodes(1)
          sol_step(ni)%inter(2)      = sol_step(ni)%val_nodes(size_quad_nodes)
        ELSE 
          sol_step(ni)%inter(1)      = eval_step(x_cell(ni),ni,   LOC= LLoc)
          sol_step(ni)%inter(2)      = eval_step(x_cell(ni+1),ni, LOC= LLoc)
        END IF

      END DO
    END DO

    ! print *,"---------------------"

    DO ni=1,nb_cell
        sol(ni)%base_poly  = sol_step(ni)%base_poly

        DO ii=1,size_quad_nodes
          sol(ni)%val_nodes(ii)  = eval_sol(x_quad(ii),ni, ii, LOC=LRef )
        END DO

        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol(ni)%inter(1)      = sol(ni)%val_nodes(1)
          sol(ni)%inter(2)      = sol(ni)%val_nodes(size_quad_nodes)
        ELSE 
          sol(ni)%inter(1)      = eval_sol(x_cell(ni),ni,  LOC=LLoc)
          sol(ni)%inter(2)      = eval_sol(x_cell(ni+1),ni,LOC=LLoc)
        END IF
    END DO


  END SUBROUTINE Time_step

  SUBROUTINE Time_step_subcell
    IMPLICIT NONE
    INTEGER :: ni, tni

    INTEGER :: ii,jj,kk
    REAL(prec) :: ti
    REAL(prec), DIMENSION(nb_subcell) :: phi_m
    REAL(prec), DIMENSION(size_base)  :: phi
    REAL(prec), DIMENSION(nb_subcell,2) :: phi_val
    REAL(prec) :: BB

    ! print *, "subcells"

    DO ni=1,nb_cell
      sol_step(ni)%base_poly  = sol(ni)%base_poly
      sol_step(ni)%val_nodes  = sol(ni)%val_nodes
      sol_step(ni)%val_subcells=sol(ni)%val_subcells
      sol_step(ni)%inter      = sol(ni)%inter
    END DO

    DO ii=1,order_t
      CALL flux_numerique
      

      DO ni = 1,nb_cell
        DO jj =1,nb_subcell
          
          sol_step(ni)%val_subcells(jj)= RK_alpha(ii,1) * sol(ni)%val_subcells(jj) &
                                    & + RK_alpha(ii,2) * sol_step(ni)%val_subcells(jj) &
                                    & - RK_beta(ii) *(2._prec *dt/(cell_size(ni)* subcell_size(jj))) * (flux_h(ni)%val_subcells(jj+1)- flux_h(ni)%val_subcells(jj))                                   
        END DO
      END DO


      DO ni = 1,nb_cell
        DO jj=1,size_base
          sol_step(ni)%base_poly(jj) = DOT_PRODUCT(Projection_VF_inv_plus(jj,:), sol_step(ni)%val_subcells)
        END DO

        DO jj=1,size_quad_nodes
          sol_step(ni)%val_nodes(jj)  = eval_step(x_quad(jj),ni, jj, LOC= LRef)
        END DO
        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol_step(ni)%inter(1)      = sol_step(ni)%val_nodes(1)
          sol_step(ni)%inter(2)      = sol_step(ni)%val_nodes(size_quad_nodes)
        ELSE 
          sol_step(ni)%inter(1)      = eval_step(x_cell(ni),ni,  LOC= LLoc)
          sol_step(ni)%inter(2)      = eval_step(x_cell(ni+1),ni,LOC= LLoc)
        END IF
      END DO

    END DO

    
    DO ni=1,nb_cell
      sol(ni)%base_poly  = sol_step(ni)%base_poly
      sol(ni)%val_nodes  = sol_step(ni)%val_nodes
      sol(ni)%val_subcells=sol_step(ni)%val_subcells
      sol(ni)%inter      = sol_step(ni)%inter
    END DO




  END SUBROUTINE Time_step_subcell


  SUBROUTINE dt_calc
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(prec) :: max_u,min_u,u_step



    IF(TRIM(flux_name) == "advection") THEN
      max_dflux = abs(vit_adv)
    ELSE !IF(TRIM(flux_name) == "burgers") THEN
      !boucle pour max de u ?
      DO i=1,nb_cell
        DO j=1,size_quad_nodes
          IF(max_u .LT. sol(i)%val_nodes(j) ) THEN
            max_u = sol(i)%val_nodes(j)
          END IF

          IF(min_u .GT. sol(i)%val_nodes(j) ) THEN
            min_u = sol(i)%val_nodes(j)
          END IF
        END DO
      END DO

      max_dflux = gamma_calc(min_u,max_u)


    END IF
    dt_old = dt
    dt = min(CFL*dx/(REAL(2*order_x-1,prec)*max_dflux),tmax-time)

    dt = min(dt, 1.05_prec * dt_old)

    IF(order_x .GT. order_t) dt = dt**(order_x/order_t)

    ! IF(dt .LT. 10._prec**(-20)) THEN
    !   write(*, fmt ='("dt trop petit : dt =",e16.6, 1x,",max dflux =",e16.6 )') dt, max_dflux
    !   CALL Emergency_stop
    ! END IF

    
    IF(sol(1)%val_nodes(1) .NE. sol(1)%val_nodes(1) ) THEN
      write(*, fmt ='("NAN found emergency stop")')
      CALL Emergency_stop
    END IF


  END SUBROUTINE dt_calc

  SUBROUTINE writout
    IMPLICIT NONE

    INTEGER :: i,j
    REAL(prec) :: out1, out2,xi
    REAL(prec) :: err1 , err2, errLi
    INTEGER :: kk

    err1 = 0._prec; err2 =0._prec; errLi = 0._prec
    IF(modulo(n_time,500) == 0)  THEN
      write(*,fmt='("---------------",i7," ",f8.6," ",e16.6, "--------------")') n_time, time, dt
    END IF

    IF(time .GE.  REAL(n_imp,prec)*t_imp-eps0)  THEN
      write(*,fmt='("---------------",i7," ",f8.6," ",e16.6, "--------------")') n_time, time, dt
      DO i=1,nb_cell
        IF(subcell_use) THEN

        DO j=1,nb_subcell
          xi = Ref_to_loc(i,x_submiddle(j))
          out1 = sol(i)%val_subcells(j)

          IF(TRIM(flux_name) == "advection") THEN 
            out2 =Q_init(xi - time*vit_adv,0)
          ELSE 
            call pied_charact(xi,time,out2)
          END IF

          write(unit=numfile_sol,  fmt='(f10.6, f16.6, f16.6)') xi,out1, out2

          errLi = max(errLi , abs(out1-out2))
          err1 = err1 + abs(out1-out2)*w_quad(j)    *cell_size(i)/2
          err2 = err2 + ((out1-out2)*w_quad(j))**2  *cell_size(i)/2
        END DO

        ELSE
          
        DO j=1,size_base
          xi = Ref_to_loc(i,x_quad(j))
          out1 = sol(i)%val_nodes(j)

          IF(TRIM(flux_name) == "advection") THEN 
            out2 =Q_init(xi - time*vit_adv,0)
          ELSE 
            call pied_charact(xi,time,out2)
          END IF

          write(unit=numfile_sol,  fmt='(f10.6, f16.6, f16.6)') xi,out1, out2

          errLi = max(errLi , abs(out1-out2))
          err1 = err1 + abs(out1-out2)*w_quad(j)    *cell_size(i)/2
          err2 = err2 + ((out1-out2)*w_quad(j))**2  *cell_size(i)/2
        END DO

        END IF
      END DO

      err2 = sqrt(err2)
      
      write(*, fmt ='("err L1 = ", e20.12)')  err1
      write(*, fmt ='("err L2 = ", e20.12)')  err2
      write(*, fmt ='("err Li = ", e20.12)')  errLi

      err_L1 = max(err1, err_L1); err_L2 = max(err2, err_L2); err_Li = max(errLi, err_Li)

      n_imp = n_imp +1
      Time_stemp(n_imp) = time
      
      write(unit=numfile_sol, fmt='("------------------------")' ) 
    END IF
  END SUBROUTINE writout
  
  subroutine pied_charact(x,t,sol)

    REAL(prec), INTENT(IN) :: x,t
    REAL(prec), intent(out):: sol
    REAL(prec) :: xd,xf

    IF(TRIM(flux_name)=="advection") THEN
      xd = xL-abs(vit_adv)*t; xf = xR + abs(vit_adv)*t
    ELSE
      xd = xL-t; xf = xR + t
    END IF
    sol = Q_init(dicho(h,xd,xf),0)

    contains
    FUNCTION h(x_)
        use precis
        REAL(prec),INTENT(IN) :: x_
        REAL(prec)  :: h
        h = flux_d(Q_init(x_,0))*t + x_ -x
        return 
    END FUNCTION h

  END subroutine pied_charact

  FUNCTION Newton_search(x,t) result(q)
      implicit none
      REAL(prec), INTENT(IN)    :: x,t
      REAL(prec)                :: xk, err, q, epsi = 1e-20
      integer             :: n=0
      
    !     n = 0
    !     err = abs(flux(Q_init(xk))*t+ xk-x)
    !     ! print *, err, epsi
    !     xk = x

    !     do while(err>epsi .and. n<50)
    !         xk = xk -   (flux(Q_init(xk))*t + xk-x)/(flux_d(Q_init(xk))*Q_init_d(xk)*t +1)
    !         err =    abs(flux(Q_init(xk))*t + xk-x)

    !         n = n+1
    !     END do
        

    !     q = Q_init(xk)
    !     return

  END FUNCTION Newton_search

  SUBROUTINE Emergency_stop
    INTEGER :: i,j
    REAL(prec) :: out1, out2,xi
    REAL(prec) :: err1 , err2, errLi
    INTEGER :: kk

    write(*,fmt='("--------------",i5," ",f8.4," ",e16.6, "--------------")') n_time, time, dt
    DO i=1,nb_cell
      DO j=1,size_base
        xi = Ref_to_loc(i,x_quad(j))
        out1 = sol(i)%val_nodes(j)

        IF(TRIM(flux_name) == "advection") THEN 
          out2 =Q_init(xi - time*vit_adv,0)
        ELSE 
          call pied_charact(xi,time,out2)
        END IF

        write(unit=numfile_sol,  fmt='(f10.6, f16.6, f16.6)') xi,out1, out2

        errLi = max(errLi , abs(out1-out2))
        err1 = err1 + abs(out1-out2)*w_quad(j)    *cell_size(i)/2
        err2 = err2 + ((out1-out2)*w_quad(j))**2  *cell_size(i)/2
      END DO
    END DO

    err2 = sqrt(err2)
    
    write(*, fmt ='("err L1 = ", e20.12)')  err1
    write(*, fmt ='("err L2 = ", e20.12)')  err2
    write(*, fmt ='("err Li = ", e20.12)')  errLi

    err_L1 = max(err1, err_L1); err_L2 = max(err2, err_L2); err_Li = max(errLi, err_Li)

    n_imp = n_imp +1
    Time_stemp(n_imp) = time
    
    write(unit=numfile_sol, fmt='("------------------------")' ) 

    write(unit= numfile_data, fmt='("nt = "i5)') n_imp
    DO i=1,n_imp
        write(unit= numfile_data, fmt='("time ",i5," = ",f16.6)') i, Time_stemp(i)
    END DO

    close(unit=numfile_sol)
    close(unit=numfile_data)

    open(unit=numfile_conv,  file=nomfile_conv, form ='formatted', status ='old', position='append')
    write(unit=numfile_conv, fmt='("=====================")') 
    write(unit=numfile_conv, fmt='("for elements P",i1," and RK SSP of order ",i1)' ) size_base-1, order_t
    write(unit=numfile_conv, fmt='("for nx = "i5" we have error :")' ) nb_cell
    write(unit=numfile_conv, fmt='("err_L1 :" e20.12 )') err_L1
    write(unit=numfile_conv, fmt='("err_L2 :" e20.12 )') err_L2
    write(unit=numfile_conv, fmt='("err_Li :" e20.12 )') err_Li
    write(unit=numfile_conv, fmt='("=====================")') 
    close(unit=numfile_conv)

    print *,"closed"
    print *, counter1, counter2
    

    CALL DEALLOCATE_all


    print *, "EMERGENCY STOP"
  END SUBROUTINE Emergency_stop
  
END MODULE