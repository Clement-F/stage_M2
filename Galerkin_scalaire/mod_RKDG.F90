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
    REAL(prec) :: theta_loc, DF
    INTEGER, DIMENSION(2) :: voi_L, voi_R

    LOGICAL :: extrema

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

    IF(smooth_extrema) CALL extrema_detect

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
          extrema = (subcells_(voi_L(1),voi_L(2))%extrema) .OR. (subcells_(voi_R(1),voi_R(2))%extrema)

          IF(.not. extrema) THEN
          f_FV = Flux_FV(ug,ud);      
          DF = flux_h(ni)%val_subcells(jj) - f_FV; IF(abs(DF) .LT. eps0) DF = 0._prec   

          theta_loc = theta(voi_L,voi_R, DF)

          flux_h(ni)%val_subcells(jj) = f_FV + theta_loc * DF 
          END IF

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
      
    ! print *,"-----------------"

    DO tni =1,order_t
      
      CALL flux_numerique
      ! print *,"========================="

      DO ni =1,nb_cell

        V_B = MATMUL(Rigid,flux_h(ni)%base_poly)
        S_B = -(g(ni+1)*sig_2 - g(ni)*sig_1)
        BB  = (V_B + S_B)
        
        L_step = MATMUL(Masse_inv, BB  )*(2._prec/(cell_size(ni))) 

        sol_step(ni)%base_poly = RK_alpha(tni,1) * sol(ni)%base_poly + RK_alpha(tni,2) * sol_step(ni)%base_poly &
                             &+  RK_beta(tni) *dt * L_step
        
        sol_step(ni)%val_nodes = eps0 
                  
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
    REAL(prec), DIMENSION(nb_cell,nb_subcell) :: max_loc, min_loc
    REAL(prec) :: L

    ! print *,"---------------------"

    DO ni=1,nb_cell
      sol_step(ni)%base_poly  = sol(ni)%base_poly
      sol_step(ni)%val_nodes  = sol(ni)%val_nodes
      sol_step(ni)%val_subcells=sol(ni)%val_subcells
      sol_step(ni)%inter      = sol(ni)%inter
    END DO

    DO ii=1,order_t
      CALL flux_numerique
      
      IF(max_check == 1) THEN
        DO ni=1,nb_cell
          DO jj=1,nb_subcell
            max_loc(ni,jj) = minmax_loc((/ni,jj/),"max")
            min_loc(ni,jj) = minmax_loc((/ni,jj/),"min")
          END DO 
        END DO
      END IF

      DO ni = 1,nb_cell
        DO jj =1,nb_subcell
          
          ! print *,flux_h(ni)%val_subcells(jj+1), flux_h(ni)%val_subcells(jj)

          L = (flux_h(ni)%val_subcells(jj+1)- flux_h(ni)%val_subcells(jj))

          sol_step(ni)%val_subcells(jj)= RK_alpha(ii,1) * sol(ni)%val_subcells(jj) &
                                    & + RK_alpha(ii,2) * sol_step(ni)%val_subcells(jj) &
                                    & - L*RK_beta(ii) *(2._prec *dt/(cell_size(ni)* subcell_size(jj)))
                                                                      
        END DO
      END DO

      IF(max_check .ne. 0) THEN
        ! print *,"--------------------------"
        DO ni=1,nb_cell
          DO jj=1,nb_subcell
            IF(max_check == 1) THEN

              IF( (.not.((ni == 1) .and. (jj == 1))) .and. (.not.((ni == nb_cell) .and. (jj == nb_subcell) ))) THEN
                IF(sol_step(ni)%val_subcells(jj) .GT. max_loc(ni,jj)+4*eps0 ) write(*,fmt="('problem max rule at :(',i2,1x,i2,'), max : ',e12.6,' ,val : ',e12.6,' diff : ',e12.6 )") ni,jj,max_loc(ni,jj),sol_step(ni)%val_subcells(jj), (sol_step(ni)%val_subcells(jj) - max_loc(ni,jj))   
                IF(sol_step(ni)%val_subcells(jj) .LT. min_loc(ni,jj)-4*eps0 ) write(*,fmt="('problem min rule at :(',i2,1x,i2,'), min : ',e12.6,' ,val : ',e12.6,' diff : ',e12.6 )") ni,jj,min_loc(ni,jj),sol_step(ni)%val_subcells(jj), (sol_step(ni)%val_subcells(jj) - min_loc(ni,jj))     
              END IF

            ELSE IF(max_check == 2) THEN
              
              IF( (.not.((ni == 1) .and. (jj == 1))) .and. (.not.((ni == nb_cell) .and. (jj == nb_subcell) ))) THEN
                IF(sol_step(ni)%val_subcells(jj) .GT. max_glob+4*eps0 ) write(*,fmt="('problem max rule at :(',i2,1x,i2,'), max : ',e12.6,' ,val : ',e12.6,' diff : ',e12.6 )") ni,jj,max_glob,sol_step(ni)%val_subcells(jj), (sol_step(ni)%val_subcells(jj) - max_glob)    
                IF(sol_step(ni)%val_subcells(jj) .LT. min_glob-4*eps0 ) write(*,fmt="('problem min rule at :(',i2,1x,i2,'), min : ',e12.6,' ,val : ',e12.6,' diff : ',e12.6 )") ni,jj,min_glob,sol_step(ni)%val_subcells(jj), (sol_step(ni)%val_subcells(jj) - min_glob)      
              END IF
              
            END IF
          END DO
        END DO
      END IF



      DO ni = 1,nb_cell
        DO jj=1,size_base
          sol_step(ni)%base_poly(jj) = DOT_PRODUCT(Projection_VF_inv(jj,:), sol_step(ni)%val_subcells)
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

    ! IF(smooth_extrema) THEN
    !   DO ni =1,nb_cell
        
    !   END DO
    ! END IF
    
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
    INTEGER, DIMENSION(2) :: next_subcell
    REAL(prec) :: max_u,min_u,  gamma_temp



    IF(TRIM(flux_name) == "advection") THEN
      max_dflux = abs(vit_adv)
    ELSE 
      max_dflux =eps0
      if(subcell_use) THEN
        DO i=1,nb_cell
          DO j=1,nb_subcell
            next_subcell = subcells_(i,j)%L
            gamma_temp = gamma_calc(sol(i)%val_subcells(j), sol(next_subcell(1))%val_subcells(next_subcell(2)))
            
            max_dflux = max(max_dflux, gamma_temp)
          END DO
        END DO
      ELSE
        DO i=1,nb_cell        
          DO j=1,nb_subcell
            IF(max_u .LT. sol(i)%val_subcells(j) ) THEN
              max_u = sol(i)%val_subcells(j)
            END IF

            IF(min_u .GT. sol(i)%val_subcells(j) ) THEN
              min_u = sol(i)%val_subcells(j)
            END IF
          END DO
        END DO

        max_dflux = gamma_calc(min_u,max_u) 
      END IF


    END IF
    dt_old = dt
    dt = min(CFL*dx/(REAL(2*order_x-1,prec)*max_dflux),tmax-time)

    IF(order_x .GT. order_t) THEN
    dt = min(   (CFL*dx/(REAL(2*order_x-1,prec)*max_dflux)) **(Real(order_x,prec) * 1._prec/Real(order_t,prec)) ,tmax-time) 
    END IF 

    dt = min(dt, 1.05_prec * dt_old)


    IF(dt .LT. 10._prec**(-20)) THEN
      write(*, fmt ='("dt trop petit : dt =",e16.6, 1x,",max dflux =",e16.6 )') dt, max_dflux
      CALL Emergency_stop
    END IF

    
    IF(ISNAN(sol(1)%val_nodes(1))) THEN
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
      write(*,fmt='("---------------",i7,2x,f10.6,2x,e16.6, "--------------")') n_time, time, dt
    END IF

    IF(time .GE.  REAL(n_imp,prec)*t_imp-eps0)  THEN
      write(*,fmt='("---------------",i7,2x,f10.6,2x,e16.6, "--------------")') n_time, time, dt
      DO i=1,nb_cell
        IF(subcell_use .and.(.not. error_calc)) THEN

          DO j=1,nb_subcell
            xi = Ref_to_loc(i,x_submiddle(j))
            out1 = sol(i)%val_subcells(j)
            write(unit=numfile_sol,  fmt='(f10.6, f16.6, f16.6, 2x, l1)') xi,out1, 0._prec, subcells_(i,j)%extrema
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
    
    write(unit=numfile_sol, fmt='("------------------------")' ) 

    write(unit= numfile_data, fmt='("nt = ",i5)') n_imp
    DO i=1,n_imp
        write(unit= numfile_data, fmt='("time ",i5," = ",f16.6)') i, Time_stemp(i)
    END DO

    close(unit=numfile_sol)
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


    print *, "EMERGENCY STOP"
  END SUBROUTINE Emergency_stop
  
END MODULE