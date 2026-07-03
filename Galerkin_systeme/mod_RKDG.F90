MODULE mod_RKDG 
   use mod_Monolith
   use mod_Divers
  IMPLICIT NONE

CONTAINS

  
  SUBROUTINE flux_numerique
    IMPLICIT NONE
    INTEGER :: ni,ii,jj
    REAL(prec), DIMENSION(nb_nodes, nb_var) :: flux_uh_val
    REAL(prec), DIMENSION(nb_var) :: ug,ud,u_temp, DF

    REAL(prec) :: x_s
    REAL(prec) :: fh_L, fh_R
    INTEGER, DIMENSION(2) :: voi_L,voi_R

    DO ni = 1,nb_cell+1
      IF (ni ==1) THEN 

        IF (TRIM(bdry_cond) == "period") THEN 
          ug = sol_step(nb_cell)  %inter(2,:)
          ud = sol_step(1)        %inter(1,:)
        ELSE IF(TRIM(bdry_cond) == "Sym") THEN
          ug = sol_step(1)        %inter(1,:)
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
        ELSE 
          print *, "boundary condition non reconnue"
        END IF

      ELSE 
        ug = sol_step(ni-1)%inter(2,:)
        ud = sol_step(ni)  %inter(1,:)
      END IF

      g(ni,:) = (flux(ug) + flux(ud) - max_dflux*(ud-ug))  * 0.5_prec


    END DO

    DO ni = 1,nb_cell
      DO ii=1,nb_nodes
        u_temp = sol_step(ni)%val_nodes(ii,:)
        flux_uh_val(ii,:) = flux(u_temp)
      END DO

      DO ii = 1,nb_var
        CALL Projection_Pk(flux_uh,flux_h(ni)%flux_DG(:,ii),LOC= LLoc,ni =ni, nvar =ii, fct_val= flux_uh_val(:,ii))
      END DO
    END DO
    
    IF(subcell_use) THEN
      
    DO ni = 1,nb_cell; DO ii = 1,nb_var
      fh_L = eval_poly(-1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)
      fh_R = eval_poly( 1._prec,ni,flux_h(ni)%flux_DG(:,ii), LOC= LRef)

      DO jj =1,nb_subcell+1
        x_s = x_subcell(jj)

        flux_h(ni)%flux_subcells(jj,ii) = eval_poly(x_s,ni=ni, base_poly=flux_h(ni)%flux_DG(:,ii),LOC= LRef) &
                                  & - C_m(jj)*(fh_L-g(ni  ,ii)) &
                                  & - C_p(jj)*(fh_R-g(ni+1,ii))
                              
      END DO

      
    END DO; END DO

    
    IF(monolithique) THEN 

    IF(smooth_extrema .GT. 0) CALL extrema_detect

    DO ni=1,nb_cell; DO jj=1,nb_subcell+1
          voi_L = Voisin_Face(ni,jj,'L'); ug = sol_step(voi_L(1))%val_subcells(voi_L(2),:)
          voi_R = Voisin_Face(ni,jj,'R'); ud = sol_step(voi_R(1))%val_subcells(voi_R(2),:)

          flux_h(ni)%flux_vf(jj,:) = (flux(ug) + flux(ud) - max_dflux*(ud-ug))  * 0.5_prec
          DF = ( flux_h(ni)%flux_subcells(jj,:)- flux_h(ni)%flux_vf(jj,:))
          theta_(ni,jj) = theta(voi_L,voi_R, DF)
                    
          flux_h(ni)%flux_subcells(jj,:) = flux_h(ni)%flux_vf(jj,:) + theta_(ni,jj)*DF
    END DO; END DO
    END IF
    END IF

  END SUBROUTINE flux_numerique

  SUBROUTINE  Time_step
    IMPLICIT NONE
    
    INTEGER :: ni, tni

    INTEGER :: ii,jj
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

      DO ni =1,nb_cell;     DO jj=1,nb_var

        V_B = MATMUL(Rigid,flux_h(ni)%flux_DG(:,jj))
        S_B = -(g(ni+1,jj)*sig_2 - g(ni,jj)*sig_1)
        BB  = (V_B + S_B)
        ! write(*,fmt='(  3( 2(f10.6) ) )')V_B,S_B,BB
        
        L_step = MATMUL(Masse_inv, BB  )*(2._prec/(cell_size(ni))) 

        sol_step(ni)%base_poly(:,jj) = RK_alpha(tni,1) * sol(ni)%base_poly(:,jj)& 
                                    &+ RK_alpha(tni,2) * sol_step(ni)%base_poly(:,jj) &
                                    &+ RK_beta(tni)    * dt * L_step
                          
        DO ii=1,nb_nodes          
          sol_step(ni)%val_nodes(ii,jj)  = eval_step(x_quad(ii),nvar= jj,ni=ni,kk= ii,LOC= LRef )
        END DO

        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol_step(ni)%inter(1,jj)      = sol_step(ni)%val_nodes(1,jj)
          sol_step(ni)%inter(2,jj)      = sol_step(ni)%val_nodes(nb_nodes,jj)
        ELSE 
          sol_step(ni)%inter(1,jj)      = eval_step(x_cell(ni),  ni=ni, LOC= LLoc, nvar=jj)
          sol_step(ni)%inter(2,jj)      = eval_step(x_cell(ni+1),ni=ni, LOC= LLoc, nvar=jj)
        END IF

      END DO; END DO
    END DO

    ! print *,"---------------------"

    DO ni=1,nb_cell; DO jj=1,nb_var
        sol(ni)%base_poly  = sol_step(ni)%base_poly

        DO ii=1,nb_nodes          
          sol(ni)%val_nodes(ii,jj)  = eval_sol(x_quad(ii),nvar= jj,ni=ni,kk= ii,LOC= LRef )
        END DO

        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol(ni)%inter(1,jj)      = sol(ni)%val_nodes(1,jj)
          sol(ni)%inter(2,jj)      = sol(ni)%val_nodes(nb_nodes,jj)
        ELSE 
          sol(ni)%inter(1,jj)      = eval_sol(x_cell(ni),  ni=ni,LOC=LLoc,nvar=jj)
          sol(ni)%inter(2,jj)      = eval_sol(x_cell(ni+1),ni=ni,LOC=LLoc,nvar=jj)
        END IF
    END DO; END DO


  END SUBROUTINE Time_step

  SUBROUTINE Time_step_subcell
    IMPLICIT NONE
    INTEGER :: ni

    INTEGER :: ii,jj,kk
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
      
      DO ni = 1,nb_cell;  DO kk=1,nb_var 

        DO jj =1,nb_subcell
          
          L = (flux_h(ni)%flux_subcells(jj+1,kk)- flux_h(ni)%flux_subcells(jj,kk))

          sol_step(ni)%val_subcells(jj,kk)=  RK_alpha(ii,1) * sol(ni)     %val_subcells(jj,kk) &
                                        &  + RK_alpha(ii,2) * sol_step(ni)%val_subcells(jj,kk) &
                                        &  - L*RK_beta(ii)  *(2._prec *dt/(cell_size(ni)* subcell_size(jj)))
                                                                      
        END DO; END DO
      END DO

      DO ni = 1,nb_cell; DO kk=1,nb_var 

        DO jj=1,size_base
          sol_step(ni)%base_poly(jj,kk) = DOT_PRODUCT(Projection_VF_inv(jj,:), sol_step(ni)%val_subcells(:,kk))
        END DO

        DO jj=1,nb_nodes
          sol_step(ni)%val_nodes(jj,kk)  = eval_step(x_quad(jj),nvar=kk,ni=ni, kk=jj, LOC= LRef)
        END DO
        IF(TRIM(quad_meth)=="Lobatto") THEN
          sol_step(ni)%inter(1,kk)      = sol_step(ni)%val_nodes(1,kk)
          sol_step(ni)%inter(2,kk)      = sol_step(ni)%val_nodes(nb_nodes,kk)
        ELSE 
          sol_step(ni)%inter(1,kk)      = eval_step(x_cell(ni),   nvar=kk,ni=ni,LOC= LLoc)
          sol_step(ni)%inter(2,kk)      = eval_step(x_cell(ni+1), nvar=kk,ni=ni,LOC= LLoc)
        END IF

      END DO; END DO

      CALL Error_check

    END DO

    DO ni=1,nb_cell
      sol(ni)%base_poly    = sol_step(ni)%base_poly
      sol(ni)%val_nodes    = sol_step(ni)%val_nodes
      sol(ni)%val_subcells = sol_step(ni)%val_subcells
      sol(ni)%inter        = sol_step(ni)%inter
    END DO




  END SUBROUTINE Time_step_subcell

  SUBROUTINE dt_calc
    IMPLICIT NONE
    INTEGER :: i,j
    INTEGER, DIMENSION(2) :: nxt_
    REAL(prec) :: gamma_temp
    REAL(prec), DIMENSION(nb_var) :: u_,v_

    IF(TRIM(flux_name) == "advection") THEN
      max_dflux = abs(vit_adv)
    ELSE 
      max_dflux = 0._prec

      IF(subcell_use .and.(.not. error_calc)) THEN

      DO i=1,nb_cell;    DO j=1,nb_subcell
        nxt_ = Voisin_Face(i,j,'R')

        u_=sol(i)%val_subcells(j,:); v_ = sol(nxt_(1))%val_subcells(nxt_(2),:)
        gamma_temp = gamma_calc(u_, v_)
        
        max_dflux = max(max_dflux, gamma_temp)
      END DO;END DO

      ELSE 
        DO i=1,nb_cell;    DO j=1,nb_nodes
          nxt_ = Voisin_quad(i,j,'R')

          u_=sol(i)%val_nodes(j,:); v_ = sol(nxt_(1))%val_nodes(nxt_(2),:)
          gamma_temp = gamma_calc(u_, v_)
                    
          max_dflux = max(max_dflux, gamma_temp)
        END DO;END DO
      END IF
    END IF
    
    dt_old = dt
    dt = min(CFL*dx/(REAL(2*order_x-1,prec)*max_dflux),tmax-time)

    ! IF(order_x .GT. order_t) THEN
    ! dt = min(   (CFL*dx/(REAL(2*order_x-1,prec)*max_dflux)) **(Real(order_x,prec) * 1._prec/Real(order_t,prec)) ,tmax-time) 
    ! END IF 

    IF(dt .LT. 10._prec**(-20) .AND. (time+dt .LT.tmax)) THEN
      write(*, fmt ='("dt trop petit : dt =",e16.6, 1x,",max dflux =",e16.6 )') dt, max_dflux
      CALL Emergency_stop
    END IF
    
    IF(ISNAN(sol(1)%val_nodes(1,1))) THEN
      write(*, fmt ='("NAN found emergency stop")')
      CALL Emergency_stop
    END IF


  END SUBROUTINE dt_calc

  SUBROUTINE writout(switch)
    IMPLICIT NONE

    INTEGER :: i,j,k
    REAL(prec), DIMENSION(nb_var) :: out, out_ex,u_
    REAL(prec) :: xi
    REAL(prec) :: err1 , err2, errLi, pression_,p_max,p_min
    LOGICAL, optional :: switch
    LOGICAL :: force
    Character(len=63) :: save_format

    IF(present(switch)) THEN; force = switch
    ELSE; force = .FALSE.
    END IF

    err1 = 0._prec; err2 =0._prec; errLi = 0._prec; p_max = 0._prec
    IF(modulo(n_time,500) == 0)  THEN
      write(*,fmt='("---------------",i7,2x,f10.6,2x,e16.6, "--------------")') n_time, time, dt
    END IF

    IF((time .GE.  REAL(n_imp,prec)*t_imp-eps0) .OR. force )  THEN
      CALL Out_The_Mesh(time)
      write(*,fmt='("---------------",i7,2x,f10.6,2x,e16.6, "--------------")') n_time, time, dt
      DO i=1,nb_cell
        IF(subcell_use .and.(.not. error_calc)) THEN

          save_format = "(f10.6"//Repeat(",f16.6",nb_var)//")"

          DO j=1,nb_subcell; 
            xi = Ref_to_loc(i,x_submiddle(j))
            out = sol(i)%val_subcells(j,:)
            write(unit=numfile_sol,  fmt=save_format,  advance="no") xi,out

            ! IF(TRIM(flux_name)== "Euler" ) THEN
            ! u_ = sol(i)%val_subcells(j,:); pression_ = pression(u_)
            ! p_max = max(pression_, p_max); p_min = min(pression_,p_min)
            ! write(unit=numfile_sol, fmt= '(f10.6)') pression_
            ! ELSE; 
              write(unit=numfile_sol, fmt= '(1x)')
            ! END IF

          END DO; 

        ELSE
          
          DO j=1,size_base
            xi = Ref_to_loc(i,x_quad(j))

            out = sol(i)%val_nodes(j,:)

            DO k=1,nb_var
              IF(TRIM(flux_name) == "advection") THEN               
                out_ex(k) =Q_init(xi - time*vit_adv,0,nvar=k)
              ELSE 
                ! call pied_charact(xi,time,out_ex(k))
                out_ex =0._prec
              END IF
            END DO

            save_format = "(f10.6"//Repeat(",f16.6",nb_var)//")"
            write(unit=numfile_sol,  fmt=save_format,  advance="no") xi,out


            errLi = max(errLi , maxval(abs(out-out_ex)))
            err1 = err1 + SUM(abs(out-out_ex))*w_quad(j)    *cell_size(i)/2
            err2 = err2 + SUM(((out-out_ex)*w_quad(j))**2)  *cell_size(i)/2

            ! IF(TRIM(flux_name)== "Euler" ) THEN
            ! u_ = sol(i)%val_nodes(j,:); pression_ = pression(u_)
            ! p_max = max(pression_, p_max); p_min = min(pression_,p_min)
            ! write(unit=numfile_sol, fmt= '(f10.6)') pression_
            ! ELSE;
               write(unit=numfile_sol, fmt= '(1x)')
            ! END IF

          END DO

        END IF
      END DO

      err2 = sqrt(err2)
      IF(TRIM(flux_name)== "Euler" ) write(*,fmt='("pression in [", f10.6,",",f10.6,"]")') p_min,p_max
      write(*, fmt ='("err L1 = ", e20.12)')  err1
      write(*, fmt ='("err L2 = ", e20.12)')  err2
      write(*, fmt ='("err Li = ", e20.12)')  errLi

      err_L1 = max(err1, err_L1); err_L2 = max(err2, err_L2); err_Li = max(errLi, err_Li)

      n_imp = n_imp +1
      Time_stemp(n_imp) = time
      
      write(unit=numfile_sol, fmt='("----------",f10.6,"--------------")' ) time
    END IF
  END SUBROUTINE writout
  
  ! subroutine pied_charact(x,t,sol)

  !   REAL(prec), INTENT(IN) :: x,t
  !   REAL(prec), intent(out):: sol
  !   REAL(prec) :: xd,xf

  !   IF(TRIM(flux_name)=="advection") THEN
  !     xd = xL-abs(vit_adv)*t; xf = xR + abs(vit_adv)*t
  !   ELSE
  !     xd = xL-t; xf = xR + t
  !   END IF
  !   sol = Q_init(dicho(h,xd,xf),0)

  !   contains
  !   FUNCTION h(x_)
  !       use precis
  !       REAL(prec),INTENT(IN) :: x_
  !       REAL(prec)  :: h
  !       h = flux_d(Q_init(x_,0))*t + x_ -x
  !       return 
  !   END FUNCTION h

  ! END subroutine pied_charact


  SUBROUTINE Out_The_Mesh(ti)
    IMPLICIT NONE
    REAL(prec), INTENT(IN) :: ti
    INTEGER :: ni,n_sub
    REAL(prec) :: xi, out1

    write(unit=numfile_meshout, fmt='("---------",f10.6,"---------------")' ) ti
    DO ni =1,nb_cell;    DO n_sub=1,nb_subcell
      xi = Ref_to_loc(ni,x_subcell(n_sub))
      out1 = theta_(ni,n_sub)
      write(unit=numfile_meshout,  fmt='(f10.6, f16.6, f16.6, 2x, l1)') xi,out1
    END DO; END DO

  END SUBROUTINE Out_The_Mesh

  SUBROUTINE Emergency_stop
    INTEGER :: i

    CALL writout(.TRUE.)

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