MODULE mod_RKDG 
   use mod_polynome
   use mod_Divers
  IMPLICIT NONE

CONTAINS

  FUNCTION flux(u)
    IMPLICIT NONE
    REAL(prec), INTENT(in) :: u
    REAL(prec) :: flux

    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      flux = vit_adv*u
    CASE("burgers")
      flux = 0.5_prec * u**2  
    CASE DEFAULT
       WRITE(*,*) "flux non reconnu ",flux_name
       flux = 0._prec
       STOP
    END SELECT


  END FUNCTION flux



  FUNCTION flux_d(u)
    IMPLICIT NONE
    REAL(prec), INTENT(in) :: u
    REAL(prec) :: flux_d

    SELECT CASE (TRIM(flux_name))
    CASE("advection")
      flux_d = vit_adv
    CASE("burgers")
      flux_d = u  
    CASE DEFAULT
       WRITE(*,*) "flux non reconnu ",flux_name
       flux_d = 0._prec
       STOP
    END SELECT


  END FUNCTION flux_d

  FUNCTION flux_uh(x,ni)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ni
    REAL(prec), INTENT(IN) :: x
    REAL(prec) :: flux_uh

    flux_uh = flux(eval_step(x,ni))

  END FUNCTION flux_uh

  SUBROUTINE flux_numerique
    IMPLICIT NONE
    INTEGER :: ni
    REAL(prec) :: ug,ud

    ! print *, "flux"
    DO ni = 1,nb_cell
      sol_step(ni)%inter(1) = eval_step(x_cell(ni),ni)
      sol_step(ni)%inter(2) = eval_step(x_cell(ni+1),ni)
    END DO

    DO ni = 1,nb_cell
      IF (ni ==1) THEN 
        ug = sol_step(nb_cell)  %inter(2)
        ud = sol_step(1)        %inter(1)
      ELSE 
        ug = sol_step(ni-1)%inter(2)
        ud = sol_step(ni)  %inter(1)
      END IF

      flux_h(ni)  %inter(1) = (flux(ug) + flux(ud) - max_dflux*(ud-ug))  * 0.5_prec

      IF (ni ==1) THEN ;  flux_h(nb_cell)%inter(2) = flux_h(1) %inter(1)    
      ELSE ;              flux_h(ni-1)   %inter(2) = flux_h(ni)%inter(1)    
      END IF

      CALL Projection_Pk(flux_uh,flux_h(ni)%base_poly,ni)
    END DO

  END SUBROUTINE flux_numerique

  SUBROUTINE  Time_step
    IMPLICIT NONE
    
    INTEGER :: ni, tni

    INTEGER :: ii,jj
    REAL(prec) :: ti
    REAL(prec), DIMENSION(size_base) :: sig_1, sig_2
    REAL(prec), DIMENSION(size_base) :: V_B, S_B, BB

    ! print *,"calc time"

    DO ni=1,nb_cell
      sol_step(ni)%base_poly = sol(ni)%base_poly
    END DO 
      

    DO tni =1,order_t
      
        ! print *,RK_alpha(tni,1),RK_alpha(tni,2), RK_beta(tni)

      CALL flux_numerique

      DO ni =1,nb_cell

        DO ii=1,size_base
          sig_1(ii) = DG_base(Loc_to_Ref(ni,x_cell(ni)),ii); 
          sig_2(ii) = DG_base(Loc_to_Ref(ni,x_cell(ni+1)),ii); 
        END DO

        V_B = MATMUL(Rigid,flux_h(ni)%base_poly)
        S_B = -(flux_h(ni)%inter(2)*sig_2 - flux_h(ni)%inter(1)*sig_1)
        BB  = (V_B + S_B)
        L_step = MATMUL(Masse_inv, BB  )*(2._prec/(cell_size(ni))) 
        ! print *, V_B,S_b,BB
        ! print *,RK_alpha(tni,1),RK_alpha(tni,2), RK_beta(tni)

        sol_step(ni)%base_poly = RK_alpha(tni,1) * sol(ni)%base_poly + RK_alpha(tni,2) * sol_step(ni)%base_poly &
                             &+  RK_beta(tni) *dt * L_step
        
      END DO
    END DO

    ! print *,"---------------------"

    DO ni=1,nb_cell
        sol(ni)%base_poly  = sol_step(ni)%base_poly
    END DO


  END SUBROUTINE Time_step

  SUBROUTINE dt_calc
    IMPLICIT NONE

    IF(TRIM(flux_name) == "advection") THEN
      max_dflux = abs(vit_adv)
    ELSE IF(TRIM(flux_name) == "burgers") THEN
      !boucle pour max de u ? 
      max_dflux = 0._prec
      DO i=1,nb_cell
        IF(max_dflux .LT. sol(i)%base_poly(1)) THEN
          max_dflux = abs(sol(i)%base_poly(1))
        END IF
      END DO
    ELSE 
      print *,"flux non reconnue"
      max_dflux = 1._prec
    END IF
    dt = min(CFL*dx/((2*(order_x-1)+1)*max_dflux),tmax-time)

  END SUBROUTINE dt_calc

  SUBROUTINE writout
    IMPLICIT NONE

    REAL(prec) :: out1, out2,xi
    REAL(prec) :: err1 , err2, errLi

    err1 = 0._prec; err2 =0._prec; errLi = 0._prec

    IF(time >=  n_imp*t_imp)  THEN
      write(*,fmt='("--------------",i5," ",f8.4," ",e16.6, "--------------")') n_time, time, dt
      DO i=1,nb_cell
          DO j=1,size_base
              xi = Ref_to_loc(i,x_quad(j))
              out1 = eval_sol(xi,i)

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

      ! print *, sol(1)%base_poly
      ! print *, sol(nb_cell)%base_poly
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
    sol = Q_init(dicho(g,xd,xf),0)

    contains
    FUNCTION g(x_)
        use precis
        REAL(prec),INTENT(IN) :: x_
        REAL(prec)  :: g
        g = flux_d(Q_init(x_,0))*t + x_ -x
        return 
    END FUNCTION g

  END subroutine pied_charact

    ! FUNCTION Newton_search(x,t, Q_init) result(q)
    !     implicit none
    !     REAL(prec), INTENT(IN)    :: x,t
    !     REAL(prec)                :: xk, err, q, epsi = 1e-20
    !     integer             :: n=0
        
    !     interface
    !         FUNCTION Q_init(x)
    !             USE precis   
    !             REAL(prec),INTENT(IN) :: x 
    !             REAL(prec) Q_init 
    !         END FUNCTION Q_init

            
    !         FUNCTION Q_init_d(x)
    !             USE precis   
    !             REAL(prec),INTENT(IN) :: x 
    !             REAL(prec) Q_init_d 
    !         END FUNCTION Q_init_d
    !     END interface

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

    ! END FUNCTION Newton_search

END MODULE