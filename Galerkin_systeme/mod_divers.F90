MODULE mod_Divers
    use mod_declaration
    IMPLICIT NONE
    
CONTAINS

    FUNCTION dicho (fct,xd,xf)

        REAL(prec), INTENT(IN) :: xd,xf
        REAL(prec) :: dicho

        interface
            FUNCTION fct(x)
                USE precis   
                REAL(prec),INTENT(IN) :: x 
                REAL(prec) fct 
            END FUNCTION fct
        END interface

        REAL(prec) :: a,b,c
        REAL(prec) :: t
        integer :: n = 0

        if(abs(fct(xf)) .LT. eps0) then 
            dicho = xf
            return
        else if(abs(fct(xd)).LT.eps0) then
            dicho = xd
            return
        END if


        a=xd; b=xf; 

        n =0
        do while(abs(fct(c))>eps0 .and. n<50 )
            c=(a+b)/2._prec
            t = sign(1._prec,fct(c)*fct(a))
            if(t <eps0) then 
                b=c
            else 
                a=c
            END if
            n = n+1
        END do

        dicho = c
    END FUNCTION dicho

    FUNCTION minmod(x,y,z)

        implicit none
        REAL(prec), INTENT(IN)    :: x,y,z
        REAL(prec)                :: minmod

        minmod = min(0._prec,max(x,y,z)) + max(0._prec,min(x,y,z))

    END FUNCTION

    ! FUNCTION Moindre_carre(A,b,n1,n2, A_inv) result(X)
        !     IMPLICIT NONE
        !     INTEGER,  INTENT(IN)  :: n1,n2
        !     REAL(prec), DIMENSION(n1,n2), INTENT(IN) :: A
        !     REAL(prec), DIMENSION(n2,n1), INTENT(IN),  optional :: A_inv
        !     REAL(prec), DIMENSION(n2),    INTENT(IN)  :: b

        !     REAL(prec), DIMENSION(n1,n2) :: mat_inv
        !     REAL(prec), DIMENSION(n1) :: X


        !     IF(n1 == n2) THEN 
        !         IF(.not.present(A_inv)) THEN; CALL inv_mat(A,mat_inv,0);    X = MATMUL(mat_inv,b)
        !         ELSE IF(present(A_inv)) THEN; X = MATMUL(A_inv,b)
        !         END IF
        !     ELSE
        !         IF(.not.present(A_inv)) THEN; 
        !             CALL inv_mat(MATMUL(Transpose(A),A),mat_inv,0);
        !             X = MATMUL(mat_inv,MATMUL(Transpose(A),b) )
        !         ELSE IF(present(A_inv)) THEN; X = MATMUL(A_inv,MATMUL(Transpose(A),b))
        !         END IF

        !     END IF

    ! END FUNCTION Moindre_carre
    
    FUNCTION Newton_search(fct, fct_d, x0) result(xk)
            implicit none

            REAL(prec), INTENT(IN), optional :: x0


            interface
                FUNCTION fct(x)
                    USE precis   
                    REAL(prec),INTENT(IN) :: x 
                    REAL(prec) fct 
                END FUNCTION fct

                
                FUNCTION fct_d(x)
                    USE precis   
                    REAL(prec),INTENT(IN) :: x 
                    REAL(prec) fct_d
                END FUNCTION fct_d
            END interface

            REAL(prec) :: xk, err
            integer    :: n=0
        
            n = 0
            IF(present(x0)) THEN;   xk =x0
            ELSE;                   xk =0._prec
            END IF

            err = abs(fct(xk))

            do while(err>eps0 .and. n<50)
                xk = xk -   (fct(xk))/(fct_d(xk))
                err = abs(fct(xk))
                n = n+1
            END do
            
    END FUNCTION Newton_search
  
    FUNCTION factoriel(n,nb)
        INTEGER, INTENT(IN) :: n,nb
        REAL(prec) :: factoriel

        INTEGER :: i

        factoriel = 1._prec
        DO i=n,(n-nb),-1
            factoriel = REAL(i,prec)*factoriel
        END DO

    END FUNCTION factoriel

    FUNCTION Voisin_Face(ni,ns,LR)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ni,ns
        CHARACTER(len =1) :: LR
        INTEGER, DIMENSION(2):: Voisin_Face

        ! print *,ni,ns,LR

        IF(LR == "L") THEN
        IF(ni == 1 .AND. ns == 1) THEN 
            IF(TRIM(bdry_cond) == "period") THEN
            Voisin_Face(1) = nb_cell; Voisin_Face(2) = nb_subcell
            ELSE IF(TRIM(bdry_cond)=="Sym" .OR. TRIM(bdry_cond)=="Solid") THEN 
            Voisin_Face(1) = 1 ; Voisin_Face(2) = 1
            END IF

        ELSE IF(ns ==1 ) THEN
            Voisin_Face(1) = ni-1 ; Voisin_Face(2) = nb_subcell
        ELSE 
            Voisin_Face(1) = ni ; Voisin_Face(2) = ns-1
        END IF

        ELSE IF(LR == "R") THEN
        IF(ni == nb_cell .AND. ns == nb_subcell+1) THEN 
            IF(TRIM(bdry_cond) == "period") THEN
            Voisin_Face(1) = 1; Voisin_Face(2) = 1
            ELSE IF(TRIM(bdry_cond)=="Sym" .OR. TRIM(bdry_cond)=="Solid") THEN 
            Voisin_Face(1) = nb_cell ; Voisin_Face(2) = nb_subcell
            END IF
            
        ELSE IF(ns == nb_subcell+1 ) THEN
            Voisin_Face(1) = ni+1 ; Voisin_Face(2) = 1
        ELSE 
            Voisin_Face(1) = ni ; Voisin_Face(2) = ns
        END IF

        ELSE 
        print *,"direction non reconnue"
        STOP
        END IF

        ! print *,Voisin_Face

    END FUNCTION Voisin_Face

    FUNCTION Voisin_cell(ni,n_sub,LR) result(V_c)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ni, n_sub
        CHARACTER(len =1) :: LR
        INTEGER, DIMENSION(2) :: V_c

        IF(LR == "L") THEN
        IF(ni == 1 .AND. n_sub ==1 ) THEN 
            IF(TRIM(bdry_cond) == "period") THEN
            V_c = [nb_cell, nb_subcell]
            ELSE IF(TRIM(bdry_cond)=="Sym" .OR. TRIM(bdry_cond)=="Solid") THEN 
            V_c = [1, 1]
            END IF
        ELSE IF (n_sub == 1) THEN 
            V_c = [ni-1, nb_subcell]
        ELSE 
            V_c = [ni, n_sub-1]
        END IF

        ELSE IF(LR == "R") THEN
        IF(ni == nb_cell .AND. n_sub == nb_subcell) THEN 
            IF(TRIM(bdry_cond) == "period") THEN
            V_c = [1,1]
            ELSE IF(TRIM(bdry_cond)=="Sym" .OR. TRIM(bdry_cond)=="Solid") THEN 
            V_c = [nb_cell, nb_subcell]
            END IF
        ELSE IF(n_sub == nb_subcell) THEN
            V_c = [ni+1, 1]
        ELSE  
            V_c = [ni, n_sub+1]
        END IF

        ELSE 
        print *,"direction non reconnue"
        STOP
        END IF


    END FUNCTION Voisin_cell

    FUNCTION Voisin_quad(ni,nq,LR)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ni,nq
        CHARACTER(len =1) :: LR
        INTEGER, DIMENSION(2):: Voisin_quad


        IF(LR == "L") THEN
        IF(ni == 1 .AND. nq == 1) THEN 
            IF(TRIM(bdry_cond) == "period") THEN
            Voisin_quad(1) = nb_cell; Voisin_quad(2) = nb_nodes
            ELSE IF(TRIM(bdry_cond)=="Sym" .OR. TRIM(bdry_cond)=="Solid") THEN 
            Voisin_quad(1) = 1 ; Voisin_quad(2) = 1
            END IF

        ELSE IF(nq ==1 ) THEN
            Voisin_quad(1) = ni-1 ; Voisin_quad(2) = nb_nodes
        ELSE 
            Voisin_quad(1) = ni ; Voisin_quad(2) = nq-1
        END IF

        ELSE IF(LR == "R") THEN
        IF(ni == nb_cell .AND. nq == nb_nodes) THEN 
            IF(TRIM(bdry_cond) == "period") THEN
            Voisin_quad(1) = 1; Voisin_quad(2) = 1
            ELSE IF(TRIM(bdry_cond)=="Sym" .OR. TRIM(bdry_cond)=="Solid") THEN 
            Voisin_quad(1) = nb_cell ; Voisin_quad(2) = nb_nodes
            END IF
            
        ELSE IF(nq == nb_nodes ) THEN
            Voisin_quad(1) = ni+1 ; Voisin_quad(2) = 1
        ELSE 
            Voisin_quad(1) = ni ; Voisin_quad(2) = nq
        END IF

        ELSE 
        print *,"direction non reconnue"
        STOP
        END IF

        ! print *,Voisin_quad

    END FUNCTION Voisin_quad


  SUBROUTINE eval_time(prd_ini,prd_max,prd_sec,time_r)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: time_r,prd_ini,prd_max,prd_sec
    INTEGER :: prd_fin,prd
    REAL(prec) :: time_cal


    CALL SYSTEM_CLOCK(COUNT=prd_fin)
    prd=prd_fin-prd_ini
    IF (prd_fin.LT.prd_ini) prd=prd+prd_max
    time_cal=REAL(prd,prec)/REAL(prd_sec,prec)

    CALL display_time(time_cal,"COMPUTATIONAL TIME :   ")

    IF (time_r==1 .AND. time>0._prec) THEN
       CALL display_time((tmax/time-1._prec)*time_cal," ~= REMAINING TIME :   ")
    END IF

  END SUBROUTINE eval_time

  SUBROUTINE display_time(time_v,mess_v)
    IMPLICIT NONE

    REAL(prec),INTENT(IN) :: time_v
    CHARACTER(LEN=23),INTENT(IN) :: mess_v
    INTEGER :: vday,vhour,vminute,vseconde, vcent
    CHARACTER(LEN=2) :: fday,fhour,fminute,fseconde, fcent

    vday=INT(time_v)/86400
    vhour=INT(time_v-REAL(vday*86400,prec))/3600
    vminute=INT(time_v-REAL(vday*86400+vhour*3600,prec))/60
    vseconde=INT(time_v-REAL(vday*86400+vhour*3600+vminute*60,prec))
    vcent= INT(100._prec*(time_v-REAL(vday*86400+vhour*3600+vminute*60+vseconde,prec) ) )
    WRITE(fday,'(i2.2)') vday
    WRITE(fhour,'(i2.2)') vhour
    WRITE(fminute,'(i2.2)') vminute
    WRITE(fseconde,'(i2.2)') vseconde
    WRITE(fcent,'(i2.2)') vcent
    IF (vday+vhour+vminute==0) THEN
       WRITE(*,*) mess_v,fseconde," SECONDES AND ", &
            & fcent, " CENTIEME "
    ELSE IF (vday+vhour==0) THEN
       WRITE(*,*) mess_v,fminute," MINUTES AND ",&
            & fseconde," SECONDES"
    ELSE IF (vday==0) THEN
       WRITE(*,*) mess_v,fhour," HOURS, ",&
            & fminute," MINUTES AND ",fseconde," SECONDES"
    ELSE
       WRITE(*,*) mess_v,fday," DAYS,",fhour,&
            & " HOURS, ",fminute," MINUTES AND ",fseconde," SECONDES"
    END IF

    print *,time_v

  END SUBROUTINE display_time
    
  SUBROUTINE Knapsack_greedy(item_pr,item_volume,volume_cstr,item_cstr,nb_item)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nb_item
    REAL(prec), DIMENSION(nb_item), INTENT(OUT) :: item_pr
    REAL(prec), DIMENSION(nb_item), INTENT(IN) :: item_volume
    REAL(prec), DIMENSION(nb_item), INTENT(IN) :: item_cstr
    REAL(prec), INTENT(IN) :: volume_cstr

    REAL(prec), DIMENSION(nb_item) :: item_in,item_vol
    INTEGER :: i_maxloc

    item_in = item_cstr; item_vol = item_volume

    IF(volume_cstr .LT. -eps0) THEN; print *,"Knapsack : constraint error"; STOP; ENDIF;

    DO WHILE(DOT_PRODUCT(item_in, item_vol) .GT. volume_cstr-eps0)
        i_maxloc = MAXLOC(item_vol,DIM =nb_item)
        item_vol(i_maxloc) = 0._prec 

        IF(DOT_PRODUCT(item_in, item_vol) .LT. volume_cstr-eps0 &
        & .AND. DOT_PRODUCT(item_in, item_vol) .GT. eps0) THEN 
            item_in(i_maxloc) = (volume_cstr-DOT_PRODUCT(item_in, item_vol))/item_volume(i_maxloc)
        ELSE;   
            item_in(i_maxloc) = 0._prec
        END IF
        
        ! print *,DOT_PRODUCT(item_in, item_vol),volume_cstr
    END DO
    item_pr = item_in
    ! write(*, fmt ='("theta = ", f10.6, f10.6, f10.6)') item_pr
  END SUBROUTINE Knapsack_greedy

END MODULE mod_Divers