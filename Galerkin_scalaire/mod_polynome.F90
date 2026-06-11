
MODULE mod_polynome
    USE mod_declaration
    IMPLICIT NONE
    
CONTAINS


    FUNCTION sinus(x,ni)
        IMPLICIT NONE
        REAL(prec), INTENT(in) :: x
        INTEGER,    INTENT(in) :: ni
        REAL(prec) :: sinus 
        sinus = sin(2._prec*pi*x)
    END FUNCTION sinus
    
    FUNCTION creneau(x,ni)
        IMPLICIT NONE
        REAL(prec), INTENT(in) :: x
        INTEGER,    INTENT(in) :: ni
        REAL(prec) :: creneau 
        creneau = 0._prec
        if(-0.5_prec<x .and. x<0._prec) creneau = 1._prec
    END FUNCTION creneau

    FUNCTION Q_init(x,ni)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: x
        INTEGER,    INTENT(IN) :: ni
        REAL(prec) :: Q_init

        IF     (TRIM(sol_ini_name)=="sinus") THEN
            Q_init = sinus(x,ni)
        ELSE IF(TRIM(sol_ini_name)=="unit") THEN
            Q_init = unit(x,ni)
        ELSE IF(TRIM(sol_ini_name)=="Riemann") THEN
            Q_init = 0._prec
            if(x<0._prec) Q_init =1._prec
        ELSE IF(TRIM(sol_ini_name)=="creneau") THEN
            Q_init = creneau(x,ni)
        ELSE IF(TRIM(sol_ini_name)=="Burgers_choc") THEN
            Q_init = 0._prec
            if(0.3_prec<x .and. x<0.7_prec) Q_init = -1._prec
            if(x>0.7_prec) Q_init = 0.5_prec
        END IF
    END FUNCTION Q_init

! ---------------------------------------------------------------

    FUNCTION eval_poly(YY,ni, base_poly, kk, LOC)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: ni
        INTEGER, INTENT(in), optional :: kk
        CHARACTER (len=8), INTENT(IN) :: LOC
        REAL(prec), INTENT(in) :: YY
        REAL(prec), DIMENSION(size_base), INTENT(IN) :: base_poly
        REAL(prec)  :: eval_poly
        INTEGER :: ii
        eval_poly= 0._prec


        IF(.not. present(kk)) THEN
            counter1 = counter1 +1
            DO ii = 1,size_base
                eval_poly = eval_poly + base_poly(ii) * DG_base(YY,ii,LOC,ni)
            END DO
        ELSE 
            counter2 = counter2 +1
            DO ii = 1,size_base
                eval_poly = eval_poly + base_poly(ii) * sig_quad(ii,kk)
            END DO
        END IF

    END FUNCTION eval_poly

    FUNCTION eval_sol(YY,ni, kk, LOC)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: ni
        INTEGER, INTENT(in), optional :: kk
        CHARACTER (len=8), INTENT(IN) :: LOC
        REAL(prec)   , INTENT(in) :: YY
        REAL(prec) :: eval_sol 
        INTEGER :: ii
        eval_sol= 0._prec

        IF(.not. present(kk)) THEN
            counter1 = counter1 +1
            DO ii = 1,size_base
                eval_sol = eval_sol + sol(ni)%base_poly(ii) * DG_base(YY,ii,LOC,ni)
            END DO
        ELSE 
            counter2 = counter2 +1
            DO ii = 1,size_base
                eval_sol = eval_sol + sol(ni)%base_poly(ii) * sig_quad(ii,kk)
            END DO
        END IF

    END FUNCTION eval_sol
    
    FUNCTION eval_step(YY,ni,kk, LOC)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: ni
        INTEGER, INTENT(in), optional :: kk
        CHARACTER (len=8), INTENT(IN) :: LOC
        REAL(prec)   , INTENT(in) :: YY
        REAL(prec) :: eval_step 
        INTEGER :: ii
        eval_step= 0._prec
        
        IF(.not. present(kk)) THEN
            DO ii = 1,size_base
                eval_step = eval_step + sol_step(ni)%base_poly(ii) * DG_base(YY,ii, LOC,ni)
            END DO
        ELSE 
            DO ii = 1,size_base
                eval_step  = eval_step + sol_step(ni)%base_poly(ii) * sig_quad(ii,kk)
            END DO
        END IF

    END FUNCTION eval_step

    FUNCTION quadrature(fct1,opt1,fct2,opt2,LOC,ni,n_sub) 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ni
        INTEGER, INTENT(IN), optional :: n_sub 
        CHARACTER (len=8), INTENT(IN) :: LOC
        INTERFACE
            FUNCTION fct1(YY,opt)
                USE precis
                REAL(prec),INTENT(IN) :: YY
                INTEGER,   INTENT(IN) :: opt
                REAL(prec) :: fct1
            END FUNCTION fct1
            
            FUNCTION fct2(YY,opt)
                USE precis
                REAL(prec),INTENT(IN) :: YY
                INTEGER,   INTENT(IN) :: opt
                REAL(prec) :: fct2
            END FUNCTION fct2
        END INTERFACE

        INTEGER,   INTENT(IN) :: opt1,opt2

        REAL(prec) :: quadrature
        REAL(prec) :: YY
        INTEGER :: kk
        
        quadrature = 0._prec

            IF(TRIM(LOC) == "Loc") THEN
            DO kk =size_quad_nodes,1,-1
                YY = Ref_to_loc(ni,x_quad(kk))
                quadrature = quadrature + fct1(YY,opt1)*fct2(YY,opt2)*w_quad(kk)*cell_size(ni)/2._prec
            END DO

        ELSEIF(TRIM(LOC) == "Ref") THEN
            DO kk =size_quad_nodes,1,-1
                quadrature = quadrature + fct1(x_quad(kk),opt1)*fct2(x_quad(kk),opt2)*w_quad(kk)
            END DO

        ELSEIF(TRIM(LOC) == "SubRef") THEN
            DO kk =size_quad_nodes,1,-1
                YY = Ref_to_Refsub(x_quad(kk),n_sub)
                quadrature = quadrature + fct1(YY,opt1)*fct2(YY,opt2)*w_quad(kk)*(subcell_size(n_sub)/2._prec)
            END DO
        END IF   

    END FUNCTION quadrature

    FUNCTION integration(fct1,opt1,fct2,opt2)
        IMPLICIT NONE

        INTERFACE
        FUNCTION fct1(YY,opt)
            USE precis
            REAL(prec),INTENT(IN) :: YY
            INTEGER,   INTENT(IN) :: opt
            REAL(prec) :: fct1
        END FUNCTION fct1
        
        FUNCTION fct2(YY,opt)
            USE precis
            REAL(prec),INTENT(IN) :: YY
            INTEGER,   INTENT(IN) :: opt
            REAL(prec) :: fct2
        END FUNCTION fct2
        END INTERFACE

        INTEGER,   INTENT(IN) :: opt1,opt2

        INTEGER :: ni
        REAL(prec) :: integration

        integration = 0._prec

        DO ni = 1,nb_cell
            integration = integration + quadrature(fct1,opt1, fct2,opt2,LLoc,ni)
        END DO


    END FUNCTION integration

    FUNCTION unit(x, ni)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: x
        INTEGER,    INTENT(IN) :: ni
        REAL(prec) :: unit 
        unit =1._prec
    END FUNCTION unit

    FUNCTION Lagrange_basis(XX,ii)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: XX 
        INTEGER,    INTENT(IN) :: ii
        INTEGER :: jj 
        REAL(prec) :: Lagrange_basis
        Lagrange_basis = 1._prec

        DO jj =1,size_base

            IF(jj .NE. ii) THEN 
                Lagrange_basis = Lagrange_basis * (XX-pts_DG(jj))/(pts_DG(ii)-pts_DG(jj))
            END IF

        END DO

    END FUNCTION Lagrange_basis

    SUBROUTINE print_mat(Mat,t1,t2)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t1,t2
        REAL(prec), DIMENSION(t1,t2) :: Mat
        INTEGER :: ii,jj

        DO ii=1,t1
            write(*,fmt ='()')
            write(*,fmt ='(a)', advance="no") '|'
            DO jj =1,t2
                write(*,fmt ='(f16.6,a)',advance="no") Mat(ii,jj), '|'
            END DO
        END DO
            write(*,fmt ='()')

    END SUBROUTINE print_mat
    
    SUBROUTINE Coeff_DG_init
        IMPLICIT NONE

        REAL(prec), DIMENSION(size_base ) :: temp
    
        IF(TRIM(DG_meth) == "Taylor") THEN
            call base_Taylor_init
            coeff_DG(:,:) = coeff_Taylor(:,:)

        ELSE IF (TRIM(DG_meth) == "Legendre") THEN
            call base_Legendre_init
            coeff_DG(:,:) = coeff_legendre(:size_base,:size_base)

        ELSE IF (TRIM(DG_meth) == "Lobatto") THEN
            CALL quad_1D_lobatto(size_base,pts_DG,temp)

        END IF

    END SUBROUTINE Coeff_DG_init

    SUBROUTINE Coeff_quad_init
        IMPLICIT NONE

        IF(TRIM(quad_meth) == "Lobatto") THEN
            CALL quad_1D_lobatto(size_quad_nodes,x_quad,w_quad)
        ELSE IF(TRIM(quad_meth) == "Legendre") THEN
            CALL quad_1D_legendre(size_quad_nodes,x_quad,w_quad)
        END IF

    END SUBROUTINE Coeff_quad_init
    
    SUBROUTINE base_Legendre_init
        IMPLICIT NONE

        IF (order_x.GT.9) THEN
        WRITE(*,*) "PROBLEM IN init_coef_legendre"
        WRITE(*,*) "THE order_x IS TOO HIGH"
        STOP
        END IF
        
        coeff_legendre=0._prec
        coeff_legendre(1,1)=1._prec
        coeff_legendre(2,2)=1._prec
        coeff_legendre(3,3)=3._prec/2._prec
        coeff_legendre(3,1)=-1._prec/2._prec
        coeff_legendre(4,4)=5._prec/2._prec
        coeff_legendre(4,2)=-3._prec/2._prec
        coeff_legendre(5,5)=35._prec/8._prec
        coeff_legendre(5,3)=-15._prec/4._prec
        coeff_legendre(5,1)=3._prec/8._prec
        coeff_legendre(6,6)=63._prec/8._prec
        coeff_legendre(6,4)=-35._prec/4._prec
        coeff_legendre(6,2)=15._prec/8._prec
        coeff_legendre(7,7)=231._prec/16._prec
        coeff_legendre(7,5)=-315._prec/16._prec
        coeff_legendre(7,3)=105._prec/16._prec
        coeff_legendre(7,1)=-5._prec/16._prec
        coeff_legendre(8,8)=429._prec/16._prec
        coeff_legendre(8,6)=-693._prec/16._prec
        coeff_legendre(8,4)=315._prec/16._prec
        coeff_legendre(8,2)=-35._prec/16._prec
        coeff_legendre(9,9)=6435._prec/128._prec
        coeff_legendre(9,7)=-3003._prec/32._prec
        coeff_legendre(9,5)=3465._prec/64._prec
        coeff_legendre(9,3)=-315._prec/32._prec
        coeff_legendre(9,1)=35._prec/128._prec
            
    END SUBROUTINE base_Legendre_init

    SUBROUTINE base_Taylor_init
        IMPLICIT NONE

        IF (order_x.GT.9) THEN
        WRITE(*,*) "PROBLEM IN init_coef_moment"
        WRITE(*,*) "THE order_x IS TOO HIGH"
        STOP
        END IF
        
        coeff_Taylor=0._prec
        coeff_Taylor(1,2)=1._prec
        coeff_Taylor(2,3)=1._prec
        coeff_Taylor(2,1)=-1._prec/12._prec
        coeff_Taylor(3,4)=1._prec
        coeff_Taylor(3,2)=-1._prec/4._prec
        coeff_Taylor(4,5)=1._prec
        coeff_Taylor(4,3)=-1._prec/2._prec
        coeff_Taylor(4,1)=7._prec/240._prec
        coeff_Taylor(5,6)=1._prec
        coeff_Taylor(5,4)=-5._prec/6._prec
        coeff_Taylor(5,2)=7._prec/48._prec
        coeff_Taylor(6,7)=1._prec
        coeff_Taylor(6,5)=-5._prec/4._prec
        coeff_Taylor(6,3)=7._prec/16._prec
        coeff_Taylor(6,1)=-31._prec/1344._prec
        coeff_Taylor(7,8)=1._prec
        coeff_Taylor(7,6)=-7._prec/4._prec
        coeff_Taylor(7,4)=49._prec/48._prec
        coeff_Taylor(7,2)=-31._prec/192._prec
        coeff_Taylor(8,9)=1._prec
        coeff_Taylor(8,7)=-7._prec/3._prec
        coeff_Taylor(8,5)=49._prec/24._prec
        coeff_Taylor(8,3)=-31._prec/48._prec
        coeff_Taylor(8,1)=127._prec/3840._prec

    END SUBROUTINE base_Taylor_init

    SUBROUTINE quad_1D_lobatto(n_quad_1D,x_quad,w_quad)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: n_quad_1D 
        REAL(prec),DIMENSION(:),INTENT(OUT) :: x_quad,w_quad
        REAL(prec) :: aa,bb

        x_quad(1)=-1._prec
        x_quad(n_quad_1D)=1._prec
        SELECT CASE (n_quad_1D)
        CASE (2)
        w_quad=1._prec
        CASE (3)
        x_quad(2)=0._prec
        w_quad(1)=1._prec/3._prec
        w_quad(2)=4._prec/3._prec
        w_quad(3)=1._prec/3._prec
        CASE (4)
        aa=SQRT(1._prec/5._prec)
        x_quad(2)=-aa
        x_quad(3)=aa
        bb=1._prec/9._prec
        w_quad(1)=1._prec/6._prec
        w_quad(2)=5._prec/6._prec
        w_quad(3)=5._prec/6._prec
        w_quad(4)=1._prec/6._prec
        CASE (5)
        aa=SQRT(3._prec/7._prec)
        x_quad(2)=-aa
        x_quad(3)=0._prec
        x_quad(4)=aa
        w_quad(1)=1._prec/10._prec
        w_quad(2)=49._prec/90._prec
        w_quad(3)=32._prec/45._prec
        w_quad(4)=49._prec/90._prec
        w_quad(5)=1._prec/10._prec
        CASE (6)
        aa=SQRT(1._prec/3._prec-2._prec*SQRT(7._prec)/21._prec)
        bb=SQRT(1._prec/3._prec+2._prec*SQRT(7._prec)/21._prec)
        x_quad(2)=-bb
        x_quad(3)=-aa
        x_quad(4)=aa
        x_quad(5)=bb
        aa=(14._prec+SQRT(7._prec))/30._prec
        bb=(14._prec-SQRT(7._prec))/30._prec
        w_quad(1)=1._prec/15._prec
        w_quad(2)=bb
        w_quad(3)=aa
        w_quad(4)=aa
        w_quad(5)=bb
        w_quad(6)=1._prec/15._prec
        CASE (7)
        aa=SQRT(5._prec/11._prec-2._prec*SQRT(5._prec/3._prec)/11._prec)
        bb=SQRT(5._prec/11._prec+2._prec*SQRT(5._prec/3._prec)/11._prec)
        x_quad(2)=-bb
        x_quad(3)=-aa
        x_quad(4)=0._prec
        x_quad(5)=aa
        x_quad(6)=bb
        aa=(124._prec+7._prec*SQRT(15._prec))/350._prec
        bb=(124._prec-7._prec*SQRT(15._prec))/350._prec
        w_quad(1)=1._prec/21._prec
        w_quad(2)=bb
        w_quad(3)=aa
        w_quad(4)=256._prec/525._prec
        w_quad(5)=aa
        w_quad(6)=bb
        w_quad(7)=1._prec/21._prec
        CASE (8)
        x_quad(5)=0.2092992179024788687686572603453513_prec
        x_quad(6)=0.5917001814331423021445107313979532_prec
        x_quad(7)=0.8717401485096066153374457612206634_prec
        x_quad(4)=-x_quad(5)
        x_quad(3)=-x_quad(6)
        x_quad(2)=-x_quad(7)
        
        w_quad(5)=0.41245879465870388156705297140221_prec
        w_quad(6)=0.341122692483504364764240677107748_prec
        w_quad(7)=0.2107042271435060393829920657757563_prec
        w_quad(8)=0.03571428571428571428571428571428571_prec
        w_quad(4)=w_quad(5)
        w_quad(3)=w_quad(6)
        w_quad(2)=w_quad(7)
        w_quad(1)=w_quad(8)
        CASE (9)
        x_quad(5)=0._prec
        x_quad(6)=0.3631174638261781587107520687086592_prec
        x_quad(7)=0.6771862795107377534458854270913425_prec
        x_quad(8)=0.899757995411460157312345244418338_prec
        x_quad(4)=-x_quad(6)
        x_quad(3)=-x_quad(7)
        x_quad(2)=-x_quad(8)
        
        w_quad(5)=0.3715192743764172335600907029478458_prec
        w_quad(6)=0.346428510973046345115131532139718_prec
        w_quad(7)=0.2745387125001617352807056185793727_prec
        w_quad(8)=0.1654953615608055250463397200292083_prec
        w_quad(9)=0.02777777777777777777777777777777778_prec
        w_quad(4)=w_quad(6)
        w_quad(3)=w_quad(7)
        w_quad(2)=w_quad(8)
        w_quad(1)=w_quad(9)
        CASE (10)
        x_quad(6)=0.1652789576663870246262197659581735_prec
        x_quad(7)=0.477924949810444495661175092731258_prec
        x_quad(8)=0.7387738651055050750031061748598307_prec
        x_quad(9)=0.9195339081664588138289326608223381_prec
        x_quad(5)=-x_quad(6)
        x_quad(4)=-x_quad(7)
        x_quad(3)=-x_quad(8)
        x_quad(2)=-x_quad(9)
        
        w_quad(6)=0.327539761183897456656510527916893_prec
        w_quad(7)=0.292042683679683757875582257374444_prec
        w_quad(8)=0.2248893420631264521194578217310478_prec
        w_quad(9)=0.1333059908510701111262271707553929_prec
        w_quad(10)=0.02222222222222222222222222222222222_prec
        w_quad(5)=w_quad(6)
        w_quad(4)=w_quad(7)
        w_quad(3)=w_quad(8)
        w_quad(2)=w_quad(9)
        w_quad(1)=w_quad(10)
        CASE DEFAULT
        WRITE(*,*) "PROBLEM IN quad_1D_lobatto",n_quad_1D
        STOP
        END SELECT

    END SUBROUTINE quad_1D_lobatto
    
    SUBROUTINE quad_1D_legendre(n_quad_1D,x_quad,w_quad)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: n_quad_1D 
        REAL(prec),DIMENSION(:),INTENT(OUT) :: x_quad,w_quad
        REAL(prec) :: aa,bb

        SELECT CASE (n_quad_1D)
        CASE (1)
        x_quad=0._prec
        w_quad=2._prec
        CASE (2)
        aa=1._prec/SQRT(3._prec)
        x_quad(1)=-aa
        x_quad(2)=aa
        w_quad=1._prec
        CASE (3)
        aa=SQRT(3._prec/5._prec)
        x_quad(1)=-aa
        x_quad(2)=0._prec
        x_quad(3)=aa
        bb=1._prec/9._prec
        w_quad(1)=5._prec*bb
        w_quad(2)=8._prec*bb
        w_quad(3)=5._prec*bb
        CASE (4)
        aa=SQRT(3._prec/7._prec-2._prec/7._prec*SQRT(6._prec/5._prec))
        bb=SQRT(3._prec/7._prec+2._prec/7._prec*SQRT(6._prec/5._prec))
        x_quad(1)=-bb
        x_quad(2)=-aa
        x_quad(3)=aa
        x_quad(4)=bb
        aa=SQRT(30._prec)
        bb=1._prec/36._prec
        w_quad(1)=(18._prec-aa)*bb
        w_quad(2)=(18._prec+aa)*bb
        w_quad(3)=(18._prec+aa)*bb
        w_quad(4)=(18._prec-aa)*bb
        CASE (5)
        aa=SQRT(5._prec-2._prec*SQRT(10._prec/7._prec))/3._prec
        bb=SQRT(5._prec+2._prec*SQRT(10._prec/7._prec))/3._prec
        x_quad(1)=-bb
        x_quad(2)=-aa
        x_quad(3)=0._prec
        x_quad(4)=aa
        x_quad(5)=bb
        aa=SQRT(70._prec)
        bb=1._prec/900._prec
        w_quad(1)=(322._prec-13._prec*aa)*bb
        w_quad(2)=(322._prec+13._prec*aa)*bb
        w_quad(3)=128._prec/225._prec
        w_quad(4)=(322._prec+13._prec*aa)*bb
        w_quad(5)=(322._prec-13._prec*aa)*bb
        CASE DEFAULT
        
        END SELECT

    END SUBROUTINE quad_1D_legendre

! ---------------------------------------------------------------
! ---------------------------------------------------------------

    FUNCTION DG_base(x,ordre_poly,LOC,ni)
        IMPLICIT NONE
        INTEGER,            INTENT(IN) :: ordre_poly
        REAL(prec),         INTENT(IN) :: x
        CHARACTER (len=8),  INTENT(IN) :: LOC
        INTEGER,            INTENT(IN) :: ni

        REAL(PREC) :: xx
        REAL(prec) :: DG_base
        INTEGER :: ii

        IF(TRIM(LOC) == "Loc")      XX = Loc_to_Ref(ni,x)
        IF(TRIM(LOC) == "Ref")      XX = x
        IF(TRIM(LOC) == "SubRef")   XX = Refsub_to_Ref(x,ni)

        DG_base = 0._prec

        IF(DG_meth == "Legendre") THEN
            DO ii =ordre_poly,1,-1
                DG_base = DG_base + coeff_DG(ordre_poly,ii) * (XX**(ii-1))
            END DO

        ELSE IF(DG_meth == "Lobatto") THEN
            DG_base = Lagrange_basis(XX,ordre_poly)

        END IF
            

    END FUNCTION DG_base
    
    FUNCTION dDG_base(x,ordre_poly, LOC,ni)
        IMPLICIT NONE
        INTEGER,            INTENT(IN) :: ordre_poly
        CHARACTER (len=8),  INTENT(IN) :: LOC
        INTEGER,            INTENT(IN) :: ni
        REAL(Prec),         INTENT(IN) :: x

        REAL(prec) :: XX
        REAL(prec) :: dDG_base
        REAL(prec) :: temp
        INTEGER :: ii,jj

        IF(TRIM(LOC) == "Loc")      XX = Loc_to_Ref(ni,x)
        IF(TRIM(LOC) == "Ref")      XX = x
        IF(TRIM(LOC) == "SubRef")   XX = Refsub_to_Ref(x,ni)

        dDG_base = 0._prec

        IF(DG_meth == "Legendre") THEN
        DO ii =ordre_poly,2,-1
            dDG_base = dDG_base + (ii-1)*coeff_DG(ordre_poly,ii) * XX**(ii-2)
        END DO
        ELSE IF(DG_meth == "Lobatto") THEN
        DO ii =1,size_base
            temp =1._prec

            DO jj =1,size_base
                IF((jj .NE. ii) .and. (jj .ne.ordre_poly)) THEN 
                temp = temp * (XX-pts_DG(jj))/(pts_DG(ordre_poly)-pts_DG(jj))
                END IF
            END DO

            IF(ii .NE. ordre_poly)  dDG_base = dDG_base + temp/(pts_DG(ordre_poly)-pts_DG(ii))
        END DO
        END IF
        

    END FUNCTION dDG_base

    FUNCTION Ref_to_Refsub(XX,n_sub)  result(ZZ)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: XX
        INTEGER,    INTENT(IN) :: n_sub
        REAL(prec)  :: ZZ

        ZZ = 2._prec * (XX-x_submiddle(n_sub))/subcell_size(n_sub)

    END FUNCTION Ref_to_Refsub

    FUNCTION Refsub_to_Ref(ZZ,n_sub)  result(XX)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: ZZ
        INTEGER,    INTENT(IN) :: n_sub
        REAL(prec)  :: XX

        XX = ZZ*subcell_size(n_sub)/2._prec + x_submiddle(n_sub)

    END FUNCTION Refsub_to_Ref

    FUNCTION Ref_to_loc(ni,XX, n_sub) result(YY)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ni
        INTEGER, INTENT(IN), optional :: n_sub
        REAL(prec), INTENT(IN) :: XX
        REAL(prec) :: xL,xR
        REAL(prec) :: YY

        IF(.not. present(n_sub)) THEN
            YY = XX*cell_size(ni)/2._prec + x_middle(ni)
        ELSE 
            YY = Refsub_to_Ref(XX,n_sub)*cell_size(ni)/2._prec + x_middle(ni)
        END IF

    END FUNCTION Ref_to_loc

    FUNCTION Loc_to_Ref(ni,YY, n_sub) result(XX)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ni
        INTEGER, INTENT(IN), optional :: n_sub
        REAL(prec), INTENT(IN) :: YY
        REAL(prec) :: XX

        IF(.not. present(n_sub)) THEN
            XX = 2._prec * (YY-x_middle(ni))/cell_size(ni)
        ELSE
            XX = 2._prec * (YY-x_middle(ni))/cell_size(ni)
            XX = Ref_to_Refsub(XX,n_sub)
        END IF

    END FUNCTION Loc_to_Ref

    SUBROUTINE Projection_Pk(fct, fct_h ,LOC ,ni ,fct_val)
        IMPLICIT NONE
        
        INTERFACE
        FUNCTION fct(YY,ni)
            USE precis
            REAL(prec),INTENT(IN) :: YY
            INTEGER,   INTENT(IN) :: ni
            REAL(prec) :: fct
        END FUNCTION fct
        END INTERFACE

        INTEGER, INTENT(IN) :: ni
        REAL(prec), DIMENSION(size_base), INTENT(OUT) :: fct_h
        CHARACTER (len=8),  INTENT(IN) :: LOC
        REAL(prec), DIMENSION(size_quad_nodes), INTENT(IN), optional :: fct_val

        REAL(prec), DIMENSION(size_base) :: f_prod
        REAL(prec) :: YY
        INTEGER :: jj, kk


        f_prod = 0._prec
        IF(.not. present(fct_val)) THEN
        DO jj =size_base,1,-1
            DO kk =1,size_quad_nodes

                IF(TRIM(LOC) == "Loc") YY = Ref_to_loc(ni,x_quad(kk))
                IF(TRIM(LOC) == "Ref") YY = x_quad(kk)

                f_prod(jj) = f_prod(jj) + fct(YY,ni)*sig_quad(jj,kk)*w_quad(kk)
            END DO
        END DO

        ELSEIF(  present(fct_val)) THEN
        DO jj =size_base,1,-1
            DO kk =1,size_quad_nodes

                f_prod(jj) = f_prod(jj) + fct_val(kk)*sig_quad(jj,kk)*w_quad(kk)

            END DO
        END DO

        END IF
        fct_h = MATMUL(Masse_inv,f_prod)

    END SUBROUTINE Projection_Pk

    SUBROUTINE Matrice_Masse_init
        IMPLICIT NONE

        INTEGER :: ii,jj,kk
        Masse =0._prec

        DO ii = 1,size_base
            DO jj = ii,size_base
                
                DO kk =size_quad_nodes,1,-1
                    Masse(ii,jj) = Masse(ii,jj) + sig_quad(ii,kk)*sig_quad(jj,kk)*w_quad(kk)
                END DO

                Masse(jj,ii) = Masse(ii,jj)
            END DO
        END DO
        
        CALL inv_mat(Masse,Masse_inv,1)


    END SUBROUTINE Matrice_Masse_init

    SUBROUTINE Matrice_Rigid_init
        IMPLICIT NONE

        INTEGER :: ii,jj

        DO ii = 1,size_base
            DO jj = 1,size_base
                Rigid(jj,ii) = quadrature(T_DG_base,ii,T_dDG_base,jj,LRef,0)
            END DO
        END DO
        
        CONTAINS
            FUNCTION T_DG_base(x,ordre_poly)
                IMPLICIT NONE
                REAL(prec), INTENT(IN) :: x
                INTEGER,    INTENT(IN) :: ordre_poly
                REAL(prec) :: T_DG_base
                T_DG_base  = DG_base(x,ordre_poly,LRef,0)
            END FUNCTION T_DG_base
            
            FUNCTION T_dDG_base(x,ordre_poly)
                IMPLICIT NONE
                REAL(prec), INTENT(IN) :: x
                INTEGER,    INTENT(IN) :: ordre_poly
                REAL(prec) :: T_dDG_base
                T_dDG_base  = dDG_base(x,ordre_poly,LRef,0)
            END FUNCTION T_dDG_base
        ! CALL inv_mat(Rigid,Rigid_inv,1)


    END SUBROUTINE Matrice_Rigid_init


! ---------------------------------------------------------------
! ---------------------------------------------------------------
    
    SUBROUTINE Projection_VF_init
        IMPLICIT NONE
        INTEGER ::i,j,m,kk
        integer :: ord
        REAL(prec) ::YY
        REAL(prec), DIMENSION(size_base)    :: phi_m

        ! Projection_VF = 0._prec
        ! DO m =1,nb_subcell
        !     CALL Projection_Pk(unit_sm,phi_m,LOC =LRef,ni= m) 
        !     write (*,fmt='(3(f10.6))')  phi_m
        !     Projection_VF(m,:) = MATMUL(Masse,phi_m)/(subcell_size(m))
        ! END DO

        Projection_VF = 0._prec
        DO j =1,nb_subcell
            DO i=1,size_base
                DO kk =size_quad_nodes,1,-1
                    YY = x_quad(kk)
                    Projection_VF(j,i) = Projection_VF(j,i) + DG_base(YY,i,LOC=LSub, ni=j)*w_quad(kk)/2._prec
                END DO
            END DO
        END DO

        CALL print_mat(Projection_VF,nb_subcell,size_base )        

        IF(nb_subcell == size_base) THEN
            CALL inv_mat(Projection_VF,Projection_VF_inv,0)
        ELSE 
            CALL inv_mat(MATMUL(Transpose(Projection_VF),Projection_VF),Projection_VF_inv,0)
            Projection_VF_inv = MATMUL(Projection_VF_inv, Transpose(Projection_VF))
        END IF

        CALL print_mat(MATMUL(Transpose(Projection_VF),Projection_VF), size_base, size_base)
        CALL print_mat(Projection_VF_inv,size_base,nb_subcell)


    END SUBROUTINE Projection_VF_init

    SUBROUTINE sub_cells_init
        IMPLICIT NONE
        INTEGER :: i,j
        REAL(prec) :: temp
        REAL(prec), DIMENSION(size_base)    :: phi
        REAL(prec), DIMENSION(nb_subcell)   :: phi_m
        REAL(prec), DIMENSION(nb_subcell,2) :: phi_val
        
        DO i =1,nb_subcell+1 
            x_subcell(i) = -1._prec + (i-1)*sub_dx
        END DO

        print *,x_subcell

        DO i =1,nb_subcell
            subcell_size(i)= x_subcell(i+1)-x_subcell(i)  
            x_submiddle(i) = x_subcell(i) + subcell_size(i)/2._prec
        END DO

        CALL Projection_VF_init
        CALL print_mat(Projection_VF,nb_subcell,size_base)


        DO i=1,nb_cell
            sol(i)%val_subcells =MATMUL(Projection_VF,sol(i)%base_poly)
        END DO

        ! phi_val = 0._prec
        ! DO i=1,nb_subcell
        !     phi_m = 0._prec; phi_m(i) = 1._prec
        !     phi = MATMUL(Projection_VF_inv,phi_m)
            
        !     print *,"phi",i,phi
            
        !     phi_val(i,1) = eval_poly(-1._prec,0, phi, LOC=LRef)
        !     phi_val(i,2) = eval_poly( 1._prec,0, phi, LOC=LRef)
        ! END DO

        phi_val(:,1) = subcell_size*MATMUL(Projection_VF,MATMUL(Masse_inv,sig_1))
        phi_val(:,2) = subcell_size*MATMUL(Projection_VF,MATMUL(Masse_inv,sig_2))

        C_m =0._prec; 
        C_p = 0._prec;
        print *,1,phi_val(:,1)
        print *,2,phi_val(:,2)
        DO i=1,nb_subcell+1
            DO j=1,i-1
                C_p(i) = C_p(i) + phi_val(j,2)
            END DO
            DO j=1,i-1
                C_m(i) = C_m(i) + phi_val(j,1)
            END DO
            C_m(i) = 1._prec -C_m(i)
        END DO

        ! C_p(1)            =0._prec; C_m(1)            =1._prec
        ! C_p(nb_subcell+1) =1._prec; C_m(nb_subcell+1) =0._prec

        print *,"Cp =",C_p
        print *,"Cm =",C_m
    
    END SUBROUTINE sub_cells_init

    FUNCTION unit_sm(x,n_sub)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: x
        INTEGER,    INTENT(IN) :: n_sub
        REAL(prec) :: xls, xrs
        REAL(prec) :: unit_sm

        xls = x_subcell(n_sub)
        xrs = x_subcell(n_sub+1)

        unit_sm = 0._prec
        if((x .LE. xrs ) .and. (x .GE. xls)) unit_sm =1._prec

    END FUNCTION unit_sm
! ---------------------------------------------------------------
! ---------------------------------------------------------------

  !----------------------------------------!
  ! INVERSION PROCEDURE OF A SQUARE MATRIX !
  !----------------------------------------!
  SUBROUTINE inv_mat(M,M1,symm)
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: symm
    REAL(prec),DIMENSION(:,:),INTENT(IN) :: M
    REAL(prec),DIMENSION(:,:),INTENT(OUT) :: M1
    INTEGER :: ndim

    !! Check if the matrix is square !!
    ndim=SIZE(M(1,:))
    IF (SIZE(M(:,1))/=ndim) THEN
       WRITE(*,*) "NON SQUARE MATRIX ",ndim,SIZE(M(:,1))
       STOP
    END IF
    
    IF (ndim==1) THEN
       M1(1,1)=1._prec/M(1,1)
    ELSE IF (ndim==2) THEN
       CALL invers_M22(M,M1)
    ELSE IF (ndim==3) THEN
       CALL invers_M33(M,M1)
    ELSE
       IF (symm==1) THEN
          CALL invers_matrix_sym(M,M1)
       ELSE
          CALL invers_matrix(M,M1)
       END IF
    END IF

  END SUBROUTINE inv_mat

  !-------------------------------------!
  ! INVERSION PROCEDURE OF A 2x2 MATRIX !
  !-------------------------------------!
  SUBROUTINE invers_M22(M,M1)
    IMPLICIT NONE
    
    REAL(prec),DIMENSION(:,:),INTENT(IN) :: M
    REAL(prec),DIMENSION(:,:),INTENT(OUT) :: M1
    REAL(prec) :: det

    det=M(1,1)*M(2,2)-M(1,2)*M(2,1)
    M1(1,1)=M(2,2);M1(1,2)=-M(1,2)
    M1(2,1)=-M(2,1);M1(2,2)=M(1,1)

    IF (ABS(det).LE.0.1_prec**prec) THEN
       WRITE(*,*) "PROBLEM IN THE 2x2 MATRIX INVERSION",ABS(det)
       WRITE(*,*) "|",M(1,1),",",M(1,2),"|"
       WRITE(*,*) "|",M(2,1),",",M(2,2),"|"
    END IF

    M1=M1/det

  END SUBROUTINE invers_M22

  !-------------------------------------!
  ! INVERSION PROCEDURE OF A 3x3 MATRIX !
  !-------------------------------------!
  SUBROUTINE invers_M33(M,M1)
    IMPLICIT NONE
    
    REAL(prec),DIMENSION(:,:),INTENT(IN) :: M
    REAL(prec),DIMENSION(:,:),INTENT(OUT) :: M1
    REAL(prec) :: det,aa,bb,cc,dd,ee,ff,gg,hh,ii

    aa=M(1,1);bb=M(1,2);cc=M(1,3)
    dd=M(2,1);ee=M(2,2);ff=M(2,3)
    gg=M(3,1);hh=M(3,2);ii=M(3,3)

    det=aa*ee*ii+dd*hh*cc+gg*bb*ff&
     & -gg*ee*cc-aa*hh*ff-ii*bb*dd

    M1(1,1)=ee*ii-hh*ff;M1(1,2)=cc*hh-bb*ii;M1(1,3)=bb*ff-ee*cc
    M1(2,1)=ff*gg-dd*ii;M1(2,2)=aa*ii-cc*gg;M1(2,3)=dd*cc-aa*ff
    M1(3,1)=dd*hh-ee*gg;M1(3,2)=gg*bb-hh*aa;M1(3,3)=aa*ee-bb*dd

    IF (ABS(det).LE.0.1_prec**prec) THEN
       WRITE(*,*) "PROBLEM IN THE 3x3 MATRIX INVERSION",ABS(det)
       WRITE(*,*) "|",M(1,1),",",M(1,2),",",M(1,3),"|"
       WRITE(*,*) "|",M(2,1),",",M(2,2),",",M(2,3),"|"
       WRITE(*,*) "|",M(3,1),",",M(3,2),",",M(3,3),"|"
    END IF

    M1=M1/det

  END SUBROUTINE invers_M33

  !-------------------------------------------!
  ! INVERSION PROCEDURE OF A SYMMETRIC MATRIX !
  !-------------------------------------------!
  SUBROUTINE invers_matrix_sym(M,M1)
    IMPLICIT NONE
    
    REAL(prec),DIMENSION(:,:),INTENT(IN) :: M
    REAL(prec),DIMENSION(:,:),INTENT(OUT) :: M1
    REAL(prec),DIMENSION(size(M(1,:)),size(M(1,:))) :: L,L1
    INTEGER :: n,i,j,k
    
    
    n=SIZE(M(1,:))
    M1=0._prec;L1=0._prec
    
    CALL cholesky(M,L)
    
    !$OMP PARALLEL NUM_THREADS(num_thrd)
    !$OMP DO PRIVATE(i,k) 
    DO j=1,n
       L1(j,j)=1._prec/L(j,j)
       DO i=j+1,n
          DO k=j,i-1
             L1(i,j)=L1(i,j)-L(i,k)*L1(k,j)
          END DO
          L1(i,j)=L1(i,j)/L(i,i)
       END DO
    END DO
    !$OMP END DO

    !$OMP DO PRIVATE(k) COLLAPSE(2)
    DO i=1,n
       DO j=1,n
          DO k=MAX(i,j),n
             M1(i,j)=M1(i,j)+L1(k,i)*L1(k,j)
          END DO
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
  END SUBROUTINE invers_matrix_sym
  
  !-------------------------------------------------------------------!
  ! L.L^t DECOMPOSITION FOR a POSITIVE SEMI-DEFINITE SYMMETRIC MATRIX !
  !-------------------------------------------------------------------!
  SUBROUTINE cholesky(M,L)
    IMPLICIT NONE
    
    REAL(prec),DIMENSION(:,:),INTENT(IN) :: M
    REAL(prec),DIMENSION(:,:),INTENT(OUT) :: L
    INTEGER :: n,i,j,k
    
    
    L=0._prec
    n=SIZE(M(1,:))
    DO i=1,n

       L(i,i)=M(i,i)
       DO k=1,i-1
          L(i,i)=L(i,i)-L(i,k)**2
       END DO
       L(i,i)=SQRT(L(i,i))

       DO j=i+1,n
          L(j,i)=M(j,i)
          DO k=1,i-1
             L(j,i)=L(j,i)-L(i,k)*L(j,k)
          END DO
          L(j,i)=L(j,i)/L(i,i)
       END DO

    END DO

  END SUBROUTINE cholesky

  !---------------------------------!
  ! INVERSION PROCEDURE OF A MATRIX !
  !---------------------------------!
  SUBROUTINE invers_matrix(M,M1)
    IMPLICIT NONE
    
    REAL(prec),DIMENSION(:,:),INTENT(IN) :: M
    REAL(prec),DIMENSION(:,:),INTENT(OUT) :: M1
    REAL(prec),DIMENSION(size(M(1,:)),size(M(1,:))) :: L,L1,U,U1
    INTEGER :: nn,i,j,k    
    
    nn=SIZE(M(1,:))
    M1=0._prec
    L1=0._prec
    U1=0._prec
    
    CALL decomp_LU(M,L,U)
    
    !$OMP PARALLEL NUM_THREADS(num_thrd)
    !$OMP DO PRIVATE(i,k)
    DO j=1,nn
       L1(j,j)=1._prec!/L(j,j)
       U1(j,j)=1._prec/U(j,j)
       DO i=j+1,nn
          DO k=j,i-1
             L1(i,j)=L1(i,j)-L(i,k)*L1(k,j)
             U1(j,i)=U1(j,i)-U1(j,k)*U(k,i)
          END DO
          !L1(i,j)=L1(i,j)/L(i,i)
          U1(j,i)=U1(j,i)/U(i,i)
       END DO
    END DO
    !$OMP END DO

    !$OMP DO PRIVATE(k) COLLAPSE(2)
    DO i=1,nn
       DO j=1,nn
          DO k=MAX(i,j),nn
             M1(i,j)=M1(i,j)+U1(i,k)*L1(k,j)
          END DO
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
  END SUBROUTINE invers_matrix

  !-------------------------------!
  ! L.U DECOMPOSITION (DoolittLE) !
  !-------------------------------!
  SUBROUTINE decomp_LU(M,L,U)
    IMPLICIT NONE
    
    REAL(prec),DIMENSION(:,:),INTENT(IN) :: M
    REAL(prec),DIMENSION(:,:),INTENT(OUT) :: L,U
    INTEGER :: nn,i,j,k
    
    
    L=0._prec;U=0._prec
    nn=size(M(1,:))
    DO i=1,nn

       L(i,i)=1._prec
       U(i,i)=M(i,i)
       DO k=1,i-1
          U(i,i)=U(i,i)-L(i,k)*U(k,i)
       END DO

       DO j=i+1,nn
          U(i,j)=M(i,j)
          L(j,i)=M(j,i)
          DO k=1,i-1
             U(i,j)=U(i,j)-L(i,k)*U(k,j)
             L(j,i)=L(j,i)-L(j,k)*U(k,i)
          END DO
          IF (ABS(U(i,i)).LT.eps0) THEN
             WRITE(*,*) "ERROR DECOMPOSITION LU"
          ELSE
             L(j,i)=L(j,i)/U(i,i)
          END IF
       END DO

    END DO

  END SUBROUTINE decomp_LU

END MODULE mod_polynome