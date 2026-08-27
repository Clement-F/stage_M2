MODULE mod_SolIni
    USE mod_polynome
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
        creneau = eps0
        if(-0.5_prec<x .and. x<eps0) creneau = 0.1_prec
    END FUNCTION creneau

    FUNCTION composite(x,ni) result(res)
        IMPLICIT NONE
        REAL(prec), INTENT(in) :: x
        INTEGER,    INTENT(in) :: ni
        REAL(prec) :: res 

        REAL(prec), parameter :: delta =0.005_prec, alpha = 10._prec ,z=-0.7_prec,a=0.5_prec
        REAL(prec), parameter :: beta =LOG(2._prec)/(36._prec*delta**2)

        res = eps0
        
            IF((-0.8_prec .LT. x )  .AND. (x .LT. -0.6_prec) ) THEN
        res = 1._prec/6._prec *(G(x,beta,z-delta) +G(x,beta,z+delta)+4*G(x,beta,z))
        
        ELSEIF((-0.4_prec .LT. x )  .AND. (x .LT. -0.2_prec) ) THEN
        res = 1._prec
        
        ELSEIF((  0._prec .LT. x )  .AND. (x .LT.  0.2_prec) ) THEN
        res = 1-abs(10*(x-0.1_prec))
        
        ELSEIF(( 0.4_prec .LT. x )  .AND. (x .LT.  0.6_prec) ) THEN
        res = 1._prec/6._prec *(F(x,alpha,a-delta) +F(x,alpha,a+delta)+4*F(x,alpha,a))

        END IF

        CONTAINS
        FUNCTION F(x,alpha,a)
            IMPLICIT NONE
            REAL(prec), INTENT(IN) :: x,alpha,a
            REAL(prec) :: F
            F = sqrt(max(1-alpha**2 *(x-a)**2, 0._prec))
        END FUNCTION

        FUNCTION G(x,beta,z)
            IMPLICIT NONE
            REAL(prec), INTENT(IN) :: x,beta,z
            REAL(prec) :: G
            G = exp(-beta* (x-z)**2)
        END FUNCTION

    END FUNCTION composite

    FUNCTION Riemann_Buckley(x,ni) result(s)
        IMPLICIT NONE
        REAL(prec), INTENT(in) :: x
        INTEGER,    INTENT(in) :: ni
        REAL(prec) :: s 
        s = eps0
        IF(x<0._prec) THEN 
            s = -3._prec
        ELSE 
            s = 3._prec
        END IF
    END FUNCTION Riemann_Buckley

    FUNCTION sod_tube(x,ni)
        IMPLICIT NONE
        REAL(prec), INTENT(IN):: x
        INTEGER, INTENT(IN) :: ni
        REAL(prec), DIMENSION(nb_var) :: sod_tube 

        sod_tube = eps0
        IF(x .LT. 0.3_prec) THEN
            sod_tube(1) = 1._prec
            sod_tube(2) = 0.75_prec
            sod_tube(3) = 1._prec/(gamma_iso-1._prec) + 0.5_prec* (0.75_prec)**2 
        ELSE
            sod_tube(1) = 0.125_prec
            sod_tube(2) = eps0
            sod_tube(3) = 0.1_prec/(gamma_iso-1._prec) 
        END IF

    END FUNCTION sod_tube

    FUNCTION smooth_isentropique(x,ni)
        IMPLICIT NONE
        REAL(prec), INTENT(IN):: x
        INTEGER, INTENT(IN) :: ni
        REAL(prec), DIMENSION(nb_var) :: smooth_isentropique 

        smooth_isentropique = eps0

        smooth_isentropique(1) = 1._prec  +  0.999999999_prec* sin(pi*x)
        smooth_isentropique(3) = ((1._prec  +  0.99999999_prec* sin(pi*x))**(gamma_iso))/(gamma_iso-1._prec)

    END FUNCTION smooth_isentropique

    FUNCTION acoustic_wave(x,ni) result(res)
        IMPLICIT NONE
        REAL(prec), INTENT(IN):: x
        INTEGER, INTENT(IN) :: ni
        REAL(prec), DIMENSION(nb_var) :: res 

        res = eps0
        IF(x .LT. -4._prec) THEN
            res(1) = 3.857143_prec
            res(2) = 2.629369_prec *  3.857143_prec
            res(3) = 10.333333_prec/(gamma_iso-1._prec)  + 0.5_prec* (2.629369_prec)**2 *  3.857143_prec
        ELSE
            res(1) = 1._prec + 0.2_prec* sin(5._prec*x)
            res(2) = eps0
            res(3) = 1._prec/(gamma_iso-1._prec)
        END IF

    END FUNCTION acoustic_wave

    FUNCTION Blast(x,ni) result(res)
        IMPLICIT NONE
        REAL(prec), INTENT(IN):: x
        INTEGER, INTENT(IN) :: ni
        REAL(prec), DIMENSION(nb_var) :: res 

        res = eps0
        
        res(1) = 1._prec
        res(2) = eps0
        
        IF(x .LT. 0.1_prec)                      res(3) = Real(1D3 ,prec)/(gamma_iso-1._prec) 
        IF(0.1_prec .LT. x .AND. x .LT. 0.9_prec)res(3) = Real(1D-2,prec)/(gamma_iso-1._prec)  
        IF(0.9_prec .LT. x .AND. x .LT. 1._prec) res(3) = Real(1D2 ,prec)/(gamma_iso-1._prec) 

    END FUNCTION Blast


    FUNCTION Q_init(x,ni,nb_var)
        IMPLICIT NONE
        REAL(prec), INTENT(IN) :: x
        INTEGER,    INTENT(IN) :: ni, nb_var
        REAL(prec), DIMENSION(nb_var) :: Q_init

        SELECT CASE(TRIM(sol_ini_name))
        CASE("sinus")
        Q_init = sinus(x,ni)

        CASE("unit")
            Q_init = unit(x,ni)

        CASE("Riemann")
            Q_init = 0._prec
            if(x<0._prec) Q_init =1._prec

        CASE("creneau")
            Q_init = creneau(x,ni)

        CASE("Riemann_Buckley")
            Q_init = Riemann_Buckley(x,ni)

        CASE("composite")
            Q_init = composite(x,ni)

        CASE("Burgers_choc")
            Q_init = 0._prec
            if(0.3_prec<x .and. x<0.7_prec) Q_init = -1._prec
            if(x>0.7_prec) Q_init = 0.5_prec

        CASE("Sod")
            Q_init = sod_tube(x,ni)
            print *, Q_init

        CASE("isentropique")
            Q_init = smooth_isentropique(x,ni)

        CASE("acoustic_wave")
            Q_init = acoustic_wave(x,ni)
        CASE("Blast")
            Q_init = Blast(x,ni)

        CASE DEFAULT
        WRITE(*,*) "Sol initiale non reconnue"
        STOP
        END SELECT

    END FUNCTION Q_init


END MODULE mod_SolIni