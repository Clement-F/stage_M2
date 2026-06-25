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

END MODULE mod_Divers