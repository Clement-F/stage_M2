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

        if(abs(fct(xf)) .LT. 1e-20) then 
            dicho = xf
            return
        else if(abs(fct(xd)).LT.1e-20) then
            dicho = xd
            return
        END if


        a=xd; b=xf; 

        n =0
        do while(abs(fct(c))>1e-20 .and. n<50 )
            c=(a+b)/2._prec
            t = sign(1._prec,fct(c)*fct(a))
            if(t <0._prec) then 
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

        minmod = min(0.,max(x,y,z)) + max(0.,min(x,y,z))

    END FUNCTION
END MODULE mod_Divers