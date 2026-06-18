program test1
    IMPLICIT NONE
    CHARACTER(len =2) c1,c2

    c1 = "ab"; c2 ="a "

    print *,c1,"/"
    print *,c2,"/"
    print *,TRIM(c1),"/"
    print *,TRIM(c2),"/"

end program test1