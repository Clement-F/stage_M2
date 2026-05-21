program test1
    implicit none
    INTEGER, DIMENSION(:), POINTER :: p1,p2
    INTEGER :: a1,a2
    
    ALLOCATE(p1(1),p2(1))

    print *,"init"
    a1 = 1; a2= 2;
    p1 = a1; p2 = a2;
    print *, p1,p2,'/', a1,a2
    a1 = 0; 
    print *, p1,p2,'/', a1,a2
    a1 = 1;
    p1 = p2;
    print *, p1,p2,'/', a1,a2
    p2 = 0;
    print *, p1,p2,'/', a1,a2

end program test1