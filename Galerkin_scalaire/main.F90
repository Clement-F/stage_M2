PROGRAM MAIN

   use mod_polynome
   implicit none

    REAL(prec) :: t
    INTEGER :: i

    xL = -1.5_prec ; xR = 2._prec
    order_x =3
    nb_cell =100

    print *, "init"

    ALLOCATE(x_quad(order_x+1), w_quad(order_x+1))
    ALLOCATE(coeff_DG(order_x+1, order_x+1))
    ALLOCATE(x_cell(nb_cell+1))
    ALLOCATE(coeff_legendre(10,10))

    dx = REAL((xR-xL)/nb_cell, prec)

    DO i =1,nb_cell+1 
        x_cell(i) = xL + (i-1)*dx
    END DO


    print * ,"alloc"

    DG_meth = "Legendre"
    quad_meth = "Legendre"
    CALL Coeff_quad_init
    CALL Coeff_DG_init

    print *, "poly"
    
    print *, x_quad,w_quad, coeff_DG

    CALL Matrice_Masse_init

    print *, "masse"

END PROGRAM