PROGRAM conv_rate
  IMPLICIT NONE

  REAL(kind=8) :: Ex_L1,Ex2_L1,odr_L1
  REAL(kind=8) :: Ex_L2,Ex2_L2,odr_L2

  PRINT*, "Enter L1 and L2 errors for Dx:"
  READ*, Ex_L1,Ex_L2

  PRINT*, "Enter L1 and L2 errors for Dx/2:"
  READ*, Ex2_L1,Ex2_L2

  odr_L1=LOG(Ex_L1/Ex2_L1)/LOG(2.D0)
  odr_L2=LOG(Ex_L2/Ex2_L2)/LOG(2.D0)

  PRINT*, "L1 Order =",odr_L1
  PRINT*, "L2 Order =",odr_L2

END PROGRAM conv_rate
