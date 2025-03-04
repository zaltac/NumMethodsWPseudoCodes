PROGRAM Test_Gauss_Legendre_Quad
! ==============================================================================
!  The main program to test SUBROUTINE Gauss_Legendre_Quad
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(90) :: x, w
DOUBLE PRECISION :: eps
INTEGER :: n, i
!
PRINT*, "Enter n"
READ *, n
!
eps= 1.0d-8
!     
CALL Gauss_Legendre_Quad(n,eps,x,w)
!
PRINT*, "========= Gauss-Legendre Quads ========="
WRITE(*,1)
DO i=1,n
   WRITE(*,2) i, x(i), w(i)
ENDDO 
1 FORMAT(5x," i",11x,"x_i",16x,"w_i",/,5x,3('-'),2x,17('-'),2x,17('-'))
2 FORMAT(5x,i2,1x,F18.12,1x,F18.12)
END PROGRAM Test_Gauss_Legendre_Quad

SUBROUTINE Gauss_Legendre_Quad(n,eps,x,w)

!  ==================================================================================
!  CODE8.6-GAUSS_LEGENDRE_QUAD.f95. A Fortran95 module implementing Pseudocode 8.6.               
!  
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
!  First Edition. (c) By Zekeriya ALTAÇ (2024).
!  ISBN: 978-1-032-75474-1 (hbk)
!  ISBN: 978-1-032-75642-4 (pbk)
!  ISBN: 978-1-003-47494-4 (ebk)
!  
!  DOI : 10.1201/9781003474944
!  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
!  
!  This free software is complimented by the author to accompany the textbook.
!  E-mail: altacz@gmail.com.

!  DESCRIPTION: A subroutine to generate N-point Gauss-Legendre quadrature                     
!    abscissas and weights on [-1,1].                                                          
!                                                                                              
!  ON ENTRY                                                                                    
!     n   :: Number of quadrature points;                                                      
!     eps :: Tolerance, i.e., desired level of numerical accuracy.                             
!                                                                                              
!  ON EXIT                                                                                     
!     x   :: Array of length N containing the abscissas;                                       
!     w   :: Array of length N containing the weights.                                         
!                                                                                              
!  USES                                                                                        
!    DABS  :: Built-in Intrinsic function returning the absolute value of a real value;         
!    DCOS  :: Built-in Intrinsic function returning trigonometric cosine value;                          
!   DFLOAT :: A built-in intrinsic function that converts an integer argument to a double float.  
!                                                                                              
!  REVISION DATE :: 03/04/2025                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN) :: eps
DOUBLE PRECISION, INTENT(OUT), DIMENSION(n) :: x, w
DOUBLE PRECISION :: u, del, P0, P1, P2, dP, pi
INTEGER :: m, i, k
DATA pi/3.1415926535897932385d0/ 
!     
m = (n+1)/2  
DO i = 1, m
    u = DCOS( pi * DFLOAT(4*i-1) / DFLOAT(4*n+2) )
!========= start Newton-Raphson iterations 
   del= 1.0d0    
   DO WHILE (DABS(del) > eps)
      P0 = 1.0d0
      P1 = u
      DO k = 2, n
         P2 = (DFLOAT(2*k-1) * u * P1 - DFLOAT(k-1) * P0) / DFLOAT(k)
         P0 = P1 
         P1 = P2              
      ENDDO
      dP = DFLOAT(n) * ( u * P2- P0 )/(u * u - 1.0d0)
      del= P2 / dP
       u = u - del
   END DO 
   x(i) =-u
   w(i) = 2.0d0/(1.0d0 - u * u) / dP**2
   x(n+1-i) = u
   w(n+1-i) = w(i)
ENDDO  
END SUBROUTINE Gauss_Legendre_Quad
