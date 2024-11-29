PROGRAM Test_Richardson
! ==============================================================================
!  The main program to test SUBROUTINE Richardson
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(0:10,0:10) :: D
DOUBLE PRECISION :: x0, h, err, eps, deriv
INTEGER :: k, m, nr 
!      
x0  = 150.0d0
 h  =  5.0d0 
err =  1.0d0
eps =  1.0D-06

CALL Richardson(x0,h,eps,D,nr,deriv)
!
DO k=0,nr
   WRITE( *,1) k,(D(k,m),m=0,k)
ENDDO 
print*, " ---------------------------------"
print*, "Derivative is=",deriv
print*, " ---------------------------------"
1 FORMAT(1x,i2,1x,5(G19.11))
END PROGRAM Test_Richardson

 
SUBROUTINE Richardson(x0,h,eps,D,nr,deriv)
!  ==================================================================================
!  CODE5-2-RICHARDSON.f95. A Fortran95 module implementing Pseudocode 5.2.                        
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

!  DESCRIPTION: A subroutine to compute the first derivative of an explicitly                      
!      defined function using Richardson's extrapolation.                                    
!                                                                                              
!  ON ENTRY                                                                                    
!      x0  :: Point at which derivative is to be computed;                                     
!      h   :: Initial interval size;                                                           
!      eps :: Tolerance desired.                                                               
!                                                                                              
!  ON EXIT                                                                                     
!      D   :: A matrix containing the Richardson's table (0..n, 0..n)                        
!      nr  :: Size of the table;                                                               
!    deriv :: Estimated derivative.                                                            
!                                                                                              
!  USES                                                                                        
!    DABS  :: Built-in Intrinsic function returning the absolute value of a DP value;        
!   DFLOAT :: A built-in intrinsic function that converts an integer argument to a DP value. 
!                                                                                              
!  ALSO REQUIRED                                                                               
!     FUNC  :: User-defined external function providing the nonlinear equation.                
!                                                                                              
!  REVISION DATE :: 06/13/2024                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(OUT) :: nr
DOUBLE PRECISION, INTENT(IN) :: x0, eps
DOUBLE PRECISION, INTENT(INOUT) :: h
DOUBLE PRECISION, INTENT(OUT), DIMENSION(0:10,0:10) :: D
DOUBLE PRECISION, INTENT(OUT) :: deriv
DOUBLE PRECISION :: err, FUNC
INTEGER :: k, m 
!       
 m = 0
 k = 0
err= 1.0d0
!
DO WHILE (err>eps)
   D(k,0) = (FUNC(x0+h)-FUNC(x0-h)) / (2.0d0*h)         ! 1st derivative
   DO m=1,k
      D(k,m)= (4**m * D(k,m-1) - D(k-1,m-1)) / DFLOAT(4**m-1)
   END DO
   IF( k>=1 ) THEN  ! Estimate diagonalwise differentiation error
       err = DABS(D(k,k) - D(k-1,k-1))
   ENDIF
   h = h/2.0d0
   k = k + 1
END DO        
   nr=k-1
deriv=D(nr,nr)
END SUBROUTINE Richardson
        
!  ==================================================================================
!  USER-DEFINED FUNCTION "FUNC" OF ONE-VARIABLE
!  ==================================================================================
DOUBLE PRECISION FUNCTION FUNC(x)
DOUBLE PRECISION, INTENT(IN) :: x
 FUNC = 25000.d0/(-57.0d0 + x) - 5.2d6/x**2
END FUNCTION FUNC
 
