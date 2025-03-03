PROGRAM TEST_LINEARIZE_REGRESS
! ==============================================================================
!  The main program to test SUBROUTINE LINEARIZE_REGRESS
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER ::  n=6
DOUBLE PRECISION :: x(n), y(n), a0, b0, E, S, r2
INTEGER :: model

x = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
y = (/ 5.3d0, 6.5d0, 6.8d0, 7.1d0, 7.8d0, 7.5d0 /)       

model=1  
CALL LINEARIZE_REGRESS(n,x,y,model,a0,b0,E,S,r2)

PRINT 1, a0, b0, E, S, r2
1 FORMAT(10x,/"******* Best-Fit Coefficients and Parameters ********",/, &
  3x,"  a =",F12.6,4x," b =",F12.6,/, &
  3x,"  E =",F12.6,4x," S =",F12.6,3x," r-squared =",F12.6/)
END PROGRAM TEST_LINEARIZE_REGRESS

SUBROUTINE LINEARIZE_REGRESS(n,x,y,model,a0,b0,E,S,r2)
!  ==================================================================================
!  CODE7.2-LINEARIZE_REGRESS.f95. A Fortran95 module implementing Pseudocode 7.2.                 
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

!  DESCRIPTION: A module to obtain least-squares-fit to the power model. Users                 
!    can likewise incorporate the other linearizable models into the module.                   
!                                                                                              
!  ON ENTRY                                                                                    
!     n   :: The number of data in the set;                                                    
!    x,y  :: Arrays of length n containing the data;                                           
!   model :: Model flag, Model = 1 corresponds to the power model, Y= a0*x^b0.                 
!                                                                                              
!  ON EXIT                                                                                     
!    a0,b0:: Model parameters;                                                                 
!     E   :: Sum of the Squares of Residuals (SSR);                                            
!     S   :: Sum of the Squares of Mean Deviation (SSMD);                                      
!    r2   :: r-squared, coefficient of determination.                                          
!                                                                                              
!  USES                                                                                        
!     EXP :: Built-in Intrinsic function returning exponential of a real value, e^x.           
!     LOG :: Built-in Intrinsic function returning the natural log of a real value.            
!                                                                                              
!  REVISION DATE :: 03/03/2025                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n, model
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n)  :: x, y
DOUBLE PRECISION, INTENT(OUT) :: a0, b0, E, S, r2
DOUBLE PRECISION :: c(2,2), b(2), xx(n), yy(n)
DOUBLE PRECISION :: DD, D1, D2, yavg, yi, xk
INTEGER :: i, j, k, p
!
DO k = 1, n
   xx(k) = x(k)
   yy(k) = y(k)    ! save original data on xx & yy to compute E, S, r2.
   IF(model==1) THEN
      x(k) = dlog(x(k))
      y(k) = dlog(y(k))
   ELSE IF( (model<1).OR.(model>1)) THEN
      print*, "Undefined model..."
      EXIT
   ENDIF
END DO
! *** Construct the 2x2 linear system 
DO i = 1, 2
   b(i) = 0.0d0 
   DO j = 1, 2
      c(i,j) = 0.0d0
      DO k = 1, n
         p = i + j - 2
         xk= 1.0d0
         IF(p/=0) xk = x(k)**p
         c(i,j) = c(i,j) + xk
      ENDDO
   ENDDO 
   DO k = 1, n
      p = i - 1
      xk= 1.0d0
      IF(p/=0) xk = x(k)**p
      b(i) = b(i) + xk * y(k)        
   ENDDO
ENDDO              
! Solve 2x2 system using the Cramer's Rule.
dd = c(1,1)*c(2,2)-c(2,1)*c(1,2)
d1 = b(1)*c(2,2)-b(2)*c(1,2)
d2 = b(2)*c(1,1)-b(1)*c(2,1)
! 
a0 = d1 / dd
b0 = d2 / dd
!
yavg =0.0d0
DO k = 1, n 
   yavg= yavg + yy(k)
ENDDO
!
yavg= yavg / dfloat(n) ! avgerage of the original data
S = 0.0d0
E = 0.0d0
!
IF(model<=2) THEN
   a0 = DEXP(a0)
ELSE 
   PRINT*, "Undefined model..."
   STOP 
ENDIF
! Compute S & E for the original data
DO k = 1, n
   S = S + (yy(k) - yavg)**2
   yi = a0 * xx(k)**b0
   E = E + (yi - yy(k))**2
ENDDO
r2 = 1.0d0 - E / S
END SUBROUTINE LINEARIZE_REGRESS
