PROGRAM Test_Simpsons_Rule_DF
! ==============================================================================
!  The main program to test SUBROUTINE Simpsons_Rule_DF
! ==============================================================================
IMPLICIT NONE
INTEGER :: n, i
DOUBLE PRECISION :: f(0:99), a, b, x, h, intg 
!
PRINT*, "Enter Number of Panels "
read*, n
! *** Construct a discrete function using f(x)=x^4
a = 0.0d0
b = 1.0d0
h = (b-a) / dfloat(n)
x = a
DO i=0,n  ! Discrete dataset generation from y=FUNC(x)
   f(i) = x**4  
      x = x + h
END DO
!
CALL Simpsons_Rule_DF(n,h,f,intg)
!
PRINT*, n, " panel Simpson's 1/3 Rule = ", intg
PRINT*, " "
!    
END PROGRAM Test_Simpsons_Rule_DF

SUBROUTINE Simpsons_Rule_DF(n,h,f,intg)
!  ==================================================================================
!  CODE8.2-Simpsons_Rule_DF.f95. A Fortran95 module implementing Pseudocode 8.2.                  
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
!
!  DESCRIPTION: A module to estimate the integral of a discrete function f on [a,b]                        
!   using the Simpson's 1/3 rule.                                                              
!                                                                                              
!  ON ENTRY                                                                                    
!     n   :: Number of panels (must be even!..);                                               
!     h   :: Interval size (uniform spacing, x_{i+1}-x_i);                                     
!     f   :: Array of length (n+1) containing the ordinates, f_0, f_1, ..., f_n.               
!                                                                                              
!  ON EXIT                                                                                     
!   intg  :: Numerical estimate for the integral.                                              
!                                                                                              
!  USES                                                                                        
!    MOD  :: Modulo function, returning the remainder after number is divided by divisor.      
!                                                                                              
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN) :: h, f(0:n)
DOUBLE PRECISION, INTENT(OUT) :: intg
INTEGER :: i, m
DOUBLE PRECISION :: odd, even
! ------- Standard Simpson's Rule           
m = MOD(n,2)
IF(m/=0) THEN
  print *, "Number of panels is not EVEN"
  STOP
END IF
odd=0.0d0
DO i=1,n-1,2
   odd = odd + f(i)
ENDDO
Even=0.0d0
DO i=2,n-2,2
   even = even + f(i)
ENDDO
intg = f(0) + f(n) + 4.0d0 * Odd + 2.0d0 * Even
intg = intg * h / 3.0d0
END SUBROUTINE Simpsons_Rule_DF
        
 
      
