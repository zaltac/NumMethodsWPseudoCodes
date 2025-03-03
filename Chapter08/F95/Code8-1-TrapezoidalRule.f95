PROGRAM Test_TrapezoidalRule
! ==============================================================================
!  The main program to test SUBROUTINE Trapezoidal_Rule_RF
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION :: a, b, intg, intgc
INTEGER :: n
!
a=0.0d0
b=1.0d0
PRINT*, "Enter Number of Panels "
READ *, n
!
CALL Trapezoidal_Rule_RF(n,a,b,intg,intgc)
!
PRINT 1, n,intg, intgc
1 FORMAT(5x," Panel no. =",i3,/, &
         5x," Estimate = ",f20.14,4x,"Trapezoidal rule",/, &
         5x," Estimate = ",f20.14,4x,"Trapezoidal rule with end-correction")
END PROGRAM Test_TrapezoidalRule

SUBROUTINE Trapezoidal_Rule_RF(n,a,b,intg,intgc)
!  ==================================================================================
!  CODE8.1-Trapezoidal_Rule_RF.f95. A Fortran95 module implementing Pseudocode 8.1.               
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

!  DESCRIPTION: A subroutine to estimate the integral of y=f(x) on [a,b]                       
!   using the Trapezoidal rule with/without end correction.                                    
!                                                                                              
!  ON ENTRY                                                                                    
!     n   :: Number of panels (i.e., n+1 integration points);                                  
!   [a,b] :: Integration interval.                                                             
!                                                                                              
!  ON EXIT                                                                                     
!   intg   :: Integral estimate using the ordinary Trapezoidal rule;                           
!   intgc  :: Integral estimate using the Trapezoidal rule with the end-point correction.      
!                                                                                              
!  ALSO REQUIRED                                                                               
!     FX   :: User-defined external function providing the function, f(x);                     
!     FU   :: User-defined external function providing the first derivative, f'(x).            
!                                                                                              
!  USES                                                                                        
!    DFLOAT :: A built-in intrinsic function that converts an integer argument to a real value. 
!                                                                                              
!  REVISION DATE :: 03/03/2025                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN) :: a, b 
DOUBLE PRECISION, INTENT(OUT) :: intg, intgc
DOUBLE PRECISION :: xi, h, corr, FX, FU
INTEGER :: i
! ------- Standard Trapezoidal Rule  
h = (b-a) / DFLOAT(n)
intg = 0.50d0 * (FX(a) + FX(b))      
xi = a
DO i = 1, n-1
   xi = xi + h
   intg = intg + FX(xi)
ENDDO
intg= h * intg 
! ------- Trapezodial Rule with End-Correction  
corr = -h*h*(FU(b) - FU(a)) / 12.0d0
intgc= intg + corr
END SUBROUTINE Trapezoidal_Rule_RF
        
DOUBLE PRECISION FUNCTION FX(x)
! ==============================================================================
! DESCRIPTION: User-defined function providing y=f(x) to be integrated. 
!
! ARGUMENTS:
!      x   :: a real input value.
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x
  FX=x**2
END FUNCTION FX    

DOUBLE PRECISION FUNCTION FU(x)
! ==============================================================================
! DESCRIPTION: User-defined function providing first derivative, f'(x), explicitly. 
!
! ARGUMENTS:
!      x   :: a real input value.
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x
  FU=2.0d0*x
END FUNCTION FU

      
