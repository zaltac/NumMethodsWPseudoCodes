PROGRAM Test_Explicit_Euler
! ==============================================================================
!  The main program to test SUBROUTINE Explicit_Euler
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION :: x0, y0, xlast, h
! 
! *** initialize problem
!   
  x0  = 0.0d0
  y0  = 2.0d0 
 xlast= 1.0d0  
    h = 0.1d0 
!
 CALL Explicit_Euler(h,x0,y0,xlast)
!
END PROGRAM Test_Explicit_Euler



SUBROUTINE Explicit_Euler(h,x0,y0,xlast)
!  ==================================================================================
!  CODE9.1-Explicit_Euler.f95. A Fortran95 module implementing Pseudocode 9.1.                    
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

!  DESCRIPTION: A subroutine to estimate the solution of a first order IVP on [x0,xlast]           
!    using the explicit Euler method. Numerical estimates are printed out, not stored.         
!    With minor modifications, this PROGRAM can also be used to solve explicit methods         
!    such as MIDPOINT RULE and MODIFIED EULER.                                                 
!                                                                                              
!  ON ENTRY                                                                                    
!   h     :: Step size (it must be uniform);                                                   
!   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
!   xlast :: End point of the solution interval.                                               
!                                                                                              
!  Other Internal Variables                                                                             
!   x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
!                                                                                              
!  USES                                                                                        
!    DABS :: Built-in Intrinsic function returning the absolute value of a real value.                                                                                           
!    FCN  :: User-defined external function providing y'=f(x,y).                               
!                                                                                              
!  REVISION DATE :: 03/05/2025                                                                 
!  ==================================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT( IN) :: h, xlast
DOUBLE PRECISION, INTENT(INOUT) :: x0, y0
DOUBLE PRECISION :: FCN, exact, x, y
!    
WRITE(*,1)
WRITE(*,2) x0, y0   
x = x0
DO WHILE (x<xlast)  
   x = x0 + h 
   y = y0 + h * FCN(x0, y0)      ! Explicit-Euler
   WRITE(*,2) x, y, dabs(exact(x)-y)
   x0 = x 
   y0 = y
END DO
1 FORMAT(8x,"x",13x,"y",11x,"Abs. Error",/,4x,42("-"))
2 FORMAT(2x,2(f12.7,2x),1PE14.3)
END SUBROUTINE Explicit_Euler      


DOUBLE PRECISION FUNCTION fcn(x,y)
! ==============================================================================
! DESCRIPTION: A function subprogram providing y'=f(x,y)
!
! ARGUMENTS:
!      x, y  :: Real input values.
!
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x, y
 FCN =-y/(x+1.0d0)
END FUNCTION fcn


DOUBLE PRECISION FUNCTION exact(x)
! ==============================================================================
! DESCRIPTION: A function subprogram providing the true solution y=(x) for 
!    testing the module. 
!
! ARGUMENTS:
!      x   :: A real input, independent variable.
!
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x
  exact = 2.0d0/(x+1.0d0)
END FUNCTION exact
