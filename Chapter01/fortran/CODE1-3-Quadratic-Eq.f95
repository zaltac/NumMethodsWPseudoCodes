PROGRAM TestQuadraticEq
! ==============================================================================
!  The main program to test SUBROUTINE QUADRATIC_EQ
! ==============================================================================
IMPLICIT NONE 
DOUBLE PRECISION :: p, q, xr(2), xi(2) 

WRITE(*,"(A26)") "Type in values for p and q" 
READ(*,*) p, q 

! Solve quadratic equation
CALL QUADRATIC_EQ(p,q,xr,xi)

! Print results
WRITE(*,"(A9,F9.4,1X,A1,F9.4,A2)") &   
   "1st Root ", xr(1), "+", xi(1), " i"
WRITE(*,"(A9,F9.4,1X,A1,F9.4,A2)") &  
   "2nd Root ", xr(2), "+", xi(2), " i"  
END PROGRAM TestQuadraticEq

SUBROUTINE QUADRATIC_EQ(p,q,re,im)
!  ==================================================================================
!  CODE1.3-Quadratic_Eq.f95. A fortran module for implementing Pseudocode 1.3.          
! 
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
!  First Edition. (c) By Zekeriya ALTAC (2024).
!  ISBN: 978-1-032-75474-1 (hbk)
!  ISBN: 978-1-032-75642-4 (pbk)
!  ISBN: 978-1-003-47494-4 (ebk)
!  
!  DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
!  
!  This free software is complimented by the author to accompany the textbook.
!  E-mail: altacz@gmail.com.
!                                                                                    
!  DESCRIPTION: A Subroutine to find the roots of a quadratic equation of the           
!    form : x**2 + p * x + q = 0.                                                 
!                                                                                    
!  ON ENTRY                                                                          
!    p, q :: Coefficients of the quadratic equation;                                 
!                                                                                    
!  ON EXIT                                                                           
!    re   :: Array of length 2 containing real parts of the roots: re1, re2;         
!    im   :: Array of length 2 containing imaginary parts of the roots: im1, im2.    
!                                                                                    
!  USES                                                                              
!    SQRT :: Built-in Intrinsic function to evaluate the square root of a real value.
!                                                                                    
!  REVISION DATE :: 03/18/2024                                                       
!                                                                                    
!  ==================================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)  :: p, q         ! Input Variables
DOUBLE PRECISION, INTENT(OUT) :: re(2), im(2) ! Output Variables
DOUBLE PRECISION :: d                         ! Internal Variables  
!
d = p*p-4.0*q 
if(d < 0.0d0) then  ! real imaginary roots
   d = sqrt(-d)
   re(1)=-p/2.0d0 ; re(2)= re(1) 
   im(1)=-d/2.0d0 ; im(2)=-im(1)
else                ! d>=0, Real Roots 
   d = sqrt(d)
   re(1)=(-p-d)/2.0d0 ; im(1)= 0.0d0
   re(2)=(-p+d)/2.0d0 ; im(2)= 0.0d0
end if
!
END SUBROUTINE QUADRATIC_EQ

 

