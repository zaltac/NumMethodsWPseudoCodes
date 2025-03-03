PROGRAM Test_Lagrange
! ==============================================================================
!  The main program to test FUNCTION LAGRANGE_EVAL
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n = 3
DOUBLE PRECISION, DIMENSION(0:n) :: x, f
DOUBLE PRECISION :: xval, LAGRANGE_EVAL
!
x(0)= 0.08d0; x(1)= 0.250d0; x(2)= 0.50d0; x(3)= 0.90d0
f(0)= 0.25d0; f(1)= 0.625d0; f(2)= 0.81d0; f(3)= 0.43d0 
!
PRINT*, "Enter x "
READ *, xval    
!
PRINT*, "fval= ",LAGRANGE_EVAL(n,xval,x,f)
!
END PROGRAM Test_Lagrange



SUBROUTINE LAGRANGE_P(n,xval,x,L)
!  ==================================================================================
!  CODE6.1-LAGRANGE-P.F95. A Fortran95 module implementing Pseudocode 6.1.                      
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

!  DESCRIPTION: A subroutine to compute the Lagrange polynomials at x=xval.                        
!                                                                                              
!  ON ENTRY                                                                                    
!    n    :: The number of data in the set minus 1;                                            
!    xval :: x value at which the dataset is to be interpolated;                               
!    x    :: Array of length n+1 containing abscissas, k=0,1,2,...,n.                          
!                                                                                              
!  ON EXIT                                                                                     
!    L    :: Array of length (n+1) containing Lagrange polynomials                             
!            that is, L(k) = L_k(xval) for k=0,1,2,...,n.                                      
!                                                                                              
!  REVISION DATE :: 02/28/2025                                                                 
!  ==================================================================================
IMPLICIT NONE 
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN)  :: x(0:n), xval
DOUBLE PRECISION, INTENT(OUT) :: L(0:n)
INTEGER :: i, k
!
DO k=0,n
   L(k)=1.0d0
   DO i=0,n
      IF(i/=k) THEN
         L(k)=L(k)*(xval-x(i))/(x(k)-x(i))
      END IF
   ENDDO
ENDDO
END SUBROUTINE LAGRANGE_P



DOUBLE PRECISION FUNCTION LAGRANGE_EVAL(n,xval,x,f)
! ==============================================================================
!  CODE6.1-LAGRANGE_EVAL.F95. A Fortran95 module implementing Pseudocode 6.1.                        
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
!  DESCRIPTION: A function to evaluate Lagrange interpolating polynomial for an arbitrary xval, 
!    f=f(xval)=fval within x0 <= xval <= xn.                                                   
!                                                                                              
!  ARGUMENTS                                                                                   
!    n    :: The number of data in the set minus 1;                                            
!    xval :: x value at which the dataset is to be interpolated;                               
!    x    :: Array of length (n+1) containing abscissas, k=0,1,2,...,n;                        
!    f    :: Array of length (n+1) containing ordinates, k=0,1,2,...,n.                        
!                                                                                              
!  USES                                                                                        
!    LAGRANGE_P :: A module generating the Lagrange polynomials.                                         
!                                                                                              
!  REVISION DATE :: 02/28/2025        
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN) :: x(0:n), f(0:n) 
DOUBLE PRECISION :: L(0:n), xval, fval
INTEGER :: k
!
CALL LAGRANGE_P(n,xval,x,L)
!
fval=0.0d0
DO k=0,n
    fval = fval + L(k)*f(k)
ENDDO        
LAGRANGE_EVAL = fval
END FUNCTION LAGRANGE_EVAL

 
