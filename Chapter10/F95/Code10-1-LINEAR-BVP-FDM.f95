PROGRAM Test_LBVP_Solve
! ==============================================================================
!  The main program to test SUBROUTINE LBVP_SOLVE.F95
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: nmax=201
DOUBLE PRECISION :: x(nmax), y(nmax), alpha(2), beta(2), gamma(2)  
DOUBLE PRECISION :: error, xa, xb, h, EXACT
INTEGER :: neq, n, i
!
PRINT*, "Enter No of intervals "
READ*, n
IF(n>nmax) THEN
   Print*, "Max grid size is ",nmax
   Print*, "Increase NMAX, and then try it again."
   STOP
ENDIF
!
! Setup the ODE, Grids & BCs       
neq = n +1 
xa = 1.0d0   ; xb = 2.0d0   
h  = (xb-xa)/dfloat(n)
DO i=1,neq
   x (i) = xa + DFLOAT(i-1)*h
END DO 
alpha(1)= 1.0d0;  beta(1)= 2.0d0;  gamma(1)= 4.0d0 
alpha(2)= 1.0d0;  beta(2)=-2.0d0;  gamma(2)= 2.0d0
!     
CALL LBVP_SOLVE(neq,x,alpha,beta,gamma,y)
!      
PRINT 1,
DO i=1,neq
   error=dabs(y(i)-exact(x(i)))
   PRINT 2,   x(i),exact(x(i)),y(i),error
ENDDO
1 FORMAT(//,8x,"x",7x,"Exact",6x,"N.Approx",6x,"Abs Error")
2 FORMAT(1x,f10.4,2F12.7,1PE16.5)        
END PROGRAM Test_LBVP_Solve
 

SUBROUTINE LBVP_SOLVE(neq,x,alpha,beta,gamma,y)
!  ==================================================================================
!  CODE10.1-LBVP_SOLVE.f95. A Fortran95 module implementing Pseudocode 10.1.                      
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

!  DESCRIPTION: A subroutine to find approximate solution of the following linear                  
!    differential equation using the Finite Difference Method:                                 
!          p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]                                   
!    subject to                                                                                
!          alpha1 * y'(a)+ beta1 * y(a) = gamma1                                               
!          alpha2 * y'(b)+ beta2 * y(b) = gamma2                                               
!                                                                                              
!  CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1.  
!                                                                                              
!  ON ENTRY                                                                                    
!    neq  :: Number of (equations) grid poinds;                                                
!     x   :: Array of length neq containing the abscissa of the grid points;                   
!    alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                   
!           the boundary conditions as stated above;                                           
!                                                                                              
!  ON EXIT                                                                                     
!     y   :: Array of length neq containing the approximate solution.                          
!                                                                                              
!  USES                                                                                        
!   COEFFS  :: A subroutine containing the coefficients and rhs of the linear ODE, i.e., p(x), q(x)
!   TRIDIAGONAL:: A subroutine solving a tridiagonal system of equations using the Thomas algorithm
!                                                                                              
!  REVISION DATE :: 03/07/2025                                                                 
!  ==================================================================================
IMPLICIT NONE           
INTEGER, INTENT(IN) :: neq
DOUBLE PRECISION, INTENT( IN) :: x(neq), alpha(2), beta(2), gamma(2)
DOUBLE PRECISION, INTENT(OUT) :: y(neq)
DOUBLE PRECISION :: b(neq), a(neq), d(neq), cy(neq)
DOUBLE PRECISION :: h, px, qx, rx, gx
INTEGER :: BC_Left, BC_Rigt, k, s1
!        
BC_Left = INT(alpha(1)) 
BC_Rigt = INT(alpha(2))
h  = x(2)-x(1)   
! Initialize the arrays of the tridiagonal system
a = 0.0d0; b=0.0d0; d=0.0d0; CY=0.0d0
!
DO k = 1, neq
   call COEFFS(x(k),px,qx,rx,gx)
   d (k)= rx   -2.0d0*px/h**2
   a (k)= (px/h + 0.5d0*qx)/h
   b (k)= (px/h - 0.5d0*qx)/h
   cy(k)= gx  
END DO
!
IF(BC_Left==0) THEN
   d(1)= 1.0d0
   a(1)= 0.0d0; b(1)=0.0d0
  cy(1)= gamma(1)/beta(1)
ELSE     
   a(1)= a(1) +b(1)
   d(1)= d(1) + 2.0d0*h*b(1)*beta(1)/alpha(1)  
  cy(1)= cy(1)+ 2.0d0*h*b(1)*gamma(1)/alpha(1)        
ENDIF        
 
IF(BC_Rigt==0) THEN
   d(neq)= 1.0d0
   a(neq)= 0.0d0; b(neq)=0.0d0
  cy(neq)= gamma(2)/beta(2)
ELSE        
   b(neq)= a(neq) +b(neq)
   d(neq)= d(neq) -2.0d0*h*a(neq)*beta(2)/alpha(2) 
  cy(neq)= cy(neq)-2.0d0*h*a(neq)*gamma(2)/alpha(2)        
ENDIF
!
s1=1     	 
CALL TRIDIAGONAL(s1, neq, b, d, a, cy, y)
!
END SUBROUTINE LBVP_SOLVE        
 
 
SUBROUTINE COEFFS(x,p,q,r,g)
! ==============================================================================
! DESCRIPTION: A user-defined module suppling the coefficients p(x), q(x), r(x) and
!    rhs g(x) of the linear ODE given in the following form: 
!         p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]
!
! ON ENTRY
!    x   :: Independent variable (a<= x<= b);
!
! ON EXIT
!   p, q, r, abd g :: Coefficients & rhs evaluated at x. 
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT( IN) :: x
DOUBLE PRECISION, INTENT(OUT) :: p, q, r, g
p = x*x 
q = -5.0d0*x
r = 8.0d0 
g = 0.0d0
END SUBROUTINE COEFFS
 

DOUBLE PRECISION FUNCTION Exact(x)
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
!   
Exact= x*x*(x*x-0.5d0)  
!   
END FUNCTION Exact        
!======================================================================
 

SUBROUTINE TRIDIAGONAL(s1, sn, b, d, a, c, x)
! ==============================================================================
! CODE2-13-TRIDIAGONAL.F95. A fortran program for implementing Pseudocode 2.13
!     
! NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
! First Edition. (c) By Zekeriya ALTAÇ (2024).
! ISBN: 9781032754741 (hbk)
! ISBN: 9781032756424 (pbk)
! ISBN: 9781003474944 (ebk)
!
! DOI : 10.1201/9781003474944
! C&H/CRC PRESS, Boca Raton & London. 
!  
! This free software is complimented by the author to accompany the textbook.
! E-mail: altacz@gmail.com
!
! DESCRIPTION: A subroutine to solve a tridiagonal system of linear equations
!   using Thomas algorithm. 
!
! ON ENTRY
!    s1 :: Subscript of the first unknown (usually 1); 
!    sn :: Subscript of the last unknown (usually No. of Eqs); 
!     b :: Array of length n containing coefficients of below diagonal elements;
!     d :: Array of length n containing coefficients of diagonal elements;
!     a :: Array of length n containing coefficients of above diagonal elements;
!     c :: Array of length n containing coefficients of rhs.
!
! ON EXIT
!     x :: An array of length n containing the solution.
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE 
INTEGER, INTENT(IN) :: s1, sn
DOUBLE PRECISION, INTENT( IN), DIMENSION(sn) :: a, b
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(sn) :: d, c
DOUBLE PRECISION, INTENT(OUT), DIMENSION(sn) :: x
INTEGER :: i 
DOUBLE PRECISION :: ratio
!
DO i = s1+1, sn
  ratio = b(i)/d(i-1)
   d(i) = d(i) - Ratio * a(i-1)
   c(i) = c(i) - Ratio * c(i-1)
ENDDO

x(sn) = c(sn) / d(sn)
DO i = sn-1, s1, (-1)
   x(i) = (c(i) - a(i) * x(i+1) )/d(i)
ENDDO
END SUBROUTINE TRIDIAGONAL
