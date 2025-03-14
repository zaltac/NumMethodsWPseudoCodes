PROGRAM Test_POWER_METHOD_S
! ==============================================================================
!  The main program to test SUBROUTINE POWER_METHOD_S.F95
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=4
INTEGER :: i, maxit
DOUBLE PRECISION   :: A(n,n), x(n), eps, Aerr, lambda
!
A = RESHAPE( (/  8.d0,  9.d0,  10.d0,   9.d0, &   
               -13.d0,-12.d0, -12.d0, -11.d0, & 
               -18.d0, -9.d0, -20.d0,  -9.d0, &        
                11.d0,  1.d0,  10.d0,   0.d0/), (/ 4, 4 /), ORDER=(/2,1/) ) 
!
! Initialize eigenvalue & eigenvector (initial guess)
lambda= 1.0d0
x=0.0d0 
x(1) = 1.0d0
maxit = 299
eps   = 1.0d-3
!  
! *** Apply Power Method  
CALL POWER_METHOD_S(n,A,X,lambda,Aerr,eps,maxit)    
!   
! *** Print out the results
IF(Aerr<eps) then
   WRITE(*,"(//20('='))")
   WRITE(*,"('Largest lambdavalue=',F12.7)") lambda
   WRITE(*,"('lambdavector       =',9F12.4)") (x(i),i=1,n)
   WRITE(*,"(20('='))")
ELSE
   WRITE(*,*)'Max number of iterations is reached.'       
ENDIF 
5 FORMAT(9(f10.3,3x))       
6 FORMAT(1x,"iter=",i3,"lambdaval=",F12.7," Abs error",1PG12.3)
END PROGRAM Test_POWER_METHOD_S



SUBROUTINE POWER_METHOD_S(n, A, x, lambda, error, eps, maxit)
!  ==================================================================================
!  CODE11.1-POWER_METHOD_S.f95. A Fortran95 module implementing Pseudocode 11.1.                  
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
!   DESCRIPTION: A module to find the dominant eigenvalue (largest in absolute value)          
!      using the Power Method with scaling technique.                                          
!                                                                                              
!   ON ENTRY                                                                                   
!      n     :: Size of the matrix;                                                            
!      A     :: A real square matrix (nxn);                                                    
!      x     :: Array of length n containing the initial guess for the eigenvector;            
!     lambda :: An initial guess for the dominant eigenvalue;                                  
!      eps   :: Convergence tolerance;                                                         
!     maxit  :: Maximum number of iterations permitted.                                        
!                                                                                              
!   ON EXIT                                                                                    
!     lambda :: Estimated dominant (largest in absolute value) eigenvalue;                     
!      x     :: Array of length n containing the estimated eigenvector;                        
!     error  :: Error, max of both L2-norm of displacement vector and relative error           
!                for eigenvalue.                                                               
!                                                                                              
!   USES                                                                                       
!     dabs    :: Built-in intrinsic function returning the absolute value of a double float;                                  
!     AX      :: A module for computing Ax matrix-vector product (see Pseudocode 2.5);          
!     MAX_SIZE:: A function module providing largest (in absolute) value of a vector.          
!                                                                                              
!   REVISION DATE :: 03/14/2025                                                                
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, maxit
DOUBLE PRECISION, INTENT( IN) :: A(n,n), eps
DOUBLE PRECISION, INTENT(INOUT) :: X(n), lambda
DOUBLE PRECISION, INTENT(OUT) :: error
INTEGER :: i, p
DOUBLE PRECISION ::  Xn(n), MAX_SIZE, err1, err2, lambdao, ENORM
!
error = MAX_SIZE(n,x)
lambdao=lambda
   p = 0
!
DO WHILE (Error>eps .AND. p<maxit) 
   p = p + 1       ! count iterations
   CALL AX(n,A,x,xn)         ! Solve [A][x]=[y]
!   print*, xn
   lambda=MAX_SIZE(n,xn)     ! Find L-infty of x^(p+1)   
   DO i=1,n                  ! Normalize x
      xn(i)=xn(i)/lambda
   END DO
   Err1 = ENORM(n,DABS(xn-x))          
   Err2 = DABS(1.0d0 - lambdao / lambda)
   print*, "Err1, Err2",Err1,Err2
   Error= DMAX1(Err1, Err2)           
   print 1, p, lambda, error
   lambdao= lambda
   x = xn
   pause
END DO  
1 FORMAT(1x,i3,1x,F12.7,1x,1PE14.6)  
END SUBROUTINE POWER_METHOD_S
       
DOUBLE PRECISION FUNCTION MAX_SIZE(n,x)
!  ==================================================================================
!  CODE11.1-MAX_SIZE.f95 of a module in Pseudocode 11.1.                  
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

!   DESCRIPTION: A function module to find the largest element (in absolute value) of an array.
!                                                                                              
!   ARGUMENTS                                                                                  
!      n   :: Length of the array;                                                          
!      x   :: Array of length n.                                                            
!                                                                                              
!   USES                                                                                       
!     DABS :: Built-in intrinsic function returning the absolute value of a real value in double precision.       
!                                                                                              
!   REVISION DATE :: 03/10/2025                                                                
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: x(n)
DOUBLE PRECISION :: xmax
INTEGER :: i
!
xmax=x(1)
DO i=2,n
   IF (DABS(x(i)) > DABS(xmax)) xmax=x(i)
END DO
MAX_SIZE=xmax
END FUNCTION MAX_SIZE

SUBROUTINE AX(n,A,x,b)
! ==============================================================================
! CODE2-5-AX.F95. A fortran program for implementing Pseudocode 2.5.
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
! DESCRIPTION: A subroutine to perform A * x = b matrix-vector multiplication.
!
! ON ENTRY
!    n   :: Dimension attributes of input/output matrices;
!    A   :: An input matrix of size n×n;
!    x   :: An input vector of length n.
!
! ON EXIT
!    b   :: The output vector of length n.  
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN), DIMENSION(n,n):: A
DOUBLE PRECISION, INTENT( IN), DIMENSION(n)  :: x
DOUBLE PRECISION, INTENT(OUT), DIMENSION(n)  :: b
INTEGER :: i, k
!
DO i=1,n
   b(i) = 0.0d0
   DO k=1,n
      b(i) = b(i) + A(i,k) * x(k)
   END DO 
END DO 
!
END SUBROUTINE AX

DOUBLE PRECISION FUNCTION ENORM(n,x)
!  ==================================================================================
!  CODE3.1-JACOBI.F95. A Fortran95 module implementing ENORM module of Pseudocode 3.1.                            
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

!  DESCRIPTION: A function module to compute Euclidean (L2) norm of a vector.                  
!                                                                                              
!  ARGUMENTS                                                                                   
!      n  :: The length of an input vector;                                                    
!      x  :: A vector (array) of length n.                                                     
!                                                                                              
!  USES                                                                                        
!    DSQRT:: Built-in Intrinsic function returning the square root of a double float.  
!  REVISION DATE :: 11/09/2024                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n 
DOUBLE PRECISION, INTENT(IN) :: x(n) 
INTEGER :: i
DOUBLE PRECISION :: delta      
!
delta = 0.0d0
DO i = 1, n
   delta = delta + x(i) * x(i)
ENDDO 
ENorm = DSQRT(delta)      
END FUNCTION ENORM




