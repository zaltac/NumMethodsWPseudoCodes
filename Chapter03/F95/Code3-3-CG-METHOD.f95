PROGRAM Test_CGM
! ==============================================================================
!  The main program to test SUBROUTINE CGM
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=10
DOUBLE PRECISION :: a(n,n), b(n), x(n), error, eps
INTEGER :: i, j, iter, maxit
!
a =RESHAPE((/ 2., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., &
        -1.,  2., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,       &
         0., -1.,  2., -1.,  0.,  0.,  0.,  0.,  0.,  0.,       & 
         0.,  0., -1.,  2., -1.,  0.,  0.,  0.,  0.,  0., &
         0.,  0.,  0., -1.,  2., -1.,  0.,  0.,  0.,  0., &
         0.,  0.,  0.,  0., -1.,  2., -1.,  0.,  0.,  0., &
         0.,  0.,  0.,  0.,  0., -1.,  2., -1.,  0.,  0., &
         0.,  0.,  0.,  0.,  0.,  0., -1.,  2., -1.,  0., &
         0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  2., -1., &
         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  2. /), (/10,10/) )

b = (/ 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.2/)

!-------- write input data ----------------------------  
 WRITE(9,"(4x,'Coefficient Matrix',4x,'RHS',4x,'Initial Guess')")

Do i=1,n
  x(i) = 0.0d0  ! ***  SET INITIAL GUESS
  print 1, (A(i,j),j=1,n),b(i),x(i)
enddo     
!
eps  = 1.0D-07
maxit= 999 
!
! ***  GO TO CGM MODULE 
CALL CGM(n,eps,a,b,x,maxit,iter,error)
!
! ***  PRINT OUT THE RESULTS
!
PRINT "(//4x,' *** Solution ***'/)",
DO i=1,n
   WRITE(*,2) i, x(i)
ENDDO          
!
WRITE(*,3) iter, error 
1 FORMAT(3x,13(F5.1),3x,f4.1)
2 FORMAT(2x,i3,2x,f11.6)
3 FORMAT(/2x," Total no of iterations = ",i4,/,  &
          2x," Maximum Error          = ",1PE10.3//)         
END PROGRAM Test_CGM
              
SUBROUTINE CGM(n,eps,A,b,x,maxit,iter,error)
! ==============================================================================
! CODE3-3-CG-METHOD.F95. A fortran program for implementing Pseudocode 3.3.
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
! DESCRIPTION: A subroutine to solve Ax=b linear system with the Conjugate Gradient Method.
!
! ON ENTRY
!    n   :: Number of equations; 
!    A   :: Coefficient matrix (n×n);  
!    b   :: Array of length n containing the right hand side;
!    x   :: Array of length n containing the initial guess;    
!   eps  :: Convergence tolerance;
!  maxit :: Maximum number of iterations.
!
! ON EXIT
!    x   :: Array of length n containing the approximate solution;    
!  iter  :: Total number of iterations performed;
!  error :: Euclidean (L2-) norm of displacement at exit.
!
! USES
!   AX   :: Subroutine to evaluate A*x matrix vector product;
!   SQRT :: Built-in Intrinsic function returning the square root of a real value;
!   XdotY:: Function to evaluate the dot product of two vectors.
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, maxit
DOUBLE PRECISION, INTENT( IN) :: A(n,n), b(n), eps
INTEGER, INTENT(OUT) :: iter
DOUBLE PRECISION, INTENT(OUT) :: x(n), error
DOUBLE PRECISION, DIMENSION(n) :: c, r, d 
DOUBLE PRECISION :: Enorm, rho0, XdotY, rho, alpha, beta
INTEGER ::  p 
!          
CALL AX(n,A,x,r)
! following statements involve "whole array" arithmetics
r = b - r          ! [r]^(0)=[b]-[A][x]^(0)
d = r              ! [d]^(0)=[r]^(0)
rho0 = XdotY(n,r,r)       ! rho^(0)=[r^(0)].[r^(0)]
PRINT "(//4x,'*** Iteration history ***')"
! -------------------------------------------------------------------
DO p=0,MAXIT                   ! BEGIN CGM ITERATION LOOP  
   Enorm= DSQRT(rho0)          ! E-norm=Sqrt([r]^(p-1).[r]^(p-1))           
   PRINT 1, p,Enorm            ! PRINTOUT ITERATION PROGRESS  
   IF (Enorm<eps) EXIT         ! CHECK FOR CONVERGENCE, EXIT if converged... 
   CALL AX(n,A,d,c)            ! [c]^(p-1)=[A][d]^(p-1)
   rho=XdotY(n,d,c)            ! rho=[d]^(p-1).[c]^(p-1)
   alpha=rho0/rho              ! alpha^(p)=[r].[r]/([d].[c])
   x = x + alpha * d           ! [x]^(p)=[x]^(p-1)+alfa^(p)*[d]^(p-1)    
   r = r - alpha * c           ! [r]^(p)=[r]^(p-1)-alfa^(p)*[d]^(p-1) 
   rho=XdotY(n,r,r)            ! rho^(p+1)=[r^(p)].[r^(p)]
   beta=rho/rho0               ! beta^(p) = rho/rho0
   d =r + beta * d             ! [d]^(p)=[r]^(p)+beta*[d]^(p-1)
   rho0=rho                    ! rho^(p)<==rho^(p+1)
END DO 
! -------------------------------------------------------------------       
iter=p                         ! Set total number of iterations
error=Enorm                    ! Set recent Enorm as error       
IF(iter==maxit) PRINT 2, maxit! PRINTout WARNING if iterations not converged 
1 FORMAT(2x,"iter=",i4,4x,"E-norm= ",1PE11.5)
2 FORMAT(//2x,"Faied to converge after",i4," iterations"//) 
END SUBROUTINE CGM

 
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



DOUBLE PRECISION FUNCTION XdotY(n,x,y)
! ==============================================================================
! CODE2-3-XdotY.F95. A fortran program for implementing Pseudocode 2.3
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
! DESCRIPTION: A function to compute the dot product of two vectors, x and y.
!
! ARGUMENTS
!    n   :: Dimension attribute of the input vectors; 
!   x, y :: The input vectors of length n.
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN), DIMENSION(n) :: x, y
DOUBLE PRECISION ::sums
INTEGER :: i
!
sums = 0.0d0
DO i=1,n
   sums = sums + x(i)*y(i)
END DO 
XdotY = sums
END FUNCTION XdotY  
