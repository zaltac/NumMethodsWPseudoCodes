PROGRAM Test_Jacobi_Method
! ==============================================================================
!  The main program to test SUBROUTINE JACOBI.F95
! ==============================================================================
IMPLICIT NONE
INTEGER, parameter :: n=10
DOUBLE PRECISION   :: a(n,n), & ! Coefficient matrix 
                        b(n), & ! rhs vector
                       xo(n), & ! initial guess at start
                        x(n)    ! solution estimate
INTEGER            :: i, j, iter, maxit, Iprint
DOUBLE PRECISION   :: Errmax, eps 
! Create and input data file

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

!-------- read input data ----------------------------  
WRITE(9,"(4x,'Coefficient Matrix',4x,'RHS',4x,'Initial Guess')")

DO i=1,n
  xo(i) = 1.0d0
  print 1, (A(i,j),j=1,n),b(i),xo(i)
ENDDO

eps  =1.0d-6
maxit=999 
Iprint=1  ! Iprint=1, print iteration history to the output file; Iprint=0 does not!  
iter=iprint
!  
! ***  CALL SUBROUTINE JACOBI TO ITERATE
!
CALL JACOBI(n,eps,A,b,xo,x,maxit,iter,Errmax)
!
! ***  PRINT OUT THE RESULTS
!
PRINT*, " ";Print*, "----- Solution --------"
DO i=1,n
   WRITE(*,2) i, x(i)
ENDDO          
PRINT*, "-----------------------"
WRITE(*,3) iter,errmax
1 FORMAT(2x,10(F6.1),3x,2f12.2)
2 FORMAT(2x,i3,2x,f11.6)
3 FORMAT(/2x,"Total no of iterations = ",i4,/, &
          2x,"Maximum Error          = ",1PE10.3//)
END PROGRAM Test_Jacobi_Method
                   
SUBROUTINE JACOBI(n,eps,A,b,xo,x,maxit,iter,error)
! ==============================================================================
! Code3-1-JACOBI-METHOD.F95. A fortran program for implementing Pseudocode 3.1.
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
! DESCRIPTION: A subroutine to iteratively solve Ax=b using the Jacobi method.
!
! ON ENTRY
!    n   :: Number of equations (size of A); 
!    A   :: Input coefficient matrix (n×n);
!    b   :: Array of length n containing the right-hand;  
!    x   :: Array of length n containing the estimate at (p+1)'th step;
!    xo  :: Array of length n containing the initial guess, or iterates at estimate at p'th step;
!   eps  :: Convergence tolerance;
!  maxit :: Maximum permitted number of iterations.
!
! ON EXIT
!    x   :: Array of length n containing the estimated solution;
!  iter  :: Total number of iterations performed;
!  error :: L2 norm of the displacement vector.
!
! USES
!   JACOBI_DRV :: Accompanying subroutine performing one step Jacobi iteration.
!
! REVISION DATE :: 11/09/2024 
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, maxit 
DOUBLE PRECISION, INTENT( IN) :: a(n,n), b(n), eps 
DOUBLE PRECISION, INTENT(INOUT) :: xo(n) 
DOUBLE PRECISION, INTENT(OUT) :: x(n) 
DOUBLE PRECISION, INTENT(OUT) :: error
INTEGER, INTENT(OUT) :: iter
DOUBLE PRECISION :: del, del0
INTEGER :: p
!
! ***   START OF ITERATION LOOP
! 
p= 0
del0=1.0d0 
DO
  p=p+1
  CALL JACOBI_DRV(n,A,b,xo,x,del)  ! APPLY ONE-STEP JACOBI ITERATION 
! ***   PRINTOUT THE ITERATION PROGRESS        
  PRINT 1,   p,del,del/del0
  xo  = x
  del0= del
  IF( del<eps .OR. p==maxit ) EXIT  ! CHECK FOR CONVERGENCE OR ITERATION BOUND 
ENDDO 
!
! ***   END OF ITERATION LOOP              
error= del
iter = p
IF(p==maxit) THEN
   PRINT 2, maxit  ! ***   IF ITERATIONS REACH THE UPPER BOUND
ENDIF
!        
1 format(2x,"iter=",i4,3x," Error= ",1PE11.5,3x,"  Ratio= ",1PE11.5)
2 format(//2x,"Jacobi method faied to converge after",i4," iterations",/, &
           2x,"within the specified EPS tolerance.",//)          
4 format(1x,8F12.6)
END SUBROUTINE JACOBI !=====================================
        
SUBROUTINE JACOBI_DRV(n,A,b,xo,x,del)
! ==============================================================================
! Code3-1-JACOBI-METHOD.F95. A fortran program for implementing JACOBI_DRV in Pseudocode 3.1
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
! DESCRIPTION: A subroutine to to perform one step Jacobi iteration and compute
!    the Euclidean norm of the displacement vector.
!
! ON ENTRY
!    n   :: Number of equations (size of A); 
!    A   :: Input coefficient matrix (n×n);
!    b   :: Array of length n containing the right-hand;  
!    x   :: Array of length n containing the estimate at (p+1)'th step;
!    xo  :: Array of length n containing the estimate at p'th step.
!
! ON EXIT
!    x   :: Array of length n containing the estimated solution;
!    del :: Maximum absolute error achieved.
!
! USES
!   ENORM:: User-defined function calculating the Euclidean vector (L2 norm) of a vector.
!
! REVISION DATE :: 11/09/2024 
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n 
DOUBLE PRECISION, INTENT( IN) :: A(n,n), b(n) 
DOUBLE PRECISION, INTENT( IN) :: xo(n) 
DOUBLE PRECISION, INTENT(OUT) :: x(n) 
DOUBLE PRECISION, INTENT(OUT) :: del
DOUBLE PRECISION :: sums, ENorm
DOUBLE PRECISION :: d(n) 
INTEGER :: i, j
!
DO i=1,n        
   sums = 0.0d0 
   DO j=1,n
      IF (i.NE.j) THEN 
         sums = sums + A(i,j) * xo(j)
      ENDIF
   ENDDO
   x(i) = (b(i) - sums) / A(i,i)   
   d(i) = x(i) - xo(i)   ! find the displacement vector  
ENDDO                       
del = Enorm(n,d) ! find L2- (E-)norm of the displacement vector
!    
END SUBROUTINE JACOBI_DRV ! ===================================

DOUBLE PRECISION FUNCTION ENorm(n,x)
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
!    SQRT :: Built-in Intrinsic function returning the square root of a real value.            
!                                                                                              
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
END FUNCTION ENorm

