PROGRAM Test_Linear_Solve 
! ==============================================================================
!  The main program to test SUBROUTINE Linear_Solve
! ==============================================================================      
IMPLICIT NONE
INTEGER, PARAMETER :: n=7
DOUBLE PRECISION   :: A(n,n), b(n), x(n)
INTEGER :: i, j, opt
!
A = RESHAPE( (/ 2.0d0,	1.0d0, -1.0d0,	0.0d0,	3.0d0,	1.0d0,	0.0d0,  &
               -1.0d0,	3.0d0,	1.0d0,	2.0d0,	4.0d0, -2.0d0,	1.0d0,  &
                4.0d0,	1.0d0,	5.0d0,	1.0d0, -3.0d0,	2.0d0,	2.0d0,  &
               -2.0d0,	1.0d0, -2.0d0, -2.0d0,	3.0d0,	3.0d0,	1.0d0,  &
                1.0d0,	1.0d0,	1.0d0,	3.0d0,	2.0d0, -3.0d0,	2.0d0,  &
                4.0d0,	1.0d0, -1.0d0,	0.0d0,	3.0d0, -5.0d0,	1.0d0,  &
                2.0d0,	7.0d0, -3.0d0, -4.0d0, -1.0d0,	2.0d0,	2.0d0 /), (/ 7, 7 /), ORDER=(/2,1/) ) 
B=(/ -10.0d0, 7.0d0, 28.0d0, -30.0d0, 16.0d0, 3.0d0, -12.0d0/)     

DO i=1,n
   WRITE(*,1) (A(i,j),j=1,n), b(i)
ENDDO

opt=1
! *** Solve the system of equations          
 CALL Linear_Solve(n,A,b,x,opt)
!   
IF(opt==0) PRINT '(/" Naive Gauss Elimination ")'
IF(opt/=0) PRINT '(/" Gauss Jordan Elimination")'     
PRINT*, " Solution Vector "
DO i=1,n
   WRITE(*,2) i,x(i)
ENDDO 

1 FORMAT(3X,8(F8.3,1x))
2 FORMAT(3X,"x(",I1,")=",F6.3)
END PROGRAM Test_Linear_Solve 
! *********END OF MAIN PROGRAM *********************************

 
SUBROUTINE Linear_Solve(n,A,b,x,opt)
!  ==================================================================================
!  CODE2.9-LINEAR_SOLVE.f95. A Fortran95 module implementing Pseudocode 2.9.                      
! 
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
!  First Edition. (c) By Zekeriya ALTA« (2024).
!  ISBN: 978-1-032-75474-1 (hbk)
!  ISBN: 978-1-032-75642-4 (pbk)
!  ISBN: 978-1-003-47494-4 (ebk)
!  
!  DOI : 10.1201/9781003474944
!  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
!  
!  This free software is complimented by the author to accompany the textbook.
!  E-mail: altacz@gmail.com.

!  DESCRIPTION: A subroutine to solve a system of linear equations using naive                 
!    Gauss Elimination (opt=0) or Gauss-Jordan Elimination (opt/=0) algorithm.                 
!                                                                                              
!  ON ENTRY                                                                                    
!      n  :: Number of unknowns;                                                               
!      A  :: Input coefficient matrix of size n√ón;                                            
!      b  :: An input array of length n containing the rhs;                                    
!     opt :: Option key (=0, Naive Gauss-Elim.; /=0, Gauss-Jordan Elimn).                      
!                                                                                              
!  ON EXIT                                                                                     
!      x  :: The output array of length n containing the solution.                             
!                                                                                              
!  USES                                                                                        
!    ABS  :: Built-in Intrinsic function returning the absolute value of a real value;         
!    Back_Substitute :: A subrotine to solve an upper-triangular system.                       
!                                                                                              
!  REVISION DATE :: 03/18/2024                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, opt
DOUBLE PRECISION, INTENT(INOUT) :: A(n,n), b(n)
DOUBLE PRECISION, INTENT(OUT) :: x(n)
INTEGER :: i, j, k, m
DOUBLE PRECISION :: eps, val, ajj, s 
!
eps=1.0D-12  ! Internally set positive small number, it may be replaced by MACHEPS
!
!  *** FORWARD ELIMINATION STEPS
!        
DO j=1,n
   ajj = A(j,j)
   IF( abs(ajj) < eps) then
     PRINT*, "Pivot is zero at j=",j
     PRINT*, "Execution is halted!"
     STOP
   ELSE
     val = 1.0d0 / ajj
     b(j) = b(j) * val
     DO m = j, n
        A(j,m) = A(j,m) * val
     ENDDO
     DO i = j+1, n
        s = A(i,j)
        A(i,j) = 0.0d0
        DO k = j+1, n
           A(i,k) = A(i,k) - s * A(j,k)
        ENDDO
        b(i) = b(i) - s * b(j)
     ENDDO  
   ENDIF          
ENDDO
!
! *** BACK SUBSTITUTION STEPS
!
IF(opt==0) THEN
    CALL Back_Substitute(n,A,b,x)
ELSE
    DO j=n,2,-1
       DO i=j-1,1,-1
          b(i) = b(i) - A(i,j) * b(j)
          A(i,j)= 0.0d0
       ENDDO
       x(j) = b(j)
    ENDDO
    x(1) = b(1)
ENDIF
END SUBROUTINE Linear_Solve
 
SUBROUTINE Back_Substitute(n,U,b,x)
! ==============================================================================
! CODE2-8-BACK_SUBSTITUTE.F95. A fortran program for implementing Pseudocode 2.8.
!     
! NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
! First Edition. (c) By Zekeriya ALTA« (2024).
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
! DESCRIPTION: A subroutine to find the solution of a upper-triangular system 
!   using back substitution algorithm.
!
! ON ENTRY
!     n  :: Number of unknowns; 
!     U  :: Input coefficient (upper-triangular) matrix (n◊n);
!     b  :: Input array of size n containing the rhs.
!
! ON EXIT
!     x  :: Output array of size n containing the solution.
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN), DIMENSION(n,n) :: U
DOUBLE PRECISION, INTENT( IN), DIMENSION(n)   :: b
DOUBLE PRECISION, INTENT(OUT), DIMENSION(n)   :: x
DOUBLE PRECISION :: sums
INTEGER :: k, j
!       
x(n) = b(n) / U(n,n)
DO k=(n-1),1,(-1)
   sums = 0.0d0
   DO j=(k+1),n
      sums = sums + U(k,j) * x(j)
   ENDDO
   x(k) = (b(k) - sums)/U(k,k)
ENDDO               
END SUBROUTINE Back_Substitute 
