PROGRAM Test_Back_Substitute
! ==============================================================================
!  The main program to test SUBROUTINE FORW_SUBSTITUTE
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=5
DOUBLE PRECISION :: a(n,n), b(n), x(n)
INTEGER :: i, j
!
A = RESHAPE( (/ 1.d0, 2.d0, 4.d0,-3.d0, 1.d0, &
                0.d0, 3.d0, 4.d0,-4.d0, 1.d0, &      
                0.d0, 0.d0, 4.d0, 3.d0, 1.d0, &     
                0.d0, 0.d0, 0.d0, 2.d0, 1.d0, &     
                0.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/ 5, 5 /), ORDER=(/2,1/) ) 

b = (/10.0d0, 7.0d0, 29.0d0, 13.0d0, 5.0d0 /) 
!
PRINT '(/" ********** Input Matrix ***********")'
DO i=1,n
   WRITE(*,1) (A(i,j),j=1,n)
END DO          
Print*, " "          
!
CALL Back_Substitute(n,a,b,x)
!
PRINT '(/" ********** Input Matrix ***********")'
DO i=1,n
   WRITE(*,1) (A(i,j),j=1,n)
END DO           

PRINT '(/" *** Solution Vector ***")' 
DO i=1,n
   WRITE(*,*) i,x(i)
END DO              
PRINT*, "-----------------------------------";
1 FORMAT(5(f10.5,3x),5x,f10.2)
END PROGRAM Test_Back_Substitute
! *********END OF MAIN PROGRAM *********************************


              
SUBROUTINE Back_Substitute(n,U,b,x)
! ==============================================================================
! CODE2-8-BACK_SUBSTITUTE.F95. A fortran program for implementing Pseudocode 2.8.
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
! DESCRIPTION: A subroutine to find the solution of a upper-triangular system 
!   using back substitution algorithm.
!
! ON ENTRY
!     n  :: Number of unknowns; 
!     U  :: Input coefficient (upper-triangular) matrix (n×n);
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
