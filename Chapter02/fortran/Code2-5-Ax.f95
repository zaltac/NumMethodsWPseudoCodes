PROGRAM Test_AX
! ==============================================================================
!  The main program to test SUBROUTINE AX
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=4 
INTEGER            :: i, j
DOUBLE PRECISION   :: A(n,n), x(n), b(n)
 
!
A = RESHAPE( (/ 1.0d0, 2.0d0, 3.0d0, -2.0d0, & 
                4.0d0, 1.0d0, 2.0d0,  3.0d0, &
                3.0d0, 2.0d0, 1.0d0,  2.0d0, &
               -2.0d0, 3.0d0, 4.0d0,  1.0d0 /), (/ 4, 4 /), ORDER=(/2,1/) ) 
x = (/ 2.0d0, 5.0d0, 3.0d0, -3.0d0/)

WRITE(*,"(/5x,'INPUT MATRIX A')")
DO i=1,n
  PRINT 1, (A(i,j),j=1,n)
ENDDO
WRITE(*,"(/5x,'INPUT VECTOR x')/")
DO i=1,n
  PRINT 2, i,x(i)
ENDDO       
!
CALL Ax(n,A,x,b)  ! performs A**(-1)*A=C matrix operation
!
PRINT*, " "
PRINT*, "------ A(n,n)X(n) product is "

!
WRITE(*,"(/5x,'OUTPUT VECTOR b')/")
DO i=1,n
  WRITE(*,1) b(i)
ENDDO           
1 FORMAT(3x,10(f10.5,2x))
2 FORMAT(5x,"x(",i2,") =",f11.7)
END PROGRAM Test_AX
! *********END OF MAIN PROGRAM *********************************              

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
