PROGRAM Test_XdotY
! ==============================================================================
!  The main program to test FUNCTION XdotY
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=4
INTEGER            :: i
DOUBLE PRECISION   :: X(n), Y(n), XdotY
!
X=(/1.0d0, 2.0d0, 3.0d0, 4.0d0/)
Y=(/4.0d0, 3.0d0, 2.0d0, 1.0d0/)
!
WRITE(*,"(/5x,'Input Vector X'/)")
DO i=1,n
  PRINT*, "x(",i,")=",x(i)  
ENDDO
WRITE(*,"(/5x,'Input Vector Y'/)")
DO i=1,n
  PRINT*, "y(",i,")=",y(i)  
ENDDO    
!
PRINT*, " "
PRINT*, "------ X*Y dot product is ==> ", XdotY(n,X,y)
!  
END PROGRAM Test_XdotY
! *********END OF MAIN PROGRAM *********************************              

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
