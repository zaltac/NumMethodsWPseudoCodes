PROGRAM Test_MAT_ADD
! ==============================================================================
!  The main program to test SUBROUTINE MAT_ADD
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=3, m=3
INTEGER            :: i, j
DOUBLE PRECISION, DIMENSION(m,n) :: A, B, C

A = reshape((/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0 /), shape(A))
B = reshape((/ 3.d0,-2.d0, 1.d0,-2.d0, 2.d0, 4.d0, 3.d0,-5.d0, 1.d0 /), shape(B))
!
WRITE(*,"(/5x,'Input Matrix A'/)")
DO i=1,m
  PRINT 1, (A(i,j),j=1,n)  
ENDDO
WRITE(*,"(/5x,'Input Matrix B'/)")
DO i=1,m
  PRINT 1, (B(i,j),j=1,n)  
ENDDO
!
CALL MAT_ADD(m,n,A,B,C)
WRITE(*,"(/5x,'Output Matrix C'/)")
DO i=1,m
  PRINT 1, (C(i,j),j=1,n)  
ENDDO
 
1 FORMAT(2x,10(F10.5))  
END PROGRAM Test_MAT_ADD
! *********END OF MAIN PROGRAM *********************************
              

SUBROUTINE MAT_ADD(m,n,A,B,C)
!  ==================================================================================
!  CODE2.1-MAT_ADD.f95. A Fortran95 module implementing Pseudocode 2.1.                           
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
!   DESCRIPTION: A Subroutine to perform C=A+B matrix addition.                                
!                                                                                              
!   ON ENTRY                                                                                   
!     m,n :: Dimension attributes of the matrices;                                             
!      A  :: An input matrix (mxn);                                                            
!      B  :: An input matrix (mxn).                                                            
!                                                                                              
!   ON EXIT                                                                                    
!      C :: The output matrix (mxn).                                                           
!                                                                                              
!   REVISION DATE :: 03/18/2024                                                                
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: m, n
DOUBLE PRECISION, INTENT( IN), DIMENSION(m,n) :: A, B
DOUBLE PRECISION, INTENT(OUT), DIMENSION(m,n) :: C
INTEGER :: i, j
!
DO i=1,m
   DO j=1,n
      C(i,j) = A(i,j) + B(i,j)
   ENDDO
ENDDO
END SUBROUTINE MAT_ADD
