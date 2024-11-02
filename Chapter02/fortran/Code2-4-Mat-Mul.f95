PROGRAM Test_MAT_MUL
! ==============================================================================
!  The main program to test SUBROUTINE MAT_MUL
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: m=2, p=4, n=3
INTEGER            :: i, j
DOUBLE PRECISION   :: A(m,p), B(p,n), C(m,n)
!
A = RESHAPE( (/1.0,-3.0, 5.0, 1.0, -2.0, 4.0, 1.0, 2.0/), (/ 2, 4 /), ORDER=(/2,1/) ) 
B = RESHAPE( (/3.0, 1.0,-2.0, 2.0,  0.0, 4.0,-1.0, 1.0, &
     -3.0, 2.0, 5.0, 3.0 /), (/4, 3/),ORDER=(/2,1/))

WRITE(*,"(/5x,'INPUT MATRIX A'/)")
DO i=1,m
  WRITE(*,1) (A(i,j),j=1,p)  ! Print matrix A to the monitor
ENDDO

WRITE(*,"(/5x,'INPUT MATRIX B')/")
DO i=1,p
  WRITE(*,1) (B(i,j),j=1,n)  ! Print matrix B to the monitor
ENDDO       
PRINT*, " ********* End of Input data *********"
!
PRINT*, " "
PRINT*, "------ A(m,p)xB(p,n) matrix product is ------"
CALL MAT_MUL(m,p,n,A,B,C)  ! performs A**(-1)*A=C matrix operation
!
WRITE(*,"(/5x,'OUTPUT MATRIX C')/")
DO i=1,m
  WRITE(*,1) (C(i,j),j=1,n)
ENDDO           
1 FORMAT(2x,10(f10.5,2x))
END PROGRAM Test_MAT_MUL
! *********END OF MAIN PROGRAM *********************************
              

SUBROUTINE MAT_MUL(m,p,n,A,B,C)
!  ==================================================================================
!  CODE2.4-MAT_MUL.f95. A Fortran95 module implementing Pseudocode 2.4.                           
! 
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
!  First Edition. (c) By Zekeriya ALTA� (2024).
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
!   DESCRIPTION: A subroutine to find A*B=C matrix multiplication.                             
!                                                                                              
!   ON ENTRY                                                                                   
!    m,p,n :: Dimension attributes of input/output matrices;                                   
!       A  :: An input matrix of size m×p;                                                    
!       B  :: An input matrix of size p×n.                                                    
!                                                                                              
!   ON EXIT                                                                                    
!       C  :: The output matrix of size m×n.                                                  
!                                                                                              
!   REVISION DATE :: 03/18/2024                                                                
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: m, p, n
DOUBLE PRECISION, INTENT( IN) :: A(m,p), B(p,n)
DOUBLE PRECISION, INTENT(OUT) :: C(m,n)
INTEGER :: i, j, k
!
DO i=1,m
   DO j=1,n
      c(i,j) = 0.0d0
      DO k=1,p
         c(i,j) = c(i,j) + A(i,k) * B(k,j)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE MAT_MUL
