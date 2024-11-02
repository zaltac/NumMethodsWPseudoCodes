PROGRAM Test_Inv_Mat
! ==============================================================================
!  The main program to test SUBROUTINE Inv_Mat
!  USES
!     CODE2-4-MAT_MUL.F95
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=5
INTEGER            :: i, j
DOUBLE PRECISION, DIMENSION(n,n) :: A, AI, C, B
!
A = RESHAPE( (/1.0d0,  2.0d0,  1.0d0,  2.0d0,  3.0d0, &
              11.0d0, -1.0d0,  1.0d0,  4.0d0,  1.0d0, &
               4.0d0, -1.0d0,  1.0d0,  1.0d0, -1.0d0, &
              -3.0d0,  1.0d0, -8.0d0, -1.0d0,  5.0d0, &
              -1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0/), (/ 5, 5 /), ORDER=(/2,1/) ) 

PRINT*, "-------- Input matrix A(-1) -------";       
DO i=1,n
  WRITE(*,1) (A(i,j),j=1,n)
ENDDO
DO i=1,n; DO j=1,n
   B(i,j)=A(i,j)  ! A will be destroyed so save A on B for checking
ENDDO ; ENDDO              
PRINT*, " "
!
CALL INV_MAT(n,A,AI)
!
PRINT*, "------ Inverse matrix A(-1) -------"
DO i=1,n
   WRITE(*,1) (AI(i,j),j=1,n)
ENDDO        
PRINT*, "-----------------------------------"
CALL MAT_MUL(n,n,n,B,AI,C)  ! performs A**(-1)*A=C matrix operation

!  Check the accuracy of the result
PRINT*, "------ matrix A*A(-1) -------------"
DO i=1,n
   WRITE(*,1) (C(i,j),j=1,n)
ENDDO              
PRINT*, "-----------------------------------"
1 FORMAT(5(f10.5,3x))
END PROGRAM Test_Inv_Mat
! *********END OF MAIN PROGRAM *********************************

              
SUBROUTINE Inv_Mat(n,A,B)
!  ==================================================================================
!  CODE2.6-INV_MAT.f95. A Fortran95 module implementing Pseudocode 2.6.                           
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
!  DESCRIPTION: A subroutine to find the inverse of a square matrix (with no pivoting).        
!                                                                                              
!  ON ENTRY                                                                                    
!     n  :: Dimension attribute of input matrix A;                                             
!     A  :: An input matrix (nxn).                                                             
!                                                                                              
!  ON EXIT                                                                                     
!     B  :: Inverse of A (nxn).                                                                
!                                                                                              
!  REVISION DATE :: 03/18/2024                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,n) :: A
DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,n) :: B
DOUBLE PRECISION ::  p, s
INTEGER :: i, j, k
!    
B = 0.0d0 
DO i=1,n
   B(i,i) = 1.0d0   ! Set B = I, i.e., construct identity matrix
ENDDO   
!     
DO j=1,n
   p = 1.0d0 / A(j,j)
   DO k=1,n      
      A(j,k) = p * A(j,k)
      B(j,k) = p * B(j,k)
   ENDDO
   DO i=1,n
      s = A(i,j)
      IF(i/=j) THEN
         DO k=1,n
           A(i,k) = A(i,k) - s * A(j,k)
           B(i,k) = B(i,k) - s * B(j,k)
         ENDDO 
      ENDIF 
   ENDDO
ENDDO       
END SUBROUTINE Inv_Mat


SUBROUTINE MAT_MUL(m,p,n,A,B,C)
! ==============================================================================
! CODE2-4-MAT_MUL.F95. A fortran program for implementing Pseudocode 2.4.
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
! DESCRIPTION: A subroutine to find A*B=C matrix multiplication. 
!
! ON ENTRY
!  m,p,n :: Dimension attributes of input/output matrices; 
!     A  :: An input matrix of size m×p;
!     B  :: An input matrix of size p×n.
!
! ON EXIT
!     C  :: The output matrix of size m×n. 
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: m, p, n
DOUBLE PRECISION, INTENT( IN) :: A(m,p), B(p,m)
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
