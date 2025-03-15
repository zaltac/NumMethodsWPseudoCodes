PROGRAM Test_InverseL
! ==============================================================================
!  The main program to test SUBROUTINE InverseL
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=5
INTEGER :: i, j
DOUBLE PRECISION :: L(n,n), eL(n,n) 
!        
L = RESHAPE( (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                4.d0, 3.d0, 0.d0, 0.d0, 0.d0, &      
                5.d0, 2.d0, 4.d0, 0.d0, 0.d0, &     
                3.d0, 1.d0, 8.d0, 2.d0, 0.d0, &     
               -1.d0, 1.d0, 1.d0, 1.d0, 1.d0/), (/ 5, 5 /), ORDER=(/2,1/) ) 
!
PRINT '(/" ----- Input Lower-T Matrix L -----")'
DO i = 1, n
   PRINT 1, (L(i,j),j=1,n)
ENDDO
!
CALL InverseL(n,L,eL)
!
PRINT '(//" ------ Output : Inverse L ------")'
DO i = 1, n
   PRINT 1, (eL(i,j),j=1,n)
ENDDO  
1 FORMAT(1x,6(f12.7,1x))
END PROGRAM Test_InverseL



SUBROUTINE InverseL(n,L,eL)
!  ==================================================================================
!  CODE11.4-InverseL.f95. A Fortran95 module implementing Pseudocode 11.4.                        
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

!   DESCRIPTION: A subroutine to invert a lower-triangular matrix.                                 
!                                                                                              
!   ON ENTRY                                                                                   
!     n    :: Dimension attribute of the matrix L;                                             
!     L    :: A lower-triangular matrix, nxn.                                                  
!                                                                                              
!   ON EXIT                                                                                    
!     eL   :: Inverse of L, also a lower-triangular matrix, nxn.                               
!                                                                                              
!   REVISION DATE :: 03/15/2025                                                                
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN) :: L(n,n) 
DOUBLE PRECISION, INTENT(OUT) :: eL(n,n) 
DOUBLE PRECISION :: sums
INTEGER :: i, j, k
!
eL = 0.d0  ! Initialize matrix eL with zero.
!
eL(1,1) = 1.0d0/L(1,1)
DO i = 2, n
   eL(i,i) =1.0d0 / L(i,i)
   DO j = i-1, 1, -1
      sums = 0.d0
      DO k = j+1, i
         sums = sums + eL(i,k) * L(k,j)
      ENDDO 
   eL(i,j) = -sums / L(j,j)
   ENDDO 
ENDDO 
END SUBROUTINE InverseL
     
