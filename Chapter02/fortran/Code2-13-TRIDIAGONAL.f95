PROGRAM Test_Tridiagonal
! ==============================================================================
!  The main program to test SUBROUTINE TRIDIAGONAL
! ==============================================================================              
IMPLICIT NONE
INTEGER, PARAMETER :: n=9
DOUBLE PRECISION, DIMENSION(n) :: b, d, a, c, x
INTEGER :: i, s1, sn
!
b =(/ 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /)
d =(/-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0 /)
a =(/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /)
c =(/-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -14.0/)
!  
PRINT '(/" ********** Input Tridiagonal Matrix & RHS ***********")'
DO i=1,n
   WRITE(*,1) b(i), d(i), a(i), c(i)
END DO        
s1=1
sn=n
CALL TRIDIAGONAL(s1,sn,b,d,a,c,x)
! 
! ***     PRINT OUT THE RESULTS
!         
PRINT*, " ";PRINT*, '*** Solution ***';PRINT*, " ";
DO i=1,n
   WRITE(*,2) i,x(i)
ENDDO
WRITE(*,*) " "
1 FORMAT(3(f10.5,3x),5x,f10.5)
2 FORMAT(3X,"x(",I1,")=",F6.3)
END PROGRAM Test_Tridiagonal


SUBROUTINE TRIDIAGONAL(s1, sn, b, d, a, c, x)
!  ==================================================================================
!  CODE2.13-TRIDIAGONAL.f95. A Fortran95 module implementing Pseudocode 2.13.                     
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
!  DESCRIPTION: A subroutine to solve a tridiagonal system of linear equations                 
!    using Thomas algorithm.                                                                   
!                                                                                              
!  ON ENTRY                                                                                    
!     s1 :: Subscript of the first unknown (usually 1);                                        
!     sn :: Subscript of the last unknown (usually no. of eqs, n);                             
!      b :: Array of length n containing coefficients of below diagonal elements;              
!      d :: Array of length n containing coefficients of diagonal elements;                    
!      a :: Array of length n containing coefficients of above diagonal elements;              
!      c :: Array of length n containing coefficients of rhs.                                  
!                                                                                              
!  ON EXIT                                                                                     
!      x :: An array of length n containing the solution.                                      
!                                                                                              
!  REVISION DATE :: 03/18/2024                                                                 
!  ==================================================================================
IMPLICIT NONE 
INTEGER, INTENT(IN) :: s1, sn
DOUBLE PRECISION, INTENT( IN), DIMENSION(sn) :: a, b
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(sn) :: d, c
DOUBLE PRECISION, INTENT(OUT), DIMENSION(sn) :: x
INTEGER :: i 
DOUBLE PRECISION :: ratio
!
DO i = s1+1, sn
  ratio = b(i)/d(i-1)
   d(i) = d(i) - Ratio * a(i-1)
   c(i) = c(i) - Ratio * c(i-1)
ENDDO

x(sn) = c(sn) / d(sn)
DO i = sn-1, s1, (-1)
   x(i) = (c(i) - a(i) * x(i+1) )/d(i)
ENDDO
END SUBROUTINE TRIDIAGONAL
