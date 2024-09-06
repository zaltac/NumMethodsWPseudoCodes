PROGRAM Test_MACHEPS
! ==============================================================================
!  The main program to test FUNCTION MACHEPS
! ==============================================================================
IMPLICIT NONE 
DOUBLE PRECISION :: x, MACHEPS 
WRITE(*,"(A28)") "Enter a real value in 0<x<=1" 
READ(*,*) x 
! Print result
WRITE(*,*) "Machine Epsilon is=", MACHEPS(x) 
END PROGRAM Test_MACHEPS
 

DOUBLE PRECISION FUNCTION MACHEPS(x)
!  ==================================================================================
!  CODE1.5-MachEps.f95. A fortran module for implementing Pseudocode 1.5.               
! 
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
!  First Edition. (c) By Zekeriya ALTAÇ (2024).
!  ISBN: 978-1-032-75474-1 (hbk)
!  ISBN: 978-1-032-75642-4 (pbk)
!  ISBN: 978-1-003-47494-4 (ebk)
!  
!  DOI : 10.1201/9781003474944
!  C&H/CRC PRESS, Boca Raton & London.
!  
!  This free software is complimented by the author to accompany the textbook.
!  E-mail: altacz@gmail.com.
!                                                                                    
!  DESCRIPTION: A function module to calculate the machine epsilon, which            
!     provides a positive machine value that is almost negligible compared to 1.     
!                                                                                    
!  INPUT ARGUMENT                                                                    
!     x   :: A real number, 0<x<=1.                                                  
!                                                                                    
!  REVISION DATE :: 03/18/2024                                                       
!                                                                                    
!  ==================================================================================
DOUBLE PRECISION, INTENT(IN) :: x
DO WHILE (1.0d0 + x/2.0d0 > 1.0d0)
   x = x/2.0d0
END DO
MACHEPS = x 
END FUNCTION MACHEPS
 
