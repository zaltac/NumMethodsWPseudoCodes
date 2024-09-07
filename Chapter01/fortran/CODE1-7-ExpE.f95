PROGRAM Test_EXPE
! ==============================================================================
!  The main program to test FUNCTION EXPE
! ==============================================================================
DOUBLE PRECISION :: x, EXPE
INTEGER :: n
!
WRITE(*,"(A7)") "Enter x, n " 
READ(*,*) x, n
WRITE(*,1) x,EXPE(x,n)
1 format(1x,"e^",f5.3," =",f10.7)
END PROGRAM Test_EXPE


DOUBLE PRECISION FUNCTION EXPE(x,n)
!  ==================================================================================
!  CODE1.7-ExpE.f95. A Fortran95 module implementing Pseudocode 1.7.                              
! 
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
!  First Edition. (c) By Zekeriya ALTAÇ (2024).
!  ISBN: 978-1-032-75474-1 (hbk)
!  ISBN: 978-1-032-75642-4 (pbk)
!  ISBN: 978-1-003-47494-4 (ebk)                                                                                                                                                             
!  DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.                           
!                                                                                              
!  This free software is complimented by the author to accompany the textbook.                 
!  E-mail: altacz@gmail.com                                                                    
!                                                                                              
!  DESCRIPTION: A function to compute e^x using the MacLaurin series with specified            
!     number of terms.                                                                         
!                                                                                              
!  ARGUMENTS                                                                                   
!     x   :: A real input (exponent) value;                                                    
!     n   :: The number of terms of the MacLauring series to be included.                      
!                                                                                              
!  USES                                                                                        
!    FLOAT:: A built-in intrinsic function that converts an integer argument to a real value.  
!                                                                                              
!  REVISION DATE :: 03/18/2024                                                                                                                                                               
!  ==================================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x          ! Input Variables
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION    :: sums, term ! Internal Variables
INTEGER :: k
!
term= 1.0d0
sums= 1.0d0 
DO k = 1, n-1
     term = term*x/dFLOAT(k)
     sums = sums + term
END DO
EXPE = sums
END FUNCTION EXPE
