PROGRAM Test_EXPA
! ==============================================================================
!  The main program to test FUNCTION EXPA
! ==============================================================================
DOUBLE PRECISION :: x, eps
WRITE(*,"(A7)") "Enter x" 
READ(*,*) x
eps=1.0d-6
WRITE(*,1) x,EXPA(x,eps)
1 format(1x,"e^",f5.3," =",f10.7)
END PROGRAM Test_EXPA


DOUBLE PRECISION FUNCTION EXPA(x,eps)
!  ==================================================================================
!  CODE1.8-ExpA.f95. A Fortran95 module implementing Pseudocode 1.8.                              
! 
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
!  First Edition. (c) By Zekeriya ALTAÇ (2024).
!  ISBN: 978-1-032-75474-1 (hbk)
!  ISBN: 978-1-032-75642-4 (pbk)
!  ISBN: 978-1-003-47494-4 (ebk)
!  
!  DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
!  
!  This free software is complimented by the author to accompany the textbook.
!  E-mail: altacz@gmail.com.
!                                                                                              
!  DESCRIPTION: A function to compute e^x adaptively using the MacLaurin series                
!     within a user-defined tolerance.                                                         
!                                                                                              
!  ARGUMENTS                                                                                   
!     x   :: A real input (exponent) value;                                                    
!    eps  :: A user-defined convergence tolerance.                                             
!                                                                                              
!  USES                                                                                        
!   DFLOAT:: A built-in intrinsic function that converts an integer argument to a double precision value.  
!   DABS  :: A built-in intrinsic function returning the absolute value of a double precision value.         
!                                                                                              
!  REVISION DATE :: 04/11/2024                                                                 
!                                                                                              
!  ==================================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x      ! Input Variables
DOUBLE PRECISION, INTENT(IN) :: eps    
DOUBLE PRECISION  :: sums, term        ! Internal Variables
INTEGER :: k
!
k=0 
term= 1.0d0
sums= 1.0d0
DO WHILE ( DABS(term) > eps * DABS(sums) )
   k = k + 1
   term = term*x/DFLOAT(k)
   sums = sums + term
END DO
EXPA = sums
!
END FUNCTION EXPA
