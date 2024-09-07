PROGRAM TestFactorial
! ==============================================================================
!  The main program to test RECURSIVE FUNCTION FACTORIAL
! ==============================================================================
INTEGER :: n, factorial

WRITE(*,"(A7)") "Enter n" 
READ(*,*) n

WRITE(*,*) n,"! = ",FACTORIAL(n)

END PROGRAM

INTEGER RECURSIVE FUNCTION FACTORIAL(n) RESULT(fact)
!  ==================================================================================
!  CODE1.4-factorial.f95. A fortran module implementing Pseudocode 1.4.          
! 
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
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
!  DESCRIPTION: A recursive function for computing n!                                
!                                                                                    
!  INPUT ARGUMENT                                                                    
!    n    :: Integer (n>=0)                                                          
!                                                                                    
!  ON EXIT                                                                           
!   Result:: n!                                                                      
!                                                                                    
!  REVISION DATE :: 03/21/2024                                                       
!                                                                                    
!  ==================================================================================
INTEGER, INTENT(IN) :: n    ! Input Variable
! 
IF(n<0) THEN
  write(*,*) "Error, Illegal Input"
  write(*,*) "Argument should be n>=0, you entered ",n
  stop
ELSE
   IF(n<=1) THEN
       fact = 1
   ELSE
       fact = n*FACTORIAL(n-1)
   ENDIF  
ENDIF  
END FUNCTION FACTORIAL 
