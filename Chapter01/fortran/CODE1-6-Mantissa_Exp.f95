PROGRAM Test_Mantissa_Exp 
! ==============================================================================
!  The main program to test SUBROUTINE MANTISSA_EXP
! ==============================================================================
IMPLICIT NONE
INTEGER :: e
DOUBLE PRECISION :: fl, m
!
PRINT*, "Enter a real number "
READ *, fl
!
CALL Mantissa_Exp(fl,m,e)
PRINT*, " Mantissa = ",m
PRINT*, " Exponent = ",e
!
END PROGRAM Test_Mantissa_Exp 


SUBROUTINE MANTISSA_EXP(fl, m, e)
! ==============================================================================
!  CODE1.6-MantissaExp.f95. A fortran module implementing Pseudocode 1.6.           
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
!  DESCRIPTION: A module to determine a floating-point number's                      
!     mantissa and exponent, i.e., fl=mx10^e                                         
!                                                                                    
!  ON ENTRY                                                                          
!     fl  :: A floating point number.                                                
!                                                                                    
!  ON EXIT                                                                           
!     M   :: Mantissa;                                                               
!     e   :: Exponent.                                                               
!                                                                                    
!  USES                                                                              
!    DABS  :: Built-in Intrinsic function returning the absolute value of a real value
!    FLOOR :: Built-in Intrinsic function returning the greatest integer less than    
!      or equal to a real value;                                                     
!    DLOG10:: Built-in Intrinsic function returning the base 10 logarithm of a real value.
!                                                                                    
!  REVISION DATE :: 03/22/2024   
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: fl
DOUBLE PRECISION, INTENT(OUT):: m
INTEGER, INTENT(OUT) :: e
!
IF( DABS(fl)>0.0d0) THEN
  e = FLOOR(DLOG10(fl)) + 1  
  m = fl * 10.0d0**(-e)
ELSE
  e = 0.0d0
  m = 0.0d0
ENDIF
!
END SUBROUTINE MANTISSA_EXP
 

 
