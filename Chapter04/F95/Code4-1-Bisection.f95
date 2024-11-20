PROGRAM Test_Bisection
! ==============================================================================
!  The main program to test SUBROUTINE BISECTION.F95
! ==============================================================================
IMPLICIT NONE
INTEGER :: Halves, maxit
DOUBLE PRECISION :: a, b, eps, root, fa, fb, FUNC
!
maxit= 99
eps  = 5.0d-5
   a = 0.0d0
   b = 4.0d0
!   
fa = FUNC(a)
fb = FUNC(b)
IF(fa*fb>0) THEN
   PRINT*, 'No root in interval (a,b). Change the interval.'
   STOP
ENDIF
!
CALL Bisection(a,b,maxit,eps,root,halves)
!
PRINT 1, root, halves       
root  = (a*fb-b*fa)/(fb-fa)
PRINT 2, root, DABS(FUNC(root))
1 FORMAT(1x,44("="),/,1x,"! Root is ",1PG14.8," after",i3, &
   " bisections !",/,1x,44("="))
2 FORMAT(//1x," *** Root (with linear interpolation)   = ",1PG14.9,/, &
         1x," *** Clossness to the root, Abs[f(root)]= ",1PG10.4//)
END PROGRAM Test_Bisection

SUBROUTINE BISECTION(a,b,maxit,eps,root,halves)
! ==============================================================================
! CODE4-1-BISECTION.F95. A fortran program for implementing Pseudocode 4.1
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
! DESCRIPTION: A fortran subroutine to find a root of a nonlinear equation in [a,b] 
!   using the Bisection method.
!
! ON ENTRY
!  [a,b] :: Initial search interval (it must bracket one root);
!  maxit :: Maximum number of iterations permitted;
!  eps   :: Convergence tolerance. 
!
! ON EXIT
!  halves:: Number of halves realized;
!  root  :: Computed approximation for the root.
!
! USES
!  ABS   :: Built-in Intrinsic function returning the absolute value of a real value;
!
! ALSO REQUIRED
!  FUNC  :: User-defined external function providing the nonlinear equation.
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT( IN) :: maxit
INTEGER, INTENT(OUT) :: halves
DOUBLE PRECISION, INTENT( IN) :: eps
DOUBLE PRECISION, INTENT(INOUT) :: a, b
DOUBLE PRECISION, INTENT(OUT) :: root
INTEGER :: p
DOUBLE PRECISION :: fa, fb, xm, fm, FUNC, interval
!   
p = 0
interval = b - a
fa = FUNC(a)
fb = FUNC(b)
!
PRINT 2,
DO 
    p = p + 1    ! **************** Start REPEAT-UNTIL loop *******************
   xm = 0.5d0*(a+b)
   fm = FUNC(xm)
   PRINT 1, p, a, b, fa, fb, xm, fm , interval
   IF(fa*fm>0.0d0) THEN
      a = xm
      fa= fm
   ELSE ! case of fa*fm<0.0d0
      b = xm
      fb= fm
   ENDIF 
   interval=0.50d0*interval
   IF( ((DABS(fm)<eps).AND.(interval<eps)).OR. p==maxit) EXIT
END DO  !    End the iteration loop  
 root  = xm
halves = p           
IF(p==maxit) THEN
  PRINT 3, p
ENDIF    
! 
1 FORMAT(1x,i3,2x,2(1x,1PG12.7),2(1x,1PE12.4),3x,1PG12.7,1x,1PE11.4,1x,1PE11.4)
2 FORMAT(3x,"p",7x,"a",12x,"b",12x," f(a)",10x,"f(b)",8x,"xm",&
  10x,"f(xm)",7x,"interval",/,97("-"))    
3 FORMAT(1x,37("!"),/,1x,"Max iteration number reached=",i3,1x,37("!"))            
END SUBROUTINE BiSection


DOUBLE PRECISION FUNCTION func(x)
! ==========================================================================          
! User-defined function providing f(x), which should be cast as Func(x)=0.
! ==========================================================================   
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x
!
func= x * x + 0.025d0 * x - 4.0d0
!
END FUNCTION func
