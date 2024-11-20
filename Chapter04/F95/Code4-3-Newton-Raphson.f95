PROGRAM test_Newton_Raphson
! ==============================================================================
!  The main program to test SUBROUTINE NEWTON_RAPHSON
! ==============================================================================
IMPLICIT NONE
INTEGER :: maxit, iter
DOUBLE PRECISION  :: eps, root
!
maxit= 99
eps  = 0.50D-4
root = 2.20d0
!
 CALL NEWTON_RAPHSON(root,maxit,eps,iter)
 print*, "-----------------------------"
 print*, "Root is",root," converged after ",iter," iterations"
 print*, "-----------------------------"        
 STOP
 END PROGRAM

 SUBROUTINE NEWTON_RAPHSON(root,maxit,eps,iter)
! =======================================================================        
! CODE4-3.NEWTON-RAPHSON.F95 A fortran program for implementing Pseudocode 4.3.
!
! DESCRIPTION: A subroutine to compute a root of a nonlinear equation using the
!   Newton-Raphson method.
!
! ON ENTRY
!   root  :: Initial guess for the root;
!   maxit :: Maximum number of iterations permitted;
!   eps   :: Convergence tolerance. 
!
! ON EXIT
!   iter  :: Number of iterations realized;
!   root  :: Computed approximation for the root.
!
! USES
!   ABS   :: Built-in Intrinsic function returning the absolute value of a real value;
!
! ALSO REQUIRED
!   FUNC  :: User-defined external function providing the nonlinear equation, f(x).
!   FUNCP :: User-defined external function providing the first derivative
!            of the nonlinear equation, f'(x).
!
! REVISION DATE :: 04/29/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT( IN) :: maxit
INTEGER, INTENT(OUT) :: iter
DOUBLE PRECISION, INTENT( IN) :: eps
DOUBLE PRECISION, INTENT(INOUT) :: root
DOUBLE PRECISION :: fn, fpn, aerr, rate, x0, xn, del, del0
DOUBLE PRECISION :: FUNC, FUNCP
INTEGER :: p 
!
PRINT 110,
!
del0 = 1.0d0 
  x0 = root   ! Initialize root 
  p  = 0      ! Initialize iteration counter
!!!!! ------- BEGIN REPEAT-UNTIL LOOP  !!!!!!!!!!!!!!!!!!!
DO
   fn  = FUNC (x0)
   fpn = FUNCP(x0)
   del =-fn/fpn   
   aerr= DABS(del)      ! absolute error
   rate= aerr/del0**2   ! estimate the convergence rate
   PRINT 111, p, x0, fn, fpn, aerr, rate
   xn  = x0 + del        ! update the estimate
   x0  = xn
   del0= DABS(del)
   p = p+1 
   IF( (DABS(fn)<eps.AND.aerr<eps).OR.(p==maxit)) EXIT ! ---- End the iteration loop
ENDDO
root = xn
iter = p
IF(p==maxit) THEN
  PRINT*, "** Max iteration number reached="
  PRINT*, "** Estimated root has NOT converged, del, f(x)= ", DABS(del),FUNC(xn)
ENDIF         
110 FORMAT(3x,"p",6x,"x^(p)",8x,"f(x^(p))",8x,"f'(x^(p))",7x,"Aerr", &
     10x,"Rate")        
111 FORMAT(1x,i3,2x,7(2x,G13.7))        
END SUBROUTINE  

DOUBLE PRECISION FUNCTION func(x)
! =======================================================================          
! User-defined function, defining f(x) which should be cast as func(x)=0.
! =======================================================================  
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x
   func= 4.0d0+x*x*(8.0d0-x*x)
END FUNCTION func

DOUBLE PRECISION FUNCTION funcp(x)
! =======================================================================          
! User-defined function, defining f'(x)=funcp(x).
! =======================================================================  
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x
    funcp= 4.0d0*x*(4.0-x*x) 
END FUNCTION funcp

   
