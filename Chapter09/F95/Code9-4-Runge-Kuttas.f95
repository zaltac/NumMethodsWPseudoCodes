PROGRAM Test_Runge_Kuttas
! ==============================================================================
!  The main program to test SUBROUTINE Runge_Kutta
! ==============================================================================
IMPLICIT NONE
INTEGER :: n
DOUBLE PRECISION :: x0, xlast, y0, h
!     
  x0 = 0.0d0; 
xlast= 5.0d0 
  y0 = 2.0d0 
  h  = 0.1d0; 
!
PRINT*, "Enter Order of RK method "
READ*, n
!    
CALL Runge_Kutta(n,h,x0,y0,xlast)
!
END PROGRAM Test_Runge_Kuttas

SUBROUTINE Runge_Kutta(n,h,x0,y0,xlast)
!  ==================================================================================
!  CODE9.4-Runge_Kutta.f95. A Fortran95 module implementing Pseudocode 9.4.                       
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

!  DESCRIPTION: A subroutine to estimate the solution of a first order IVP on [x0,xlast]       
!    using the 2nd to 4th order Runge-Kutta scheme. Numerical estimates are printed out, not stored.
!                                                                                              
!  ON ENTRY                                                                                    
!    n     :: Order of the Runge-Kutta scheme;                                                 
!    h     :: Step size (it must be uniform);                                                  
!    x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps
!   xlast  :: End point of the solution interval.                                              
!                                                                                              
!  Other Internal Variables                                                                    
!    x,y   :: Current estimates, x^(p+1) and y^(p+1).                                          
!                                                                                              
!   USES                                                                                       
!     ABS  :: Built-in Intrinsic function returning the absolute value of a real value.        
!    DRV_RK:: A driver subprogram performing one-step RK scheme.                               
!                                                                                              
!  REVISION DATE :: 03/05/2025                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN) :: h, xlast
DOUBLE PRECISION, INTENT(INOUT) :: x0, y0
DOUBLE PRECISION :: x, y, yt, aerr, Exact
!
WRITE(*,1)
WRITE(*,2) x0, y0 
x=x0
DO WHILE (x<xlast)  
   CALL DRV_RK(n,h,x0,y0,x,y)
   yt = exact(x)
   aerr= DABS(y-yt)
   WRITE(*,2) x, y, aerr
   x0=x; y0=y      
ENDDO
1 FORMAT(8x,"x",13x,"y",11x,"Abs. Error",/,4x,42("-"))
2 FORMAT(2x,2(f12.7,2x),1PE14.3)
END SUBROUTINE Runge_Kutta

SUBROUTINE DRV_RK(n,h,x0,y0,x,y)
!  ==================================================================================
!  CODE9.4-DRV_RK.f95. A Fortran95 module implementing Pseudocode 9.4.                            
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

!  DESCRIPTION: A driver sub-program employing one-step RK2, RK3, or RK4 scheme.               
!                                                                                              
!  ON ENTRY                                                                                    
!   n     :: Order of Runge-Kutta scheme;                                                      
!   h     :: Step size (it must be uniform);                                                   
!   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
!                                                                                              
!  ON EXIT                                                                                     
!   x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
!                                                                                              
!  USES                                                                                               
!    FCN  :: User-defined external function providing y'=f(x,y).                               
!                                                                                              
!  REVISION DATE :: 03/05/2025                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT( IN) :: h 
DOUBLE PRECISION, INTENT( IN) :: x0, y0
DOUBLE PRECISION, INTENT(OUT) :: x, y
DOUBLE PRECISION :: x1, xh, ym, xk, xk1, xk2, xk3, xk4, hlf, FCN
!
hlf= 0.50d0
xh = x0 + 0.5d0*h
x1 = x0 + h

SELECT CASE (n)
   CASE (2)  ! Case of RK2
       xk1 = h* FCN(x0, y0)
        ym = y0 + xk1 
       xk2 = h* FCN(x1, ym)
       xk  = hlf*(xk1+ xk2)

   CASE (3)  ! Case of RK3
       xk1 = h* FCN(x0, y0)    
        ym = y0 + hlf*xk1
       xk2 = h* FCN(xh, ym)
        ym = y0- xk1+ 2.0d0*xk2
       xk3 = h* FCN(x1, ym)      
       xk  = (xk1+ 4.0d0*xk2+ xk3 )/6.0d0      

   CASE(4)   ! Case of RK4
       
       xk1 = h* FCN(x0, y0)
       print 11, "x0, y0, k1 ",x0,y0,xk1
        ym = y0 + hlf*xk1    
       xk2 = h* FCN(xh, ym)
       print 11, "xh, ym, k2 ",xh,ym,xk2       
        ym = y0 + hlf*xk2 
       xk3 = h* FCN(xh, ym)  
       print 11, "xh, ym, k3 ",xh,ym,xk3   
        ym = y0 + xk3
       xk4 = h* FCN(x1, ym)
       print 11, "x1, ym, k4 ",x1,ym,xk4 
       xk  = (xk1+ 2.0d0*xk2+ 2.0d0*xk3+ xk4)/6.0d0 
       pause
11     format(A11,F6.3,2x,f10.7,3x,f10.7)       
    CASE (:1) ! Case of n<=1
       WRITE(*,99) n
       STOP

    CASE (5:) ! Case of n>=5
       WRITE(*,99) n
       STOP

END SELECT 
y = y0 + xk      
x = x1
99 FORMAT(' PROGRAM DOES NOT HANDLE CASE OF N=',i3)
END SUBROUTINE DRV_RK



DOUBLE PRECISION FUNCTION FCN(x,y)
! ==============================================================================
! DESCRIPTION: A function subprogram providing y'=f(x,y)
!
! ARGUMENTS:
!      x, y  :: Real input values.
!
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x, y
  fcn = x/(y+1.0d0)
END FUNCTION fcn
    
   
DOUBLE PRECISION FUNCTION exact(x)
! ==============================================================================
! DESCRIPTION: A function subprogram providing the true solution y=(x) for 
!    testing the module. 
!
! ARGUMENTS:
!      x   :: A real input, independent variable.
!
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x 
   exact = -1.0d0 + dsqrt(x*x+9.0d0) 
END FUNCTION exact
