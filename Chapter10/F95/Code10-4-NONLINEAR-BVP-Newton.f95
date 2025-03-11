PROGRAM Test_NonLinear_Newton_Method
! ==============================================================================
!  The main program to test SUBROUTINE NONLINEAR_NEWTON.F95
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: nmax=201
DOUBLE PRECISION :: x(nmax), y(nmax), yo(nmax), alpha(2), beta(2), gamma(2)  
DOUBLE PRECISION :: xa, xb, h, eps, guess, error, EXACT
INTEGER :: neq, n, i, maxit, bc_left, bc_rigt 
!  
!
PRINT*, "Enter No of intervals "
READ*, n
IF(n>nmax) THEN
   PRINT*, "Max grid size is ",nmax
   PRINT*, "Increase NMAX, and then try it again."
   STOP
ENDIF
!
! Setup the ODE, Grids & BCs 
eps  = 1.0D-6
guess= 0.4d0
maxit= 99
neq  = n + 1 
xa = 0.0d0   ; xb = 1.0d0   
h  = (xb-xa)/dfloat(n)
DO i=1,neq
   x (i) = xa + DFLOAT(i-1)*h
ENDDO 
alpha(1)= 0.0d0;  beta(1)= 1.0d0;  gamma(1)= 1.0d0 
alpha(2)= 1.0d0;  beta(2)= 0.0d0;  gamma(2)= 0.0d0
!  
! Prep initial guess for the solution
bc_left= INT(dabs(alpha(1))) 
bc_rigt= INT(dabs(alpha(2)))
DO i=1,neq
   x (i) = xa + DFLOAT(i-1)*h
   yo(i) = guess
   IF((bc_left.eq.0).AND.(bc_rigt.eq.0)) THEN
        yo(i)=gamma(1) +(gamma(2)-gamma(1))*(x(i)-xa)/(xb-xa)
   ENDIF 
ENDDO  
! ***
CALL NONLINEAR_NEWTON(neq,eps,x,yo,y,alpha,beta,gamma,maxit)
! 
! PRINT OUT THE RESULTS 
PRINT 1,   
DO i=1,neq
   error=dabs(Exact(x(i))-y(i))
   PRINT 2, x(i),Exact(x(i)),y(i),error 
ENDDO
PRINT*,  
PRINT*, "*** D O N E ***"
PRINT*,         
3 FORMAT(5x,"Iter=",i4,3x,"E-norm=",1PE14.5)
1 FORMAT(//,8x,"x",8x,"Exact",6x,"N.Approx",3x,"Abs Error")
2 FORMAT(2x,f10.5,2f12.7,1PE12.3)  
END PROGRAM Test_NonLinear_Newton_Method



SUBROUTINE NONLINEAR_NEWTON(M,eps,x,yo,y,alpha,beta,gamma,maxit)
!  ==================================================================================
!  CODE10.4-NONLINEAR_NEWTON.f95. A Fortran95 module implementing Pseudocode 10.4.                
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
!   DESCRIPTION: A module to find approximate solution of a two-point nonlinear differential   
!     equation using the Newton's Method. The nonlinear equations is cast in the following form
!           y'' = f(x,y,y')  on [a,b]                                                          
!     subject to                                                                               
!           alpha1 * y'(a)+ beta1 * y(a) = gamma1                                              
!           alpha2 * y'(b)+ beta2 * y(b) = gamma2                                              
!                                                                                              
!   CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1. 
!                                                                                              
!   ON ENTRY                                                                                   
!      M   :: Number of (equations) grid poinds;                                               
!      x   :: Array of length M containing the abscissas of the grid points;                   
!      yo  :: Array of length M containing the initial guess for the solution;                 
!     alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                  
!            the boundary conditions as stated above;                                          
!     eps  :: Convergence tolerance;                                                           
!     maxit:: Maximum number of iterations permitted.                                          
!                                                                                              
!   ON EXIT                                                                                    
!      y   :: Array of length M containing the approximate solution.                           
!                                                                                              
!   USES                                                                                       
!     FUNCS :: A user-defined external function module providing the coefficients of the 
!             nonlinear two-point BVP; 
!     ENORM:: A function module to calculate the Euclidean vector (L2 norm) of a vector;       
!     TRIDIAGONAL :: A module to solve a tridiagonal system of equations with Thomas algorithm.
!                                                                                              
!   REVISION DATE :: 03/10/2025                                                                
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: m, maxit
DOUBLE PRECISION, INTENT( IN) :: x(m), alpha(2), beta(2), gamma(2), eps
DOUBLE PRECISION, INTENT(INOUT) :: yo(m), y(m)
DOUBLE PRECISION :: b(m), a(m), d(m), c(m)
DOUBLE PRECISION :: h, error, bb2h, hb2, c1d, c1r, c2d, c2r  
DOUBLE PRECISION :: hsqr, yx, yxx, f, fy, fp, ENORM  
INTEGER ::   BC_left, BC_rigt, p, k
!     
bc_left= INT(dabs(alpha(1))) 
bc_rigt= INT(dabs(alpha(2)))
!    
 h  = x(2) -x(1)  ; hsqr = h*h 
bb2h= 0.5d0/h     ; hb2  = 0.50d0*h
!      
error= 1.0d0
p  = 0
y = yo  ! initialize numerical solution
DO 
   DO k = 1, m
      IF(k.eq.1) THEN
         IF(BC_left==0) THEN
             d(k)= 1.0d0
             yo(k)= gamma(1)/beta(1)
             a(k)= 0.0d0
             b(k)= 0.0d0
             c(k)= 0.0d0
         ELSE    
             c1d = 2.0d0*h*beta (1)/alpha(1)
             c1r = 2.0d0*h*gamma(1)/alpha(1)
             yx  = (gamma(1)-beta(1)*yo(k))/alpha(1)  
             CALL FUNCS(x(k),yo(k),yx,f,fy,fp)
             a(k) = 2.0d0
             b(k) = 0.0d0 
             d(k) =-2.0d0+c1d- hsqr*(fy-fp*beta(1)*yo(k)/alpha(1)) 
             c(k) =-c1r + (-2.0d0+c1d)*yo(k) +2.0d0*yo(k+1) -hsqr*f 
         ENDIF 
       ELSE IF(k==m) THEN
         IF(BC_rigt==0) THEN
             d(k)= 1.0d0
             yo(k)= gamma(2)/beta(2)
             a(k)= 0.0d0
             b(k)= 0.0d0
             c(k)= 0.0d0
         ELSE     
             c2d = 2.0d0*h*beta (2)/alpha(2)
             c2r = 2.0d0*h*gamma(2)/alpha(2) 
             yx  = (gamma(2)-beta(2)*yo(k))/alpha(2)  
             CALL FUNCS(x(k),yo(k),yx,f,fy,fp) 
             b(k)= 2.0d0  
             d(k)= -2.0d0-c2d-hsqr*(fy-fp*beta(2)*yo(k)/alpha(2))
             a(k)= 0.0d0
             c(k)= c2r -(2.0d0+c2d)*yo(k) +2.0d0*yo(k-1) -hsqr*f   
          ENDIF 
       ELSE
              yx  = (yo(k+1)-yo(k-1))/(2.0d0*h)
             yxx  = yo(k+1) -2.0d0*yo(k) + yo(k-1)             
             CALL FUNCS(x(k),yo(k),yx,f,fy,fp)
             b(k) = 1.0d0 + hb2 * fp
             d(k) =-2.0d0 - hsqr* fy
             a(k) = 1.0d0 - hb2 * fp
             c(k) = yxx - hsqr*f          
       ENDIF       
   ENDDO
! 
! === SOLVE TRIDIAGONAL SYSTEM OF EQUATIONS FOR CORRECTION    
   CALL TRIDIAGONAL(1,m,b,d,a,c,c)
   error = ENORM(m,c)
   DO k = 1, m
      y(k) = yo(k) - c(k)    ! Apply correction to prev estimates
   ENDDO
   p = p + 1           
   PRINT 3, p,error
   DO k = 1, m
      yo(k) = y(k)
   ENDDO
   IF (error<eps .OR. p==maxit) EXIT
ENDDO
!       
3 FORMAT(5x,"Iter=",i3,3x,"E-norm=",1PE14.5)
END SUBROUTINE 

 
 
SUBROUTINE FUNCS(x,y,yp,f,dfdy,dfdp)
! ==============================================================================
! DESCRIPTION: A user-defined function supplying the coefficients of the nonlinear 
!    two-point BVP given in the form: 
!              y'' = f(x,y,y')  on [a,b]
!
! ON ENTRY
!   x    :: Independent variable (a<= x<= b);
!   y    :: Dependent variable y=y(x).
!
! ON EXIT
!   yp   :: First derivative of the dependent variable y'=y'(x);
!   f    :: The nonlinear ode defines as f(x,y,y') at (x,y);
!   dfdy :: Partial derivative of f wrt y at (x,y), df/dy;
!   dfdp :: Partial derivative of f wrt y' at (x,y), df/dy';
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT( IN) :: x, y, yp
DOUBLE PRECISION, INTENT(OUT) :: f, dfdy, dfdp
f    = y - yp * yp / y
dfdy = 1.0d0 + (yp / y) ** 2
dfdp = -2.0d0 * yp / y 
END SUBROUTINE FUNCS  

SUBROUTINE TRIDIAGONAL(s1, sn, b, d, a, c, x)
! ==============================================================================
! CODE2-13-TRIDIAGONAL.F95. A fortran program for implementing Pseudocode 2.13
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
! DESCRIPTION: A subroutine to solve a tridiagonal system of linear equations
!   using Thomas algorithm. 
!
! ON ENTRY
!    s1 :: Subscript of the first unknown (usually 1); 
!    sn :: Subscript of the last unknown (usually No. of Eqs); 
!     b :: Array of length n containing coefficients of below diagonal elements;
!     d :: Array of length n containing coefficients of diagonal elements;
!     a :: Array of length n containing coefficients of above diagonal elements;
!     c :: Array of length n containing coefficients of rhs.
!
! ON EXIT
!     x :: An array of length n containing the solution.
!
! REVISION DATE :: 03/18/2024
! ==============================================================================
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

DOUBLE PRECISION FUNCTION ENORM(n,x)
!  ==================================================================================
!  CODE3.1-JACOBI.F95. A Fortran95 module implementing ENORM module of Pseudocode 3.1.                            
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
!  DESCRIPTION: A function module to compute Euclidean (L2) norm of a vector.                  
!                                                                                              
!  ARGUMENTS                                                                                   
!      n  :: The length of an input vector;                                                    
!      x  :: A vector (array) of length n.                                                     
!                                                                                              
!  USES                                                                                        
!    SQRT :: Built-in Intrinsic function returning the square root of a real value.            
!                                                                                              
!  REVISION DATE :: 11/09/2024                                                                 
!  ==================================================================================
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: x
DOUBLE PRECISION :: delta  
INTEGER :: i
!    
delta = 0.0d0
DO i = 1, n
   delta = delta + x(i) * x(i)
ENDDO 
ENORM = DSQRT(delta)      
END FUNCTION ENORM

 
 
DOUBLE PRECISION FUNCTION Exact(x)
! ==============================================================================
! DESCRIPTION: A function subprogram providing the true solution y=f(x) for 
!    testing the sample problem. 
!
! ARGUMENTS:
!     x   :: A double input, independent variable.
!
! USES                                                                                        
!   dsqrt :: Built-in intrinsic function returning the square root of a real value in Double;
!   dcosh :: Built-in Intrinsic function returning the hyperbolic cosine of a real value in Double. 
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x
!
Exact = dCOSH(dSQRT(2.0d0) * (1.0d0 - x))
Exact = dSQRT(Exact / DCosh(dSQRT(2.0d0)))
!   
END FUNCTION Exact  
 
 
