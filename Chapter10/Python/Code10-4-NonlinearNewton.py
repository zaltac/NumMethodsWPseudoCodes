import numpy as np

def FCNS(x, y, yp):
# ==============================================================================
# DESCRIPTION: A user-defined function supplying the coefficients of the nonlinear 
#    two-point BVP given in the form: 
#              y'' = f(x,y,y')  on [a,b]
#
# ON ENTRY
#   x    :: Independent variable (a <=x<= b);
#   y    :: Dependent variable y=y(x).
#
# ON EXIT
#   f    :: The nonlinear two-point BVP f(x,y,y') evaluated at (x,y);
#   yp   :: First derivative of the dependent variable y';
#   dfdy :: Partial derivative of f wrt y evaluated at (x,y), df/dy;
#   dfdp :: Partial derivative of f wrt y' evaluated at (x,y), df/dy'.
#
# REVISION DATE :: 03/10/2025
# ==============================================================================
    fun = y - yp*yp/y
    dfdy = 1.0 + (yp/y)*(yp/y)
    dfdp = -2.0*yp/y
    return (fun, dfdy, dfdp)
    
def ENORM(n, c):
# ==================================================================================
# CODE3.1-ENORM.PY. A C++ module implementing ENORM module of Pseudocode 3.1.                                      
#
# NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
# First Edition. (c) By Zekeriya ALTAÇ (2024).
# ISBN: 978-1-032-75474-1 (hbk)
# ISBN: 978-1-032-75642-4 (pbk)
# ISBN: 978-1-003-47494-4 (ebk)
# 
# DOI : 10.1201/9781003474944
# C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
# 
# This free software is complimented by the author to accompany the textbook.
# E-mail: altacz@gmail.com.
# 
# DESCRIPTION: A function module to compute Euclidean (L2) norm of a vector.                  
#                                                                                             
# ARGUMENTS                                                                                   
#     n  :: The length of an input vector;                                                    
#     x  :: A vector (array) of length n.                                                     
#                                                                                             
# USES                                                                                        
#   sqrt :: Built-in NumPy library function returning the square root of a real value.            
#                                                                                             
# REVISION DATE :: 11/09/2024                                                                  
# ==================================================================================
    ENORM = 0.0
    for k in range(1, n+1):
        ENORM += c[k-1]**2
    return np.sqrt(ENORM)

def Exact(x):
# ==============================================================================
#  DESCRIPTION: A function subprogram providing the true solution y=f(x) for 
#     testing the sample problem. 
#
#  ARGUMENTS:
#      x   :: A double input, independent variable.
#
#  USES                                                                                        
#     sqrt :: Built-in NumPy library function returning the square root of a real value;
#     cosh :: Built-in NumPy library function returning the hyperbolic cosine of a real value.  
# ==============================================================================
    term = np.cosh(np.sqrt(2.0) * (1.0 - x))
    term = np.sqrt(term / np.cosh(np.sqrt(2.0)))
    return term

def TRIDIAGONAL(s1, sn, b, d, a, c):
#  ==================================================================================
#  CODE2.13-TRIDIAGONAL.py. A Python module implementing Pseudocode 2.13.                         
# 
#  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
#  First Edition. (c) By Zekeriya ALTAÇ (2024).
#  ISBN: 978-1-032-75474-1 (hbk)
#  ISBN: 978-1-032-75642-4 (pbk)
#  ISBN: 978-1-003-47494-4 (ebk)
#  
#  DOI : 10.1201/9781003474944
#  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
#  
#  This free software is complimented by the author to accompany the textbook.
#  E-mail: altacz@gmail.com.
#  
#  DESCRIPTION: A subroutine to solve a tridiagonal system of linear equations                 
#    using Thomas algorithm.                                                                   
#                                                                                              
#  ON ENTRY                                                                                    
#     s1 :: Subscript of the first unknown (usually 1);                                        
#     sn :: Subscript of the last unknown (usually no. of eqs, n);                             
#      b :: Array of length n containing coefficients of below diagonal elements;              
#      d :: Array of length n containing coefficients of diagonal elements;                    
#      a :: Array of length n containing coefficients of above diagonal elements;              
#      c :: Array of length n containing coefficients of rhs.                                  
#                                                                                              
#  ON RETURN                                                                                     
#      x :: An array of length n containing the solution.                                      
#     
#  USES
#     NumPy modules, ZEROS and RANGE. 
#                                                                                              
#  REVISION DATE :: 03/18/2024                                                                 
#  ==================================================================================
    x = np.zeros(sn-s1+1)
    for i in range(s1+1, sn+1):
        ratio = b[i] / d[i-1]
        d[i] = d[i] - ratio * a[i-1]
        c[i] = c[i] - ratio * c[i-1]

    x[sn] = c[sn] / d[sn]
    for i in range(sn-1, s1-1, -1):
        x[i] = (c[i] - a[i] * x[i+1]) / d[i]
    
    return x # End of module TRIDIAGONAL

def NONLINEAR_NEWTON(m, x, yo, y, tol, alpha, beta, gamma, maxit):
#  ==================================================================================
#  CODE10.4-NONLINEAR_NEWTON.py. A Python module implementing Pseudocode 10.4.                    
#  
#  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
#  First Edition. (c) By Zekeriya ALTAÇ (2024).
#  ISBN: 978-1-032-75474-1 (hbk)
#  ISBN: 978-1-032-75642-4 (pbk)
#  ISBN: 978-1-003-47494-4 (ebk)
#  
#  DOI : 10.1201/9781003474944
#  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
#  
#  This free software is complimented by the author to accompany the textbook.
#  E-mail: altacz@gmail.com.
#  
#   DESCRIPTION: A module to find approximate solution of a two-point nonlinear differential   
#     equation using the Newton's Method. The nonlinear equations is cast in the following form
#           y'' = f(x,y,y')  on [a,b]                                                          
#     subject to                                                                               
#           alpha1 * y'(a)+ beta1 * y(a) = gamma1                                              
#           alpha2 * y'(b)+ beta2 * y(b) = gamma2                                              
#                                                                                              
#   CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1. 
#                                                                                              
#   ON ENTRY                                                                                   
#      M   :: Number of (equations) grid poinds;                                               
#      x   :: Array of length M containing the abscissas of the grid points;                   
#      yo  :: Array of length M containing the initial guess for the solution;                 
#     alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                  
#            the boundary conditions as stated above;                                          
#     eps  :: Convergence tolerance;                                                           
#     maxit:: Maximum number of iterations permitted.                                          
#                                                                                              
#   ON EXIT                                                                                    
#      y   :: Array of length M containing the approximate solution.                           
#                                                                                              
#   USES                                                                                       
#     FCNS :: A user-defined external function module providing the coefficients of tridiagonal
#             system of linear equations;                                                      
#     ENORM:: A function module to calculate the Euclidean vector (L2 norm) of a vector;       
#     TRIDIAGONAL:: A module to solve a tridiagonal system of equations with Thomas algorithm.
#                                                                                              
#   REVISION DATE :: 03/10/2025                                                                
#  ==================================================================================
    A = np.zeros(m); B = np.zeros(m)
    D = np.zeros(m); C = np.zeros(m)
    X = np.zeros(m); y = np.zeros(m)

    h = x[1] - x[0]
    hsqr = h*h; bb2h = 0.5/h;  hb2 = 0.5*h

    if int(abs(alpha[0]))==0:
        bc_left = 0
    else:
        bc_left = 1
        c1d = 2.0*h* beta(1) / alpha(1)
        c1r = 2.0*h* gamma(1) / alpha(1)
 
    if int(abs(alpha[1]))==0:
        bc_rigt = 0
    else:
        bc_rigt = 1
        c2d = 2.0*h*beta[1]/alpha[1]
        c2r = 2.0*h*gamma[1]/alpha[1]

    error = 1.0
    p = 0
    y = yo.copy()
    
    while error > tol and p < maxit:
        for k in range(m):
            if k == 0:
                if bc_left == 0:
                    D[k] = 1.0
                    yo[k] = gamma[0]/beta[0]
                    A[k] = 0.0
                    B[k] = 0.0
                    C[k] = 0.0
                else:
                    yx = (gamma[0] - beta[0]*yo[k])/alpha[0]
                    fun, dfdy, dfdp = FCNS(x[k], yo[k], yx)
                    A[k] = 2.0
                    B[k] = 0.0
                    D[k] = -2.0 + c1d - hsqr*(dfdy - dfdp*beta[0]*yo[k]/alpha[0])
                    C[k] = -c1r + (-2.0+c1d)*yo[k] + 2.0*yo[k+1] - hsqr*fun
            elif k == m-1:
                if bc_rigt == 0:
                    D[k] = 1.0
                    yo[k] = gamma[1]/beta[1]
                    A[k] = 0.0
                    B[k] = 0.0
                    C[k] = 0.0
                else:
                    yx = (gamma[1] - beta[1]*yo[k])/alpha[1]
                    fun, dfdy, dfdp = FCNS(x[k], yo[k], yx)
                    B[k] = 2.0
                    D[k] = -2.0 - c2d - hsqr*(dfdy - dfdp*beta[1]*yo[k]/alpha[1])
                    A[k] = 0.0
                    C[k] = c2r - (2.0+c2d)*yo[k] + 2.0*yo[k-1] - hsqr*fun
            else:
                yx = bb2h*(yo[k+1] - yo[k-1])
                yxx = yo[k+1] - 2.0*yo[k] + yo[k-1]
                fun, dfdy, dfdp = FCNS(x[k], yo[k], yx)
                B[k] = 1.0 + hb2 * dfdp
                D[k] = -2.0 - hsqr* dfdy
                A[k] = 1.0 - hb2 * dfdp
                C[k] = yxx - hsqr*fun

        s1=0; sn=m-1
        C = TRIDIAGONAL(s1, sn, B, D, A, C)
        error = ENORM(m, C)
        print(f"Iter={p:3d}, E-norm={error:14.5e}")
        y[:] = yo[:] - C[:]
        p += 1
        yo[:] = y[:]
    
    return y


def test_NONLINEAR_NEWTON():
# ==============================================================================
#  The main program to test the module NONLINEAR_NEWTON.PY
# ==============================================================================
    n = int(input(" Enter no of intervals "))
    neq = n + 1 # no. of equations

    # Initialize grid & discrete solution    
    x = np.zeros(neq)
    y = np.zeros(neq)
    yo = np.zeros(neq)
    
    # Define the solution interval and grids
    xa = 0.0;     xb = 1.0
    h = (xb - xa)/n
    for i in range(neq):
        x[i] = xa + i * h
        
    # Set BCs
    alpha = np.array([0.0, 1.0])
    beta = np.array([1.0, 0.0])
    gamma = np.array([1.0, 0.0])
    bc_left = int(abs(alpha[0]))
    bc_rigt = int(abs(alpha[1]))
    
    # Initialize initial guess & iteration parameters
    tol = 1.0e-6
    guess = 0.4
    maxit = 99

    # IUse linear interpolation to generate initial guess when both BCs are Dirichlet type 
    for i in range(neq):
        yo[i] = guess
        if (bc_left == 0) and (bc_rigt == 0):
            yo[i] = gamma[0] + (gamma[1]-gamma[0])*(x[i]-xa)/(xb-xa)
    # Go to the solver to solve nonlinear ODE
    y = NONLINEAR_NEWTON(neq, x, yo, y, tol, alpha, beta, gamma, maxit)
    
    # print out the results 
    print(f"{'      x':15s}{'Exact':10s}{'N.Approx':14s}{'Abs Error':16s}")
    for i in range(neq):
        yex  = Exact(x[i])
        aerr = abs(yex - y[i])
        print(f"{x[i]:10.4f} {yex:11.7f} {y[i]:11.7f} {aerr:14.5e}")

if __name__ == "__main__":
    test_NONLINEAR_NEWTON()