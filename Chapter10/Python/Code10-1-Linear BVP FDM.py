import numpy as np

def LBVP_SOLVE(neq, x, alpha, beta, gamma):
#  ==================================================================================
#  CODE10.1-LBVP_SOLVE.py. A Python module implementing Pseudocode 10.1.                          
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
#  DESCRIPTION: A module to find approximate solution of the following linear                  
#    differential equation using the Finite Difference Method:                                 
#          p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]                                   
#    subject to                                                                                
#          alpha1 * y'(a)+ beta1 * y(a) = gamma1                                               
#          alpha2 * y'(b)+ beta2 * y(b) = gamma2                                               
#                                                                                              
#  CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1.  
#                                                                                              
#  ON ENTRY                                                                                    
#    neq  :: Number of (equations) grid poinds;                                                
#     x   :: Array of length neq containing the abscissa of the grid points;                   
#    alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                   
#           the boundary conditions as stated above;                                           
#                                                                                              
#  ON RETURN                                                                                     
#     y   :: Array of length neq containing the approximate solution.                          
#                                                                                              
#  USES                                                                                        
#  COEFFS  :: A module containing the coefficients and rhs of the linear ODE, i.e., p(x), q(x), r(x), g(x).
#   TRIDIAGONAL:: A module solving a tridiagonal system of equations using the Thomas algorithm
#                                                                                              
#  REVISION DATE :: 03/08/2025                                                                 
#  ==================================================================================    
    nbc1 = int(alpha[0]); nbc2 = int(alpha[1])
    bb = np.zeros(neq); aa = np.zeros(neq)
    dd = np.zeros(neq); cy = np.zeros(neq)
    h = x[1] - x[0]
    h2 = h * h

    for k in range(neq):
        pxk, qxk, rxk, fxk = coeffs(x[k])
        dd[k] = rxk - 2.0 * pxk / h2
        aa[k] = pxk / h2 + 0.5 * qxk / h
        bb[k] = pxk / h2 - 0.5 * qxk / h
        cy[k] = fxk

    if nbc1 == 0:
        dd[0] = 1.0
        aa[0] = 0.0
        bb[0] = 0.0
        cy[0] = gamma[0] / beta[0]
    else:
        aa[0] = aa[0] + bb[0]
        dd[0] = dd[0] + 2.0 * h * bb[0] * beta[0] / alpha[0]
        cy[0] = cy[0] + 2.0 * h * bb[0] * gamma[0] / alpha[0]

    if nbc2 == 0:
        dd[neq-1] = 1.0
        aa[neq-1] = 0.0
        bb[neq-1] = 0.0
        cy[neq-1] = gamma[1] / beta[1]
    else:
        bb[neq-1] = aa[neq-1] + bb[neq - 1]
        dd[neq-1] = dd[neq-1] - 2.0 * h * aa[neq-1] * beta[1] / alpha[1]
        cy[neq-1] = cy[neq-1] - 2.0 * h * aa[neq-1] * gamma[1] / alpha[1]

    s1 = 0; sn = neq - 1
    y = tridiagonal(s1, sn, bb, dd, aa, cy)
    return y    # End of LBVP_SOLVE

def tridiagonal(s1, sn, b, d, a, c):
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
#     NumPy modules 
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
    
def coeffs(x):
# ==============================================================================
# DESCRIPTION: A user-defined module suppling the coefficients p(x), q(x), r(x) and
#    rhs g(x) of the linear ODE given in the following form: 
#         p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]
#
# ON ENTRY
#    x   :: Independent variable (a <=x<= b);
#
# ON RETURN
#   p, q, r, abd g :: Coefficients & rhs evaluated at x. 
#
# REVISION DATE :: 03/08/2025
# ==============================================================================
    p = x*x
    q = -5.0 * x
    r = 8.0
    f = 0.0
    return p, q, r, f

def exact(x):
# ==============================================================================
# DESCRIPTION: A function subprogram providing the true solution y=f(x) for testing the module. 
#
# ARGUMENTS:
#      x   :: A real input, independent variable.
# ==============================================================================
    return x ** 2 *( x**2 - 0.5 )



def test_lbvp_solve():
# ==============================================================================
#  The main program to test LBVP_SOLVE.PY
# ==============================================================================    
    n = int(input(" Enter no of intervals "))
    neq = n + 1 # no. of equations
    x = np.zeros(neq)
    y = np.zeros(neq)
    # Define the solution interval
    xa = 1.0;     xb = 2.0
    h = (xb - xa)/n
    for i in range(neq):
        x[i] = xa + i * h
        
    # SET THE BCs
    alpha = np.array([1.0, 1.0])
    beta = np.array([2.0,-2.0])
    gamma = np.array([4.0, 2.0])
 
    y = LBVP_SOLVE(neq, x, alpha, beta, gamma)

    print(f"{'      x':15s}{'Exact':10s}{'N.Approx':14s}{'Abs Error':16s}")
    for i in range(neq):
        error = abs(y[i] - exact(x[i]))
        print(f"{x[i]:10.4f}{exact(x[i]):12.7f}{y[i]:12.7f}{error:16.5e}")

if __name__ == "__main__":
    test_lbvp_solve()