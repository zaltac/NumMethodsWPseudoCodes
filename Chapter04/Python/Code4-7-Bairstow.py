import numpy as np


#  ==================================================================================
#  CODE4.7-BAIRSTOW.PY. A Python module implementing Pseudocode 4.7.                              
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
#  DESCRIPTION: A module to find all real and/or imaginary roots of a polynomial           
#    of the n'th degree using the BAIRSTOW's method.                                         
#                                                                                              
#  ON ENTRY                                                                                    
#    n    :: Degree of the polynomial;                                                         
#   p0,q0 :: Initial guesses for a quadratic equation; i.e., for p and q;                      
#     a   :: Array of length (n+1) containing the coefficients of polynomial defined as        
#                 a0 x^n + a1 x^(n-1) + ... + an = 0                                           
#    eps  :: Convergence tolerance;                                                            
#   maxit :: Maximum number of iterations permitted;                                           
#   iprnt :: printing key, =0 do not print intermediate results, <> 0 print intermediates.     
#                                                                                              
#  ON EXIT                                                                                     
#    xre  :: Array of length n containing real parts of the roots;                             
#    xim  :: Array of length n containing imaginary parts of the roots.                        
#                                                                                              
#  OTHER VARIABLES                                                                             
#     b   :: Array of length [n] containing coefficients of quotient polynomial (0<=k<=n-2);   
#     c   :: Array of length [n] containing coefficients of partial derivatives.               
#                                                                                              
#  USES                                                                                        
#    abs  :: Built-in NumPy function returning the absolute value of a real value;         
#    QUADRATIC :: A module that solves a quadratic equation of the form x2 + p x + q = 0. (see CODE1-3)   
#                                                                                              
#  REVISION DATE :: 11/20/2024                                                                 
#  ==================================================================================
def bairstow(n, p0, q0, a, eps, maxit, iprnt):
    """
    DESCRIPTION: A module to compute real and/or imaginary roots of a polynomial of nth degree using the BAIRSTOW’s method.

    On ENTRY
        n   :: Degree of polynomial;
    p0,q0   :: Initial guesses for quadratic coefficients: p and q,
        a   :: Array of length [n+1] containing coefficients of polynomial defined as
                a0 x^n + a1 x^(n-1) + ... + an = 0 
        eps :: Convergence tolerance;
      maxit :: Maximum number of iterations permitted;
      iprnt :: Output control key, 
                = 0 does not print intermediate results, 
                = 1 prints a short intermediate results, 
                = 2 prints all possible results.  
    On RETURN 
        xre  :: Array of length [n] containing real part of the roots 
        xim  :: Array of length [n] containing imaginary part of the roots 

    OTHER VARIABLES
        b   :: Array of length [n] containing coefficients of quotient polynomial (0<=k<=n-2)
        c   :: Array of length [n] containing coefficients of partial derivatives 

    USES
        quadratic_eq, abs
   Last-modified: April 29, 2024. 
    """
    xr = np.zeros(2);  xi = np.zeros(2)
    xre= np.zeros(n);  xim= np.zeros(n)
    a = a / a[0]  # Normalize a's by a[0]
    m = n
    kount = 0
    while n > 1:
        k = 0
        p = p0
        q = q0
        delM = 1.0
        while delM > eps and k <= maxit:
            k += 1
            b = np.zeros(n+1)
            c = np.zeros(n+1)
            b[0] = 1.0
            c[0] = 1.0
            b[1] = a[1] - p
            c[1] = b[1] - p
            for i in range(2, n+1):
                b[i] = a[i] - p*b[i-1] - q*b[i-2]
                c[i] = b[i] - p*c[i-1] - q*c[i-2]
            cbar = c[n-1] - b[n-1]
            del1 = b[n-1]*c[n-2] - b[n]*c[n-3]
            del2 = b[n]*c[n-2] - b[n-1]*cbar
            delt = c[n-2]*c[n-2] - cbar*c[n-3]
            delp = del1 / delt
            delq = del2 / delt
            p += delp
            q += delq
            delM = abs(delp) + abs(delq)
            if iprnt == 1:
                print(f" Iter= {k:2d}, delM= {delM:10.4g}, p= {p:12.5g}, q= {q:12.5g}")
            elif iprnt == 2:
                print(f"\n iter={k:4d}")
                print(f"---------")
                print(f" dp ={delp:14.6g}, dq ={delq:14.6g}, delM={delM:14.6g}")
                print(f"  p ={p:14.6g},  q ={q:14.6g}")
                print(f"{'-'*47}")
                print(f" k{' '*10}a(k){' '*11}b(k){' '*11}c(k)")
                for i in range(n+1):
                    print(f"{i:5d} {a[i]:12.6f} {b[i]:12.6f} {c[i]:12.6f}")
                print(f"{'-'*47}")
        if k-1 == maxit:
            print("Quadratic factor did not converge after", k-1, "iterations")
            print("Recent values of p, q, delM are", p, q, delM)
            print("Corresponding roots may be questionable ...")

        (xr, xi) = quadratic_eq(p, q)
        xre[kount] = xr[0]
        xim[kount] = xi[0]
        kount += 1
        xre[kount] = xr[1]
        xim[kount] = xi[1]
        kount += 1
        print()
        print(f"======== FOUND A QUADRATIC FACTOR ========")
        print(f"x*x + ({p:12.6f})*x + ({q:12.6f})")
        print(f"==========================================")
        n -= 2
        a[:n+1] = b[:n+1]
        if n == 1:
            xre[kount] = -a[1]
            xim[kount] = 0.0
            print(f"======== FOUND A LINEAR FACTOR ========")
            print(f"     x  + ({a[1]:12.6f})")
            print(f"========================================")
    n = m
    return (xre, xim)


# ==============================================================================
#  Main program to test BAIRSTOW module
# ==============================================================================
def test_bairstow():

    n = 5
    a = np.array([1.0, -5.0, -15.0, 85.0, -26.0, -120.0])
    xre = np.zeros(n)
    xim = np.zeros(n)
    iprnt = 2
    maxit = 99
    p0 = 0.0
    q0 = 0.0
    eps = 0.5e-4
    
    (xre, xim) = bairstow(n, p0, q0, a, eps, maxit, iprnt)
    
    print("    ======== All the Roots are =========")
    for i in range(n):
        print(f"Root({i+1:2d}) = {xre[i]:8.5f} + ({xim[i]:8.5f}) i")
    print("    ====================================")

if __name__ == "__main__":
    test_bairstow()
