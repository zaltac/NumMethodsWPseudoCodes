import numpy as np
# ==================================================================================
# Main program to test BISECTION.PY
# ==================================================================================  
def test_bisection_method():
    
    a = 0.0
    b = 4.0
    maxit = 99
    eps = 0.50e-4
 
    if func(a) * func(b) > 0:
        print('No root in interval (a,b). Change the interval.')
        exit()
    else:
        (root,iter,tol)=bisection(a, b, maxit, eps)
        print("\n The root is ",f"{root:10.6f}")
        print(" Tolerance achieved ",f"{tol:10.4E}")
        print(" No. of bisections performed is",f"{iter:4d}")          
 

def bisection(a, b, maxit, eps):
    """
  ==================================================================================
  CODE4.1-BISECTION.py. A Python module implementing Pseudocode 4.1.                             
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTAÇ (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
  DESCRIPTION: A module to find a root of a nonlinear equation in [a,b]                       
    using the Bisection method.                                                               
                                                                                              
  ON ENTRY                                                                                    
    [a,b] :: Initial search interval (it must bracket one root);                               
    maxit :: Maximum number of iterations permitted;                                           
    eps   :: Convergence tolerance.                                                            
                                                                                              
  ON RETURN                                                                                     
    halves:: Number of halves realized;                                                        
    root  :: Computed approximation for the root;                                              
    tol   :: accuracy tolerance achieved.   
                                                                                              
  USES                                                                                        
    abs   :: Built-in Intrinsic function returning the absolute value of a real value;         
                                                                                              
  ALSO REQUIRED                                                                               
    func  :: User-defined external function providing the nonlinear equation.                  
                                                                                              
  REVISION DATE :: 11/20/2024                                                                 
  ==================================================================================
    """

    fa = func(a)
    fb = func(b)
    interval = b - a

    print(f"{'p':10} {'a':12} {'b':10} {'f(a)':12} {'f(b)':12} {'xm':11} {'f(xm)':10} {'interval':12}")
    print("-"*97)
    
    for p in range(1, maxit+1):
        xm = 0.5 * (a + b)
        fm = func(xm)
        print(f"{p:3} {a:12.7f} {b:12.7f} {fa:12.4e} {fb:12.4e} {xm:12.7f} {fm:12.4e} {interval:12.4e}")

        if fa * fm > 0.0:
            a = xm
            fa = fm
        else:
            b = xm
            fb = fm

        interval = 0.5 * interval
        if (abs(fm) < eps) and (interval < eps):
            break
    
    root = xm
    halves=p
    tol =interval
    
    if interval<abs(fm):
        tol = abs(fm)
        
    if  p==maxit: 
        print(f"\n{'!'*37} Max iteration number reached={p} {'!'*37}")
    
    return (root, halves, tol)
 
 

# ==========================================================================          
# User-defined function providing f(x), which should be cast as func(x)=0.
# ==========================================================================  
def func(x):
    return x * x + 0.025 * x - 4.0

if __name__ == "__main__":
    test_bisection_method()