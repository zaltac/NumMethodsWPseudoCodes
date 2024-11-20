import numpy as np

# ==============================================================================
#  Main program to test module NEWTON_RAPHSON.C
# ==============================================================================
def test_newton_raphson():
    
    root = 3.0
    maxit =199
    eps = 0.50e-5
 
    (root,iter) = newton_raphson(root,maxit,eps)
    print("\n The root is ",f"{root:10.8f}")
    print(" No. of iterations performed is",f"{iter:4d}")          
 

def newton_raphson(root, maxit, eps):
    """
  ==================================================================================
  CODE4.3-NEWTON_RAPHSON.py. A Python module implementing Pseudocode 4.3.                        
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTAÇ (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
  DESCRIPTION: A module to compute a root of a nonlinear equation using the                   
    Newton-Raphson method.                                                                    
                                                                                              
  ON ENTRY                                                                                    
    root  :: Initial guess for the root;                                                      
    maxit :: Maximum number of iterations permitted;                                          
    eps   :: Convergence tolerance.                                                           
                                                                                              
  ON RETURN                                                                                     
    iter  :: Number of iterations realized;                                                   
    root  :: Computed approximation for the root.                                             
                                                                                              
  USES                                                                                        
    ABS   :: Built-in Intrinsic function returning the absolute value of a real value;        
                                                                                              
  ALSO REQUIRED                                                                               
    FUNC  :: User-defined external function providing the nonlinear equation, f(x).           
    FUNCP :: User-defined external function providing the first derivative                    
             of the nonlinear equation, f'(x).                                                
                                                                                              
  REVISION DATE :: 11/20/2024                                                                 
  ==================================================================================    
    """

    p = 0
    x0 = root
    del0 = 1.0
    print(f"{'  p':8} {'x^(p)':12} {'f(x^(p))':12} {'fx(x^(p))':14} {'Error':12} {'Conv.Rate':11}")
 
    while True:
        fn = func(x0)
        fpn = funcx(x0)
        deln = -fn / fpn
        aerr = abs(deln)
        rate = aerr / del0 ** 2
        print(f"{p:3d} {x0:13.7f} {fn:13.4E} {fpn:13.7f} {aerr:13.4E} {rate:13.7f}")
        xn = x0 + deln
        x0 = xn
        del0 = abs(deln)
        p += 1
        if (abs(fn) < eps and aerr < eps) or (p == maxit):
            break

    root = xn
    iter = p

    if p == maxit:
        print("** Max iteration number reached=")
        print(f"** Estimated root has NOT converged, del, f(x)= {abs(deln):.7f}, {func(xn):.7f}")

    return (root, iter) # End newton_raphson module
 
 

# ==========================================================================          
# User-defined function providing f(x), which should be cast as func(x)=0.
# ========================================================================== 
def func(x):
    return x**3 - 3.8


# ==========================================================================          
# User-defined function providing f'(x)=funcp(x).
# ========================================================================== 
def funcx(x):
    return 3.0*x**2
    
    
if __name__ == "__main__":
    test_newton_raphson() 