import numpy as np

#  ==================================================================================
#  CODE3.1-JACOBI.PY. A Python module implementing JACOBI_DRV module of Pseudocode 3.1.                                
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
#  DESCRIPTION: A subroutine to to perform one step Jacobi iteration and compute               
#     the Euclidean norm of the displacement vector.                                           
#                                                                                              
#  ON ENTRY                                                                                    
#     n   :: Number of equations (size of A);                                                  
#     A   :: Input coefficient matrix (nxn);                                                  
#     b   :: Array of length n containing the right-hand;                                      
#     x   :: Array of length n containing the estimate at (p+1)'th step;                       
#     xo  :: Array of length n containing the estimate at p'th step.                           
#                                                                                              
#  ON EXIT                                                                                     
#     x   :: Array of length n containing the estimated solution;                              
#     del :: Maximum absolute error achieved.                                                  
#                                                                                              
#  USES                                                                                        
#   ENORM:: User-defined function calculating the Euclidean vector (L2 norm) of a vector.      
#                                                                                              
#  REVISION DATE :: 11/09/2024                                                                 
#  ==================================================================================
def jacobi_drv(n, a, b, xo):  
    d = np.zeros(n)  
    x = np.zeros(n)  
    for i in range(n):
        sums = 0.0
        for j in range(n):
            if(i!=j):
                sums += a[i, j] * xo[j]
    
        x[i] =  (b[i] - sums) / a[i, i]
        d[i] = x[i] - xo[i]
    
    delta = ENORM(d)
    return  (delta, x) # End of module JACOBI_DRV

#  ==================================================================================
#  CODE3.1-JACOBI.PY. A Python module implementing Pseudocode 3.1.                                
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
#  DESCRIPTION: A module to iteratively solve Ax=b using the Jacobi method.                
#                                                                                              
#  ON ENTRY                                                                                    
#     n   :: Number of equations (size of A);                                                  
#     A   :: Input coefficient matrix (nxn);                                                  
#     b   :: Array of length n containing the right-hand;                                      
#     x   :: Array of length n containing the estimate at (p+1)'th step;                       
#     xo  :: Array of length n containing the initial guess, or iterates at estimate at p'th st
#    eps  :: Convergence tolerance;                                                            
#   maxit :: Maximum permitted number of iterations.                                           
#                                                                                              
#  ON EXIT                                                                                     
#     x   :: Array of length n containing the estimated solution;                              
#   iter  :: Total number of iterations performed;                                             
#   error :: L2 norm of the displacement vector.                                               
#                                                                                              
#  USES                                                                                        
#    JACOBI_DRV :: Accompanying module performing one step Jacobi iteration.               
#                                                                                              
#  REVISION DATE :: 11/09/2024                                                                 
#  ==================================================================================
def jacobi(n, eps, a, b, xo, maxit):
    del0 = 1.0
    p = 0
    x = xo.copy()
    while (del0>eps and p<maxit):
        p += 1
        (del1,x)=jacobi_drv(n, a, b, xo)
        # printout the iteration progress 
        print(f" iter={p:5}  delta = {del1:8.2e}")
        xo[:] = x
        del0 = del1
    
    if p == maxit:
        print(f"\nJacobi method failed to converge after {maxit} iterations\nwithin the specified EPS tolerance.")
    error=del1
    iter=p
    return (iter, error, x) # End of module jacobi


def ENORM(x):
    # ==================================================================================
    # CODE3.1-ENORM.PY of a module in Pseudocode 3.1.                                      
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
    # DESCRIPTION: A Python function module to compute Euclidean (L2) norm of a vector.                  
    #                                                                                             
    # ARGUMENTS                                                                                   
    #     n  :: The length of an input vector;                                                    
    #     x  :: A vector (array) of length n.                                                     
    #                                                                                             
    # USES                                                                                        
    #   sqrt :: Built-in NumPy library function returning the square root of a real value.      
    #   len, range  :: Built-in NumPy library functions    
    #                                                                                             
    # REVISION DATE :: 11/09/2024                                                                  
    # ==================================================================================
    n=len(x)
    delta = 0.0
    for i in range(n):
        delta += x[i] * x[i]
    
    return np.sqrt(delta)  # End of E-norm module


def test_jacobi_method():
    # ==============================================================================
    #  The main program to test JACOBI.PY
    # ==============================================================================
    n = 3
    a = np.zeros((n, n))  # Coefficient matrix
    b = np.zeros(n)  # rhs vector
    xo = np.zeros(n)  # initial guess at start
    x = np.zeros(n)  # solution
    
    # Prepare input data (linear system)
    b[0] = 5.0; b[1] = 2.0 ; b[2] = 6.0
    a[0, 0] = 4.0; a[0, 1] = 1.0; a[0, 2] = 1.0
    a[1, 0] = 2.0; a[1, 1] = 4.0; a[1, 2] = 2.0
    a[2, 0] = 1.0; a[2, 1] = 3.0; a[2, 2] = 4.0

    print("* Coefficient Matrix *", " * RHS *", " * Initial Guess *")
    for i in range(n):
         print(f"{a[i, :]}, {b[i]}, {xo[i]}")

    eps = .5e-4
    maxit = 99
  
    
    (iter, error, x)=jacobi(n, eps, a, b, xo, maxit)

    print("\n------- Solution ---------------\n")
    for i in range(n):
        print(" x(",i+1,") = ", f"{x[i]:10.6f}")
    print("\n-----------------------")
    print("Total no of iterations = ",iter,"Maximum Error = ",error)

if __name__ == "__main__":
    test_jacobi_method()  
