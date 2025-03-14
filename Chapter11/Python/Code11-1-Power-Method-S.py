import numpy as np

def POWER_METHOD_S(n, a, x, eigeno, eps, maxit):
#  ==================================================================================
#  CODE11.1-POWER_METHOD_S.py. A Python module implementing Pseudocode 11.1.                      
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
#   DESCRIPTION: A module to find the dominant eigenvalue (largest in absolute value)          
#      using the Power Method with scaling technique.                                          
#                                                                                              
#   ON ENTRY                                                                                   
#      n     :: Size of the matrix;                                                            
#      A     :: A real square matrix (nxn);                                                    
#      x     :: Array of length n containing the initial guess for the eigenvector;            
#     lambda :: An initial guess for the dominant eigenvalue;                                  
#      eps   :: Convergence tolerance;                                                         
#     maxit  :: Maximum number of iterations permitted.                                        
#                                                                                              
#   ON EXIT                                                                                    
#     lambda :: Estimated dominant (largest in absolute value) eigenvalue;                     
#      x     :: Array of length n containing the estimated eigenvector;                        
#     error  :: Error, max of both L2-norm of displacement vector and relative error           
#                for eigenvalue.                                                               
#                                                                                              
#   USES                                                                                       
#     abs     :: Built-in NumPy library function returning the absolute value of a real value;
#     dot     :: Built-in NumPy library function returning the A.x=b;
#     AX      :: A module for computing Ax matrix-vector product (see Pseudocode 2.5);          
#     MAX_SIZE:: A function module providing largest (in absolute) value of a vector.          
#                                                                                              
#   REVISION DATE :: 03/14/2025                                                                
#  ==================================================================================
    x = np.array(x)
    a = np.array(a)
    xn = np.zeros_like(x)
    eigen = eigeno
    error = MAX_SIZE(np.abs(x))
    
    p = 0
    print(f"iter  eigenvalue   Abs error")
    while error>eps and p<maxit:
        p += 1
        xn = np.dot(a, x)
        eigen = MAX_SIZE(xn)
        xn /= eigen
        err1 = ENORM(np.abs(xn - x))
        err2 = np.abs(1.0 - eigen / eigeno)
        error = max(err1, err2)
        print(f"{p:3d} {eigen:12.7f} {error:12.4e}")
        eigeno = eigen
        x = xn.copy()

    return eigen, x, error
                               

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
    

def MAX_SIZE(x):
    #  ==================================================================================
    #  CODE11.1-POWER_METHOD_S.py of a module in Pseudocode 11.1.                      
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
    #   DESCRIPTION: A function module to find the largest element (in absolute value) of an array.
    #                                                                                              
    #   ARGUMENTS                                                                                  
    #      n   :: Length of the array;                                                          
    #      x   :: Array of length n.                                                            
    #                                                                                              
    #   USES                                                                                        
    #     abs  :: Built-in NumPy library function returning the absolute value of a real value;       
    #                                                                                              
    #   REVISION DATE :: 03/10/2025                                                                
    #  ==================================================================================
    xmax = x[0];
    for i in range(1, len(x)):    
        if(np.abs(x[i]) > np.abs(xmax)):    
            xmax = x[i]; 
    return xmax



def main():
    # ==============================================================================
    #  The main program to test the module POWER_METHOD_S.PY
    # ==============================================================================
    n = 4
    #a = np.zeros((n, n))
    eps = 1e-3; maxit = 999; eigen0 = 1.0
    a = np.array([[  8.0,   9.0,   10.0,   9.0],
                  [-13.0, -12.0,  -12.0, -11.0],
                  [-18.0,  -9.0,  -20.0,  -9.0],
                  [ 11.0,   1.0,   10.0,   0.0]])
    x0 = np.array([ 1.0, 0.0, 0.0, .0])      

    (eigen, x, error) = POWER_METHOD_S(n, a, x0, eigen0, eps, maxit)

    if error < eps:
        print("=" * 35)
        print(f" Largest eigenvalue = {eigen:12.7f}")
        print(f" and corresponding eigenvector\n")
        for i in range(n):
            print(f"  x({i+1:2d})= {x[i]:10.7f}")
        print("=" * 35)
    else:
        print("\n *** Max number of iterations is reached!!! ***")

if __name__ == "__main__":
    main()
