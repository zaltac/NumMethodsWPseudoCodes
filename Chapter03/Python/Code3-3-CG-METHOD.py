import numpy as np

def cgm(n, eps, a, b, xo, maxit):
    """
  ==================================================================================
  CODE3.3-CGM.py. A Python module implementing Pseudocode 3.3.                                   
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTA« (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
  DESCRIPTION: A subroutine to solve Ax=b linear system with the Conjugate Gradient Method.   
                                                                                              
  ON ENTRY                                                                                    
     n   :: Number of equations;                                                              
     A   :: Coefficient matrix (n√ón);                                                        
     b   :: Array of length n containing the right hand side;                                 
     x   :: Array of length n containing the initial guess;                                   
    eps  :: Convergence tolerance;                                                            
   maxit :: Maximum number of iterations.                                                     
                                                                                              
  ON EXIT                                                                                     
     x   :: Array of length n containing the approximate solution;                            
   iter  :: Total number of iterations performed;                                             
   error :: Euclidean (L2-) norm of displacement at exit.                                     
                                                                                              
  USES                                                                                           
    AX   :: Subroutine to evaluate A*x matrix vector product;                                 
    SQRT :: Built-in Intrinsic function returning the square root of a real value;            
    XdotY:: Function to evaluate the dot product of two vectors.                              
                                                                                              
  REVISION DATE :: 12/11/2024                                                                 
  ==================================================================================
    """
    c = np.zeros(n+1)
    r = np.zeros(n+1)
    d = np.zeros(n+1)
    x = xo.copy() #np.zeros(n+1)

    # Compute initial residual
    r =  Ax(a, x)
    for i in range(n+1):
        r[i] = b[i] - r[i]              
        d[i] = r[i]    

    rho0 = XdotY(r, r)
    print("\n *** ITERATION HISTORY ***")
    print(f"{' iter':10s} {' E-norm '}")
    for p in range(maxit+1):
        enorm = np.sqrt(rho0)
        print(f"{p:4d} {enorm:16.5E}")    
        if enorm < eps:
            break
        # Compute a search direction
        c = Ax(a, d)
        rden = XdotY(d, c)
        alpha = rho0 / rden
        for k in range(n+1):
            x[k] += alpha * d[k]
            r[k] -= alpha * c[k]

        rho = XdotY(r, r)
        beta = rho / rho0
        for k in range(n+1):
            d[k] = r[k] + beta * d[k]
        rho0 = rho

    if p == maxit:
        print(f"\n Failed to converge after {maxit} iterations")

    return p, enorm, x

def Ax(A, x):
    """
  ==================================================================================
  CODE2.5-Ax.py. A Python module implementing Pseudocode 2.5.                                    
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTA« (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
  DESCRIPTION: A subroutine to perform A * x = b matrix-vector multiplication.                
                                                                                              
  ON ENTRY                                                                                    
     n   :: Dimension attributes of input/output matrices;                                    
     A   :: An input matrix of size n√ón;                                                     
     x   :: An input vector of length n.                                                      
                                                                                              
  ON EXIT                                                                                     
     b   :: The output vector of length n.                                                    

  USES
     NumPy modules 
                                                                                                 
  REVISION DATE :: 03/18/2024                                                                 
  ==================================================================================
    """    
    n=len(x) # Assuming SIZE of A and LEN of X are the same!!!    
    b = np.zeros(n)
    for i in range(n):
        b[i] = 0.0
        for j in range(n):
            b[i] += A[i, j] * x[j]
    
    return b

def XdotY(x, y):
	"""
  ==================================================================================
  CODE2.3-XdotY.PY. A Python module implementing Pseudocode 2.3.                                 
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTA« (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.

   DESCRIPTION: A function to compute the dot product of two vectors, x and y.                
                                                                                              
   ARGUMENTS                                                                                  
      n   :: Dimension attribute of the input vectors;                                        
     x, y :: The input vectors of length n.                                                   
 
   USES
     NumPy modules 
                                                                                                
   REVISION DATE :: 03/18/2024                                                                
  ==================================================================================
	"""
    XdotY = 0.0
    n=len(x) # Assuming LEN of X and Y are the same!!!
    for i in range(n):
        XdotY += x[i] * y[i]

    return XdotY  # End of XdotY

# ==================================================================================
# Main program to test CGM.PY
# ==================================================================================
def test_cgm():
    neqs = 10
    n = neqs - 1
    a = np.zeros((n+1, n+1))
    b = np.zeros(n+1)
    x = np.zeros(n+1)
    error = 0.0
    eps = 1e-7
    maxit = 999
    
    # Input the linear system: matrix A and the rhs vector B
    a = np.array([[ 2., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                  [-1.,  2., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                  [ 0., -1.,  2., -1.,  0.,  0.,  0.,  0.,  0.,  0.],
                  [ 0.,  0., -1.,  2., -1.,  0.,  0.,  0.,  0.,  0.],
                  [ 0.,  0.,  0., -1.,  2., -1.,  0.,  0.,  0.,  0.],
                  [ 0.,  0.,  0.,  0., -1.,  2., -1.,  0.,  0.,  0.],
                  [ 0.,  0.,  0.,  0.,  0., -1.,  2., -1.,  0.,  0.],
                  [ 0.,  0.,  0.,  0.,  0.,  0., -1.,  2., -1.,  0.],
                  [ 0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  2., -1.],
                  [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  2.]])
 
    b = np.array([ 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.2])
    
    # Set initial guess
    xo = np.array([ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])


    # Call CGM module to iterate
    iter, error, x = cgm(n, eps, a, b, xo, maxit)

    # Print out the results
    print("\n  **** SOLUTION ****")

    print(f"{' ':3s} {'i':10s} {'x'}")
    for i in range(n+1):
        print(f" {i:4d} {x[i]:14.8f}")
    print(f"\n Total no of iterations = {iter:4d}")
    print(f" Maximum Error          = {error:1.3e}")
    
if __name__ == "__main__":
    test_cgm()