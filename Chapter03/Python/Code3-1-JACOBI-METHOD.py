import numpy as np

# Program to test Jacobi method
def test_jacobi_method():
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

def jacobi(n, eps, a, b, xo, maxit):
    """
    DESCRIPTION: A module to solve Ax=b iteratively using the Jacobi method.

    On ENTRY
        n   :: Number of equations (size of matrix A); 
        a   :: Coefficient matrix (n×n);  
        b   :: Array of lenght n containing the RHS;
        xo  :: Array of lenght n containing the initial guess;   
        eps :: Convergence tolerance;
        maxit:: Maximum number of iterations allowed.
        
    On EXIT:
        x   :: Array of lenght n containing the final estimate for the solution; 
        iter:: Total number of iterations performed;
        error:: Euclidean norm of the error.
    
    USES
        jacobi_drv:: Driver module for a sor iteration.
    """
    del0 = 1.0
    p = 0
    x = xo.copy()
    while (del0>eps and p<maxit):
        p += 1
        (del1,x)=jacobi_drv(n, a, b, xo)
        # printout the ietartion progress 
        print(f" iter={p:5}  delta = {del1:8.2e}")
        xo[:] = x
        del0 = del1
    
    if p == maxit:
        print(f"\nJacobi method failed to converge after {maxit} iterations\nwithin the specified EPS tolerance.")
    error=del1
    iter=p
    return (iter, error, x) # End of module jacobi

    
def jacobi_drv(n, a, b, xo):
    """
    DESCRIPTION: A module to perform one step jacobi iteration.

    On ENTRY
        n   :: Number of equations (size of matrix A); 
        a   :: Coefficient matrix (n×n);  
        b   :: Array of lenght n containing the RHS;
        xo  :: Array of lenght n containing the initial guess;   
        eps :: Convergence tolerance;
        maxit:: Maximum number of iterations allowed.
        
    On EXIT:
        x   :: Array of lenght n containing the final estimate for the solution; 
        iter:: Total number of iterations performed;
        error:: Euclidean norm of the error
    
    USES
        Enorm  :: Euclidean norm of a vector,
    """    
    d = np.zeros(n)  
    x = np.zeros(n)  
    for i in range(n):
        sums = 0.0
        for j in range(n):
            if(i!=j):
                sums += a[i, j] * xo[j]
    
        x[i] =  (b[i] - sums) / a[i, i]
        d[i] = x[i] - xo[i]
    
    delta = ENorm(d)
    return  (delta, x) # End of module JACOBI_DRV

def ENorm(x):
    """
    DESCRIPTION: A funcion module to compute Euclidean (L2) norm of a vector.

    On ENTRY
        n  :: Length of vector; 
        x  :: An array (vector)) of length n.

    On EXIT
      ENorm:: Euclidean Norm of the vector
    """
    n=len(x)
    delta = 0.0
    for i in range(n):
        delta += x[i] * x[i]
    
    return np.sqrt(delta)  # End of E-norm module
