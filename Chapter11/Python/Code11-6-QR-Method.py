import numpy as np

def Basic_QR(n, d, e, eps, maxit):
#  ==================================================================================
#  CODE11.6-Basic_QR.py. A Python module implementing Pseudocode 11.6.                            
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
#   DESCRIPTION: A Python module implementing the QR Factorization algorithm to a symmetric       
#      tridiagonal matrix to find its eigenvalues and eigenvectors.                            
#                                                                                              
#   ON ENTRY                                                                                   
#      n   :: Dimension attribute of the tridiagonal matrix (nxn);                             
#      d   :: An array of length n containing the main diagonal, d(1) ... d(n);                
#      e   :: An array of length n containing the subdiagonal, e(1) ... e(n-1).                
#                                                                                              
#   ON EXIT                                                                                    
#      d   :: An array of length n containing the eigenvalues;                                 
#      V   :: A square matrix containing the eigenvector(nxn).                                 
#                                                                                              
#   USES                                                                                       
#     sqrt :: Built-in NumPy library function returning the square root of a real value;       
#     zeros, eye, range are also NumPy library modules.
#                                                                                              
#   REVISION DATE :: 03/15/2025                                                                
#  ==================================================================================
    c = np.zeros(n)
    s = np.zeros(n)
    V = np.eye(n)
    
    err = 1.0
    p = 0
    print(" iter     Error")
    while (err > eps) and (p < maxit):
        t = e[0]
        for k in range(n-1):
            rho = np.sqrt(d[k]**2 + t**2)
            c[k] = d[k] / rho
            s[k] = t / rho
            d[k] = rho
            t = e[k]
            e[k] = t*c[k] + d[k+1]*s[k]
            d[k+1] = -t*s[k] + d[k+1]*c[k]
            if k < n-2:
                t = e[k+1]
                e[k+1] = t*c[k]
            for i in range(n):
                q = V[i, k]
                r = V[i, k+1]
                V[i, k] = c[k]*q + s[k]*r
                V[i, k+1] = -s[k]*q + c[k]*r
        
        for k in range(n-1):
            d[k] = d[k]*c[k] + e[k]*s[k]
            t = d[k+1]
            e[k] = t*s[k]
            d[k+1] = t*c[k]
        
        err = 0.0
        for i in range(n):
            err += e[i]**2
        err = np.sqrt(err)
        p += 1
        print(f"{p:3d}    {err:12.4E}")
    
    return d, V


def Test_QRIteration():
# ==============================================================================
#  The main program to test the module Basic_QR.
# ==============================================================================
    n = 4
    d = np.array([ 3.,  8., 6., 9.])
    e = np.array([ 4.,  2., 1., 0.])
    
    maxit = 500
    eps = 1.0e-3
    
    d, V = Basic_QR(n, d, e, eps, maxit)

    print("\n *** Eigenvalues ***")
    for i in range(n):
        print(" x(", i+1,")=",f"{d[i]:10.7f}")
        #print(" x(", i+1," )=",f"{d[i]:10.7}")

    print("\n *** Eigenvectors *** ")
    for i in range(n):
        print(" ".join(f"{v:10.7f}" for v in V[i]))

if __name__ == "__main__":
    Test_QRIteration()