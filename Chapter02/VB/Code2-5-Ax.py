import numpy as np

def Ax(A, x):
    """
    DESCRIPTION: Performs A*x = b matrix multiplication.

    On ENTRY
        n   :: Dimension attributes of input/output matrices;
        A   :: Input matrix of size n×n;
        x   :: Input vector of length n;
    
    On EXIT
        b   :: Product vector of length n.
    """
    n=len(x) # Assuming SIZE of A and LEN of X are the same!!!  
    b = np.zeros(n)
    for i in range(n):
        b[i] = 0.0
        for j in range(n):
            b[i] += A[i, j] * x[j]
    
    return b

def ax(a, x):
    """
  ==================================================================================
  CODE2.5-Ax.py. A Python module implementing Pseudocode 2.5.                                    
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTAÇ (2024).
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
     A   :: An input matrix of size nÃ—n;                                                     
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
        b[i] = np.dot(a[i], x)
    return b

# ==================================================================================    
# Test program for module Ax.py
# ==================================================================================  
def test_AX():
    n = 4
    A = np.array([[ 1.0, 4.0, 3.0,-2.0], 
                [ 2.0, 1.0, 2.0, 3.0],
                [ 3.0, 2.0, 1.0, 4.0],
                [-2.0, 3.0, 2.0, 1.0]])
    x = np.array([2.0, 5.0, 3.0, -3.0])
    #b = np.zeros(n)
    
    print(" Input Matrix A")
    print(A)
    print("\n Input Vector x")
    print(x)

    b = Ax(A, x)

    print("\n ------ A(n,n)X(n) product is ")
    print("\n Output Vector b")
    print(b) 

    b = ax(A, x)

    print("\n ------ A(n,n)X(n) product is ")
    print("\n Output Vector b")
    print(b) 

if __name__ == "__main__":
    test_AX()
