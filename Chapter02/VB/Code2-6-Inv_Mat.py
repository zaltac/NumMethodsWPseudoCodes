import numpy as np

# ==============================================================================
#  The main program to test inv_mat.py
# ==============================================================================
def test_inv_mat():
    
    n = 5
    a = np.array([[1, 2, 1, 2, 3], 
                    [11., -1, 1, 4, 1],
                    [4, -1, 1, 1, -1],
                    [-3, 1, -8, -1, 5],
                    [-1, 1, 1, 1, 1]])
    b = a.copy()

    print("\n ********* Input Matrix *********")
    print('\n'.join([' '.join([f'{a[i,j]:10.5f}' for j in range(n)]) for i in range(n)]))
    
    ainv = inv_mat(a)

    print("\n -------- matrix A(-1) --------")
    print('\n'.join([' '.join([f'{ainv[i,j]:10.5f}' for j in range(n)]) for i in range(n)]))
    print("-------------------------------")

    c = np.dot(b, ainv)

    print("\n -------- matrix A*A(-1) -------")
    print('\n'.join([' '.join([f'{c[i,j]:10.5f}' for j in range(n)]) for i in range(n)]))
    print("-------------------------------")

def inv_mat(a):
    """ 
  ==================================================================================
  CODE2.6-INV_MAT.py. A Python module implementing Pseudocode 2.6.                               
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTAÇ (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
  DESCRIPTION: A subroutine to find the inverse of a square matrix (with no pivoting).        
                                                                                              
  ON ENTRY                                                                                    
     n  :: Dimension attribute of input matrix A;                                             
     A  :: An input matrix (nxn).                                                             
                                                                                              
  ON EXIT                                                                                     
     B  :: Inverse of A (nxn).                                                                
                                        
  USES
     NumPy modules 
                                                         
  REVISION DATE :: 03/18/2024                                                                 
  ==================================================================================   
    """  
    n = a.shape[0]
    ainv = np.eye(n)
    for j in range(n):
        p = 1 / a[j,j] 
        for k in range(n):
            ainv[j,k] = p*ainv[j,k]
            a[j,k] = p*a[j,k]

        for i in range(n):
            s = a[i,j]
            if i != j:
                for k in range(n):
                    ainv[i,k] -= s * ainv[j,k]
                    a[i,k] -= s * a[j,k]
    return ainv # End of INV_MAT module


if __name__ == "__main__":
    test_inv_mat()