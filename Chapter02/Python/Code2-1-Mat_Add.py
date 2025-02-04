import numpy as np

def mat_add(m, n, A, B):
	"""
  ==================================================================================
  CODE2.1-MAT_ADD.py. A Python module implementing Pseudocode 2.1.                               
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTAÇ (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.

   DESCRIPTION: A module to perform C=A+B matrix addition.                                
                                                                                              
   ON ENTRY                                                                                   
     m,n :: Dimension attributes of the matrices;                                             
      A  :: An input matrix (mxn);                                                            
      B  :: An input matrix (mxn).                                                            
                                                                                              
   ON RETURN                                                                                  
      C :: The output matrix (mxn).                                                           
                                                                                              
   USES
     NumPy modules 
   
   REVISION DATE :: 03/18/2024                                                                
  ==================================================================================
   The matrix addition operation for A and B matrices of the same size 
   can be carried out as follows:
       c = a + b
   provided that a, b, and c are defined as arrays of the same size
	"""
    C = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            C[i, j] = A[i, j] + B[i, j]
    return C


# ==============================================================================
#  The main program to test mat_add.py
# ==============================================================================
n = 3
m = 3

A = np.reshape(np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]), (m, n))
B = np.reshape(np.array([3.0, -2.0, 1.0, -2.0, 2.0, 4.0, 3.0, -5.0, 1.0]), (m, n))

print("\n     Input Matrix A")
for i in range(m):
    print(" ".join(f"{A[i, j]:10.5f}" for j in range(n)))

print("\n     Input Matrix B")
for i in range(m):
    print(" ".join(f"{B[i, j]:10.5f}" for j in range(n)))

C = mat_add(m, n, A, B)

print("\n     Output Matrix C")
for i in range(m):
    print(" ".join(f"{C[i, j]:10.5f}" for j in range(n)))