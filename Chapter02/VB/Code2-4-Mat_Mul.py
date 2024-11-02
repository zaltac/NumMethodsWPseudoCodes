import numpy as np 

# ==============================================================================
#  The main program to test mat_muld.py
# ==============================================================================
def test_mat_mul():

    A = np.array([[1.0, -3.0, 5.0, 1.0], 
                  [-2.0, 4.0, 1.0, 2.0]])

    B = np.array([[3.0, 1.0, -2.0], 
                  [2.0, 0.0, 4.0], 
                  [-1.0, 1.0, -3.0], 
                  [2.0, 5.0, 3.0]])

    print("Input Matrix A \n")
    print(A)
    print("\n Input Matrix B \n")
    print(B)

    C = mat_mul(A, B)
    print("\n ------ A(m,p)xB(p,n) matrix product is ------\n")
    print(C)

# ---------------------------------------------------------------------
#  Using numpy's MATMUL module
# ---------------------------------------------------------------------
    print("\n ------ Matrix product with np.matmul ------ \n")
    C = np.matmul(A, B)
    print(C)

def mat_mul(A,B):
    """
  ==================================================================================
  CODE2.4-MAT_MUL.py. A Python module implementing Pseudocode 2.4.                               
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTA� (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
   DESCRIPTION: A subroutine to find A*B=C matrix multiplication.                             
                                                                                              
   ON ENTRY                                                                                   
    m,p,n :: Dimension attributes of input/output matrices;                                   
       A  :: An input matrix of size m×p;                                                    
       B  :: An input matrix of size p×n.                                                    
                                                                                              
   ON EXIT                                                                                    
       C  :: The output matrix of size m×n.                                                  
   USES
     NumPy modules 
                                                                                                 
   REVISION DATE :: 03/18/2024                                                                
  ==================================================================================
    """
    global C
    m=A.shape[0] ; pA=A.shape[1]
    pB=B.shape[0];  n=B.shape[1]
    if  pA == pB:
        C = np.zeros((m,n))
        for i in range(m): 
            for j in range(m):
                for k in range(pA):
                    C[i, j] += A[i, k] * B[k, j]
    else:
        return "Sorry, cannot multiply A and B."
    
    return C  # End of mat_mul
 

if __name__ == "__main__":
    test_mat_mul()