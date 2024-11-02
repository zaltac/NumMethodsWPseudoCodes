import numpy as np 

# ==============================================================================
#  The main program to test linear_solve.py
# ==============================================================================
def test_linear_solve():
    
    n = 3
    a = np.array([[1., 4., -6.], [-1., 6., -4.], [4., -1., -1]])
    b = np.array([0., 60., 0.])

    opt = 0
    x = linear_solve(n, a, b, opt)

    print("\n Method applied here is")
    if opt == 0:
        print(" Naive Gauss Elimination")
    else:
        print(" Gauss Jordan Elimination")

    print("\n------ Matrix A & b ------------\n")
    print(" A = ",a,",  b = ",b)
    print("\n------- Solution ---------------\n")
    for i in range(n):
        print(" x(",i+1,") = ", f"{x[i]:10.6f}")


def linear_solve(n, a, b, opt):
    """
  ==================================================================================
  CODE2.9-LINEAR_SOLVE.py. A Python module implementing Pseudocode 2.9.                          
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTA« (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
  DESCRIPTION: A subroutine to solve a system of linear equations using naive                 
    Gauss Elimination (opt=0) or Gauss-Jordan Elimination (opt/=0) algorithm.                 
                                                                                              
  ON ENTRY                                                                                    
      n  :: Number of unknowns;                                                               
      A  :: Input coefficient matrix of size n√ón;                                            
      b  :: An input array of length n containing the rhs;                                    
     opt :: Option key (=0, Naive Gauss-Elim.; /=0, Gauss-Jordan Elimn).                      
                                                                                              
  ON EXIT                                                                                     
      x  :: The output array of length n containing the solution.                             
                                                                                              
  USES                                                                                        
    abs  :: NumPy function returning the absolute value of a real value;         
    Back_Substitute :: A subrotine to solve an upper-triangular system.                       
                                                                                              
  REVISION DATE :: 03/18/2024                                                                 
  ==================================================================================
    """
    eps = 1e-12

    # Forward elimination steps
    for j in range(n):
        ajj = a[j, j]
        if abs(ajj) < eps:
            print(f"Pivot is zero at j={j}")
            print("Execution is halted!")
            break # return None
        else:
            val = 1.0 / ajj
            b[j] *= val
            a[j, :] *= val
            for i in range(j + 1, n):
                s = a[i, j]
                a[i, j] = 0.0
                a[i, j + 1:] -= s * a[j, j + 1:]
                b[i] -= s * b[j]

    # Back substitution steps
    if opt == 0:
        x = back_substitute(n, a, b)
    else:
        x = np.zeros(n)
        for j in range(n - 1, 0, -1):
            for i in range(j - 1, -1, -1):
                b[i] -= a[i, j] * b[j]
                a[i, j] = 0.0
            x[j] = b[j]
        x[0] = b[0]

    return x

def back_substitute(n, a, b):
    """
  ==================================================================================
    DESCRIPTION: Solves an upper triangular system using forward substitution.

    On ENTRY
        n   :: Number of equations;
        a   :: Upper triangular matrix of size n◊n.
        b   :: An array of size n containing the rhs.

    On RETURN
        x   :: An array of size n containing the solution.

    USES                                                                                        
      zeros  :: NumPy function creating real arrays of zeros. 
  ==================================================================================        
    """
    x = np.zeros(n)
    x[-1] = b[-1] / a[-1, -1]
    for k in range(n - 2, -1, -1):
        sums = 0.0
        for j in range(k + 1, n):
            sums += a[k, j] * x[j]
        x[k] = (b[k] - sums) / a[k, k]
    return x # End of BACK_SUBSTITUTE
 


if __name__ == "__main__":
    test_linear_solve()