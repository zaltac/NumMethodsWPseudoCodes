import numpy as np

def tridiagonal(s1, sn, b, d, a, c):
    """
  ==================================================================================
  CODE2.13-TRIDIAGONAL.py. A Python module implementing Pseudocode 2.13.                         
 
  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
  First Edition. (c) By Zekeriya ALTAÇ (2024).
  ISBN: 978-1-032-75474-1 (hbk)
  ISBN: 978-1-032-75642-4 (pbk)
  ISBN: 978-1-003-47494-4 (ebk)
  
  DOI : 10.1201/9781003474944
  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
  
  This free software is complimented by the author to accompany the textbook.
  E-mail: altacz@gmail.com.
  
  DESCRIPTION: A subroutine to solve a tridiagonal system of linear equations                 
    using Thomas algorithm.                                                                   
                                                                                              
  ON ENTRY                                                                                    
     s1 :: Subscript of the first unknown (usually 1);                                        
     sn :: Subscript of the last unknown (usually no. of eqs, n);                             
      b :: Array of length n containing coefficients of below diagonal elements;              
      d :: Array of length n containing coefficients of diagonal elements;                    
      a :: Array of length n containing coefficients of above diagonal elements;              
      c :: Array of length n containing coefficients of rhs.                                  
                                                                                              
  ON RETURN                                                                                     
      x :: An array of length n containing the solution.                                      
     
  USES
     NumPy modules 
                                                                                              
  REVISION DATE :: 03/18/2024                                                                 
  ==================================================================================
    """
    x = np.zeros(sn-s1+1)
    for i in range(s1+1, sn+1):
        ratio = b[i] / d[i-1]
        d[i] = d[i] - ratio * a[i-1]
        c[i] = c[i] - ratio * c[i-1]

    x[sn] = c[sn] / d[sn]
    for i in range(sn-1, s1-1, -1):
        x[i] = (c[i] - a[i] * x[i+1]) / d[i]
    
    return x # End of module TRIDIAGONAL


# ==============================================================================
#  The main program to test tridiagonal.py
# ==============================================================================
def test_tridiagonal():
    # n = 9
    b = np.array([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    d = np.array([-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0])
    a = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    c = np.array([-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -14.0])

    s1 = 0; sn = 8

    x=tridiagonal(s1, sn, b, d, a, c)
 
    print(" *** Solution ***")
    for i in range(s1,sn+1):
        print(f" x({i}) = {x[i]:10.8f}")
        

if __name__ == "__main__":
    test_tridiagonal()
