import numpy as np

def lagrange_p(n, xval, x):
#  ==================================================================================
#  CODE6.1-LAGRANGE-P.py. A Python module implementing Pseudocode 6.1.                            
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
#  DESCRIPTION: A Python module to compute the Lagrange polynomials at x=xval.                        
#                                                                                              
#  ON ENTRY                                                                                    
#    n    :: The number of data in the set minus 1;                                            
#    xval :: x value at which the dataset is to be interpolated;                               
#    x    :: Array of length n+1 containing abscissas, k=0,1,2,...,n.                          
#                                                                                              
#  ON EXIT                                                                                     
#    L    :: Array of length (n+1) containing Lagrange polynomials                             
#            that is, L(k) = L_k(xval) for k=0,1,2,...,n.                                      
#                                                                                              
#  REVISION DATE :: 02/28/2025                                                                 
#  ==================================================================================
    L = np.ones(n+1)
    for k in range(n+1):
        for i in range(n+1):
            if i != k:
                L[k] *= (xval - x[i]) / (x[k] - x[i])
    return L   # End of Lagrange_P 


def lagrange_eval(n, xval, x, f):
#  ==================================================================================
#  CODE6.1-LAGRANGE_EVAL.py. A Python module implementing Pseudocode 6.1.                            
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
#  DESCRIPTION: A function to evaluate an interpolation value f=f(xval)=fval using             
#    the Lagrange interpolating polynomials.                                                   
#                                                                                              
#  ARGUMENTS                                                                                   
#    n    :: The number of data in the set minus 1;                                            
#    xval :: x value at which the dataset is to be interpolated;                               
#    x    :: Array of length (n+1) containing abscissas, k=0,1,2,...,n;                        
#    f    :: Array of length (n+1) containing ordinates, k=0,1,2,...,n.                        
#                                                                                              
#  USES                                                                                        
#    LAGRANGE_P :: A subroutine generating the Lagrange polynomials.                                          
#                                                                                              
#  REVISION DATE :: 02/28/2025                                                                 
#  ==================================================================================
    L = lagrange_p(n, xval, x)
    fval = 0.0
    for k in range(n+1):
        fval += L[k] * f[k]
    return fval  # End of Lagrange_Eval 


# ==============================================================================
#  The main program to test subprograms lagrange_eval.PY 
# ==============================================================================

def test_lagrange_interp():

    n = 3
    x = np.array([0.08, 0.25, 0.5, 0.9])
    f = np.array([0.25, 0.625, 0.81, 0.43])

    xval = float(input(" Enter xval : "))
    fval = lagrange_eval(n, xval, x, f)

    print(" f(",xval,")=",f"{fval:10.8f}") 

test_lagrange_interp()  # End of Test_Lagrange