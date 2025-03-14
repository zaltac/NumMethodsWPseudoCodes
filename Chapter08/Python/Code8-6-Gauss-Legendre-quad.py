import numpy as np
import math

def Gauss_Legendre_Quad(n, eps):
#  ==================================================================================
#  CODE8.6-Gauss_Legendre_Quad.py. A Python module implementing Pseudocode 8.6.                   
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
#  DESCRIPTION: A subroutine to generate N-point Gauss-Legendre quadrature                     
#    abscissas and weights on [-1,1].                                                          
#                                                                                              
#  ON ENTRY                                                                                    
#     n   :: Number of quadrature points;                                                      
#     eps :: Tolerance, i.e., desired level of numerical accuracy.                             
#                                                                                              
#  ON EXIT                                                                                     
#     x   :: Array of length N containing the abscissas;                                       
#     w   :: Array of length N containing the weights.                                         
#                                                                                              
#  USES                                                                                        
#    abs  :: Built-in NumPy library function returning the absolute value of a real value;         
#    cos  :: Built-in NumPy library function returning trigonometric cosine value.  
#                                                                                              
#  REVISION DATE :: 03/04/2025                                                                 
#  ==================================================================================
    x = np.zeros(n)
    w = np.zeros(n)
    m = int( (n+1) / 2)

    pi = 3.1415926535897932385
    for i in range(m):  # formula is modified since python subscripts start from zero.
        u = np.cos(pi * (4. * i + 3.) / (4. * n + 2.))
        delta = 1.0
        while abs(delta) > eps:
            P0 = 1.0
            P1 = u
            for k in range(2, n+1):
                P2 = ((2. * k - 1.) * u * P1 - (k - 1.) * P0 )/ (1.*k)
                P0 = P1
                P1 = P2
            
            PP = n * (u * P1 - P0) / (u * u - 1.0)
            delta = P1 / PP
            u -= delta
        x[i] = -u
        w[i] = 2.0 / (1.0 - u * u) / PP ** 2
        x[n - i-1] = u
        w[n - i-1] = w[i]
    return x, w # End of Gauss_Legendre_Quad


# ==============================================================================
#  The main program to test module Gauss_Legendre_Quad.py
# ==============================================================================
def Test_Gauss_Legendre_Quad():
    n = int(input("Enter n: "))

    eps = 1.0e-6
    (x, w) = Gauss_Legendre_Quad(n, eps)
    print(" i         x_i             w_i")
    print("-- ------------------- -------------------")
    for i in range(n):
        print(f"{i+1:2d} {x[i]:18.12f} {w[i]:18.12f}")

Test_Gauss_Legendre_Quad()  # End of Test_Gauss_Legendre_Quad
