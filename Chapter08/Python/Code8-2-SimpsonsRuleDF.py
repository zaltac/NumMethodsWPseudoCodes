import numpy as np


#  ==================================================================================
#  CODE8.2-Simpsons_Rule_DF.py. A Python module implementing Pseudocode 8.2.                      
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
#  DESCRIPTION: A module to estimate the integral of a discrete function f on [a,b]                         
#   using the Simpson's 1/3 rule.                                                              
#                                                                                              
#  ON ENTRY                                                                                    
#     n   :: Number of panels (must be even!..);                                               
#     h   :: Interval size (uniform spacing, x_{i+1}-x_i);                                     
#     f   :: Array of length (n+1) containing the ordinates, f_0, f_1, ..., f_n.               
#                                                                                              
#  ON RETURN                                                                                     
#   intg  :: Numerical estimate for the integral.                                               
#                                                                                              
#  REVISION DATE :: 03/04/2025                                                                 
#  ==================================================================================
def simpsons_rule_df(n, h, f):
    m = n % 2
    if m != 0:
        print("Number of panels is not EVEN")
        exit()
    
    odd = 0.0
    for i in range(1, n, 2):
        odd += f[i]
    
    even = 0.0
    for i in range(2, n-1, 2):
        even += f[i]
    
    intg = f[0] + f[n] + 4.0 * odd + 2.0 * even
    intg = intg * h / 3.0
    return intg


! ==============================================================================
!  The main program to test module Simpsons_Rule_DF.PY
! ==============================================================================
def main():

    n = int(input("Enter Number of Panels: "))    
    a = 0.0
    b = 1.0
    h = (b - a) / float(n)
    f = np.zeros(n + 1)
    
    x = a
    for i in range(n + 1):  # Discrete dataset generation from y=x**4
        f[i] = x**4  
        x += h
    
    intg = simpsons_rule_df(n, h, f)
    
    print(n, " panel Simpson's 1/3 Rule = ", intg)
    print(" ")

if __name__ == "__main__":
    main()