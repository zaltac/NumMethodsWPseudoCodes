import numpy as np

#  ==============================================================================
#   The main program to test Module richardson
#  ==============================================================================
def test_richardson():
    x0 = 150.0
    h = 5.0
    eps = 1.0e-6
    
    (n, deriv, D) = richardson(x0, h, eps)

    print("  *** Richardson Table *** ")
    for k in range(n + 1):
        print(f" Level={k} {' '.join(f'{d:.11f}' for d in D[k, :k + 1])}")

    error =np.abs(D[n,n]-D[n-1,n-1])
    print("\n ------------------------------------")
    print(" Derivative f'(",x0,")=",f"{deriv:11.8f}")
    print(" Estimated error =",f"{error:8.3E}")
    print(" ------------------------------------\n")    



#  ==================================================================================
#  CODE5.2-RICHARDSON.py. A Python module implementing Pseudocode 5.2.                            
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
#  DESCRIPTION: A module to compute the first derivative of an explicitly                      
#      defined function using Richardson's extrapolation.                                    
#                                                                                              
#  ON ENTRY                                                                                    
#      x0  :: Point at which derivative is to be computed;                                     
#      h   :: Initial interval size;                                                           
#      eps :: Tolerance desired.                                                               
#                                                                                              
#  ON EXIT                                                                                     
#      D   :: A matrix containing the Richardson's table (0..n, 0..n)                        
#      nr  :: Size of the table;                                                               
#    deriv :: Estimated derivative.                                                            
#                                                                                              
#  USES                                                                                        
#     abs  :: Built-in NumPy library function returning the absolute value of a real value.                                                                              
#     func :: A user-defined external function providing the nonlinear equation.                
#                                                                                              
#  REVISION DATE :: 06/13/2024                                                                 
#  ==================================================================================  
def richardson(x0, h, eps):
    D = np.zeros((11, 11))
    m = 0
    k = 0
    err = 1.0

    while err > eps:
        D[k, 0] = (func(x0 + h) - func(x0 - h)) / (2.0 * h)  # 1st derivative
        for m in range(1, k + 1):
            D[k, m] = (4 ** m * D[k, m - 1] - D[k - 1, m - 1]) / (4 ** m - 1)
        if k >= 1:  # Estimate diagonalwise differentiation error
            err = abs(D[k, k] - D[k - 1, k - 1])
        h = h / 2
        k = k + 1

    nr = k - 1
    deriv = D[nr, nr]
    return (nr, deriv, D)

#  ==================================================================================
#  USER-DEFINED FUNCTION "FUNC" OF ONE-VARIABLE
# ==================================================================================
def func(x):
    return 25000.0 / (-57.0 + x) - 5.2e6 / x ** 2



if __name__ == "__main__":
    test_richardson()
