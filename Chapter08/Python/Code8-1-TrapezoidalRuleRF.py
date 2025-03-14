import numpy as np

# ==============================================================================
#  The main program to test pyhton modeul Trapezoidal_Rule_RF.py
# ==============================================================================
def test_Trapezoidal_Rule_RF():
 
    a = float(input(" Enter a : "))
    b = float(input(" Enter b : "))
    n = int(input(" Enter no. of panels : "))

    (int1,intc)=Trapezoidal_Rule_RF(FX,FU,a,b,n)

    print("\n Trapezoidal Rule    = ",f"{int1:11.8f}")
    print(" With end-correction = ",f"{intc:11.8f}")   
    
 
def Trapezoidal_Rule_RF(FX, FU, a, b, n):
#  ==================================================================================
#  CODE8.1-Trapezoidal_Rule_RF.py. A Python module implementing Pseudocode 8.1.                   
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
#  DESCRIPTION: A subroutine to estimate the integral of y=f(x) on [a,b]                       
#   using the Trapezoidal rule with/without end correction.                                    
#                                                                                              
#  ON ENTRY                                                                                    
#     n   :: Number of panels (i.e., n+1 integration points);                                  
#   [a,b] :: Integration interval.                                                             
#                                                                                              
#  ON RETURN                                                                                     
#   intg   :: Integral estimate using the ordinary Trapezoidal rule;                           
#   intgc  :: Integral estimate using the Trapezoidal rule with the end-point correction.      
#                                                                                              
#  USES                                                                              
#     FX   :: User-defined external function providing the function, f(x);                     
#     FU   :: User-defined external function providing the first derivative, f'(x);                                      
#    FLOAT :: A built-in NumPy library function that converts an integer argument to a real value. 
#                                                                                              
#  REVISION DATE :: 03/03/2025                                                                 
#  ==================================================================================
    h = (b-a) / float(n)
    intg = 0.5 * (FX(a) + FX(b))
    xi=a
    for i in range(1, int(n)):
        xi = xi + h
        intg = intg + FX(xi)
  
    intg = h*intg 
    corr = -h*h*(FU(b)-FU(a))/12   
    intgc = intg + corr
    return (intg,intgc) # End of Trapezoidal_Rule_RF 


def FX(x):
 # ==============================================================================
 # DESCRIPTION: User-defined function providing y=f(x) to be integrated. 
 #
 # ARGUMENTS:
 #      x   :: a real input value.
 # ==============================================================================
    return x**4
    
def FU(x):
 # ==============================================================================
 # DESCRIPTION: User-defined function providing first derivative, f'(x), explicitly. 
 #
 # ARGUMENTS:
 #      x   :: a real input value.
 # ==============================================================================
    return 4 * x ** 3   


if __name__ == "__main__":
    test_Trapezoidal_Rule_RF()







