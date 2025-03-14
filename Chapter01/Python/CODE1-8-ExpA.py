import numpy as np 


# ==============================================================================
#  The main program to test function EXPA
# ==============================================================================
def test_expa():

    x   = float(input("Enter  x : "))
    eps = float(input("Enter Tol: "))

    val=expa(x,eps)
    print("e^",x," = ",f"{val:10.8f}")  
    err=np.exp(x)-val
    print("True error is ",f"{err:.4E}")
    
 
def expa(x,eps):
#  ==================================================================================
#  CODE1.8-ExpA.py. A Python module implementing Pseudocode 1.8.                                  
# 
#  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
#  First Edition. (c) By Zekeriya ALTAÇ (2024).
#  ISBN: 978-1-032-75474-1 (hbk)
#  ISBN: 978-1-032-75642-4 (pbk)
#  ISBN: 978-1-003-47494-4 (ebk)
#  
#  DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
#  
#  This free software is complimented by the author to accompany the textbook.
#  E-mail: altacz@gmail.com.
#                                                                                              
#  DESCRIPTION: A function to compute e^x adaptively using the MacLaurin series                
#     within a user-defined tolerance.                                                         
#                                                                                              
#  ARGUMENTS                                                                                   
#     x   :: A real input (exponent) value;                                                    
#    eps  :: A user-defined convergence tolerance.                                             
#                                                                                              
#  USES                                                                                        
#   float :: Built-in NumPy library function converting an integer argument to a real value.  
#    abs  :: Built-in NumPy library function returning the absolute value of a real value.         
#                                                                                              
#  REVISION DATE :: 04/11/2024                                                                                                                                           
#  ==================================================================================
    term= 1.0
    expa= 1.0 
    k=0
    
    while (np.abs(term)>eps*np.abs(expa)):
        k += 1
        expa *= x/float(k)
        expa += term
    
    return (k, expa)

if __name__ == "__main__":
    test_expa()
