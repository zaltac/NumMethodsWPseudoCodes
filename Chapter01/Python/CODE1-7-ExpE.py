import numpy as np 

# Test module ExpE
def test_expe():

    n = int(input("Enter n : "))
    x = float(input("Enter x : "))

    val=expe(n,x)
    print("e^",x," = ",f"{val:10.8f}")  
    err=np.exp(x)-val     #  exp is a built-in NumPy function, evaluating e^x 
    print("True error is ",f"{err:.4E}")
    
 
def expe(n,x):
#  ==================================================================================
#  CODE1.7-ExpE.py. A Python module implementing Pseudocode 1.7.                                  
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
#  DESCRIPTION: A function to compute e^x using the MacLaurin series with specified            
#     number of terms.                                                                         
#                                                                                              
#  ARGUMENTS                                                                                   
#     x   :: A real input (exponent) value;                                                    
#     n   :: The number of terms of the MacLauring series to be included.                      
#                                                                                              
#  USES                                                                                        
#    float:: A built-in NumPy function that converts an integer argument to a real value.  
#                                                                                              
#  REVISION DATE :: 03/18/2024                                                                                                                     
#  ==================================================================================
    term= 1.0
    expe= 1.0 
    for k in range(1,n):
        term *= x/float(k)
        expe += term
    
    return expe

if __name__ == "__main__":
    test_expe()
