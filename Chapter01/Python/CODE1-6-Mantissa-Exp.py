import numpy as np

# A module to test mantissa_exp module
def text_mantissa_exp():

    x = float(input("Enter x : "))

    (mf,ef) = mantissa_exp(x)
    print(" Mantissa is ",f"{mf:8.6f}")  
    print(" Exponent is ",ef)

def mantissa_exp(fl):
#  ==================================================================================
#  CODE1.6-MantissaExp.py. A python module for implementing Pseudocode 1.6.           
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
#  DESCRIPTION: A module to determine a floating-point number's                      
#     mantissa and exponent, i.e., fl=Mx10^e                                         
#                                                                                    
#  ON ENTRY                                                                          
#     fl  :: A floating point number.                                                
#                                                                                    
#  ON EXIT                                                                           
#     M   :: Mantissa;                                                               
#     e   :: Exponent.                                                               
#                                                                                    
#  USES                                                                              
#    ans  :: Built-in NumPy function returning the absolute value of a real value
#    floor:: Built-in NumPy function returning the greatest integer less than    
#      or equal to a real value;                                                     
#    log10:: Built-in NumPy function returning the base 10 logarithm of a real value.
#                                                                                    
#  REVISION DATE :: 03/22/2024                                                                                                                                     
#  ==================================================================================
    if np.abs(fl)>0.0: 
        e = int(np.floor(np.log10(np.abs(fl))) + 1)  
        m = fl * 10.0**(-e)
    else:
        e = 0
        m = 0.0
        
    return (m, e)  # End of FACTORIAL

if __name__ == "__main__":
    text_mantissa_exp()