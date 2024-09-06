import numpy as np

def test_quadratic_eq():

    a = float(input("Insert coefficient a: "))
    b = float(input("Insert coefficient b: "))
    c = float(input("Insert coefficient c: "))
	
    (re,im)=quadratic_eq(b/a, c/a)
    
    print("**** The roots are **** ")
    print("x1 = ",f"{re[0]:10.6f}"," + (",f"{im[0]:10.6f}",") i")   
    print("x2 = ",f"{re[1]:10.6f}"," + (",f"{im[1]:10.6f}",") i")       
    
    
    
def quadratic_eq(p, q):
#  ==================================================================================
#  CODE1.3-Quadratic_Eq.py. A python module for implementing Pseudocode 1.3 in          
# 
#  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
#  First Edition. (c) By Zekeriya ALTAC (2024).
#  ISBN: 978-1-032-75474-1 (hbk)
#  ISBN: 978-1-032-75642-4 (pbk)
#  ISBN: 978-1-003-47494-4 (ebk)
#  
#  DOI : 10.1201/9781003474944
#  C&H/CRC PRESS, Boca Raton & London.
#  
#  This free software is complimented by the author to accompany the textbook.
#  E-mail: altacz@gmail.com.
#                                                                                    
#  DESCRIPTION: A python module to find the roots of a quadratic equation of            
#    the form : x*x + p * x + q = 0.                                                 
#                                                                                    
#  ON ENTRY                                                                          
#    p, q :: Coefficients of the quadratic equation;                                 
#                                                                                    
#  ON EXIT                                                                           
#    re   :: Array of length 2 containing real parts of the roots: re1, re2;         
#    im   :: Array of length 2 containing imaginary parts of the roots: im1, im2.    
#                                                                                    
#  USES                                                                              
#    SQRT :: Built-in Intrinsic function to evaluate the square root of a real value.
#                                                                                    
#  REVISION DATE :: 03/18/2024                                                       
#                                                                                    
#  ==================================================================================

    d = p * p - 4.0 * q
    re = np.zeros(2)
    im = np.zeros(2)
    if d<0:
        d = np.sqrt(-d)
        re[0]=-p/2.0 ; re[1]= re[0] 
        im[0]=-d/2.0 ; im[1]=-im[0]
    else:
        d = np.sqrt(d)
        re[0]=(-p-d)/2.0  
        re[1]=(-p+d)/2.0  

    return (re, im)

if __name__ == "__main__":
    test_quadratic_eq()

