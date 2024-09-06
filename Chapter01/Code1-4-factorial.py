import numpy as np

# A module to test recursive factorial function
def test_factorial():

    n = int(input("Enter n : "))
    val=factorial(n)
    print(n,"! = ",val) 

def factorial(n):
#  ==================================================================================
#  CODE1.4-factorial.py. A python module for implementing Pseudocode 1.4.          
# 
#  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
#  First Edition. (c) By Zekeriya ALTAÇ (2024).
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
#  DESCRIPTION: A recursive function for computing n!                                
#                                                                                    
#  INPUT ARGUMENT                                                                    
#    n    :: Integer (n>=0)                                                          
#                                                                                    
#  ON EXIT                                                                           
#   Result:: n!                                                                      
#                                                                                    
#  REVISION DATE :: 03/21/2024                                                       
#                                                                                    
#  ==================================================================================
    if n<0:
        print("Error, Illegal Input")
        print("Argument should be n>=0, you entered",n)
        exit()
    elif  n == 1:
       fact = 1
    else:
       fact = n * factorial(n-1)
       
    return fact # End of FACTORIAL

if __name__ == "__main__":
    test_factorial()