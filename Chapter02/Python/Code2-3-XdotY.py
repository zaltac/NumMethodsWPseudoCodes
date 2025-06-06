import numpy as np

def XdotY(x, y):
#  ==================================================================================
#  CODE2.3-XdotY.py. A Python module implementing Pseudocode 2.3.                                 
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
#   DESCRIPTION: A function to compute the dot product of two vectors, x and y.                
#                                                                                              
#   ARGUMENTS                                                                                  
#      n   :: Dimension attribute of the input vectors;                                        
#     x, y :: The input vectors of length n.                                                   
# 
#   USES
#     NumPy modules 
#                                                                                                
#   REVISION DATE :: 03/18/2024                                                                
#  ==================================================================================
    n=len(x) # Assuming LEN of X and Y are the same!!!
    XdotY = 0.0
    for i in range(n):
        XdotY += x[i] * y[i]
    return XdotY  # End of XdotY

    """
    THE MODULE CAN BE IMPROVED AS FOLLOWS:
    
def XdotY(x, y):
    n=len(x) # Assuming LEN of X and Y are the same!!!
    m=len(y)
    if m != n:
        print(" Size of two vectors are unequal. ")
        print(" Execution Halted! ")
        exit()
    else:
        XdotY = 0.0
        for i in range(n):
            XdotY += x[i] * y[i]
    
    return XdotY  # End of XdotY
    """
    
# ==============================================================================
#  The main program to test XdotY.py
# ==============================================================================

def test_XdotY():
    n = 5
    x = np.array([ 1.0, 3.0, 2.0, -1.0, -1.0])
    y = np.array([-1.0, 1.0, 0.5,  1.0,  1.0])

    print("\n Vector x ");    print(x)
    print("\n Vector y ");    print(y)
    
    dotp = XdotY(x,y)
    print("\n Enorm is ",f"{dotp:10.6f}")
    
if __name__ == "__main__":
    test_XdotY()