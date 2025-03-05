import numpy as np

def explicit_euler(h, x0, y0, xlast):
#  ==================================================================================
#  CODE9.1-Explicit_Euler.py. A Python module implementing Pseudocode 9.1.                        
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
#  DESCRIPTION: A Python module to estimate the solution of a first order IVP on [x0,xlast]           
#    using the Explicit Euler method. Numerical estimates are printed out, not stored.         
#    With minor modifications, this PROGRAM can also be used to solve explicit methods         
#    such as MIDPOINT RULE and MODIFIED EULER.                                                 
#                                                                                              
#  ON ENTRY                                                                                    
#   h     :: Step size (it must be uniform);                                                   
#   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
#   xlast :: End point of the solution interval.                                               
#                                                                                              
#  Other Internal Variables                                                                             
#   x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
#                                                                                              
#  USES                                                                                        
#    abs  :: Built-in NP library function returning the absolute value of a real value;
#    fcn  :: User-defined external function providing y'=f(x,y).                               
#                                                                                              
#  REVISION DATE :: 03/05/2025                                                                 
#  ==================================================================================
    x , y = x0, y0
    print(f'{"x":>6}  {"y":>9} {"True Error":>14}')
    print(f"{x0:10.7f} {y0:10.7f}")
    while x < xlast:
        x = x0 + h
        y = y0 + h * fcn(x0, y0)  # Explicit-Euler
        err= abs(y - exact(x))
        print(f"{x:10.7f} {y:10.7f} {err:10.3E}")
        x0 = x
        y0 = y

def fcn(x, y):
# ==============================================================================
# DESCRIPTION: A function subprogram providing y'=f(x,y)
#
# ARGUMENTS:
#      x, y  :: Real input values.
#
# ==============================================================================
    return -y/(x + 1.0)

def exact(x):
# ==============================================================================
# DESCRIPTION: A function subprogram providing the true solution y=(x) for 
#    testing the module. 
#
# ARGUMENTS:
#      x   :: A real input, independent variable.
#
# ==============================================================================
    return 2.0 / ( x + 1.0 )


 
def test_explicit_euler():  
# ==============================================================================
#  The main program to test Explicit_Euler.PY
# ==============================================================================
    # initialize problem
    x0 = 0.0
    y0 = 2.0
    xlast = 1.0
    h = 0.1
    
    explicit_euler(h, x0, y0, xlast)

if __name__ == "__main__":
    test_explicit_euler()
