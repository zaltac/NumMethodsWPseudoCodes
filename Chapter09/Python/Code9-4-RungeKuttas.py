import numpy as np

def fcn(x, y):
# ==============================================================================
# DESCRIPTION: A function subprogram providing y'=f(x,y)
#
# ARGUMENTS:
#      x, y  :: Real input values.
#
# ==============================================================================
    return -y/ (x + 1.0) 

def exact(x):
# ==============================================================================
# DESCRIPTION: A module providing the true solution y=(x) for 
#    testing the module. 
#
# ARGUMENTS:
#      x   :: A real input, independent variable.
#
# ==============================================================================
    return 2.0 / ( x + 1.0 )


def DRV_RK(n,h,x0,y0):
#  ==================================================================================
#  CODE9.4-DRV_RK.py. A Python module implementing Pseudocode 9.4.                                
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
#  DESCRIPTION: A driver module employing one-step RK2, RK3, or RK4 scheme.               
#                                                                                              
#  ON ENTRY                                                                                    
#   n     :: Order of Runge-Kutta scheme;                                                      
#   h     :: Step size (it must be uniform);                                                   
#   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
#                                                                                              
#  ON EXIT                                                                                     
#   x, y  :: Current estimates, x^(p+1) and y^(p+1).                                           
#                                                                                              
#  USES                                                                                             
#    FCN  :: User-defined external function providing y'=f(x,y).                               
#                                                                                              
#  REVISION DATE :: 03/05/2025                                                                 
#  ==================================================================================
    hlf= 0.50
    xh = x0 + 0.50*h
    x1 = x0 + h

    match n:

        case 2:  # Case of RK2
            xk1 = h* fcn(x0, y0)
            ym = y0 + xk1 
            xk2 = h* fcn(x1, ym)
            xk  = hlf*(xk1+ xk2)
        case 3:  # Case of RK3
            xk1 = h* fcn(x0, y0)    
            ym = y0 + hlf*xk1
            xk2 = h* fcn(xh, ym)
            ym = y0- xk1+ 2.0*xk2
            xk3 = h* fcn(x1, ym)      
            xk  = (xk1+ 4.0*xk2+ xk3 )/6.0      
        case 4:  # Case of RK4
            xk1 = h* fcn(x0, y0)
            ym = y0 + hlf*xk1    
            xk2 = h* fcn(xh, ym)
            ym = y0 + hlf*xk2 
            xk3 = h* fcn(xh, ym)     
            ym = y0 + xk3
            xk4 = h* fcn(x1, ym)
            xk  = (xk1+ 2.0*xk2+ 2.0*xk3+ xk4)/6.0 
        case _:  # Case of n<=1 or Case of n>=5
            print(' PROGRAM DOES NOT HANDLE CASE OF n= ',n)
            end
    y = y0 + xk      
    x = x1
    return (x,y)

def Runge_Kutta(order, h, x0, y0, xlast):
#  ==================================================================================
#  CODE9.4-Runge_Kutta.py. A Python module implementing Pseudocode 9.4.                           
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
#  DESCRIPTION: A module to estimate the solution of a first order IVP on [x0,xlast]       
#    using the 2nd to 4th order Runge-Kutta scheme. Numerical estimates are printed out, not stored.           
#                                                                                              
#  ON ENTRY                                                                                    
#     n    :: Order of the Runge-Kutta scheme;                                                 
#     h    :: Step size (it must be uniform);                                                  
#    x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps
#   xlast  :: End point of the solution interval.                                              
#                                                                                              
#  Other Internal Variables                                                                    
#    x, y  :: Current estimates, x^(p+1) and y^(p+1).                                          
#                                                                                              
#   USES                                                                                       
#     abs  :: Built-in NumPy library function returning the absolute value of a real value.        
#    DRV_RK:: A driver subprogram performing one-step RK scheme.                               
#                                                                                              
#  REVISION DATE :: 03/05/2025                                                                 
#  ==================================================================================
    print(f" *** {order}th order Runge-Kutta Method ***")
    print(f'{"x":>8}  {"y":>11} {"y_true":>14} {"True Error":>17}')
    print(f"{x0:12.7f} {y0:12.7f}")
    
    x = x0
    while x < xlast:
        (x ,y) = DRV_RK(order,h,x0,y0)
        yt= exact(x)
        aerr = abs(y-yt)
        
        print(f"{x:12.7f} {y:12.7f} {yt:12.7f} {aerr:14.3E}")
        x0 = x
        y0 = y



def test_Runge_Kutta():
# ==============================================================================
#  The main program to test Runge_Kutta.py
# ==============================================================================
    x0 = 0.0
    y0 = 2.0
    xlast = 1.0
    
    h = float(input("Enter h : "))
    order = int(input("Enter order of RK :"))

    Runge_Kutta(order, h, x0, y0, xlast)

test_Runge_Kutta()
