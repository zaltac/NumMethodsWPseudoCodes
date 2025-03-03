import math
import numpy as np

def Linearize_Regress(ndata, x, y, model):
#  ==================================================================================
#  CODE7.2-LINEARIZE_REGRESS.py. A Python module implementing Pseudocode 7.2.                     
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
#  DESCRIPTION: A Python module to obtain least-squares-fit to the "power model" only. Users
#    can likewise incorporate the other linearizable models into the module.          
#                                                                                              
#  ON ENTRY                                                                                    
#     n   :: The number of data in the set;                                                    
#    x,y  :: Arrays of length n containing the data;                                           
#   model :: Model flag, Model = 1 corresponds to the power model, Y= a0*x^b0.                 
#                                                                                              
#  ON EXIT                                                                                     
#    a0,b0:: Model parameters;                                                                 
#     E   :: Sum of the Squares of Residuals (SSR);                                            
#     S   :: Sum of the Squares of Mean Deviation (SSMD);                                      
#    r2   :: r-squared, coefficient of determination.                                          
#                                                                                              
#  USES                                                                                        
#     Exp :: Built-in Math Library function returning exponential of a real value, e^x.    
#     Pow :: Built-in Math Library function raising a number to the power of another number.          
#     Log :: Built-in Math Library function returning the natural log of a real value.            
#                                                                                              
#  REVISION DATE :: 03/03/2025                                                                 
#  ==================================================================================
    xx = np.zeros(ndata)
    yy = np.zeros(ndata)
    b = np.zeros(2)
    c = np.zeros((2, 2))
    DD, D1, D2, yavg, yi, xk = 0, 0, 0, 0, 0, 0

    for k in range(ndata):
        if model == 1:
            xx[k] = math.log(x[k])
            yy[k] = math.log(y[k])
        else:
            print("Undefined model...")
            return None, None, None, None, None

    for i in range(2):
        for j in range(2):
            p = i + j
            for k in range(ndata):
                xk = math.pow(xx[k], p) if p != 0 else 1.0
                c[i][j] += xk
        for k in range(ndata):
            p = i
            xk = math.pow(xx[k], p) if p != 0 else 1.0
            b[i] += xk * yy[k]

    DD = c[0][0] * c[1][1] - c[1][0] * c[0][1]
    D1 = b[0] * c[1][1] - b[1] * c[0][1]
    D2 = b[1] * c[0][0] - b[0] * c[1][0]

    a0 = D1 / DD
    b0 = D2 / DD

    yavg = np.mean(y)

    S = 0.0
    E = 0.0

    if model <= 2:
        a0 = math.exp(a0)
    else:
        print("Undefined model...")
        exit()

    for k in range(ndata):
        S += math.pow(y[k] - yavg, 2)
        yi = a0 * math.pow(x[k], b0)
        E += math.pow(yi - y[k], 2)

    r2 = 1.0 - E / S

    return a0, b0, E, S, r2

# ==============================================================================
#  The main program to test Linearize_Regress.py
# ==============================================================================
def main():
    ndata = 6
    x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    y = [5.3, 6.5, 6.8, 7.1, 7.8, 7.5]
    model = 1

    a0, b0, E, S, r2 = Linearize_Regress(ndata, x, y, model)

    print("******* Best-Fit Coefficients and Parameters ********")
    print(f"  a = {a0:.6f}   b = {b0}")
    print(f"  E = {E}   S = {S}   r-squared = {r2}")

if __name__ == "__main__":
    main()