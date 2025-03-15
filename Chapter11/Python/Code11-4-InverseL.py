import numpy as np

def InverseL(n, L):
    #  ==================================================================================
    #  CODE11.4-InverseL.py. A Python module implementing Pseudocode 11.4.                            
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
    #   DESCRIPTION: A Python module to invert a lower-triangular matrix.                                 
                                                                                              
    #   ON ENTRY                                                                                   
    #     n    :: Dimension attribute of the matrix L;                                             
    #     L    :: A lower-triangular matrix, nxn.                                                  
    #                                                                                              
    #   ON EXIT                                                                                    
    #     eL   :: Inverse of L, also a lower-triangular matrix, nxn.                               
    #     
    #   USES
    #     zeros, range modules of NumPy      
    #                                                                                    
    #   REVISION DATE :: 03/15/2025                                                                
    #  ==================================================================================
    eL = np.zeros((n, n))
    
    eL[0][0] = 1.0 / L[0][0]
    for i in range(1, n):
        eL[i][i] = 1.0 / L[i][i]
        for j in range(i - 1, -1, -1):
            sum = 0.0
            for k in range(j + 1, i + 1):
                sum += eL[i][k] * L[k][j]
            eL[i][j] = -sum / L[j][j]
    
    return eL


def main():
    # ==============================================================================
    #   The main program to test module InverseL.py
    # ==============================================================================
    n = 5
    L = np.array([[1.0, 0.0, 0.0, 0.0, 0.0], 
                  [4.0, 3.0, 0.0, 0.0, 0.0],       
                  [5.0, 2.0, 4.0, 0.0, 0.0],      
                  [3.0, 1.0, 8.0, 2.0, 0.0],      
                  [-1.0, 1.0, 1.0, 1.0, 1.0]])
    
    print(" Input Lower Matrix, L ")
    for row in L:
        print(" ".join(f"{value:.7f}" for value in row))

    eL = InverseL(n, L)

    print(" ------ Output : Inverse L ------ ")
    for row in eL:
        print(" ".join(f"{value:.7f}" for value in row))

if __name__ == "__main__":
    main()