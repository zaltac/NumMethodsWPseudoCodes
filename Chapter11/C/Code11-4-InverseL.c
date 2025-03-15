#include <stdio.h>

void InverseL(int n, double L[5][5], double eL[5][5]) {
    // ==================================================================================
    // CODE11.4-InverseL.C. A C module implementing Pseudocode 11.4.                                  
    // 
    // NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
    // First Edition. (c) By Zekeriya ALTAÇ (2024).
    // ISBN: 978-1-032-75474-1 (hbk)
    // ISBN: 978-1-032-75642-4 (pbk)
    // ISBN: 978-1-003-47494-4 (ebk)
    // 
    // DOI : 10.1201/9781003474944
    // C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
    // 
    // This free software is complimented by the author to accompany the textbook.
    // E-mail: altacz@gmail.com.
    // 
    //  DESCRIPTION: A C module to invert a lower-triangular matrix.                                 
    //                                                                                             
    //  ON ENTRY                                                                                   
    //    n    :: Dimension attribute of the matrix L;                                             
    //    L    :: A lower-triangular matrix, nxn.                                                  
    //                                                                                             
    //  ON EXIT                                                                                    
    //    eL   :: Inverse of L, also a lower-triangular matrix, nxn.                               
    //                                                                                             
    //  REVISION DATE :: 03/15/2025                                                                
    // ==================================================================================
    int i, j, k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            eL[i][j] = 0.0;
        }
    }

    eL[0][0] = 1.0 / L[0][0];
    for (i = 1; i < n; i++) {
        eL[i][i] = 1.0 / L[i][i];
        for (j = i - 1; j >= 0; j--) {
            sum = 0.0;
            for (k = j + 1; k <= i; k++) {
                sum += eL[i][k] * L[k][j];
            }
            eL[i][j] = -sum / L[j][j];
        }
    }
}



int main() {
    // ==============================================================================
    //  The main program to test the module InverseL.C
    // ==============================================================================
    const int n = 5;
    int i, j;
    double L[5][5] = {{ 1.0, 0.0, 0.0, 0.0, 0.0}, 
                      { 4.0, 3.0, 0.0, 0.0, 0.0},       
                      { 5.0, 2.0, 4.0, 0.0, 0.0},      
                      { 3.0, 1.0, 8.0, 2.0, 0.0},      
                      {-1.0, 1.0, 1.0, 1.0, 1.0}};
    double eL[5][5];
    
    printf(" Input Lower Matrix, L \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%0.7f ", L[i][j]);
        }
        printf("\n");
    }

    InverseL(n, L, eL);

    printf(" ------ Output : Inverse L ------ \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%0.7f ", eL[i][j]);
        }
        printf("\n");
    }
    return 0;
}