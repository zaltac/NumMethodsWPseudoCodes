#include <iomanip>
#include <iostream>
#include <cmath>

// ==================================================================================
// CODE5.1-EXAMPLE5-4.CPP. A C++ module implementing Pseudocode 5.1.                                  
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
// DESCRIPTION: A program that calculate the first and second derivatives of                
//   a set of position (discrete) data.                                                                   
//                                                                                             
// INPUT VARIABLES                                                                             
//    dt  :: Time increment (in seconds);                                                      
//    s   :: An array of length n, providing the distance traveled in meters.                  
//                                                                                             
// OUTPUT VARIABLES                                                                            
//    v   :: Array of length n containing velocities (m/s) at discrete points;                 
//    a   :: Array of length n containing acelerations (m2/s) at discrete points.              
//                                                                                             
// REVISION DATE :: 06/13/2024                                                                 
// ==================================================================================
int main() {
    const int n = 6;
    double s[n], v[n], a[n], dt, t;

    dt = 5.0;
    s[0] = 0.0; s[1] = 5.45; s[2] = 21.3; s[3] = 82.84; s[4] = 212.86; s[5] = 473.6;

    // Apply 2nd order forward difference formulas
    v[0] = (-s[2] + 4.0 * s[1] - 3.0 * s[0]) / (2.0 * dt);
    a[0] = (-s[3] + 4.0 * s[2] - 5.0 * s[1] + 2.0 * s[0]) / (dt * dt);

    // Apply central difference formulas
    for (int i = 1; i < n - 1; i++) {
        v[i] = (s[i + 1] - s[i - 1]) / (2.0 * dt);
        a[i] = (s[i + 1] - 2.0 * s[i] + s[i - 1]) / (dt * dt);
    }

    // Apply 2nd order backward difference formulas
    v[n - 1] = (s[n - 3] - 4.0 * s[n - 2] + 3.0 * s[n - 1]) / (2.0 * dt);
    a[n - 1] = (-s[n - 4] + 4.0 * s[n - 3] - 5.0 * s[n - 2] + 2.0 * s[n - 1]) / (dt * dt);

    // Print out the results
    t = 0.0;
    std::cout << "Time (s)     s (m)         v (m/s)      a (m/s^2)" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << std::fixed << std::setprecision(3) << t << "   " << s[i] << "   " << v[i] << "   " << a[i] << std::endl;
        t += dt;
    }

    return 0;
}