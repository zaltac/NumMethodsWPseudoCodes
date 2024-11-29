import numpy as np


#  ==================================================================================
#  CODE5.1-EXAMPLE5-4.PY. A Python module implementing Pseudocode 5.1.                            
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
#  DESCRIPTION: A program that calculates the first and second derivatives of                
#    a set of position (discrete) data.                                                                   
#                                                                                              
#  INPUT VARIABLES                                                                             
#     dt  :: Time increment (in seconds);                                                      
#     s   :: An array of length n, providing the distance traveled in meters.                  
#                                                                                              
#  OUTPUT VARIABLES                                                                            
#     v   :: Array of length n containing velocities (m/s) at discrete points;                 
#     a   :: Array of length n containing acelerations (m2/s) at discrete points.              
#                                                                                              
#  REVISION DATE :: 06/13/2024                                                                 
#  ==================================================================================
def Example5_4():
    n = 6
    s = np.array([0.0, 5.45, 21.3, 82.84, 212.86, 473.6])
    dt = 5.0

    # Apply 2nd order forward difference formulas
    v = np.zeros(n)
    a = np.zeros(n)

    v[0] = (-s[2] + 4.0 * s[1] - 3.0 * s[0]) / (2.0 * dt)
    a[0] = (-s[3] + 4.0 * s[2] - 5.0 * s[1] + 2.0 * s[0]) / dt**2

    # Apply central difference formulas
    for i in range(1, n-1):
        v[i] = (s[i+1] - s[i-1]) / (2.0 * dt)
        a[i] = (s[i+1] - 2.0 * s[i] + s[i-1]) / dt**2

    # Apply 2nd order backward difference formulas
    v[-1] = (s[-3] - 4.0 * s[-2] + 3.0 * s[-1]) / (2.0 * dt)
    a[-1] = (-s[-4] + 4.0 * s[-3] - 5.0 * s[-2] + 2.0 * s[-1]) / dt**2

    # Print out the results
    t = np.arange(0, n * dt, dt)
    print(f"{' Time (s)':>7s} {'s (m)':>11s} {'v (m/s)':>14s} {'a (m/s^2)':>12s}")
    for i in range(n):
        print(f"{t[i]:8.3f} {s[i]:12.4f} {v[i]:12.4f} {a[i]:12.4f}")

Example5_4()