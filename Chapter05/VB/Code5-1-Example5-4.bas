' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Example5_4
	'  ==================================================================================
	'  CODE5.1-EXAMPLE5-4.BAS. A Basic (VB) Sub implementing Pseudocode 5.1.                          
	' 
	'  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
	'  First Edition. (c) By Zekeriya ALTAÇ (2024).
	'  ISBN: 978-1-032-75474-1 (hbk)
	'  ISBN: 978-1-032-75642-4 (pbk)
	'  ISBN: 978-1-003-47494-4 (ebk)
	'  
	'  DOI : 10.1201/9781003474944
	'  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
	'  
	'  This free software is complimented by the author to accompany the textbook.
	'  E-mail: altacz@gmail.com.
	'  
	'  DESCRIPTION: A program that calculates the first and second derivatives of                
	'    a set of position (discrete) data.                                                                   
	'                                                                                              
	'  INPUT VARIABLES                                                                             
	'     dt  :: Time increment (in seconds);                                                      
	'     s   :: An array of length n, providing the distance traveled in meters.                  
	'                                                                                              
	'  OUTPUT VARIABLES                                                                            
	'     v   :: Array of length n containing velocities (m/s) at discrete points;                 
	'     a   :: Array of length n containing acelerations (m2/s) at discrete points.              
	'                                                                                              
	'  REVISION DATE :: 06/13/2024                                                                 
	'  ==================================================================================
	Public Sub Main()
        Dim n As Integer = 6
        Dim s(n), v(n), a(n) As Double
        Dim dt As Double
        Dim t As Double
        Dim i, m As Integer

		m = n - 1
        dt = 5.0
        s = {0.0, 5.45, 21.3, 82.84, 212.86, 473.60}

        ' *** Apply 2nd order forward difference formulas
        v(0) = (-s(2) + 4.0 * s(1) - 3.0 * s(0)) / (2.0 * dt)
        a(0) = (-s(3) + 4.0 * s(2) - 5.0 * s(1) + 2.0 * s(0)) / dt ^ 2

        ' *** Apply central difference formulas
        For i = 1 To (m - 1)
            v(i) = (s(i + 1) - s(i - 1)) / (2.0 * dt)
            a(i) = (s(i + 1) - 2.0 * s(i) + s(i - 1)) / dt ^ 2
        Next

        ' *** Apply 2nd order backward difference formulas
        v(m) = ( s(m - 2) - 4.0 * s(m - 1) + 3.0 * s(m)) / (2.0 * dt)
        a(m) = (-s(m - 3) + 4.0 * s(m - 2) - 5.0 * s(m - 1) + 2.0 * s(m)) / dt ^ 2

        ' *** Print out the results
        Console.WriteLine("Time (s)       s (m)          v (m/s)       a (m/s^2)")
        t = 0.0
		For i = 0 To m
			Console.WriteLine("{0,4:F0} {1,15:F7} {2,15:F7} {3,15:F7}  ",t,s(i),v(i),a(i))
            t += dt
        Next
    End Sub
End Module