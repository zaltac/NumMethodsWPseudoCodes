' ==============================================================================
'  The main program to test "QuadraticEq"
' ==============================================================================
Imports System

Public Module TestQuadraticEq
	Public Sub Main()
        Dim p, q As Double
        Dim xr(1), xi(1) As Double

        Console.WriteLine("Type in values for p and q")
        p = Console.ReadLine()
        q = Console.ReadLine()

        ' Solve quadratic equation
        Call QuadraticEq(p, q, xr, xi)

        ' Print results
        Console.WriteLine("1st Root {0:F4} + {1:F4} i", xr(0), xi(0))
        Console.WriteLine("2nd Root {0:F4} + {1:F4} i", xr(1), xi(1))
    End Sub

	Public Sub QuadraticEq(ByVal p As Double, ByVal q As Double, ByRef re() As Double, ByRef im() As Double)
'  ==================================================================================
'  CODE1.3-Quadratic_Eq.bas. A Basic sub for implementing Pseudocode 1.3.          
' 
'  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
'  First Edition. (c) By Zekeriya ALTAC (2024).
'  ISBN: 978-1-032-75474-1 (hbk)
'  ISBN: 978-1-032-75642-4 (pbk)
'  ISBN: 978-1-003-47494-4 (ebk)
'  
'  DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
'  
'  This free software is complimented by the author to accompany the textbook.
'  E-mail: altacz@gmail.com.
'                                                                                    
'  DESCRIPTION: A Subroutine to find the roots of a quadratic equation of            
'    the form : x*x + p * x + q = 0.                                                 
'                                                                                    
'  ON ENTRY                                                                          
'    p, q :: Coefficients of the quadratic equation;                                 
'                                                                                    
'  ON EXIT                                                                           
'    re   :: Array of length 2 containing real parts of the roots: re1, re2;         
'    im   :: Array of length 2 containing imaginary parts of the roots: im1, im2.    
'                                                                                    
'  USES                                                                              
'    Sqrt :: Built-in Intrinsic function to evaluate the square root of a real value.
'                                                                                    
'  REVISION DATE :: 03/18/2024                                                       
'                                                                                    
'  ==================================================================================
        Dim d As Double
        d = p * p - 4.0 * q
        If d < 0.0 Then ' real imaginary roots
            d = Math.Sqrt(-d)
            re(0) = -p / 2.0 : re(1) = re(0)
            im(0) = -d / 2.0 : im(1) = -im(0)
        Else ' d>=0, Real Roots 
            d = Math.Sqrt(d)
            re(0) = (-p - d) / 2.0 : im(0) = 0.0
            re(1) = (-p + d) / 2.0 : im(1) = 0.0
        End If
    End Sub
End Module