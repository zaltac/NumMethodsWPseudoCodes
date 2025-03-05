' runs on    https://dotnetfiddle.net/8K7fkp
' ==============================================================================
'  The main program to test Sub program Explicit_Euler.BAS
' ==============================================================================
Imports System

Public Module Test_Explicit_Euler
' ==============================================================================
'  The main program to test Explicit_Euler.BAS
' ==============================================================================
	Public Sub Main()
        Dim x0 As Double = 0.0
        Dim y0 As Double = 2.0
        Dim xlast As Double = 1.0
        Dim h As Double = 0.1

        Call Explicit_Euler(h, x0, y0, xlast)
    End Sub


	Public Sub Explicit_Euler(h As Double, x0 As Double, y0 As Double, xlast As Double)
	'  ==================================================================================
	'  CODE9.1-Explicit_Euler.bas. A Basic (VB) Sub implementing Pseudocode 9.1.                      
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
	'  DESCRIPTION: A module to estimate the solution of a first order IVP on [x0,xlast]           
	'    using the Explicit Euler method. Numerical estimates are printed out, not stored.         
	'    With minor modifications, this PROGRAM can also be used to solve explicit methods         
	'    such as MIDPOINT RULE and MODIFIED EULER.                                                 
	'                                                                                              
	'  ON ENTRY                                                                                    
	'   h     :: Step size (it must be uniform);                                                   
	'   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
	'   xlast :: End point of the solution interval.                                               
	'                                                                                              
	'  Other Internal Variables                                                                             
	'   x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
	'                                                                                              
	'  USES                                                                                        
	'    Abs  :: Built-in Math library function returning the absolute value of a double float value;   
	'    FCN  :: User-defined external function providing y'=f(x,y).                               
	'                                                                                              
	'  REVISION DATE :: 03/05/2025                                                                 
	'  ==================================================================================
        Dim x As Double = x0
        Dim y As Double = y0

        Console.WriteLine(New String("-"c, 42))
	Console.WriteLine("{0,5} {1,10} {2,18}", "x", "y", "Abs. Error")
        Console.WriteLine(New String("-"c, 42))

        While x < xlast
            x = x0 + h
            y = y0 + h * FCN(x0, y0)      ' Explicit-Euler
            Console.WriteLine("{0,8:F5} {1,12:F8} {2,14:E4}", x, y, Math.Abs(y - Exact(x)))
            x0 = x
            y0 = y
        End While
    End Sub

	Public Function FCN(x As Double, y As Double) As Double
	' ==============================================================================
	' DESCRIPTION: A function subprogram providing y'=f(x,y)
	'
	' ARGUMENTS:
	'      x, y  :: Real input values.
	'
	' ==============================================================================
        Return -y/(x + 1.0)
    End Function

	Public Function Exact(x As Double) As Double
	' ==============================================================================
	' DESCRIPTION: A function subprogram providing the true solution y=(x) for 
	'    testing the module. 
	'
	' ARGUMENTS:
	'      x   :: A real input, independent variable.
	'
	' ==============================================================================
        Return 2.0 /( x + 1.0 )
    End Function
End Module