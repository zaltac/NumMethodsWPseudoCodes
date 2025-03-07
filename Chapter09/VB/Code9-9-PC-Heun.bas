' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================

Imports System
				
Public Module Test_PC_Heuns_Method
' ==============================================================================
' The main program to test PC_Heun.BAS
' ==============================================================================
	Public Sub Main()
        Dim x0 As Double = 0.0
        Dim y0 As Double = 2.0
        Dim h As Double = 0.1
        Dim xlast As Double = 2.0

        Call PC_Heun(h, x0, y0, xlast)
    End Sub

	Public Sub PC_Heun(h As Double, x0 As Double, y0 As Double, xlast As Double)
	'  ==================================================================================
	'  CODE9.9-PC_HEUNS_METHOD.bas. A Basic (VB) Sub implementing Pseudocode 9.9.                     
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
	'  DESCRIPTION: A SUB program to estimate the solution of a first order IVP on [x0,xlast]           
	'    using the Heun's Predictor-Correcter method.                                              
	'                                                                                              
	'   ON ENTRY                                                                                   
	'    h     :: Step size (it must be uniform);                                                  
	'    x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps
	'    xlast :: End point of the solution interval.                                              
	'                                                                                              
	'   Other Internal Variables                                                                   
	'    x,y   :: Current estimates, x^(p+1) and y^(p+1).                                          
	'                                                                                              
	'  USES                                                                                        
	'    ABS  :: Built-in Intrinsic function returning the absolute value of a real value.         
	'    FCN  :: User-defined external function providing y'=f(x,y).                               
	'                                                                                              
	'  REVISION DATE :: 03/07/2025                                                                 
	'  ==================================================================================
        Dim x As Double
        Dim y As Double
        Dim ys As Double
        Dim k1 As Double
        Dim k2 As Double
        Dim err As Double

	Console.WriteLine("{0,5} {1,10} {2,18}", "x", "y", "Abs. Error")
        Console.WriteLine(New String("-"c, 42))
	Console.WriteLine("{0,8:F5} {1,12:F8}", x0, y0)
 
        While x < xlast
            ' ==== PREDICTOR STEP
            k1 = h * FCN(x0, y0)
            ys = y0 + k1

            ' ==== CORRECTOR STEP
            x = x0 + h
            k2 = h * FCN(x, ys)
            y = y0 + 0.5 * (k1 + k2)
            err = Math.Abs(y - Exact(x))
			Console.WriteLine("{0,8:F5} {1,12:F8} {2,14:E4}", x, y, err)
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