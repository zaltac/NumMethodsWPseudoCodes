' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System
				
Public Module Test_Simpsons_Rule_DF
    ' ==============================================================================
    ' The main program to test Sub program Simpsons_Rule_DF.BAS
    ' ==============================================================================
	Public Sub Main()
		Dim n, i As Integer
		Dim a, b, x, h, intg As Double
		Dim f() As Double

        Console.WriteLine("Enter Number of Panels ")
        n = Console.ReadLine()

		ReDim f(n)
		
        ' Construct a discrete function using Function FX
        a = 0.0
        b = 1.0
        h = (b - a) / CDbl(n) ' conver n to double
        x = a
        For i = 0 To n  ' Discrete dataset generation from y=FUNC(x)
            f(i) = x ^ 4
            x = x + h
        Next

	Call Simpsons_Rule_DF(n, h, f, intg)

        Console.WriteLine("=== With Simpson's 1/3 Rule ===")
	Console.WriteLine("Panel no. = {0,3:F0} Estimate = {1,10:F8} with Simpson's rule",n-1,intg)
    End Sub

	Public Sub Simpsons_Rule_DF(ByVal n As Integer, ByVal h As Double, ByRef f() As Double, ByRef intg As Double)
	'  ==================================================================================
	'  CODE8.2-Simpsons_Rule_DF.bas. A Basic (VB) Sub implementing Pseudocode 8.2.                    
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
	'  DESCRIPTION: A module to estimate the integral of a discrete function f on [a,b]                      
	'     using the Simpson's 1/3 rule.                                                              
	'                                                                                              
	'  ON ENTRY                                                                                    
	'     n   :: Number of panels (must be even!..);                                               
	'     h   :: Interval size (uniform spacing, x_{i+1}-x_i);                                     
	'     f   :: Array of length (n+1) containing the ordinates, f_0, f_1, ..., f_n.               
	'                                                                                              
	'  ON EXIT                                                                                     
	'   intg  :: Numerical estimate for the integral.                                              
	'                                                                                              
	'  USES                                                                                        
	'    Mod  :: Modulo function, returning the remainder after number is divided by divisor.      
	'                                                                                              
	'  REVISION DATE :: 03/04/2025                                                                 
	'  ==================================================================================
        Dim m, i As Integer
	Dim odd, even As Double

        ' Standard Simpson's Rule           
        m = n Mod 2
        If m <> 0 Then
		Console.WriteLine(" Number of panels is not EVEN")
		Return
        End If

        odd = 0.0
        For i = 1 To n - 1 Step 2
            odd = odd + f(i)
        Next
        even = 0.0
        For i = 2 To n - 2 Step 2
            even = even + f(i)
        Next

        intg = f(0) + f(n) + 4.0 * odd + 2.0 * even
        intg = intg * h / 3.0
    End Sub

End Module