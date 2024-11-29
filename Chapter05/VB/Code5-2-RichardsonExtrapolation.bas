' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System


' ==============================================================================
'  The main program to test Sub program RICHARDSON
' ==============================================================================
Public Module Test_Richardson
	Public Sub Main()
        Dim D(0 To 10, 0 To 10) As Double
        Dim x0, h, err, eps, deriv As Double
        Dim nr As Integer

        x0 = 1.0
        h = 0.1
        err = 1.0
        eps = 1.0E-04

        Call Richardson(x0, h, eps, D, nr, deriv)

        'For k = 0 To nr
        '    Console.WriteLine(k.ToString() & " " & String.Join(" ", D(k, 0 To k).Select(Function(x) x.ToString("G19.11"))))
        'Next
        Console.WriteLine("---------------------------------")
		Console.WriteLine("Derivative is = {0,12:F8} ", deriv)
        Console.WriteLine("---------------------------------")
    End Sub


	'  ==================================================================================
	'  CODE5.2-RICHARDSON.BAS. A Basic (VB) Sub implementing Pseudocode 5.2.                          
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
	'  DESCRIPTION: A module to compute the first derivative of an explicitly                      
	'      defined function using Richardson's extrapolation.                                    
	'                                                                                              
	'  ON ENTRY                                                                                    
	'      x0  :: Point at which derivative is to be computed;                                     
	'      h   :: Initial interval size;                                                           
	'      eps :: Tolerance desired.                                                               
	'                                                                                              
	'  ON EXIT                                                                                     
	'      D   :: A matrix containing the Richardson's table (0..n, 0..n)                        
	'      nr  :: Size of the table;                                                               
	'    deriv :: Estimated derivative.                                                            
	'                                                                                              
	'  USES                                                                                        
	'     ABS  :: Built-in function in MATH returning the absolute value of a real value. 
	'                                                                                              
	'  ALSO REQUIRED                                                                               
	'     FUNC  :: User-defined external function providing the nonlinear equation.                
	'                                                                                              
	'  REVISION DATE :: 06/13/2024                                                                 
	'  ==================================================================================
	Public Sub Richardson(ByVal x0 As Double, ByRef h As Double, ByRef eps As Double, ByRef D(,) As Double, ByRef nr As Integer, ByRef deriv As Double)
        Dim err As Double 
        Dim k, m As Integer

        m = 0
        k = 0
        err = 1.0

        While err > eps
            D(k, 0) = (FUNC(x0 + h) - FUNC(x0 - h)) / (2.0 * h) ' 1st derivative
			Console.Write(" {0,6:F4} {1,12:F9} ", h, D(k,0))
            For m = 1 To k
                D(k, m) = (4 ^ m * D(k, m - 1) - D(k - 1, m - 1)) / CDbl(4 ^ m - 1)
				Console.Write(" {0,12:F9} ", D(k,m))
            Next
			Console.WriteLine(" ")
            If k >= 1 Then ' Estimate diagonalwise differentiation error
                err = Math.Abs(D(k, k) - D(k - 1, k - 1))
            End If
            h = h / 2.0
            k = k + 1
		'Console.WriteLine(" {0,10:F6}",x0)
		'Console.WriteLine("Derivative is {0,12:F8} ", D(k,m))
        End While
        nr = k - 1
        deriv = D(nr, nr)
    End Sub

	'  ==================================================================================
	'  USER-DEFINED FUNCTION "FUNC" OF ONE-VARIABLE
	'  ==================================================================================
	Public Function FUNC(ByVal x As Double) As Double
        Return 25000.0 / (-57.0 + x) - 5.2e6 / (x * x)
    End Function
End Module