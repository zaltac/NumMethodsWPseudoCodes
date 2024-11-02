' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_XdotY
    ' ==============================================================================
    ' The main program to test FUNCTION XdotY
    ' ==============================================================================
	Public Sub Main()
        Const n As Integer = 4
        Dim i As Integer
        Dim X(n) As Double
        Dim Y(n) As Double
        Dim result As Double
		
        X = {1.0, 2.0, 3.0, 4.0}
        Y = {4.0, 3.0, 2.0, 1.0}

		Console.WriteLine("Input Vector X")
		
        For i = 0 To n-1
			Console.WriteLine("x( {0:F1} )= {1:F3} ", i, X(i))
        Next

	    Console.WriteLine()
        Console.WriteLine("Input Vector Y")
        For i = 0 To n-1
            Console.WriteLine("y( {0:F1} )= {1:F3} ", i, Y(i))
        Next

        Console.WriteLine()
        result = XdotY(n, X, Y)
        Console.WriteLine(" X*Y dot product is ==> " & result)
    End Sub

	Public Function XdotY(ByVal n As Integer, ByVal x() As Double, ByVal y() As Double) As Double
'  ==================================================================================
'  CODE2.3-XdotY.BAS. A Basic (VB) Sub implementing Pseudocode 2.3.                               
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
'   DESCRIPTION: A function to compute the dot product of two vectors, x and y.                
'                                                                                              
'   ARGUMENTS                                                                                  
'      n   :: Dimension attribute of the input vectors;                                        
'     x, y :: The input vectors of length n.                                                   
'                                                                                              
'   REVISION DATE :: 03/18/2024                                                                
'  ==================================================================================
        Dim sums As Double = 0.0
        Dim i As Integer

        For i = 0 To n-1
            sums += x(i) * y(i)
        Next

        XdotY = sums
    End Function
End Module