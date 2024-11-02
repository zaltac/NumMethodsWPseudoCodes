' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_AX

    Public Sub Main()
        Dim n As Integer = 3
        Dim i, j As Integer
        Dim A(n, n) As Double, x(n) As Double, b(n) As Double
        
        A = {{ 1.0,  2.0,  3.0, -2.0}, _
             { 4.0,  1.0,  2.0,  3.0}, _
             { 3.0,  2.0,  1.0,  2.0}, _
             {-2.0,  3.0,  4.0,  1.0}}
		
        x = {2.0, 5.0, 3.0, -3.0}
        
		Console.WriteLine("Input Matrix A")
        For i = 0 To n
            For j = 0 To n
				Console.Write("{0,8:F4} ", A(i, j))
            Next j
            Console.WriteLine()
        Next i
        
	    Console.WriteLine()
		Console.WriteLine("Input Vector  x")
        For i = 0 To n
			Console.WriteLine("  x({0}) = {1,8:F4}", i, x(i))
        Next i
        
        Call Ax(n, A, x, b) ' performs A*x = b matrix operation
        
        Console.WriteLine()
		Console.WriteLine(" Output vector is b(n) = A(n,n)X(n) ")
        For i = 0 To n -1
			Console.WriteLine("  b({0}) = {1,8:F4}", i, b(i))
        Next i
    End Sub

    Public Sub Ax(ByVal n As Integer, ByVal A(,) As Double, ByVal x() As Double, ByRef b() As Double)
        Dim i, k As Integer
'  ==================================================================================
'  CODE2.5-Ax.bas. A Basic (VB) Sub implementing Pseudocode 2.5.                                  
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
'  DESCRIPTION: A sub program to perform A * x = b matrix-vector multiplication.                
'                                                                                              
'  ON ENTRY                                                                                    
'     n   :: Dimension attributes of input/output matrices;                                    
'     A   :: An input matrix of size nÃ—n;                                                     
'     x   :: An input vector of length n.                                                      
'                                                                                              
'  ON EXIT                                                                                     
'     b   :: The output vector of length n.                                                    
'                                                                                              
'  REVISION DATE :: 03/18/2024                                                                 
'  ==================================================================================        
        For i = 0 To n
            b(i) = 0.0
            For k = 0 To n
                b(i) += A(i, k) * x(k)
            Next k
        Next i
    End Sub
End Module