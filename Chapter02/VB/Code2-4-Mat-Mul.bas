' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_MAT_MUL
	Public Sub Main()
        Dim m As Integer = 2, p As Integer = 4, n As Integer = 3
        Dim i As Integer
        Dim A(1, 3) As Double, B(3, 2) As Double, C(1, 2) As Double

        A = {{1.0, -3.0, 5.0, 1.0}, {-2.0, 4.0, 1.0, 2.0}}
        B = {{3.0, 1.0, -2.0}, {2.0, 0.0, 4.0}, {-1.0, 1.0, -3.0}, {2.0, 5.0, 3.0}}

		Console.WriteLine("Input Matrix A")
		For i = 0 To m - 1
		For j = 0 To p - 1
				Console.Write("{0,8:F4} ", A(i, j))
            Next
            Console.WriteLine()
        Next
        
		Console.WriteLine()
		Console.WriteLine("Input Matrix B")
        
		For i = 0 To p - 1
            For j = 0 To n - 1
				Console.Write("{0,8:F4} ", B(i, j))
            Next
            Console.WriteLine()
        Next

        Call MAT_MUL(m, p, n, A, B, C)

        Console.WriteLine()
		Console.WriteLine("Output Matrix C(m,n) = A(m,p)xB(p,n) matrix product")
		For i = 0 To m - 1
            For j = 0 To n - 1
				Console.Write("{0,8:F4} ", C(i, j))
            Next
            Console.WriteLine()
        Next
    End Sub

	Public Sub MAT_MUL(ByVal m As Integer, ByVal p As Integer, ByVal n As Integer, ByVal A(,) As Double, ByVal B(,) As Double, C(,) As Double)
        Dim i As Integer, j As Integer, k As Integer
'  ==================================================================================
'  CODE2.4-MAT_MUL.BAS. A Basic (VB) Sub implementing Pseudocode 2.4.                             
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
'   DESCRIPTION: A sub program to find A*B=C matrix multiplication.                             
'                                                                                              
'   ON ENTRY                                                                                   
'    m,p,n :: Dimension attributes of input/output matrices;                                   
'       A  :: An input matrix of size mÃ—p;                                                    
'       B  :: An input matrix of size pÃ—n.                                                    
'                                                                                              
'   ON EXIT                                                                                    
'       C  :: The output matrix of size mÃ—n.                                                  
'                                                                                              
'   REVISION DATE :: 03/18/2024                                                                
'  ==================================================================================
        For i = 0 To m - 1
            For j = 0 To n - 1
                C(i, j) = 0.0
                For k = 0 To p - 1
                    C(i, j) = C(i, j) + A(i, k) * B(k, j)
                Next
            Next
        Next
    End Sub
End Module
