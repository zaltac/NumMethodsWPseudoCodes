' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_MAT_ADD
	Public Sub Main()
        Dim n As Integer = 3, m As Integer = 3
        Dim i As Integer, j As Integer
		Dim A(m - 1, n - 1) As Double, B(m - 1, n - 1) As Double, C(m - 1, n - 1) As Double, D(m - 1, n - 1) As Double

        A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}
        B = {{3.0, -2.0, 1.0}, {-2.0, 2.0, 4.0}, {3.0, -5.0, 1.0}}

        Console.WriteLine("Input Matrix A")
        For i = 0 To m - 1
            Console.WriteLine(String.Format("{0,10:F5} {1,10:F5} {2,10:F5}", A(i, 0), A(i, 1), A(i, 2)))
        Next

		Console.WriteLine()
        Console.WriteLine("Input Matrix B")
        For i = 0 To m - 1
            Console.WriteLine(String.Format("{0,10:F5} {1,10:F5} {2,10:F5}", B(i, 0), B(i, 1), B(i, 2)))
        Next

        Call MAT_ADD(m, n, A, B, C)

		Console.WriteLine()
        Console.WriteLine("Output Matrix C")
        For i = 0 To m - 1
            Console.WriteLine(String.Format("{0,10:F5} {1,10:F5} {2,10:F5}", C(i, 0), C(i, 1), C(i, 2)))
        Next

    End Sub


'  ==================================================================================
'  CODE2.1-MAT_ADD.BAS. A Basic (VB) Sub implementing Pseudocode 2.1.                             
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
'   DESCRIPTION: A Sub program to perform C=A+B matrix addition.                                
'                                                                                              
'   ON ENTRY                                                                                   
'     m,n :: Dimension attributes of the matrices;                                             
'      A  :: An input matrix (mxn);                                                            
'      B  :: An input matrix (mxn).                                                            
'                                                                                              
'   ON EXIT                                                                                    
'      C :: The output matrix (mxn).                                                           
'                                                                                              
'   REVISION DATE :: 03/18/2024                                                                
'  ==================================================================================
	Public Sub MAT_ADD(ByVal m As Integer, ByVal n As Integer, ByVal A(,) As Double, ByVal B(,) As Double, C(,) As Double)
        Dim i As Integer, j As Integer

        For i = 0 To m - 1
            For j = 0 To n - 1
                C(i, j) = A(i, j) + B(i, j)
            Next
        Next
    End Sub
End Module