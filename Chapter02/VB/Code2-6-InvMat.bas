' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_Inv_Mat
	Public Sub Main()
        Dim n As Integer = 5
        Dim i, j As Integer
		Dim A(n - 1, n - 1), AI(n - 1, n - 1), C(n - 1, n - 1), B(n - 1, n - 1) As Double

        A = {{ 1.0,  2.0,  1.0,  2.0,  3.0}, _
             {11.0, -1.0,  1.0,  4.0,  1.0}, _
             { 4.0, -1.0,  1.0,  1.0, -1.0}, _
             {-3.0,  1.0, -8.0, -1.0,  5.0}, _
             {-1.0,  1.0,  1.0,  1.0,  1.0}}

        Console.WriteLine("-------- Input matrix A(-1) -------")
        For i = 0 To n - 1
            For j = 0 To n - 1
				Console.Write("{0,8:F2} ", A(i, j))
            Next
            Console.WriteLine()
        Next
        Console.WriteLine()

        For i = 0 To n - 1
            For j = 0 To n - 1
                B(i, j) = A(i, j)  ' A will be destroyed so save A on B for checking
            Next
        Next

        Console.WriteLine()

        Call Inv_Mat(n, B, AI)

        Console.WriteLine("------ Inverse matrix A(-1) -------")
        For i = 0 To n - 1
            For j = 0 To n - 1
				Console.Write("{0,8:F2} ", AI(i, j))
            Next
            Console.WriteLine()
        Next
        Console.WriteLine(" ")

        call MAT_MUL(n, n, n, A, AI, C)  ' performs A**(-1)*A=C matrix operation

        ' Check the accuracy of the result
        Console.WriteLine("------ matrix A*A(-1) = C-------------")
        For i = 0 To n - 1
            For j = 0 To n - 1
				Console.Write("{0,1:F0} ", C(i, j))
            Next
            Console.WriteLine()
        Next
        Console.WriteLine("-----------------------------------")
    End Sub
    ' *********END OF MAIN PROGRAM *********************************

'  ==================================================================================
'  CODE2.6-INV_MAT.bas. A Basic (VB) Sub implementing Pseudocode 2.6.                             
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
'  DESCRIPTION: A sub program to find the inverse of a square matrix (with no pivoting).        
'                                                                                              
'  ON ENTRY                                                                                    
'     n  :: Dimension attribute of input matrix A;                                             
'     A  :: An input matrix (nxn).                                                             
'                                                                                              
'  ON EXIT                                                                                     
'     B  :: Inverse of A (nxn).                                                                
'                                                                                              
'  REVISION DATE :: 03/18/2024                                                                 
'  ==================================================================================
	Public Sub Inv_Mat(n As Integer, ByRef A(,) As Double, ByRef B(,) As Double)
        Dim p, s As Double
        Dim i, j, k As Integer

        B = New Double(n, n) {}

        For i = 0 To n - 1
            B(i, i) = 1.0  ' Set B = I, i.e., construct identity matrix
        Next

        For j = 0 To n - 1
            p = 1.0 / A(j, j)
            For k = 0 To n - 1
                A(j, k) = p * A(j, k)
                B(j, k) = p * B(j, k)
            Next
            For i = 0 To n - 1
                s = A(i, j)
                If i <> j Then
                    For k = 0 To n - 1
                        A(i, k) = A(i, k) - s * A(j, k)
                        B(i, k) = B(i, k) - s * B(j, k)
                    Next
                End If
            Next
        Next
    End Sub

	Public Sub MAT_MUL(m As Integer, p As Integer, n As Integer, A(,) As Double, B(,) As Double, C(,) As Double)
        Dim i As Integer, j As Integer, k As Integer

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