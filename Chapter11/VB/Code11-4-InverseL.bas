' runs on    https://dotnetfiddle.net/
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================

Imports System

Public Module Test_InverseL
	' ==============================================================================
	'  The main program to test Sub program InverseL.BAS
	' ==============================================================================
	Public Sub Main()
        Const n As Integer = 5
        Dim i, j As Integer
        Dim L(n, n), eL(n, n) As Double
        
        L = {{ 1.0, 0.0, 0.0, 0.0, 0.0},
             { 4.0, 3.0, 0.0, 0.0, 0.0},
             { 5.0, 2.0, 4.0, 0.0, 0.0},
             { 3.0, 1.0, 8.0, 2.0, 0.0},
             {-1.0, 1.0, 1.0, 1.0, 1.0} }
        
	Console.WriteLine(" ----- Input Lower Tridiagonal Matrix L -----")
        For i = 0 To n - 1
            For j = 0 To n - 1
		Console.Write("{0,10:F7} ", L(i, j))
            Next
            Console.WriteLine()
        Next
        
        Call InverseL(n, L, eL)

        Console.WriteLine()
        Console.WriteLine(" ------ Output : Inverse L ------")
        For i = 0 To n - 1
            For j = 0 To n - 1
				Console.Write("{0,10:F7} ", eL(i, j))
            Next
            Console.WriteLine()
        Next
    End Sub

	Public Sub InverseL(n As Integer, ByRef L(,) As Double, ByRef eL(,) As Double)
	'  ==================================================================================
	'  CODE11.4-InverseL.bas. A Basic (VB) Sub implementing Pseudocode 11.4.                          
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
	'   DESCRIPTION: A VB Sub program to invert a lower-triangular matrix.                                 
	'                                                                                              
	'   ON ENTRY                                                                                   
	'     n    :: Dimension attribute of the matrix L;                                             
	'     L    :: A lower-triangular matrix, nxn.                                                  
	'                                                                                              
	'   ON EXIT                                                                                    
	'     eL   :: Inverse of L, also a lower-triangular matrix, nxn.                               
	'                                                                                              
	'   REVISION DATE :: 03/15/2025                                                                
	'  ==================================================================================
        Dim sums As Double
        Dim i, j, k As Integer
        
        eL = New Double(n, n) {} ' Initialize matrix eL with zeros.
        
        eL(0, 0) = 1.0 / L(0, 0)
        For i = 1 To n - 1
            eL(i, i) = 1.0 / L(i, i)
            For j = i - 1 To 0 Step -1
                sums = 0.0
                For k = j + 1 To i
                    sums += eL(i, k) * L(k, j)
                Next
                eL(i, j) = -sums / L(j, j)
            Next
        Next
    End Sub

End Module 