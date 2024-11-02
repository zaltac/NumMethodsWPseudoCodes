' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_Back_Substitute
	Public Sub Main()
        Dim n As Integer = 5
        Dim A(n-1, n-1) As Double, b(n-1) As Double, x(n-1) As Double
        Dim i, j As Integer

        A = {{1.0,  2.0,  4.0, -3.0,  1.0}, _
             {0.0,  3.0,  4.0, -4.0,  1.0}, _
             {0.0,  0.0,  4.0,  3.0,  1.0}, _
             {0.0,  0.0,  0.0,  2.0,  1.0}, _
             {0.0,  0.0,  0.0,  0.0,  1.0}}
        
        ' Initialize vector b
        b = {10.0, 7.0, 29.0, 13.0, 5.0}
        
        Console.WriteLine(" ********** Input Matrix ***********")
		For i = 0 To n - 1
            For j = 0 To n - 1
				Console.Write("{0,8:F2} ", A(i, j))
            Next
            Console.WriteLine()
        Next
        Console.WriteLine()
        
        Call Back_Substitute(n, A, b, x)
        
        Console.WriteLine(" *** RHS Vector ***")
        For i = 0 To n - 1
			Console.WriteLine("  b({0:F0}) = {1,8:F3} ",i, b(i))
        Next i
        Console.WriteLine(" ")
        Console.WriteLine(" *** Solution Vector ***")
        For i = 0 To n - 1
			Console.WriteLine("  x({0:F0}) = {1,8:F3} ",i, x(i))
        Next i
        
        Console.WriteLine("-----------------------------------")
    End Sub

	Public Sub Back_Substitute(n As Integer, U(,) As Double, b() As Double, x() As Double)
'  ==================================================================================
'  CODE2.8-BACK_SUBSTITUTE.bas. A Basic (VB) Sub implementing Pseudocode 2.8.                     
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
'  DESCRIPTION: A sub programe to find the solution of a upper-triangular system                 
'    using back substitution algorithm.                                                        
'                                                                                              
'  ON ENTRY                                                                                    
'      n  :: Number of unknowns;                                                               
'      U  :: Input coefficient (upper-triangular) matrix (nxn);                               
'      b  :: Input array of size n containing the rhs.                                         
'                                                                                              
'  ON EXIT                                                                                     
'      x  :: Output array of size n containing the solution.                                   
'                                                                                              
'  REVISION DATE :: 03/18/2024                                                                 
'  ==================================================================================
        Dim sums As Double
        Dim k, j As Integer
        
        x(n - 1) = b(n - 1) / U(n - 1, n - 1)
        For k = n - 2 To 0 Step -1
            sums = 0.0
            For j = k + 1 To n - 1
                sums += U(k, j) * x(j)
            Next
            x(k) = (b(k) - sums) / U(k, k)
        Next
    End Sub
End Module