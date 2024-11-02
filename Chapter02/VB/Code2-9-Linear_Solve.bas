' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_Linear_Solve
	Public Sub Main()
        Dim n As Integer = 7
        Dim A(n - 1, n - 1) As Double
        Dim b(n - 1) As Double
        Dim x(n - 1) As Double
        Dim i As Integer, j As Integer, opt As Integer
 
        A = {{  2.0,  1.0, -1.0,  0.0,  3.0,  1.0,  0.0}, _
             { -1.0,  3.0,  1.0,  2.0,  4.0, -2.0,  1.0}, _
             {  4.0,  1.0,  5.0,  1.0, -3.0,  2.0,  2.0}, _
             { -2.0,  1.0, -2.0, -2.0,  3.0,  3.0,  1.0}, _
             {  1.0,  1.0,  1.0,  3.0,  2.0, -3.0,  2.0}, _
             {  4.0,  1.0, -1.0,  0.0,  3.0, -5.0,  1.0}, _
             {  2.0,  7.0, -3.0, -4.0, -1.0,  2.0,  2.0} }

        b = {-10.0, 7.0, 28.0, -30.0, 16.0, 3.0, -12.0}

		Console.WriteLine(" ********** Coefficient Matrix ***********")
		For i = 0 To n - 1
            For j = 0 To n - 1
				Console.Write("{0,8:F2} ", A(i, j))
            Next
            Console.WriteLine()
        Next

        Console.WriteLine(" ")
        Console.WriteLine(" *** RHS Vector ***")
        For i = 0 To n - 1
			Console.WriteLine("b({0:F0}) = {1,8:F3} ",i, b(i))
        Next i

        opt = 1
        Call Linear_Solve(n, A, b, x, opt)
        
		Console.WriteLine(" ")
		If opt = 0 Then
			Console.WriteLine("Naive Gauss Elimination was used.")
        Else
			Console.WriteLine("Gauss Jordan Elimination was used.")
        End If

        Console.WriteLine(" ")
        Console.WriteLine(" *** Solution Vector ***")
        For i = 0 To n - 1
			Console.WriteLine("x({0:F0}) = {1,8:F3} ",i, x(i))
        Next i
	End Sub

'  ==================================================================================
'  CODE2.9-LINEAR_SOLVE.bas. A Basic (VB) Sub implementing Pseudocode 2.9.                        
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
'  DESCRIPTION: A sub program to solve a system of linear equations using naive                 
'    Gauss Elimination (opt=0) or Gauss-Jordan Elimination (opt/=0) algorithm.                 
'                                                                                              
'  ON ENTRY                                                                                    
'      n  :: Number of unknowns;                                                               
'      A  :: Input coefficient matrix of size nÃ—n;                                            
'      b  :: An input array of length n containing the rhs;                                    
'     opt :: Option key (=0, Naive Gauss-Elim.; /=0, Gauss-Jordan Elimn).                      
'                                                                                              
'  ON EXIT                                                                                     
'      x  :: The output array of length n containing the solution.                             
'                                                                                              
'  USES                                                                                        
'    Math.Abs  :: Built-in Math library function returning the absolute value of a real value;         
'    Back_Substitute :: A subrotine to solve an upper-triangular system.                       
'                                                                                              
'  REVISION DATE :: 03/18/2024                                                                 
'  ==================================================================================
    Sub Linear_Solve(n As Integer, A As Double(,), b As Double(), x As Double(), opt As Integer)
        Dim i As Integer, j As Integer, k As Integer, m As Integer
        Dim eps As Double = 1.0E-12
        Dim val As Double, ajj As Double, s As Double

        ' Forward Elimination Steps
        For j = 0 To n - 1
            ajj = A(j, j)
            If Math.Abs(ajj) < eps Then
                Console.WriteLine("Pivot is zero at j=" & (j + 1))
                Console.WriteLine("Execution is halted!")
                Stop
            Else
                val = 1.0 / ajj
                b(j) *= val
                For m = j To n - 1
                    A(j, m) *= val
                Next
                For i = j + 1 To n - 1
                    s = A(i, j)
                    A(i, j) = 0.0
                    For k = j + 1 To n - 1
                        A(i, k) = A(i, k) - s * A(j, k)
                    Next
                    b(i) = b(i) - s * b(j)
                Next
            End If
        Next

        ' Back Substitution Steps
        If opt = 0 Then
            Call Back_Substitute(n, A, b, x)
        Else
            For j = n - 1 To 1 Step -1
                For i = j - 1 To 0 Step -1
                    b(i) = b(i) - A(i, j) * b(j)
                    A(i, j) = 0.0
                Next
                x(j) = b(j)
            Next
            x(0) = b(0)
        End If
    End Sub


	Public Sub Back_Substitute(ByVal n As Integer, ByVal U(,) As Double, ByVal b() As Double, x() As Double)
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