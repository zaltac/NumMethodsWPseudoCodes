' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_SOR_Method

' ==============================================================================
'  The main program to test JACOBI.BAS
' ==============================================================================
	Public Sub Main()
	Const n As Integer = 10
	Dim a(n - 1, n - 1) As Double ' Coefficient matrix
    	Dim b(n - 1) As Double ' rhs vector
	Dim xo(n - 1) As Double ' initial guess at start
    	Dim x(n - 1) As Double ' solution
	Dim i, j, iter, maxit As Integer
	Dim Errmax, eps, omg, guess As Double
        
	a = {{ 2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             {-1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0} }

        b = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.2}
		
	guess = 1.0
	
        Console.WriteLine("Coefficient Matrix    RHS    Initial Guess")
		For i = 0 To n - 1
			xo(i) = guess
			For j = 0 To n - 1
				Console.Write("{0,8:F2} ", A(i, j))
         		Next
	        	Console.Write("{0,8:F2} ", b(i))
        		Console.WriteLine("{0,8:F2} ", xo(i))
		Next

        eps = 1.0E-6
        maxit = 999
        Call JACOBI(n, eps, a, b, xo, x, maxit, iter, Errmax)

        Console.WriteLine()
        Console.WriteLine("----- Solution --------")
        For i = 0 To n - 1
            Console.WriteLine(" x({0}) ={1,11:F6}", i, x(i))
        Next
        Console.WriteLine("-----------------------")
        Console.WriteLine("Total no of iterations = {0}", iter)
	Console.WriteLine("Maximum Error          = {0,10:E4}", Errmax)
    End Sub


'  ==================================================================================
'  CODE3.1-JACOBI.BAS. A Basic (VB) Sub implementing Pseudocode 3.1.                              
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
'  DESCRIPTION: A sub program to iteratively solve Ax=b using the Jacobi method.                
'                                                                                              
'  ON ENTRY                                                                                    
'     n   :: Number of equations (size of A);                                                  
'     A   :: Input coefficient matrix (nxn);                                                  
'     b   :: Array of length n containing the right-hand;                                      
'     x   :: Array of length n containing the estimate at (p+1)'th step;                       
'     xo  :: Array of length n containing the initial guess, or iterates at estimate at p'th st
'    eps  :: Convergence tolerance;                                                            
'   maxit :: Maximum permitted number of iterations.                                           
'                                                                                              
'  ON EXIT                                                                                     
'     x   :: Array of length n containing the estimated solution;                              
'   iter  :: Total number of iterations performed;                                             
'   error :: L2 norm of the displacement vector.                                               
'                                                                                              
'  USES                                                                                        
'    JACOBI_DRV :: Accompanying sub program performing one step Jacobi iteration.               
'                                                                                              
'  REVISION DATE :: 11/09/2024                                                                 
'  ==================================================================================
     Sub JACOBI(n As Integer, eps As Double, A(,) As Double, b() As Double, xo() As Double, x() As Double, maxit As Integer, ByRef iter As Integer, ByRef err As Double)
        Dim p As Integer
        Dim del, del0 As Double

        p = 0
        del0 = 1.0	
	Console.WriteLine(" ")
        Do
            p += 1
            Call JACOBI_DRV(n, A, b, xo, x, del)
	    Console.WriteLine("iter= {0}   Error=  {1,11:E4}   Ratio= {2,8:F5}", p, del, del / del0)
            
	    For i = 0 To n - 1
	        xo (i) = x (i)
	    Next
            del0 = del
            If (del < eps) OrElse (p = maxit) Then Exit Do
        Loop

        err = del
        iter = p
	Console.WriteLine("  Error= {0,11:E4} ", err)
        If p = maxit Then
            Console.WriteLine("Jacobi method faied to converge after {0} iterations within the specified EPS tolerance.", maxit)
        End If
    End Sub

'  ==================================================================================
'  CODE3.1-JACOBI.BAS. A Basic (VB) Sub implementing JACOBI_DRV sub of Pseudocode 3.1.                              
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
'  DESCRIPTION: A sub program to perform one step Jacobi iteration and compute               
'     the Euclidean norm of the displacement vector.                                           
'                                                                                              
'  ON ENTRY                                                                                    
'     n   :: Number of equations (size of A);                                                  
'     A   :: Input coefficient matrix (nxn);                                                  
'     b   :: Array of length n containing the right-hand;                                      
'     x   :: Array of length n containing the estimate at (p+1)'th step;                       
'     xo  :: Array of length n containing the estimate at p'th step.                           
'                                                                                              
'  ON EXIT                                                                                     
'     x   :: Array of length n containing the estimated solution;                              
'     del :: Maximum absolute error achieved.                                                  
'                                                                                              
'  USES                                                                                        
'   ENORM:: User-defined function calculating the Euclidean vector (L2 norm) of a vector.      
'                                                                                              
'  REVISION DATE :: 11/09/2024                                                                 
'  ==================================================================================
    Sub JACOBI_DRV(n As Integer, A(,) As Double, b() As Double, xo() As Double, ByRef x() As Double, ByRef del As Double)
	Dim d(n - 1) As Double
	Dim sums As Double 
        Dim i, j As Integer

        For i = 0 To n - 1
            sums = 0.0
            For j = 0 To n - 1
                If i <> j Then
                    sums += A(i, j) * xo(j)
                End If
            Next
            x(i) = (b(i) - sums) / A(i, i)
            d(i) = x(i) - xo(i)
        Next

        del = Enorm(n, d)
    End Sub

'  ==================================================================================
'  CODE3.1-JACOBI.BAS. A Basic (VB) Sub implementing ENORM function of Pseudocode 3.1.                              
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
'  DESCRIPTION: A function module to compute Euclidean (L2) norm of a vector.                  
'                                                                                              
'  ARGUMENTS                                                                                   
'      n  :: The length of an input vector;                                                    
'      x  :: A vector (array) of length n.                                                     
'                                                                                              
'  USES                                                                                        
'    SQRT :: Built-in Intrinsic function returning the square root of a real value.            
'                                                                                              
'  REVISION DATE :: 11/09/2024                                                                  
'  ==================================================================================
	Public Function Enorm(n As Integer, x() As Double) As Double
        Dim sum As Double = 0.0
        For i = 0 To n - 1
            sum += x(i) ^ 2
        Next
        Return Math.Sqrt(sum)
    End Function

End Module