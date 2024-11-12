' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_CGM
    ' The main program to test SUBROUTINE CGM

' ==================================================================================
' Main program to test CGM.BAS
' ==================================================================================
	Public Sub Main()
        Dim n As Integer = 10
        Dim a(n - 1, n - 1) As Double
        Dim b(n - 1) As Double
        Dim x(n - 1) As Double
		Dim error_val, eps, guess As Double
        Dim i, j, iter, maxit As Integer

        a = {{ 2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             {-1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0}, _
             { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0}}

        b = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.2}

		guess = 1.0
        '-------- write input data ----------------------------

        Console.WriteLine("Coefficient Matrix    RHS    Initial Guess")
		For i = 0 To n - 1
			x(i) = guess
			For j = 0 To n - 1
				Console.Write("{0,8:F2} ", A(i, j))
            Next
        	Console.Write("{0,8:F2} ", b(i))
        	Console.WriteLine("{0,8:F2} ", x(i))
		Next

        eps = 1.0E-07
        maxit = 999

        ' ***  GO TO CGM MODULE 
        Call CGM(n, eps, a, b, x, maxit, iter, error_val)

        ' ***  PRINT OUT THE RESULTS
        Console.WriteLine()
        Console.WriteLine("----- Solution --------")
        For i = 0 To n - 1
            Console.WriteLine(" x({0}) ={1,12:F7}", i, x(i))
        Next
        Console.WriteLine("-----------------------")
        Console.WriteLine("Total no of iterations = {0}", iter)
		Console.WriteLine("Maximum Error          = {0,12:E4}", error_val)
    End Sub


	Public Sub CGM(ByVal n As Integer, ByVal eps As Double, ByRef A(,) As Double, ByRef b() As Double, ByRef x() As Double, ByVal maxit As Integer, ByRef iter As Integer, ByRef err As Double)
	'  ==================================================================================
	'  CODE3.3-CGM.bas. A Basic (VB) Sub implementing Pseudocode 3.3.                                 
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
	'  DESCRIPTION: A subroutine to solve Ax=b linear system with the Conjugate Gradient Method.   
	'                                                                                              
	'  ON ENTRY                                                                                    
	'     n   :: Number of equations;                                                              
	'     A   :: Coefficient matrix (nÃ—n);                                                        
	'     b   :: Array of length n containing the right hand side;                                 
	'     x   :: Array of length n containing the initial guess;                                   
	'    eps  :: Convergence tolerance;                                                            
	'   maxit :: Maximum number of iterations.                                                     
	'                                                                                              
	'  ON EXIT                                                                                     
	'     x   :: Array of length n containing the approximate solution;                            
	'   iter  :: Total number of iterations performed;                                             
	'   error :: Euclidean (L2-) norm of displacement at exit.                                     
	'                                                                                              
	'  USES                                                                                           
	'    AX   :: Subroutine to evaluate A*x matrix vector product;                                 
	'    SQRT :: Built-in Intrinsic function returning the square root of a real value;            
	'    XdotY:: Function to evaluate the dot product of two vectors.                              
	'                                                                                              
	'  REVISION DATE :: 12/11/2024                                                                 
	'  ==================================================================================
    	Dim c(n), r(n), d(n) As Double
		Dim rho0, rho, alpha, beta, enorm As Double
    	Dim i, p As Integer

    	Call AX(n, A, x, r)
		For i = 0 To n - 1
    		r(i) = r(i) - b(i)       ' [r]^(0) = [b] - [A][x]^(0)
    		d(i) = r(i)              ' [d]^(0) = [r]^(0)
        Next
		rho0 = XdotY(n, r, r)      ' rho^(0) = [r^(0)] · [r^(0)]

		Console.WriteLine(" ")
		Console.WriteLine("    *** Iteration history ***")
    
    	For p = 0 To maxit              ' BEGIN CGM ITERATION LOOP
			enorm = Math.Sqrt(rho0)     ' E-norm = Sqrt([r]^(p-1) · [r]^(p-1))
			Console.WriteLine("iter = {0,2:F0}    E-norm = {1,12:E4} ",p, Enorm) ' PRINTOUT ITERATION PROGRESS
    	    If Enorm < eps Then Exit For ' CHECK FOR CONVERGENCE, EXIT if converged
        	Call AX(n, A, d, c)         ' [c]^(p-1) = [A][d]^(p-1)
			rho = XdotY(n, d, c)        ' rho = [d]^(p-1) · [c]^(p-1)
    	    alpha = rho0 / rho          ' alpha^(p) = [r] · [r] / ([d] · [c])
			For k = 0 To n - 1
				x(k) = x(k) - alpha * d(k)       ' [x]^(p) = [x]^(p-1) + alfa^(p) · [d]^(p-1)
    			r(k) = r(k) - alpha * c(k)       ' [r]^(p) = [r]^(p-1) - alfa^(p) · [d]^(p-1)
        	Next
			rho = XdotY(n, r, r)        ' rho^(p+1) = [r^(p)] · [r^(p)]
        	beta = rho / rho0           ' beta^(p) = rho / rho0
			For k = 0 To n - 1
				d(k) = r(k) + beta * d(k) ' [d]^(p) = [r]^(p) + beta · [d]^(p-1)
			Next
			rho0 = rho                 ' rho^(p) <== rho^(p+1)
    	Next

    	iter = p                    ' Set total number of iterations
		err = enorm                 ' Set recent Enorm as error
    	If iter = maxit Then
	        Console.WriteLine("Failed to converge after {0} iterations ", maxit)
    	End If
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
	Public Function Enorm(ByVal n As Integer, ByVal x() As Double) As Double
        Dim sum As Double = 0.0
        For i = 0 To n - 1
            sum += x(i) ^ 2
        Next
        Return Math.Sqrt(sum)
    End Function

	Public Function XdotY(ByVal n As Integer, ByVal x() As Double, ByVal y() As Double) As Double
        ' ==============================================================================
        ' CODE2-3-XdotY.BAS. A Basic program for implementing Pseudocode 2.3
        '     
        ' NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
        ' First Edition. (c) By Zekeriya ALTAÇ (2024).
        ' ISBN: 9781032754741 (hbk)
        ' ISBN: 9781032756424 (pbk)
        ' ISBN: 9781003474944 (ebk)
        '
        ' DOI : 10.1201/9781003474944
        ' C&H/CRC PRESS, Boca Raton & London. 
        '  
        ' This free software is complimented by the author to accompany the textbook.
        ' E-mail: altacz@gmail.com
        '
        ' DESCRIPTION: A sub function to compute the dot product of two vectors, x and y.
        '
        ' ARGUMENTS
        '    n   :: Dimension attribute of the input vectors; 
        '   x, y :: The input vectors of length n.
        '
        ' REVISION DATE :: 03/18/2024
        ' ==============================================================================
        Dim sums As Double = 0.0
        Dim i As Integer

        For i = 0 To n-1
            sums += x(i) * y(i)
        Next

        XdotY = sums
    End Function

    Public Sub Ax(ByVal n As Integer, ByRef A(,) As Double, ByRef x() As Double, ByRef b() As Double)
    '  ==================================================================================
    '  CODE2.5-Ax.BAS. A Basic (VB) Sub implementing Pseudocode 2.5.                                  
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
        Dim i, k As Integer
        
        For i = 0 To n - 1
            b(i) = 0.0
            For k = 0 To n - 1
                b(i) += A(i, k) * x(k)
            Next k
        Next i
    End Sub
End Module