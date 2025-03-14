' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================

Imports System
				
Public Module Test_POWER_METHOD_S

	Public Sub Main()
	' ==============================================================================
	'  The main program to test the module POWER_METHOD_S.BAS
	' ==============================================================================
        Dim n As Integer = 4
        Dim maxit As Integer
        Dim A(n, n) As Double
        Dim x(n) As Double
        Dim eps, Aerr, lambda As Double
        
        ' Initialize matrix A
        A = {{  8.0,   9.0,  10.0,   9.0},
             {-13.0, -12.0, -12.0, -11.0},
             {-18.0,  -9.0, -20.0,  -9.0},
             { 11.0,   1.0,  10.0,   0.0}}
        
        ' Initialize eigenvalue & eigenvector (initial guess)
        lambda = 1.0
        x = {0.0, 0.0, 0.0, 0.0}
        x(0) = 1.0
        maxit = 999
        eps = 0.001
        
        ' Apply Power Method
        Call POWER_METHOD_S(n, A, x, lambda, Aerr, eps, maxit)
        
        ' Print out the results
        If Aerr < eps Then
            'Console.WriteLine(StrDup(20, "="c))
	     Console.WriteLine("Largest lambda value = {0,12:F7} ",lambda)
            'Console.WriteLine("Lambda vector= " & String.Join(", ", x.Select(Function(num) num.ToString("F12.4"))))
            'Console.WriteLine(StrDup(20, "="c))
        Else
            Console.WriteLine("Max number of iterations is reached.")
        End If
    End Sub

	Public Sub POWER_METHOD_S(n As Integer, A(,) As Double, x() As Double, ByRef lambda As Double, ByRef err As Double, eps As Double, maxit As Integer)
	'  ==================================================================================
	'  CODE11.1-POWER_METHOD_S.bas. A Basic (VB) Sub implementing Pseudocode 11.1.                    
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
	'   DESCRIPTION: A module to find the dominant eigenvalue (largest in absolute value)          
	'      using the Power Method with scaling technique.                                          
	'                                                                                              	
	'   ON ENTRY                                                                                   
	'      n     :: Size of the matrix;                                                            
	'      A     :: A real square matrix (nxn);                                                    
	'      x     :: Array of length n containing the initial guess for the eigenvector;            
	'     lambda :: An initial guess for the dominant eigenvalue;                                  
	'      eps   :: Convergence tolerance;                                                         
	'     maxit  :: Maximum number of iterations permitted.                                        
	'                                                                                              
	'   ON EXIT                                                                                    
	'     lambda :: Estimated dominant (largest in absolute value) eigenvalue;                     
	'      x     :: Array of length n containing the estimated eigenvector;                        
	'      err   :: Error, max of both L2-norm of displacement vector and relative error for eigenvalue.
	'                                                                                              
	'   USES                                                                                       
	'     abs     :: Built-in Math library function returning the absolute value of a real value;
	'     AX      :: A module for computing Ax matrix-vector product (see Pseudocode 2.5);          
	'     MAX_SIZE:: A function module providing largest (in absolute) value of a vector.          
	'                                                                                              
	'   REVISION DATE :: 03/10/2025                                                                
	'  ==================================================================================
        Dim i, p As Integer
	Dim xn(n), d(n) As Double
        Dim lambdao, Err1, Err2 As Double
        
        err = MAX_SIZE(n, x)
        lambdao = lambda
        p = 0
        
        Do While (err > eps And p < maxit)
            p += 1 ' count iterations
            Call AX(n, A, x, xn) ' Solve [A][x]=[y]
            lambda = MAX_SIZE(n, xn) ' Find L-infinity of x^(p+1)
	
            For i = 0 To n - 1 ' Normalize x
                xn(i) = xn(i) / lambda
            Next

			'Console.WriteLine(String.Join(", ", xn))	
			For i = 0 To n - 1
				d(i) = Math.Abs( xn(i) - x(i) )
        	Next
		    Err1 = ENORM(n, d)
	        'Console.WriteLine(" Enorm= {0,12:F7} ", Err1)
            Err2 = Math.Abs(1.0 - lambdao / lambda)
            err = Math.Max(Err1, Err2)
			Console.WriteLine(" iter= {0,3} lambda= {1,12:F7} Error {2,13:E4}", p, lambda, err)
			'Console.WriteLine(" Eigenvector=")
			'Console.WriteLine(String.Join(", ", xn))
	        Console.Write(" Eigenvector =")
        	For i = 0 To n - 1
	            Console.Write("{0,12:F7} ", X(i))  
	        Next
			Console.WriteLine()
			Console.WriteLine()
            lambdao = lambda
            For i = 0 To n - 1 
                x(i) = xn(i)
            Next
        Loop
    End Sub

	Public Function MAX_SIZE(n As Integer, x() As Double) As Double
	'  ==================================================================================
	'  CODE11.1-POWER_METHOD_S.BAS of a module in Pseudocode 11.1.                    
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
	'   DESCRIPTION: A function module to find the largest element (in absolute value) of an array.
	'                                                                                              
	'   ARGUMENTS                                                                                  
	'      n   :: Length of the array;                                                          
	'      x   :: Array of length n.                                                            
	'                                                                                              
	'   USES                                                                                       
	'     abs  :: Built-in Math library function returning the absolute value of a real value;       
	'                                                                                              
	'   REVISION DATE :: 03/10/2025                                                                
	'  ==================================================================================
        Dim xmax As Double = x(0)
        For i = 1 To n - 1
            If Math.Abs(x(i)) > Math.Abs(xmax) Then
                xmax = x(i)
            End If
        Next
        Return xmax
    End Function

	Public Sub AX(n As Integer, A(,) As Double, x() As Double, b() As Double)
	'  ==================================================================================
	'  CODE2.5-Ax.bas of a module in Pseudocode 2.5.                                 
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
            Next
        Next
    End Sub

	Public Function ENORM(n As Integer, x() As Double) As Double
	'  ==================================================================================
	'  CODE3.1-JACOBI.BAS of a module in Pseudocode 3.1.                              
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
	'    SQRT :: Built-in Math library function returning the square root of a real value.            
	'                                                                                              
	'  REVISION DATE :: 11/09/2024                                                                  
	'  ==================================================================================

        Dim sum As Double = 0.0
        For i = 0 To n - 1
            sum += x(i) ^ 2
        Next
        Return Math.Sqrt(sum)
    End Function

End Module