' runs on    https://dotnetfiddle.net/
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System
				
Public Module Test_NonLinear_Newtons_Method

	Public Sub Main()
	' ==============================================================================
	'  The main program to test NONLINEAR_NEWTON.BAS
	' ==============================================================================
        Dim x(), y(), yo(), alpha(1), beta(1), gamma(1) As Double
        Dim xa, xb, h, eps, guess, err As Double
        Dim neq, n, i, maxit, bc_left, bc_rigt As Integer

        Console.WriteLine("Enter No of intervals ")
        n = Console.ReadLine()
        neq = n + 1
		
		ReDim x(neq), y(neq), yo(neq)
		
        ' Setup the ODE, Grids & BCs 
        eps = 1.0E-6
        guess = 0.4
        maxit = 99

        xa = 0.0 
		xb = 1.0
        h = (xb - xa) / n

        alpha(0) = 0.0 : beta(0) = 1.0 : gamma(0) = 1.0
        alpha(1) = 1.0 : beta(1) = 0.0 : gamma(1) = 0.0

        ' Prep initial guess for the solution
        bc_left = Math.Abs(CInt(alpha(0)))
        bc_rigt = Math.Abs(CInt(alpha(1)))

        For i = 0 To neq
            x(i) = xa + i * h
            yo(i) = guess
            If (bc_left = 0 AndAlso bc_rigt = 0) Then
                yo(i) = gamma(0) + (gamma(1) - gamma(0)) * (x(i) - xa) / (xb - xa)
            End If
        Next

        ' ***
        Call NONLINEAR_NEWTON(neq, eps, x, yo, y, alpha, beta, gamma, maxit)

        ' PRINT OUT THE RESULTS 
        Console.WriteLine()
	    Console.WriteLine("{0,8}  {1,12}  {2,12}  {3,12}", "x", "Exact", "N.Approx", "Abs Error")
        For i = 0 To neq - 1
            err = Math.Abs(y(i) - Exact(x(i)))
            Console.WriteLine("{0,10:F4}  {1,11:F7}  {2,11:F7}  {3,12:E3}", x(i), Exact(x(i)), y(i), err)
        Next
    End Sub

	Public Sub NONLINEAR_NEWTON(m As Integer, eps As Double, x() As Double, yo() As Double, y() As Double, alpha() As Double, beta() As Double, gamma() As Double, maxit As Integer)
	'  ==================================================================================
	'  CODE10.4-NONLINEAR_NEWTON.bas. A Basic (VB) Sub implementing Pseudocode 10.4.                  
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
	'   DESCRIPTION: A module to find approximate solution of a two-point nonlinear differential   
	'     equation using the Newton's Method. The nonlinear equations is cast in the following form
	'           y'' = f(x,y,y')  on [a,b]                                                          
	'     subject to                                                                               
	'           alpha1 * y'(a)+ beta1 * y(a) = gamma1                                              
	'           alpha2 * y'(b)+ beta2 * y(b) = gamma2                                              	
	'                                                                                              
	'   CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1. 
	'                                                                                              
	'   ON ENTRY                                                                                   
	'      M   :: Number of (equations) grid poinds;                                               
	'      x   :: Array of length M containing the abscissas of the grid points;                   
	'      yo  :: Array of length M containing the initial guess for the solution;                 
	'     alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                  
	'            the boundary conditions as stated above;                                          
	'     eps  :: Convergence tolerance;                                                           
	'     maxit:: Maximum number of iterations permitted.                                          
	'                                                                                              
	'   ON EXIT                                                                                    
	'      y   :: Array of length M containing the approximate solution.                           
	'                                                                                              
	'   USES                                                                                       
	'     FUNCS :: A user-defined external function module providing the coefficients of tridiagonal
	'             system of linear equations;                                                      
	'     ENORM:: A function module to calculate the Euclidean vector (L2 norm) of a vector;       
	'     TRIDIAGONAL :: A module to solve a tridiagonal system of equations with Thomas algorithm.
	'                                                                                              
	'   REVISION DATE :: 03/10/2025                                                                
	'  ==================================================================================
		Dim a(m-1), b(m-1), d(m-1), c(m-1) As Double
        Dim h, bb2h, hb2, c1d, c1r, c2d, c2r, err, yx, yxx, f, fy, fp, hsqr As Double
        Dim bc_left, bc_rigt, p, k, s1, sn As Integer

        bc_left = Math.Abs(CInt(alpha(0)))
        bc_rigt = Math.Abs(CInt(alpha(1)))

        h = x(1) - x(0) 
		hsqr = h * h
        bb2h = 0.5 / h : hb2 = 0.5 * h
        c1d = 2.0 * h * beta(0) / alpha(0)
        c1r = 2.0 * h * gamma(0) / alpha(0)
        c2d = 2.0 * h * beta(1) / alpha(1)
        c2r = 2.0 * h * gamma(1) / alpha(1)

        err = 1.0
        p = 0
        ' initialize numerical solution
        For i = 0 To m - 1
			y(i) = yo(i)
        Next
		Do
            For k = 0 To m - 1
                If k = 0 Then
                    If bc_left = 0 Then
                        d(k) = 1.0
                        yo(k) = gamma(0) / beta(0)
                        a(k) = 0.0
                        b(k) = 0.0
                        c(k) = 0.0
                    Else
                        yx = (gamma(0) - beta(0) * yo(k)) / alpha(0)
                        Call FUNCS(x(k), yo(k), yx, f, fy, fp)
                        a(k) = 2.0
                        b(k) = 0.0
                        d(k) = -2.0 + c1d - hsqr * (fy - fp * beta(0) * yo(k) / alpha(0))
                        c(k) = -c1r + (-2.0 + c1d) * yo(k) + 2.0 * yo(k + 1) - hsqr * f
                    End If
                ElseIf k = m - 1 Then
                    If bc_rigt = 0 Then
                        d(k) = 1.0
                        yo(k) = gamma(1) / beta(1)
                        a(k) = 0.0
                        b(k) = 0.0
                        c(k) = 0.0
                    Else
                        yx = (gamma(1) - beta(1) * yo(k)) / alpha(1)
                        Call FUNCS(x(k), yo(k), yx, f, fy, fp)
                        b(k) = 2.0
                        d(k) = -2.0 - c2d - hsqr * (fy - fp * beta(1) * yo(k) / alpha(1))
                        a(k) = 0.0
                        c(k) = c2r - (2.0 + c2d) * yo(k) + 2.0 * yo(k - 1) - hsqr * f
                    End If
                Else
                    yx = (yo(k + 1) - yo(k - 1)) / (2.0 * h)
                    yxx = yo(k + 1) - 2.0 * yo(k) + yo(k - 1)
                    Call FUNCS(x(k), yo(k), yx, f, fy, fp)
                    b(k) = 1.0 + hb2 * fp
                    d(k) = -2.0 - hsqr * fy
                    a(k) = 1.0 - hb2 * fp
                    c(k) = yxx - hsqr * f
                End If
            Next

            ' === SOLVE TRIDIAGONAL SYSTEM OF EQUATIONS FOR CORRECTION    
			s1 = 0: sn = m - 1
			Call TRIDIAGONAL(s1, sn, b, d, a, c, c)
            err = Enorm(m, c)
            For k = 0 To m - 1
                y(k) = yo(k) - c(k) ' Apply correction to prev estimates
            Next
            p = p + 1
            Console.WriteLine("Iter={0,3}  E-norm={1,14:E5}", p, err)
            For k = 0 To m - 1
                yo(k) = y(k)
            Next
            If err < eps OrElse p = maxit Then Exit Do
        Loop
    End Sub


	Public Sub TRIDIAGONAL(ByVal s1 As Integer, ByVal sn As Integer, ByRef b() As Double, ByRef d() As Double, ByRef a() As Double, ByRef c() As Double, ByRef x() As Double)
	'  ==================================================================================
	'  CODE2.13-TRIDIAGONAL.bas. A Basic (VB) Sub implementing Pseudocode 2.13.                       
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
	'  DESCRIPTION: A sub program to solve a tridiagonal system of linear equations                 
	'    using Thomas algorithm.                                                                   
	'                                                                                              
	'  ON ENTRY                                                                                    
	'     s1 :: Subscript of the first unknown (usually 1);                                        
	'     sn :: Subscript of the last unknown (usually no. of eqs, n);                             
	'      b :: Array of length n containing coefficients of below diagonal elements;              
	'      d :: Array of length n containing coefficients of diagonal elements;                    
	'      a :: Array of length n containing coefficients of above diagonal elements;              
	'      c :: Array of length n containing coefficients of rhs.                                  
	'                                                                                              
	'  ON EXIT                                                                                     
	'      x :: An array of length n containing the solution.                                      
	'                                                                                              
	'  REVISION DATE :: 03/18/2024                                                                 
	'  ==================================================================================
        Dim i As Integer
        Dim ratio As Double
		
		For i = s1 + 1 To sn
            ratio = b(i) / d(i - 1)
            d(i) = d(i) - ratio * a(i - 1)
            c(i) = c(i) - ratio * c(i - 1)
        Next

        x(sn) = c(sn) / d(sn)
        For i = sn - 1 To s1 Step -1
            x(i) = (c(i) - a(i) * x(i + 1)) / d(i)
        Next
    End Sub

	Public Function Enorm(n As Integer, x() As Double) As Double
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
	'    sqrt :: Built-in Intrinsic function returning the square root of a real value.            
	'                                                                                              
	'  REVISION DATE :: 11/09/2024                                                                  
	'  ==================================================================================
        Dim sum As Double = 0.0
        For i = 0 To n - 1
            sum += x(i) ^ 2
        Next
        Return Math.Sqrt(sum)
    End Function

	Public Sub FUNCS(ByVal x As Double, ByVal y As Double, ByVal yp As Double, ByRef fun As Double, ByRef dfdy As Double, ByRef dfdp As Double) 
	' ==============================================================================
	' DESCRIPTION: A user-defined function supplying the coefficients of the nonlinear 
	'    two-point BVP given in the form: 
	'              y'' = f(x,y,y')  on [a,b]
	'
	' ON ENTRY
	'   x    :: Independent variable (a <=x<= b);
	'   y    :: Dependent variable y=y(x).
	'
	' ON EXIT
	'   f    :: The nonlinear two-point BVP f(x,y,y') evaluated at (x,y);
	'   yp   :: First derivative of the dependent variable y';
	'   dfdy :: Partial derivative of f wrt y evaluated at (x,y), df/dy;
	'   dfdp :: Partial derivative of f wrt y' evaluated at (x,y), df/dy'.
	'
	' REVISION DATE :: 03/10/2025
	' ==============================================================================
		fun = y - yp ^2 / y
		dfdy = 1.0 + ( yp / y ) ^ 2
		dfdp = -2.0 * yp / y
	End Sub

	Public Function Exact(x As Double) As Double
    	' -------------------------------------------------------------------------------
	' DESCRIPTION: A function subprogram providing the true solution y=f(x)
	' -------------------------------------------------------------------------------
    	Exact = Math.Cosh(Math.Sqrt(2.0) * (1.0 - x))
    	Exact = Math.Sqrt(Exact / Math.Cosh(Math.Sqrt(2.0)))
    End Function
End Module