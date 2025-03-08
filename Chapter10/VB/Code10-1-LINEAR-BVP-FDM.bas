' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_LBVP_SOLVE
	Public Sub Main()
	' ==============================================================================
	'  The main program to test LBVP_SOLVE.BAS
	' ==============================================================================
        Dim nmax As Integer = 201
        Dim x(nmax) As Double, y(nmax) As Double, alpha(2) As Double, beta(2) As Double, gamma(2) As Double
        Dim err As Double, xa As Double, xb As Double, h As Double
        Dim neq As Integer, n As Integer, i As Integer

        Console.WriteLine("Enter No of intervals ")
        n = Console.ReadLine()
        If n > nmax Then
            Console.WriteLine("Max grid size is {0}", nmax)
            Console.WriteLine("Increase NMAX, and then try it again.")
            Return
        End If

        ' Setup the ODE, Grids & BCs       
        neq = n + 1
        xa = 1.0 : xb = 2.0
        h = (xb - xa) / CDbl(n)
        For i = 1 To neq
            x(i) = xa + CDbl(i - 1) * h
        Next
        alpha(1) = 1.0 : beta(1) = 2.0 : gamma(1) = 4.0
        alpha(2) = 1.0 : beta(2) = -2.0 : gamma(2) = 2.0

        Call LBVP_SOLVE(neq, x, alpha, beta, gamma, y)

        Console.WriteLine()
        Console.WriteLine("   x        Exact      N.Approx     Abs Error")
        For i = 1 To neq
            err = Math.Abs(y(i) - Exact(x(i)))
            Console.WriteLine("{0,6:F3}  {1,11:F7}  {2,11:F7}  {3,12:E4}", x(i), Exact(x(i)), y(i), err)
        Next
    End Sub

	Public Sub LBVP_SOLVE(neq As Integer, x() As Double, alpha() As Double, beta() As Double, gamma() As Double, y() As Double)
	'  ==================================================================================
	'  CODE10.1-LBVP_SOLVE.bas. A Basic (VB) Sub implementing Pseudocode 10.1.                        
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
	'  DESCRIPTION: A Sub program to find approximate solution of the following linear                  
	'    differential equation using the Finite Difference Method:                                 
	'          p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]                                   
	'    subject to                                                                                
	'          alpha1 * y'(a)+ beta1 * y(a) = gamma1                                               
	'          alpha2 * y'(b)+ beta2 * y(b) = gamma2                                               
	'                                                                                              
	'  CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1.  
	'                                                                                              
	'  ON ENTRY                                                                                    
	'    neq  :: Number of (equations) grid poinds;                                                
	'     x   :: Array of length neq containing the abscissa of the grid points;                   
	'    alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                   
	'           the boundary conditions as stated above;                                           
	'                                                                                              
	'  ON EXIT                                                                                     
	'     y   :: Array of length neq containing the approximate solution.                          
	'                                                                                              
	'  USES                                                                                        
	'   COEFFS  :: A Sub program containing the coefficients and rhs of the linear ODE, i.e., p(x), q(x), r(x), g(x).
	'   TRIDIAGONAL:: A Sub program solving a tridiagonal system of equations using the Thomas algorithm
	'                                                                                              
	'  REVISION DATE :: 03/08/2025                                                                 
	'  ==================================================================================
        Dim b(neq) As Double, a(neq) As Double, d(neq) As Double, cy(neq) As Double
        Dim h As Double, px As Double, qx As Double, rx As Double, gx As Double
        Dim BC_Left As Integer, BC_Rigt As Integer, k As Integer, s1 As Integer

        BC_Left = CInt(alpha(1))
        BC_Rigt = CInt(alpha(2))
        h = x(2) - x(1)

        ' Initialize the arrays of the tridiagonal system
        For k = 1 To neq
            a(k) = 0.0: b(k) = 0.0: d(k) = 0.0: cy(k) = 0.0
        Next
        
        For k = 1 To neq
            Call Coeffs(x(k), px, qx, rx, gx)
            d(k) = rx - 2.0 * px / (h * h)
            a(k) = (px / h + 0.5 * qx) / h
            b(k) = (px / h - 0.5 * qx) / h
            cy(k) = gx
        Next

        If BC_Left = 0 Then
            d(1) = 1.0
            a(1) = 0.0 : b(1) = 0.0
            cy(1) = gamma(1) / beta(1)
        Else
            a(1) = a(1) + b(1)
            d(1) = d(1) + 2.0 * h * b(1) * beta(1) / alpha(1)
            cy(1) = cy(1) + 2.0 * h * b(1) * gamma(1) / alpha(1)
        End If

        If BC_Rigt = 0 Then
            d(neq) = 1.0
            a(neq) = 0.0 : b(neq) = 0.0
            cy(neq) = gamma(2) / beta(2)
        Else
            b(neq) = a(neq) + b(neq)
            d(neq) = d(neq) - 2.0 * h * a(neq) * beta(2) / alpha(2)
            cy(neq) = cy(neq) - 2.0 * h * a(neq) * gamma(2) / alpha(2)
        End If

        s1 = 1
        Call TriDiagonal(s1, neq, b, d, a, cy, y)
    End Sub

	Public Sub Coeffs(x As Double, ByRef p As Double, ByRef q As Double, ByRef r As Double, ByRef g As Double)
	' ==============================================================================
	' DESCRIPTION: A user-defined sub program suppling the coefficients p(x), q(x), r(x) and
	'    rhs g(x) of the linear ODE given in the following form: 
	'         p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]
	'
	' ON ENTRY
	'    x   :: Independent variable (a<= x<= b);
	'
	' ON EXIT
	'   p, q, r, abd g :: Coefficients & rhs evaluated at x. 
	'
	' REVISION DATE :: 03/08/2025
	' ==============================================================================
        p = x * x
        q = -5.0 * x
        r = 8.0
        g = 0.0
    End Sub

	Public Function Exact(x As Double) As Double
	' ==============================================================================
	' DESCRIPTION: A user-defined function providing the true solution y=(x) for testing the module. 
	'
	' ARGUMENTS:
	'      x   :: A real input, independent variable.
	' ==============================================================================
        Exact = x * x * (x * x - 0.5)
    End Function

	Public Sub TriDiagonal(s1 As Integer, sn As Integer, b() As Double, d() As Double, a() As Double, c() As Double, x() As Double)
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
End Module