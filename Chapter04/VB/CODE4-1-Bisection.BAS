' runs on    https://dotnetfiddle.net/8K7fkp
 
Imports System


' ==================================================================================
' Main program to test BISECTION.BAS
' ==================================================================================
Public Module Test_Bisection
	Public Sub Main()
        Dim Halves As Integer, maxit As Integer
        Dim a As Double, b As Double, eps As Double, root As Double, fa As Double, fb As Double

        maxit = 99
        eps = 0.0005
        a = 0.0
        b = 4.0

        fa = FUNC(a)
        fb = FUNC(b)
        If fa * fb > 0 Then
            Console.WriteLine("No root in interval (a,b). Change the interval.")
            Return
        End If

        Call Bisection(a, b, maxit, eps, root, Halves)

        Console.WriteLine(New String("=", 44))
		Console.WriteLine("! Root is {0,14:F8} after {1} bisections !",root, Halves)
        Console.WriteLine(New String("=", 44))

        root = (a * fb - b * fa) / (fb - fa)
        Console.WriteLine()
		Console.WriteLine(" *** Root (with linear interpolation)   = {0,14:F7}", root)
		Console.WriteLine(" *** Clossness to the root, Abs[f(root)]= {0,14:E4}", Math.Abs(FUNC(root)))
        Console.WriteLine()
    End Sub

	Public Sub Bisection(ByRef a As Double, ByRef b As Double, ByVal maxit As Integer, ByVal eps As Double, ByRef root As Double, ByRef Halves As Integer)
	'  ==================================================================================
	'  CODE4.1-BISECTION.BAS. A Basic (VB) Sub implementing Pseudocode 4.1.                           
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
	'  DESCRIPTION: A Visual Basic module to find a root of a nonlinear equation in [a,b]                       
	'    using the Bisection method.                                                               
	'                                                                                              
	'  ON ENTRY                                                                                    
	'   [a,b] :: Initial search interval (it must bracket one root);                               
	'   maxit :: Maximum number of iterations permitted;                                           
	'   eps   :: Convergence tolerance.                                                            
	'                                                                                              
	'  ON EXIT                                                                                     
	'   halves:: Number of halves realized;                                                        
	'   root  :: Computed approximation for the root.                                              
	'                                                                                              
	'  USES                                                                                        
	'   ABS   :: Built-in Intrinsic function returning the absolute value of a real value;         
	'                                                                                              
	'  ALSO REQUIRED                                                                               
	'   FUNC  :: User-defined external function providing the nonlinear equation.                  
	'                                                                                              
	'  REVISION DATE :: 11/20/2024                                                                 
	'  ==================================================================================
        Dim p As Integer, interval As Double, fa As Double, fb As Double, xm As Double, fm As Double

        p = 0
        interval = b - a
        fa = FUNC(a)
        fb = FUNC(b)

        Console.WriteLine(New String("-", 97))
        Console.WriteLine("{0,3} {1,10} {2,12} {3,12} {4,12} {5,12} {6,12} {7,12}", "p", "a", "b", "f(a)", "f(b)", "xm", "f(xm)", "  interval size")
        Console.WriteLine(New String("-", 97))

        Do
            p += 1 ' **************** Start REPEAT-UNTIL loop *******************
            xm = 0.50 * (a + b)
            fm = FUNC(xm)
			Console.WriteLine("{0,3:F0} {1,12:F7} {2,12:F7} {3,12:E3} {4,12:E3} {5,12:F7} {6,12:E3} {7,12:E4}", p, a, b, fa, fb, xm, fm, interval)
            If fa * fm > 0.0 Then
                a = xm
                fa = fm
            Else ' case of fa*fm<0.0D0
                b = xm
                fb = fm
            End If
            interval = 0.50 * interval
            If ((Math.Abs(fm) < eps) And (interval < eps)) Or p = maxit Then Exit Do
        Loop ' End the iteration loop

        root = xm
        Halves = p
        If p = maxit Then
            Console.WriteLine(New String("!", 37))
            Console.WriteLine("Max iteration number reached={0}", p)
            Console.WriteLine(New String("!", 37))
        End If
    End Sub

	Public Function FUNC(ByVal x As Double) As Double
        ' =======================================================================          
        ' User-defined Function providing f(x) which should be cast as f(x)=0.
        ' =======================================================================  
        Return x * x + 0.025 * x - 4.0
    End Function
End Module