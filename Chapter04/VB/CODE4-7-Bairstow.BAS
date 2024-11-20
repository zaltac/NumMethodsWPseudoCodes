' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

' ==============================================================================
'  The main program to test Sub BAIRSTOW
' ==============================================================================
Public Module Test_Bairstow
	Public Sub Main()
        Dim A(0 To 19) As Double
        Dim XRE(19), XIM(19) As Double
		Dim n, i, maxit, iprnt As Integer
		Dim p0, q0, eps As Double

        n = 5
		A = {1.0, -5.0, -15.0, 85.0, -26.0, -120.0}
		
        iprnt = 2
        maxit = 99
        p0 = 0.0
        q0 = 0.0
		eps = 0.5e-4

        Call Bairstow(n, p0, q0, A, eps, maxit, XRE, XIM)

        Console.WriteLine("    ======== All the Roots are =========")
        For i = 0 To n - 1
			Console.WriteLine(" {0} {1,12:F7} {2,12:F7}", i, XRE(i), XIM(i))
        Next
        Console.WriteLine("    ====================================")
    End Sub



'  ==================================================================================
'  CODE4.7-BAIRSTOW.BAS. A Basic (VB) Sub implementing Pseudocode 4.7.                            
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
'  DESCRIPTION: A VB subprogram to find all real and/or imaginary roots of a polynomial           
'    of the n'th degree using the BAIRSTOW's method.                                         
'                                                                                              
'  ON ENTRY                                                                                    
'    n    :: Degree of the polynomial;                                                         
'   p0,q0 :: Initial guesses for a quadratic equation; i.e., for p and q;                      
'     a   :: Array of length (n+1) containing the coefficients of polynomial defined as        
'                 a0 x^n + a1 x^(n-1) + ... + an = 0                                           
'    eps  :: Convergence tolerance;                                                            
'   maxit :: Maximum number of iterations permitted;                                           
'   iprnt :: printing key, =0 do not print intermediate results, <> 0 print intermediates.     
'                                                                                              
'  ON EXIT                                                                                     
'    xre  :: Array of length n containing real parts of the roots;                             
'    xim  :: Array of length n containing imaginary parts of the roots.                        
'                                                                                              
'  OTHER VARIABLES                                                                             
'     b   :: Array of length [n] containing coefficients of quotient polynomial (0<=k<=n-2);   
'     c   :: Array of length [n] containing coefficients of partial derivatives.               
'                                                                                              
'  USES                                                                                        
'    MATH.ABS  :: Built-in Intrinsic function returning the absolute value of a real value;         
'    QUADRATIC :: Subroutine that solves a quadratic equation of the form x2 + p x + q = 0 (see CODE1-3).    
'                                                                                              
'  REVISION DATE :: 04/29/2024                                                                 
'  ==================================================================================
	Public Sub Bairstow(ByVal n As Integer, ByRef p0 As Double, ByRef q0 As Double, ByRef A() As Double, _
		ByVal eps As Double, ByVal maxit As Integer, ByRef XRE() As Double, ByRef XIM() As Double)
        Dim B(0 To n) As Double, C(0 To n) As Double, XR(2) As Double, XI(2) As Double
        Dim p As Double, delp As Double, q As Double, delq As Double, delM As Double, cbar As Double, del As Double, del1 As Double, del2 As Double
        Dim i As Integer, k As Integer, m As Integer, kount As Integer

        For k = n To 0 Step -1
            A(k) = A(k) / A(0)
        Next

        m = n
        kount = -1
        While n > 1
            p = p0 : q = q0
            k = 0
            delM = 1.0R
			Console.WriteLine("{0} {1,12} {2,12} {3,12} ","Iter","delM", "p","q")
            While delM > eps AndAlso k <= maxit
                k += 1
                B(0) = 1.0R : C(0) = 1.0R
                B(1) = A(1) - p : C(1) = B(1) - p
                For i = 2 To n
                    B(i) = A(i) - p * B(i - 1) - q * B(i - 2)
                    C(i) = B(i) - p * C(i - 1) - q * C(i - 2)
                Next
                cbar = C(n - 1) - B(n - 1)
                del = C(n - 2) * C(n - 2) - cbar * C(n - 3)
                del1 = B(n - 1) * C(n - 2) - B(n) * C(n - 3)
                del2 = B(n) * C(n - 2) - B(n - 1) * cbar
                delp = del1 / del : delq = del2 / del
                p = p + delp : q = q + delq
                delM = Math.Abs(delp) + Math.Abs(delq)
				Console.WriteLine("{0,3:F0}  {1,14:E3}  {2,12:F7} {3,12:F7}",k, delM, p, q)		
            End While
            If k - 1 = maxit Then
                Console.WriteLine("Quadratic factor did not converge after {0} iterations", k - 1)
                Console.WriteLine("Recent values of p, q, delM are {0}, {1}, {2}", p, q, delM)
                Console.WriteLine("Corresponding roots may be questionable ...")
            End If
            Call QuadraticEq(p, q, XR, XI)
            kount += 1
            XRE(kount) = XR(0) : XIM(kount) = XI(0)
            kount += 1
            XRE(kount) = XR(1) : XIM(kount) = XI(1)
			Console.WriteLine("======== Found a Quadratic Factor ========")
			Console.WriteLine("  x*x + ({0,10:F6}) * x + ({1,10:F6})", p, q)
            Console.WriteLine("==========================================")
			Console.WriteLine(" ")
            n = n - 2
            For i = 0 To n
                A(i) = B(i)
            Next
            If n = 1 Then
                kount += 1
                XRE(kount) = -A(1)
                XIM(kount) = 0.0
				Console.WriteLine("======== Found a Linear Factor ========")
                Console.WriteLine(" ( x + ({0,10:F6})", A(1) )
                Console.WriteLine("========================================")
            End If
        End While
        n = m
    End Sub

 
End Module
