' runs on    https://dotnetfiddle.net/8K7fkp
Imports System


' ==============================================================================
'  Main program to test Sub Program NEWTON_RAPHSON 
' ==============================================================================
Public Module test_NewtonRaphson
	Public Sub Main()
        Dim maxit As Integer, iter As Integer
        Dim eps As Double, root As Double
        
        maxit = 99
        eps = 0.50E-4
		
		Console.Write("Enter an estimate for root ")
        root = Console.ReadLine()
        
        Call NEWTON_RAPHSON(root, maxit, eps, iter)
		Console.WriteLine(New String("-", 53))
		Console.WriteLine("Root is {0,12:F7} converged after {1} iterations",root,iter)
		Console.WriteLine(New String("-", 53))
        Console.ReadLine()
    End Sub


'  ==================================================================================
'  CODE4.3-NEWTON_RAPHSON.BAS. A Basic (VB) Subprogram implementing Pseudocode 4.3.                      
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
'  DESCRIPTION: A Visual Basic module to compute a root of a nonlinear equation using the                   
'    Newton-Raphson method.                                                                    
'                                                                                              
'  ON ENTRY                                                                                    
'    root  :: Initial guess for the root;                                                      
'    maxit :: Maximum number of iterations permitted;                                          
'    eps   :: Convergence tolerance.                                                           
'                                                                                              
'  ON EXIT                                                                                     
'    iter  :: Number of iterations realized;                                                   
'    root  :: Computed approximation for the root.                                             
'                                                                                              
'  USES                                                                                        
'    ABS   :: Built-in Intrinsic function returning the absolute value of a real value;        
'                                                                                              
'  ALSO REQUIRED                                                                               
'    FUNC  :: User-defined external function providing the nonlinear equation, f(x).           
'    FUNCP :: User-defined external function providing the first derivative                    
'             of the nonlinear equation, f'(x).                                                
'                                                                                              
'  REVISION DATE :: 11/20/2024                                                                 
'  ==================================================================================    
	Public Sub NEWTON_RAPHSON(ByRef root As Double, ByVal maxit As Integer, ByVal eps As Double, ByRef iter As Integer)
        Dim p As Integer
		Dim fn, fpn, aerr, rate, x0, xn, del, del0 As Double

        del0 = 1.0
        x0 = root
        p = 0
		
		Console.WriteLine(New String("-", 76))
		Console.WriteLine("{0,3} {1,12} {2,17} {3,12} {4,10} {5,12}", "p", "x^(p)", "f(x^(p))","f'(x^(p))", "error", "C.Rate")
        Console.WriteLine(New String("-", 76))        
        Do
			fn = FUNC(x0)
			fpn = FUNCP(x0)
            del = -fn / fpn
            aerr = Math.Abs(del)
            rate = aerr / del0 ^ 2
			Console.WriteLine("{0,3:F0} {1,14:F7} {2,14:F7} {3,14:E4} {4,12:E4} {5,12:F7}", p, x0, fn, fpn, aerr,rate)
            xn = x0 + del
            x0 = xn
            del0 = Math.Abs(del)
            p = p + 1
            If ((Math.Abs(fn) < eps And aerr < eps) Or (p = maxit)) Then
                Exit Do
            End If        
        Loop
        
        root = xn
        iter = p
        Console.WriteLine(" p=",p,"  iter ",iter)
        
        If p = maxit Then
            Console.WriteLine("** Max iteration number reached=")
            Console.WriteLine("** Estimated root has NOT converged, del, f(x)= " & Math.Abs(del) & ", " & func(xn))
        End If
    End Sub
    
' ==========================================================================          
' User-defined function providing f(x), which should be cast as FUNC(x)=0.
' ========================================================================== 
		Public Function FUNC(x As Double) As Double
        Return 4.0 + x * x * (8.0 - x * x)
    End Function

' ==========================================================================          
' User-defined function providing f'(x)=FUNCP(X) 
' ==========================================================================     
		Public Function FUNCP(x As Double) As Double
        Return 4.0 * x * (4.0 - x * x)
    End Function

End Module