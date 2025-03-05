 ' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================

Imports System

Public Module Test_Runge_Kuttas
' ==============================================================================
'  The main program to test ERunge_Kutta.BAS
' ==============================================================================
	Public Sub Main()
        Dim x0 As Double = 0.0
        Dim y0 As Double = 2.0
        Dim xlast As Double = 1.0
        Dim h As Double = 0.1
        Dim order As Integer = 4
	Dim maxit As Integer = 99
	Dim eps As Double = 0.5E-6	
	Dim n As Integer
		
	Console.WriteLine("Enter order of RK method")
        n = Console.ReadLine()
		
        Call Runge_Kutta(n, h, x0, y0, xlast)
    End Sub

	Public Sub Runge_Kutta(ByVal n As Integer, ByVal h As Double ,ByRef x0 As Double, ByRef y0 As Double, ByVal xlast As Double)
	'  ==================================================================================
	'  CODE9.4-Runge_Kutta.bas. A Basic (VB) Sub implementing Pseudocode 9.4.                         
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
	'  DESCRIPTION: A Visual Basic sub program to estimate the solution of a first order IVP on [x0,xlast]       
	'    using the 2nd to 4th order Runge-Kutta scheme. Numerical estimates are printed out, not stored.           
	'                                                                                              
	'  ON ENTRY                                                                                    
	'    n     :: Order of the Runge-Kutta scheme;                                                 
	'    h     :: Step size (it must be uniform);                                                  
	'    x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps	
	'   xlast  :: End point of the solution interval.                                              
	'                                                                                              
	'  Other Internal Variables                                                                    
	'    x,y   :: Current estimates, x^(p+1) and y^(p+1).                                          
	'                                                                                              
	'   USES                                                                                       
	'     Abs  :: Built-in Math libarary function returning the absolute value of a real value.        
	'    DRV_RK:: A driver subprogram performing one-step RK scheme.                               
	'                                                                                              
	'  REVISION DATE :: 03/05/2025                                                                 
	'  ==================================================================================
    	Dim yt, aerr, x, y As Double
    
	Console.WriteLine("{0,5} {1,10} {2,18}", "x", "y", "Abs. Error")
        Console.WriteLine(New String("-"c, 42))
	Console.WriteLine("{0,8:F5} {1,12:F8}", x0, y0)
    	x = x0
    
    	Do While x < xlast
        	Call DRV_RK(n, h, x0, y0, x, y)
        	yt = Exact(x)
        	aerr = Math.Abs(y - yt)
		Console.WriteLine("{0,8:F5} {1,12:F8} {2,14:E4}", x, y, aerr)			
        	x0 = x
        	y0 = y
    	Loop
    
	End Sub

	Public Sub DRV_RK(ByVal n As Integer, ByVal h As Double, ByRef x0 As Double, ByRef y0 As Double, ByRef x As Double, ByRef y As Double)
	'  ==================================================================================
	'  CODE9.4-DRV_RK.bas. A Basic (VB) Sub implementing Pseudocode 9.4.                              
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
	'  DESCRIPTION: A driver sub-program employing one-step RK2, RK3, or RK4 scheme.               
	'                                                                                              
	'  ON ENTRY                                                                                    
	'   n     :: Order of Runge-Kutta scheme;                                                      
	'   h     :: Step size (it must be uniform);                                                   
	'   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
	'                                                                                              
	'  ON EXIT                                                                                     
	'   x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
	'                                                                                              
	'  USES                                                                                               	
	'    FCN  :: User-defined external function providing y'=f(x,y).                               
	'                                                                                              
	'  REVISION DATE :: 03/05/2025                                                                 
	'  ==================================================================================
    	Dim x1, xh, ym, xk, xk1, xk2, xk3, xk4, hlf As Double
    
    	hlf = 0.5
    	xh = x0 + 0.5 * h
    	x1 = x0 + h
    
    	Select Case n
        	Case 2 ' Case of RK2
            	xk1 = h * FCN(x0, y0)
            	ym = y0 + xk1
            	xk2 = h * FCN(x1, ym)
            	xk = hlf * (xk1 + xk2)
        
        	Case 3 ' Case of RK3
            	xk1 = h * FCN(x0, y0)
            	ym = y0 + hlf * xk1
            	xk2 = h * FCN(xh, ym)
            	ym = y0 - xk1 + 2.0 * xk2
            	xk3 = h * FCN(x1, ym)
            	xk = (xk1 + 4.0 * xk2 + xk3) / 6.0
        
        	Case 4 ' Case of RK4
            	xk1 = h * FCN(x0, y0)
            	ym = y0 + hlf * xk1
            	xk2 = h * FCN(xh, ym)
            	ym = y0 + hlf * xk2
            	xk3 = h * FCN(xh, ym)
            	ym = y0 + xk3
            	xk4 = h * FCN(x1, ym)
            	xk = (xk1 + 2.0 * xk2 + 2.0 * xk3 + xk4) / 6.0
        
        	Case Is <= 1 ' Case of n<=1
            	Console.WriteLine(" PROGRAM DOES NOT HANDLE CASE OF N=" & n)
            	Environment.Exit(0)
        
	        Case Is >= 5 ' Case of n>=5
    	        Console.WriteLine(" PROGRAM DOES NOT HANDLE CASE OF N=" & n)
        	    Environment.Exit(0)
    
    	End Select
    
    	y = y0 + xk
    	x = x1
    
	End Sub

	Public Function FCN(x As Double, y As Double) As Double
	' ==============================================================================
	' DESCRIPTION: A function subprogram providing y'=f(x,y)
	'
	' ARGUMENTS:
	'      x, y  :: Real input values.
	'
	' ==============================================================================
		Return -y / ( x + 1.0 )
    End Function

	Public Function Exact(x As Double) As Double
	' ==============================================================================
	' DESCRIPTION: A function subprogram providing the true solution y=(x) for 
	'    testing the module. 
	'
	' ARGUMENTS:
	'      x   :: A real input, independent variable.
	'
	' ==============================================================================
        Return 2.0 / ( x + 1.0 )
    End Function
End Module