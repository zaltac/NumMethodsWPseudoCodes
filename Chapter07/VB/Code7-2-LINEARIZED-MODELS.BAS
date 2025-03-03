' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System



Public Module TEST_LINEARIZE_REGRESS
' ==============================================================================
'  The main program to test LINEARIZE_REGRESS.C
' ==============================================================================
	Public Sub Main()
        Dim n As Integer = 6
        Dim x(n - 1), y(n - 1) As Double
		Dim a0, b0, E, S, r2 As Double
        Dim model As Integer

        x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}
        y = {5.3, 6.5, 6.8, 7.1, 7.8, 7.5}

        model = 1
        LINEARIZE_REGRESS(n, x, y, model, a0, b0, E, S, r2)

	Console.WriteLine(" ")
        Console.WriteLine("******* Best-Fit Coefficients and Parameters ********")
        Console.Write("  a = {0,12:F8}", a0)
        Console.WriteLine("  b = {0,12:F8}", b0)
		Console.Write("  E = {0,12:F8}", E)
        Console.WriteLine("  S = {0,12:F8}", S)
        Console.WriteLine("  r-squared = {0,7:F6} ",r2)
    End Sub


'  ==================================================================================
'  CODE7.2-LINEARIZE_REGRESS.BAS. A Basic (VB) Sub implementing Pseudocode 7.2.                   
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
'  DESCRIPTION: A module to obtain least-squares-fit to the 'power model' only. Users                 
'    can likewise incorporate the other linearizable models into the module.                   
'                                                                                              
'  ON ENTRY                                                                                    
'     n   :: The number of data in the set;                                                    
'    x,y  :: Arrays of length n containing the data;                                           
'   model :: Model flag, Model = 1 corresponds to the power model, Y= a0*x^b0.                 
'                                                                                              
'  ON EXIT                                                                                     
'    a0,b0:: Model parameters;                                                                 
'     E   :: Sum of the Squares of Residuals (SSR);                                            
'     S   :: Sum of the Squares of Mean Deviation (SSMD);                                      
'    r2   :: r-squared, coefficient of determination.                                          
'                                                                                              
'  USES                                                                                        
'     Exp :: Built-in Intrinsic Math.function returning exponential of a real value, e^x.           
'     Log :: Built-in Intrinsic Math.function returning the natural log of a real value.            
'                                                                                              
'  REVISION DATE :: 03/03/2025                                                                 
'  ==================================================================================
 
	Public Sub LINEARIZE_REGRESS(n As Integer, ByRef x() As Double, ByRef y() As Double, model As Integer, ByRef a0 As Double, ByRef b0 As Double, ByRef E As Double, ByRef S As Double, ByRef r2 As Double)
        Dim c(1, 1) As Double
        Dim b(1) As Double
        Dim xx(n - 1), yy(n - 1) As Double
		Dim DD, D1, D2, yavg, yi, xk As Double
		Dim i, j, k, p As Integer

		For k = 0 To n - 1
            xx(k) = x(k)
            yy(k) = y(k)
	    Console.WriteLine(" {0}  x({0})= {1,4:F2}  y({0})= {2,4:F2} ", k, xx(k), yy(k))
            If model = 1 Then
                x(k) = Math.Log(x(k))
                y(k) = Math.Log(y(k))
            ElseIf model < 1 Or model > 1 Then
                Console.WriteLine("Undefined model...")
                Exit For
            End If
        Next

        For i = 0 To 1
            b(i) = 0.0
            For j = 0 To 1
                c(i, j) = 0.0
                For k = 0 To n - 1
                    p = i + j
                    xk = 1.0
                    If p <> 0 Then xk = x(k)^p
                    c(i, j) = c(i, j) + xk
                Next
            Next
            For k = 0 To n - 1
                p = i
                xk = 1.0
                If p <> 0 Then xk = x(k)^p 
                b(i) = b(i) + xk * y(k)
            Next
        Next

        DD = c(0,0) * c(1, 1) - c(1, 0) * c(0, 1)
        D1 = b(0) * c(1, 1) - b(1) * c(0, 1)
        D2 = b(1) * c(0, 0) - b(0) * c(1, 0)

        a0 = D1 / DD
        b0 = D2 / DD

        yavg = 0.0
        For k = 0 To n - 1
            yavg = yavg + yy(k)
        Next
        yavg = yavg / CDbl(n)

        S = 0.0
        E = 0.0

        If model <= 2 Then
            a0 = Math.Exp(a0)
        Else
            Console.WriteLine("Undefined model...")
            Return
        End If

        For k = 0 To n - 1
            S = S + (yy(k) - yavg) ^ 2
            yi = a0 * xx(k) ^ b0
            E = E + (yi - yy(k)) ^ 2
        Next

        r2 = 1.0 - E / S
    End Sub
End Module