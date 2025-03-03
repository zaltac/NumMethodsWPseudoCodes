' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System



' ==============================================================================
'  The main program to test the sub program Trapezoidal_Rule_RF
' ==============================================================================
Public Module Test_TrapezoidalRule
    Sub Main()
        Dim a As Double = 0.0
        Dim b As Double = 1.0
        Dim n As Integer
        Dim intg As Double
        Dim intgc As Double

        Console.WriteLine("Enter Number of Panels ")
        n = Convert.ToInt32(Console.ReadLine())

        Trapezoidal_Rule_RF(n, a, b, intg, intgc)

        Console.WriteLine(" Panel no. = {0}", n)
        Console.WriteLine(" Estimate = {0:F14} Trapezoidal rule", intg)
        Console.WriteLine(" Estimate = {0:F14} Trapezoidal rule with end-correction", intgc)
    End Sub



Public Sub Trapezoidal_Rule_RF(n As Integer, a As Double, b As Double, ByRef intg As Double, ByRef intgc As Double)
'  ==================================================================================
'  CODE8.1-Trapezoidal_Rule_RF.bas. A Basic (VB) Sub implementing Pseudocode 8.1.                 
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
'  DESCRIPTION: A subroutine to estimate the integral of y=f(x) on [a,b]                       
'   using the Trapezoidal rule with/without end correction.                                    
'                                                                                              
'  ON ENTRY                                                                                    
'     n   :: Number of panels (i.e., n+1 integration points);                                  
'   [a,b] :: Integration interval.                                                             
'                                                                                              
'  ON EXIT                                                                                     
'   intg   :: Integral estimate using the ordinary Trapezoidal rule;                           
'   intgc  :: Integral estimate using the Trapezoidal rule with the end-point correction.      
'                                                                                              
'  USES                                                                                
'     FX   :: User-defined external function providing the function, f(x);                     
'     FU   :: User-defined external function providing the first derivative, f'(x).                                      '                                                                                              
'  REVISION DATE :: 03/03/2025                                                                 
'  ==================================================================================
        Dim h As Double = (b - a) / n
        intg = 0.5 * (FX(a) + FX(b))
        Dim xi As Double = a

        For i As Integer = 1 To n - 1
            xi += h
            intg += FX(xi)
        Next

        intg *= h

        Dim corr As Double = -h * h * (FU(b) - FU(a)) / 12.0
        intgc = intg + corr
    End Sub


' ==============================================================================
' DESCRIPTION: User-defined function providing y=f(x) to be integrated. 
'
' ARGUMENTS:
'      x   :: a real input value.
' ==============================================================================
    Function FX(x As Double) As Double
        FX = x ^ 4
    End Function


' ==============================================================================
' DESCRIPTION: User-defined function providing first derivative, f'(x), explicitly. 
'
' ARGUMENTS:
'      x   :: a real input value.
' ==============================================================================
    Function FU(x As Double) As Double
        FU = 4.0 * x ^ 3
    End Function
End Module


