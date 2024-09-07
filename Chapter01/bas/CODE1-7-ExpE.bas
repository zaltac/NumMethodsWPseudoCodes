' runs Ok on    https://dotnetfiddle.net/ 
Imports System


Public Module Test_EXPE

' ==============================================================================
'  The main program to test function ExpE
' ==============================================================================
    Sub Main()
        Dim x As Double
        Dim n As Integer

        Console.Write("Enter x: ")
        x = Console.ReadLine()
		Console.Write("Enter n: ")
        n = Console.ReadLine()

        Console.WriteLine("e^{0:F3} = {1:F7}", x, EXPE(x, n))
    End Sub


    Public Function EXPE(ByVal x As Double, ByVal n As Integer) As Double
'  ==================================================================================
'  CODE1.7-ExpE.bas. A Basic (VB) Sub implementing Pseudocode 1.7.                                
' 
'  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
'  First Edition. (c) By Zekeriya ALTAÇ (2024).
'  ISBN: 978-1-032-75474-1 (hbk)
'  ISBN: 978-1-032-75642-4 (pbk)
'  ISBN: 978-1-003-47494-4 (ebk)
'  
'  DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
'  
'  This free software is complimented by the author to accompany the textbook.
'  E-mail: altacz@gmail.com.
'                                                                                              
'  DESCRIPTION: A function to compute e^x using the MacLaurin series with specified            
'     number of terms.                                                                         
'                                                                                              
'  ARGUMENTS                                                                                   
'     x   :: A real input (exponent) value;                                                    
'     n   :: The number of terms of the MacLauring series to be included.                      
'                                                                                              
'  USES                                                                                        
'    Convert.ToDouble:: A built-in intrinsic function that converts an integer argument to 
'     a double precision value.  
'                                                                                              
'  REVISION DATE :: 03/18/2024                                                                                                                                                        
'  ==================================================================================
        Dim sums As Double = 1.0
        Dim term As Double = 1.0

        For k As Integer = 1 To n - 1
            term *= x / Convert.ToDouble(k)
            sums += term
        Next

        Return sums
    End Function

End Module