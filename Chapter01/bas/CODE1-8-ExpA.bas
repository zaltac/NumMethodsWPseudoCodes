' runs Ok on    https://dotnetfiddle.net/ 

Imports System

Public Module Test_EXPA
' ==============================================================================
'  The main program to test function EXPA
' ==============================================================================
	Public Sub Main()
        Dim x, y As Double
        Dim eps As Double

        Console.Write("Enter x: ")
        x = Console.ReadLine()
        eps = 1.0E-6
        Console.WriteLine("x ={0:F4} and e^x = {1:F8}", x, EXPA(x, eps))
	End Sub


	Public Function EXPA(x As Double, eps As Double) As Double
'  ==================================================================================
'  CODE1.8-ExpA.bas. A Basic (VB) Sub implementing Pseudocode 1.8.                                
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
'  DESCRIPTION: A function to compute e^x adaptively using the MacLaurin series                
'     within a user-defined tolerance.                                                         
'                                                                                              
'  ARGUMENTS                                                                                   
'     x   :: A real input (exponent) value;                                                    
'    eps  :: A user-defined convergence tolerance.                                             
'                                                                                              
'  USES                                                                                        
'   Float :: A built-in Math function that converts an integer argument to a real value.  
'    Abs  :: A built-in Math function returning the absolute value of a real value.         
'                                                                                              
'  REVISION DATE :: 04/11/2024                                                                                                                                                          
'  ==================================================================================
        Dim sums As Double
        Dim term As Double
        Dim k As Integer

        k = 0
        term = 1.0
        sums = 1.0

        While Math.Abs(term) > eps * Math.Abs(sums)
            k += 1
            term = term * x / k
            sums += term
        End While

        Return sums

    End Function

End Module