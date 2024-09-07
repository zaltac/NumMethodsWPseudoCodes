Imports System

Public Module TestFactorial
	Public Sub Main()
        Dim n As Integer

        Console.Write("Enter n ")
        n = Console.ReadLine()
        
        Console.WriteLine("{0}! = {1}", n, Factorial(n))
    End Sub

	Public Function Factorial(n As Integer) As Integer
'  ==================================================================================
'  CODE1.4-factorial.bas. A Basic (VB) Sub for implementing Pseudocode 1.4.          
' 
'  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
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
'  DESCRIPTION: A recursive function for computing n!                                
'                                                                                    
'  INPUT ARGUMENT                                                                    
'    n    :: Integer (n>=0)                                                          
'                                                                                    
'  ON EXIT                                                                           
'   Result:: n!                                                                      
'                                                                                    
'  REVISION DATE :: 03/21/2024                                                       
'                                                                                    
'  ==================================================================================

        If n < 0 Then
            Console.WriteLine("Error, Illegal Input")
            Console.WriteLine("Argument should be n>=0, you entered {0}", n)
            Return -1
        Else
            If n <= 1 Then
                Return 1
            Else
                Return n * Factorial(n - 1)
            End If
        End If
    End Function
End Module