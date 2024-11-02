' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_Tridiagonal
	Public Sub Main()
        Const n As Integer = 9
        Dim b(n) As Double
        Dim d(n) As Double
        Dim a(n) As Double
        Dim c(n) As Double
        Dim x(n) As Double
        Dim i As Integer
        Dim s1 As Integer
        Dim sn As Integer

        b = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
        d = {-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0}
        a = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
        c = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -14.0}

        Console.WriteLine(" ********** Input Tridiagonal Matrix & RHS ***********")
        For i = 0 To n-1
            Console.WriteLine("{0,10:F3} {1,10:F3} {2,10:F3} {3,10:F3}", b(i), d(i), a(i), c(i))
        Next i

        s1 = 0
        sn = n-1
        Call TRIDIAGONAL(s1, sn, b, d, a, c, x)

        Console.WriteLine()
        Console.WriteLine(" *** Solution ***")
        Console.WriteLine()
        For i = 0 To n-1
            Console.WriteLine("  x({0})={1,8:F4}", i, x(i))
        Next i
    End Sub


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
Public Sub TRIDIAGONAL(ByVal s1 As Integer, ByVal sn As Integer, ByRef b() As Double, ByRef d() As Double, ByRef a() As Double, ByRef c() As Double, ByRef x() As Double)
        Dim i As Integer
        Dim ratio As Double

        For i = s1 + 1 To sn
            ratio = b(i) / d(i - 1)
            d(i) = d(i) - ratio * a(i - 1)
            c(i) = c(i) - ratio * c(i - 1)
        Next i

        x(sn) = c(sn) / d(sn)
        For i = sn - 1 To s1 Step -1
            x(i) = (c(i) - a(i) * x(i + 1)) / d(i)
        Next i
    End Sub
End Module