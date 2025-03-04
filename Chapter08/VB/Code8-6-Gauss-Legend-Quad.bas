' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System


! ==============================================================================
!  The main program to test Gauss_Legendre_Quad
! ==============================================================================
Public Module Test_Gauss_Legendre_Quad
	Public Sub Main()
        Dim x(), w() As Double
        Dim eps As Double = 1.0E-8
        Dim n, i As Integer

        Console.WriteLine("Enter n")
        n = Convert.ToInt32(Console.ReadLine())

	ReDim x(0 To n)
	ReDim w(0 To n)

	Call GAUSS_LEGENDRE_QUAD(n, eps, x, w)

        Console.WriteLine(" ========= Gauss-Legendre Quads =======--==")
        Console.WriteLine(" i           x_i                w_i")
        Console.WriteLine("     -----------------    -----------------")

        For i = 1 To n
            Console.WriteLine("{0,2}    {1,16:F12}    {2,16:F12}", i, x(i), w(i))
        Next
    End Sub


'  ==================================================================================
'  CODE8.6-GAUSS_LEGENDRE_QUAD.bas. A Basic (VB) Sub implementing Pseudocode 8.6.                 
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
'  DESCRIPTION: A subroutine to generate N-point Gauss-Legendre quadrature                     
'    abscissas and weights on [-1,1].                                                          
'                                                                                              
'  ON ENTRY                                                                                    
'     n   :: Number of quadrature points;                                                      
'     eps :: Tolerance, i.e., desired level of numerical accuracy.                             
'                                                                                              
'  ON EXIT                                                                                     
'     x   :: Array of length N containing the abscissas;                                       
'     w   :: Array of length N containing the weights.                                         
'                                                                                              
'  USES                                                                                        
'    Abs  :: Built-in Math function returning the absolute value of a real value;         
'    Cos  :: Built-in Math function returning trig cosine value.  
'                                                                                              
'  REVISION DATE :: 03/04/2025                                                                 
'  ==================================================================================
	Public Sub GAUSS_LEGENDRE_QUAD(ByVal n As Integer, ByVal eps As Double, ByRef x() As Double, ByRef w() As Double)
	Dim u, del As Double
	Dim P0, P1, P2, dP As Double
        Dim pi As Double = 3.1415926535897932385
        Dim m, i, k As Integer

        m = (n + 1) \ 2

        For i = 1 To m
            u = Math.Cos(pi * (4 * i - 1) / (4 * n + 2))

            del = 1.0

            Do While Math.Abs(del) > eps
                P0 = 1.0
                P1 = u
                For k = 2 To n
                    P2 = ((2 * k - 1) * u * P1 - (k - 1) * P0) / k
                    P0 = P1
                    P1 = P2
                Next
                dP = n * (u * P2 - P0) / (u * u - 1.0)
                del = P2 / dP
                u = u - del
            Loop

            x(i) = -u
            w(i) = 2.0 / ((1.0 - u * u) * dP * dP)
            x(n + 1 - i) = u
            w(n + 1 - i) = w(i)
        Next
    End Sub
End Module