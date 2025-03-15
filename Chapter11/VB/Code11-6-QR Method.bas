' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System
				
Public Module Test_Basic_QR
' ==============================================================================
'  The main program to test the module Basic_QR.BAS
' ==============================================================================
	Public Sub Main()
        Dim n As Integer = 4
        Dim d(n - 1) As Double, e(n - 1) As Double, V(n - 1, n - 1) As Double
        Dim eps As Double, i As Integer, j As Integer, maxit As Integer

        d = {3.0, 8.0, 6.0, 9.0}
        e = {4.0, 2.0, 1.0, 0.0}

        eps = 1.0E-5
        maxit = 299

        Call Basic_QR(n, d, e, eps, maxit, V)

		        ' print out the results
        Console.WriteLine(" ")
        Console.WriteLine(" Eigenvalues ")
        For i = 0 To n - 1
			Console.Write("{0,12:F6}", d(i))
        Next
        
		Console.WriteLine(" ")
        Console.WriteLine(" ")
        Console.WriteLine(" Eigenvectors")
		For i = 0 To n - 1
			For j = 0 To n - 1
				Console.Write("{0,12:F6}", V(i, j))
            Next
        	Console.WriteLine(" ")
		Next
    End Sub

	Public Sub Basic_QR(ByVal n As Integer, ByRef d() As Double, ByRef e() As Double, ByVal eps As Double, ByVal maxit As Integer, ByRef V(,) As Double)
	'  ==================================================================================
	'  CODE11.6-Basic_QR.bas. A Basic (VB) Sub implementing Pseudocode 11.6.                          
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
	'   DESCRIPTION: A Visual Basic Sub program implementing the QR Factorization algorithm to a symmetric       
	'      tridiagonal matrix to find its eigenvalues and eigenvectors.                            
	'                                                                                              
	'   ON ENTRY                                                                                   
	'      n   :: Dimension attribute of the tridiagonal matrix (nxn);                             
	'      d   :: An array of length n containing the main diagonal, d(1) ... d(n);                
	'      e   :: An array of length n containing the subdiagonal, e(1) ... e(n-1).                
	'                                                                                              
	'   ON EXIT                                                                                    
	'      d   :: An array of length n containing the eigenvalues;                                 
	'      V   :: A square matrix containing the eigenvector(nxn).                                 
	'                                                                                              
	'   USES                                                                                       
'     sqrt :: Built-in Math library function returning the square root of a real value;           
	'                                                                                              
	'   REVISION DATE :: 03/15/2025                                                                
	'  ==================================================================================
        Dim c(n - 1) As Double, s(n - 1) As Double, t As Double, err As Double, rho As Double, q As Double, r As Double
        Dim p As Integer, k As Integer, i As Integer

        err = 1.0
        p = 0
        For i = 0 To n - 1
            V(i, i) = 1.0
        Next

        Console.WriteLine("*** Iteration History ***")
        While (err > eps) And (p < maxit)
            t = e(0)
            For k = 0 To n - 2
                rho = Math.Sqrt(d(k) * d(k) + t * t)
                c(k) = d(k) / rho
                s(k) = t / rho
                d(k) = rho
                t = e(k)
                e(k) = t * c(k) + d(k + 1) * s(k)
                d(k + 1) = -t * s(k) + d(k + 1) * c(k)
                If k <> n - 2 Then
                    t = e(k + 1)
                    e(k + 1) = t * c(k)
                End If
                For i = 0 To n - 1
                    q = V(i, k)
                    r = V(i, k + 1)
                    V(i, k) = c(k) * q + s(k) * r
                    V(i, k + 1) = -s(k) * q + c(k) * r
                Next
            Next

            For k = 0 To n - 2
                d(k) = d(k) * c(k) + e(k) * s(k)
                t = d(k + 1)
                e(k) = t * s(k)
                d(k + 1) = t * c(k)
            Next

            err = 0.0
            For i = 0 To n - 1
                err = err + e(i) * e(i)
            Next
            err = Math.Sqrt(err)
            p = p + 1
            Console.WriteLine("iter={0,3}  {1,14:E4}", p, err)
        End While

        If p = maxit Then
            Console.WriteLine("Convergence within tolerance was not achieved. Error is {0}", err)
        End If
    End Sub
End Module