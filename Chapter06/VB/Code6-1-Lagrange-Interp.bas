' runs on    https://dotnetfiddle.net/8K7fkp
' ====================================================================================================
' NOTE: Since array indexes in VB start with zero, pseudocodes prepared for indexes starting with "1"
' have been changed to suit this feature.
' ====================================================================================================
Imports System

Public Module Test_Lagrange
	Sub Main()
    	Dim n As Integer = 3
	    Dim x() As Double = {0.08, 0.250, 0.50, 0.90}
    	Dim f() As Double = {0.25, 0.625, 0.81, 0.43}
		Dim xval, fval As Double

    	Console.Write("Enter x: ")
    	xval = Double.Parse(Console.ReadLine())
		
		fval = LAGRANGE_EVAL(n, xval, x, f)

		Console.WriteLine("fval = {0}", fval)
	End Sub
	
'  ==================================================================================
'  CODE6.1-LAGRANGE-P.bas. A Basic (VB) Sub implementing Pseudocode 6.1.                        
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
'  DESCRIPTION: A sub module to compute the Lagrange polynomials at x=xval.                        
'                                                                                              
'  ON ENTRY                                                                                    
'    n    :: The number of data in the set minus 1;                                            
'    xval :: x value at which the dataset is to be interpolated;                               
'    x    :: Array of length n+1 containing abscissas, k=0,1,2,...,n.                          
'                                                                                              
'  ON EXIT                                                                                     
'    L    :: Array of length (n+1) containing Lagrange polynomials                             
'            that is, L(k) = L_k(xval) for k=0,1,2,...,n.                                      
'                                                                                              
'  REVISION DATE :: 02/28/2025                                                                 
'  ==================================================================================
	Public Function LAGRANGE_P(ByVal n As Integer, ByVal xval As Double, ByVal x() As Double, ByRef L() As Double) 
		Dim i, k As Integer
		
        For k = 0 To n 
            L(k) = 1.0
            For i = 0 To n 
                If i <> k Then
                    L(k) *= (xval - x(i)) / (x(k) - x(i))
                End If
            Next
        Next
	End Function

'  ==================================================================================
'  CODE6.1-LAGRANGE_EVAL.bas. A Basic (VB) Sub implementing Pseudocode 6.1.                          
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
'  DESCRIPTION: A function to evaluate Lagrange interpolating polynomial for an arbitrary xval
'   within x0 <= xval <= xn, f=f(xval)=fval.                                                   
'                                                                                              
'  ARGUMENTS                                                                                   
'    n    :: The number of data in the set minus 1;                                            
'    xval :: x value at which the dataset is to be interpolated;                               
'    x    :: Array of length (n+1) containing abscissas, k=0,1,2,...,n;                        
'    f    :: Array of length (n+1) containing ordinates, k=0,1,2,...,n.                        
'                                                                                              
'  USES                                                                                        
'    LAGRANGE_P :: A subroutine generating the Lagrange polynomials.                                           
'                                                                                              
'  REVISION DATE :: 02/28/2025                                                                 
'  ==================================================================================
	Public Function LAGRANGE_EVAL(ByVal n As Integer, ByVal xval As Double, ByVal x() As Double, ByVal f() As Double) As Double
	    Dim L(n) As Double
	    Dim fval As Double = 0.0

	    Call LAGRANGE_P(n, xval, x, L)

	    For k As Integer = 0 To n
	        fval += L(k) * f(k)
	    Next
		Return fval
	End Function
End Module