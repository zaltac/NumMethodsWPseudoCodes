Imports System

' ==============================================================================
'  The main program to test MANTISSA_EXP
' ==============================================================================
Public Module Test_Mantissa_Exp
	Public Sub Main()
        Dim e As Integer
        Dim fl, m As Double

        Console.WriteLine("Enter a real number ")
        fl = Double.Parse(Console.ReadLine())

        Call Mantissa_Exp(fl, m, e)
        Console.WriteLine(" Mantissa = " & m)
        Console.WriteLine(" Exponent = " & e)
    End Sub

'  ==================================================================================
'  CODE1.6-MantissaExp.bas. A Basic (VB) Sub implementing Pseudocode 1.6.           
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
'  DESCRIPTION: A module to determine a floating-point number's                      
'     mantissa and exponent, i.e., fl=Mx10^e                                         
'                                                                                    
'  ON ENTRY                                                                          
'     fl  :: A floating point number.                                                
'                                                                                    
'  ON EXIT                                                                           
'     M   :: Mantissa;                                                               
'     e   :: Exponent.                                                               
'                                                                                    
'  USES                                                                              
'    Abs  :: Built-in Math Intrinsic function returning the absolute value of a real value
'    Floor:: Built-in Math Intrinsic function returning the greatest integer less than    
'      or equal to a real value;                                                     
'    Log10:: Built-in Math Intrinsic function returning the base 10 logarithm of a real value
'                                                                                    
'  REVISION DATE :: 03/22/2024                                                       
'                                                                                    
'  ==================================================================================
	Public Sub Mantissa_Exp(ByVal fl As Double, ByRef m As Double, ByRef e As Integer)
        If Math.Abs(fl) > 0.0 Then
            e = Math.Floor(Math.Log10(fl)) + 1
            m = fl * 10.0 ^ (-e)
        Else
            e = 0
            m = 0.0
        End If
    End Sub
End Module