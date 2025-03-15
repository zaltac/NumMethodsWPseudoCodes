PROGRAM Test_QRIteration
! ==============================================================================
!  The main program to test the module Basic_QR.
! ==============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=4
DOUBLE PRECISION :: d(n), e(n), V(n,n), eps 
INTEGER :: i, j, maxit
!

d=(/ 3.0d0, 8.0d0, 6.0d0, 9.0d0 /)
e=(/ 4.0d0, 2.0d0, 1.0d0, 0.0d0 /)
!     
eps  = 1.0d-3
maxit= 299
! c
CALL Basic_QR(n,d,e,eps,maxit,V)
!
PRINT "(/'*** Eigenvalues ***'/)"
PRINT 2, (d(i),i=1,n)

PRINT "(/'*** Eigenvectors ***'/)"
DO i=1,n
  PRINT 2, (V(i,j),j=1,n)
ENDDO 	
2 FORMAT(8(1x,f12.7))
END PROGRAM Test_QRIteration
        
 
SUBROUTINE Basic_QR(n,d,e,eps,maxit,V)
!  ==================================================================================
!  CODE11.6-Basic_QR.f95. A Fortran95 module implementing Pseudocode 11.6.                        
!  
!  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
!  First Edition. (c) By Zekeriya ALTAÇ (2024).
!  ISBN: 978-1-032-75474-1 (hbk)
!  ISBN: 978-1-032-75642-4 (pbk)
!  ISBN: 978-1-003-47494-4 (ebk)
!  
!  DOI : 10.1201/9781003474944
!  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
!  
!  This free software is complimented by the author to accompany the textbook.
!  E-mail: altacz@gmail.com.
!
!   DESCRIPTION: A fortran95 subroutine implementing the QR Factorization algorithm to a symmetric       
!      tridiagonal matrix to find its eigenvalues and eigenvectors.                            
!                                                                                              
!   ON ENTRY                                                                                   
!      n   :: Dimension attribute of the tridiagonal matrix (nxn);                             
!      d   :: An array of length n containing the main diagonal, d(1) ... d(n);                
!      e   :: An array of length n containing the subdiagonal, e(1) ... e(n-1).                
!                                                                                              
!   ON EXIT                                                                                    
!      d   :: An array of length n containing the eigenvalues;                                 
!      V   :: A square matrix containing the eigenvector(nxn).                                 
!                                                                                              
!   USES                                                                                       
!     SQRT :: Built-in Intrinsic function returning the square root of a real value;           
!                                                                                              
!   REVISION DATE :: 03/15/2025                                                                
!  ==================================================================================
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, maxit
DOUBLE PRECISION, INTENT(INOUT) :: d(n), e(n)
DOUBLE PRECISION, INTENT(OUT) :: V(n,n)
DOUBLE PRECISION :: c(n), s(n), t, eps, err, rho, q, r
INTEGER ::  p, k, i 
!
err= 1.0d0
 p = 0
 V =0.0d0  
DO i=1,n   ! initialize matrix of eigenvectors with the identity matrix
   V(i,i)=1.0d0
ENDDO
PRINT "('*** Iteration History ***')",
DO WHILE ((err>eps).AND.(p<maxit)) ! === Start QR Iteration loop 
   t = e(1)
   DO k = 1, n-1
      rho  = dSQRT( d(k)**2 + t*t )
      c(k) = d(k) / rho
      s(k) = t / rho
      d(k) = rho
        t  = e(k)
      e(k) = t * c(k) + d(k+1) * s(k)
     d(k+1)=-t * s(k) + d(k+1) * c(k)
     IF(k/=n-1) THEN
       t  = e(k+1)
       e(k+1)= t * c(k)
     ENDIF
     DO i=1,n    
        q = V(i, k )     
        r = V(i,k+1)          
        V(i, k )= c(k) * q + s(k) * r
        V(i,k+1)=-s(k) * q + c(k) * r             
     END DO      
   ENDDO
!       
! *** Construct RQ Matrix   	   
   DO k=1,n-1
       d(k) = d(k) * c(k) + e(k) * s(k)
          t = d(k+1)
       e(k) = t * s(k)
      d(k+1)= t * c(k)
   ENDDO 
  do i=1,n
     print*, (v(i,k), k=1,n)
  end do 
  pause "pozed"
!	 
   err=0.0d0
   DO i=1,n   !  Calculate L2 norm of the vector of superdiagonal elements
      err = err + e(i)**2
   ENDDO 
   err = dSQRT(err)
   p = p + 1
   PRINT 1, p, err
END DO  ! === End of While loop
IF (p==maxit) THEN
    Print*,"Convergence within tolerance was not achived. Error is ", err
    ENDIF
! ********** END OF ITERATIONS ****************
1 FORMAT(3x,"iter=",i3,2x,E14.4)
END SUBROUTINE Basic_QR 
