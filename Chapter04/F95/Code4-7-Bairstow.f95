PROGRAM Test_Bairstow
! ==============================================================================
!  Main program to test SUBROUTINE BAIRSTOW
! ==============================================================================
IMPLICIT NONE
DOUBLE PRECISION :: A(0:19), XRE(19), XIM(19)
INTEGER:: n, i, maxit, iprnt 
DOUBLE PRECISION :: p0, q0, eps

!
n=5
a(0)=  1.0D0; a(1)= -5.0D0; a(2)=-15.0D0; &
a(3)= 85.0D0; a(4)=-26.0D0; a(5)=-120.D0
! -------------------------------------------------------------
iprnt = 1  ! output control key  
! iprnt = 0  does not print iteration details, 
!       = 1  prints a short iteration history 
!       = 2  prints all iteration history  
! -------------------------------------------------------------
 maxit = 99
 p0  = 0.0D0 
 q0  = 0.0D0
 eps = 0.5D-4
!
CALL BAIRSTOW(n,p0,q0,a,eps,maxit,iprnt,xre,xim)
PRINT*, "    ======== All the Roots are ========="
DO i=1,n
   PRINT 1, i,xre(i),xim(i)
END DO
PRINT*, "    ===================================="        
1 FORMAT(5x,"Root(",i2,") = ",f8.5," + ( ",f8.5," ) i")
END PROGRAM Test_Bairstow


SUBROUTINE BAIRSTOW(n,p0,q0,a,eps,maxit,iprnt,xre,xim)
! ==============================================================================
! CODE4-7-BAIRSTOW.F95. A fortran program for implementing Pseudocode 4.7.
!     
! NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
! First Edition. (c) By Zekeriya ALTAÇ (2024).
! ISBN: 9781032754741 (hbk)
! ISBN: 9781032756424 (pbk)
! ISBN: 9781003474944 (ebk)
!
! DOI : 10.1201/9781003474944
! C&H/CRC PRESS, Boca Raton & London. 
!  
! This free software is complimented by the author to accompany the textbook.
! E-mail: altacz@gmail.com
!
! DESCRIPTION: A subroutine to find all real and/or imaginary roots of a polynomial  
!   of the n'th degree using the BAIRSTOW’s method.
!
! ON ENTRY
!   n    :: Degree of the polynomial;
!  p0,q0 :: Initial guesses for a quadratic equation; i.e., for p and q;
!    a   :: Array of length (n+1) containing the coefficients of polynomial defined as 
!                a0 x^n + a1 x^(n-1) + ... + an = 0 
!   eps  :: Convergence tolerance;
!  maxit :: Maximum number of iterations permitted;
!  iprnt :: printing key, =0 do not print intermediate results, <> 0 print intermediates.
!
! ON EXIT
!   xre  :: Array of length n containing real parts of the roots;
!   xim  :: Array of length n containing imaginary parts of the roots.
!
! OTHER VARIABLES
!    b   :: Array of length [n] containing coefficients of quotient polynomial (0<=k<=n-2);
!    c   :: Array of length [n] containing coefficients of partial derivatives. 
!
! USES
!   ABS  :: Built-in Intrinsic function returning the absolute value of a real value;
!   QUADRATIC_EQ :: Subroutine that solves a quadratic equation of the form x2 + p x + q = 0 (see CODE1.3).
!        
!
! REVISION DATE :: 04/29/2024
! ==============================================================================
IMPLICIT NONE
INTEGER, INTENT( IN)  :: maxit
INTEGER, INTENT(INOUT):: n
DOUBLE PRECISION, INTENT(INOUT) :: p0, q0, a(0:n)
DOUBLE PRECISION, INTENT(OUT) :: xre(n), xim(n), eps
DOUBLE PRECISION :: b(0:n), c(0:n), xr(2), xi(2)
DOUBLE PRECISION :: p, delp, q, delq, delM, cbar, del, del1, del2
INTEGER :: i, k, m, kount, iprnt
!
DO k=n,0,-1
   a(k)=a(k)/a(0)  ! Normalize a's by a(0) 
END DO  
!
m=n ! Save n for later use 
kount=0
DO WHILE (n>1)      
   p = p0  ;  q = q0   ! Initialize
   k = 0
   delM = 1.0d0 
   DO WHILE (delM>eps .AND. k<=maxit) !**** inner loop ****
      k = k + 1    
      b(0) = 1.0d0 ; c(0) = 1.0d0
      b(1) = a(1)-p; c(1) = b(1)-p
      DO i=2,n
         b(i)=a(i)-p*b(i-1)-q*b(i-2)
         c(i)=b(i)-p*c(i-1)-q*c(i-2)
      END DO          
      cbar= c(n-1)-b(n-1)     ! Construct and solve the 2x2 linear system
      del = c(n-2)*c(n-2)-cbar*c(n-3)
      del1= b(n-1)*c(n-2)-b(n)*c(n-3)
      del2= b(n)*c(n-2)-b(n-1)*cbar
      delp= del1/del;  delq= del2/del   
      p = p + delp ;  q = q + delq   ! Find new estimates
      delM= DABS(delp) + DABS(delq)  ! Calculate L1 norm
      IF(iprnt==1) THEN
          PRINT 1, k, delM , p , q           
      ELSE IF(iprnt==2) THEN 
          PRINT 2, k,delp,delq,delM,p,q
          PRINT 3
          PRINT 4
          DO i=0,n
             PRINT 5, i,a(i),b(i),c(i)
          END DO
          PRINT 4
      ENDIF
   END DO  !******************* end of inner loop ***********
   IF(k-1==maxit) THEN
     PRINT*, "Quadratic factor did not converge after",k-1," iterations"
     PRINT*, "Recent values of p, q, delM are ",p,q,delM
     PRINT*, "Corresponding roots may be questionable ..."
   ENDIF 
   CALL QUADRATIC_EQ(p,q,xr,xi)
   kount=kount+1
   xre(kount)=xr(1); xim(kount)=xi(1)
   kount=kount+1
   xre(kount)=xr(2); xim(kount)=xi(2)    
!   
   PRINT 6, p,q
   n=n-2
   DO i=0,n
      a(i)=b(i)
   ENDDO
   IF(n==1) THEN
      kount=kount+1
      xre(kount)= -a(1)   
      xim(kount)= 0.0d0 
      PRINT 7, a(1)  
   ENDIF
   END DO
   n=m
  
1 FORMAT(1x,"Iter= ",i2,2x,"delM= ",1PG10.4,3x,"p= ",1PG12.5,3x,"q= ",1PG12.5)
2 FORMAT(/5x,"iter=",i4,/,5x,"---------",/,&
   5x," dp =",G14.6,6x,"  dq =",G14.6,6x,"  delM=",G14.6,/,&
   5x,"  p =",G14.6,6x,"   q =",G14.6/)    
3 FORMAT(5x," k",10x,"a(k)",11x,"b(k)",11x,"c(k)")
4 FORMAT(5x,47('-'))
5 FORMAT(5x,i2,3(3x,f12.6))   
6 FORMAT(8x," ======== FOUND A QUADRATIC FACTOR  ======== ",/, &
         10x," x*x + (",F12.6,")*x + (",F12.6,") ",/, &
          8x," =========================================== "///)
7 FORMAT(8x," ======== FOUND A LINEAR FACTOR  ======== ",/, &
         13x,"      x  + (",F12.6,") ",/, &
          8x," ======================================== "///)
END SUBROUTINE BAIRSTOW
