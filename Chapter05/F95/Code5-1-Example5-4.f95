PROGRAM Example5_4
!  ==================================================================================
!  CODE5.1-EXAMPLE5-4.f95. A Fortran95 module implementing Pseudocode 5.1.                        
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

!  DESCRIPTION: A program that calculates the first and second derivatives of                
!    a set of position (discrete) data.                                                                   
!                                                                                              
!  INPUT VARIABLES                                                                             
!     dt  :: Time increment (in seconds);                                                      
!     s   :: An array of length n, providing the distance traveled in meters.                  
!                                                                                              
!  OUTPUT VARIABLES                                                                            
!     v   :: Array of length n containing velocities (m/s) at discrete points;                 
!     a   :: Array of length n containing acelerations (m2/s) at discrete points.              
!                                                                                              
!  REVISION DATE :: 06/13/2024                                                                 
!  ==================================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: n=6
DOUBLE PRECISION, DIMENSION(n) :: s, v, a
DOUBLE PRECISION :: dt, t
INTEGER :: i
!
dt = 5.0
s  = (/ 0.0d0, 5.45d0, 21.3d0, 82.84d0, 212.86d0, 473.60d0 /)
!
! *** Apply 2nd order forward difference formulas
v(1) = (-s(3)+4.0d0*s(2)-3.0d0*s(1)) / (2.0d0*dt)
a(1) = (-s(4)+4.0d0*s(3)-5.0d0*s(2)+2.0d0*s(1)) / dt**2
!
! *** Apply central difference formulas
DO i=2,(n-1)
   v(i) = (s(i+1)-s(i-1)) / (2.0d0*dt)
   a(i) = (s(i+1)-2.0*s(i)+s(i-1)) / dt**2
ENDDO
!
! *** Apply 2nd order backward difference formulas
v(n) = ( s(n-2)-4.0d0*s(n-1)+3.0d0*s(n)) / (2.0d0*dt)
a(n) = (-s(n-3)+4.0d0*s(n-2)-5.0d0*s(n-1)+2.0d0*s(n)) / dt**2
!
! *** Print out the results
t=0.0
WRITE(*,1)
DO i=1,n
   WRITE(*,2) t, s(i), v(i), a(i)
   t = t + dt
END DO
1 FORMAT(3x,"Time (s)",7x,"s (m)",10x,"v (m/s)",7x,"a (m/s^2)")
2 FORMAT(1x,F8.3,3(3x,F12.4))       
END PROGRAM Example5_4      
