program Yukawa 
! Rupak Mukherjee, Sobhan Sounda, arXiv:1705.02444, INJP: 92 (2), 197-203 (2018).
implicit none
integer n,i
real*8 x,y(2),yp(2),xst,xmax,h

open(unit=10,file='Time_Evolution.dat',status='unknown')
open(unit=20,file='Phase.dat',status='unknown')

n = 2
y(1) = 0.0d0 ! u
y(2) = 1.0d0 ! y

xst = 0.0d0
xmax = 33.0d0
h = 0.010d0

do x=xst,xmax,h
  call derivative(x,y,yp,n)
  call rk4(x,h,y,yp,n)
  
  write(10,*) x+h,(y(i),i=1,n)
  write(20,*) y(1)*dcos(x+h), y(1)*dsin(x+h)
enddo

end program yukawa

!================ DIFFERENTIAL EQUATION =============================================

subroutine derivative(x,y,yp,n)
real*8 y(2),yp(2),x,c,f,tn,sn
integer n,i

 c = -9.910000455280d0
 
  yp(1) = y(2)
! yp(2) = -y(1) ! SHM

f = c*(1.0d0+1.0d0/y(1))*dexp(-1.0d0/y(1))
  
if (dabs(y(1)) .le. 6.0d-1) then                  ! Avoiding the numerical singularity
  f = 0.0d0
  yp(2) = -y(1) - f
else
  yp(2) = -y(1) - c*(1.0d0+1.0d0/y(1))*dexp(-1.0d0/y(1))
endif

  return
end 

!==================== RK 4 SUBROUTINE ===============================================

subroutine rk4(x,h,y,yp,n)
implicit none
real*8 x,h,y(2),yp(2),k1(2),k2(2),k3(2),k4(2),dum(2)
integer n,i

 do i = 1,n,1
  k1(i) = yp(i)
  dum(i) = y(i) + k1(i)*h/2.0d0
 enddo
   call derivative(x+h/2.0d0,dum,yp,n)
 do i = 1,n,1
  k2(i) = yp(i)
  dum(i) = y(i) + k2(i)*h/2.0d0
 enddo
   call derivative(x+h/2.0d0,dum,yp,n)
 do i = 1,n,1
  k3(i) = yp(i)
  dum(i) = y(i) + k3(i)*h
 enddo
   call derivative(x+h,dum,yp,n)
 do i = 1,n,1
  k4(i) = yp(i)
  y(i) = y(i) + h/6.0d0*(k1(i) + 2.0d0*k2(i) + 2.0d0*k3(i) + k4(i))
 enddo
 
 return
end

