program pendulum 
implicit none
integer n,i
real*8 t,x(2),xp(2),y(2),yp(2),z(2),zp(2),tst,tmax,h

open(unit=10,file='Time_Evolution.dat',status='unknown')
open(unit=20,file='Phase.dat',status='unknown')

n = 2
x(1) = 1.0d0
x(2) = 0.0d0
y(1) = 0.30d0 ! 0.50d0
y(2) = 0.0d0 
z(1) = 0.0d0
z(2) = 0.0d0

tst = 0.0d0
tmax = 10.0d0
h = 0.0010d0

do t=tst,tmax,h
  call derivative(t,x,xp,y,yp,z,zp,n)
  call rk4(t,h,x,xp,y,yp,z,zp,n)
  
  !write(10,*) t+h,(x(i),i=1,n),(y(i),i=1,n),(z(i),i=1,n)
  !write(20,*) y(1)*dcos(x+h), y(1)*dsin(x+h)
   write(10,*) t+h,x(1),y(1),z(1)
   write(20,*) x(1)*dcos(t+h),x(1)*dsin(t+h),y(1)*dcos(t+h),y(1)*dsin(t+h),z(1)*dcos(t+h),z(1)*dsin(t+h)
   
enddo

end program pendulum

!================ DIFFERENTIAL EQUATION =============================================

subroutine derivative(t,x,xp,y,yp,z,zp,n)
implicit none
real*8 x(2),xp(2),y(2),yp(2),z(2),zp(2),t,wx,wy,wz
integer n,i

  wx = 1.0d0
  wy = wx
  wz = wx
  
  xp(1) = + x(2)
  xp(2) = - wx*x(1) + 0.5d0*wx*y(1) + 0.5d0*wx*z(1)
  
  yp(1) = + y(2)
  yp(2) = - wy*y(1) + 0.5d0*wy*z(1) + 0.5d0*wy*x(1)
  
  zp(1) = + z(2)
  zp(2) = - wz*z(1) + 0.5d0*wz*x(1) + 0.5d0*wz*y(1)
  
  return
end 

!==================== RK 4 SUBROUTINE ===============================================

subroutine rk4(t,h,x,xp,y,yp,z,zp,n)
implicit none
real*8 t,h,x(2),xp(2),y(2),yp(2),z(2),zp(2)
real*8 k1x(2),k2x(2),k3x(2),k4x(2),k1y(2),k2y(2),k3y(2),k4y(2),k1z(2),k2z(2),k3z(2),k4z(2)
real*8 dumx(2),dumy(2),dumz(2)
integer n,i

 do i = 1,n,1
  k1x(i) = xp(i)
  dumx(i) = x(i) + k1x(i)*h/2.0d0
  k1y(i) = yp(i)
  dumy(i) = y(i) + k1y(i)*h/2.0d0
  k1z(i) = zp(i)
  dumz(i) = z(i) + k1z(i)*h/2.0d0
 enddo
   call derivative(t+h/2.0d0,dumx,xp,dumy,yp,dumz,zp,n)
 do i = 1,n,1
  k2x(i) = xp(i)
  dumx(i) = x(i) + k2x(i)*h/2.0d0
  k2y(i) = yp(i)
  dumy(i) = y(i) + k2y(i)*h/2.0d0
  k2z(i) = zp(i)
  dumz(i) = z(i) + k2z(i)*h/2.0d0
 enddo
   call derivative(t+h/2.0d0,dumx,xp,dumy,yp,dumz,zp,n)
 do i = 1,n,1
  k3x(i) = xp(i)
  dumx(i) = x(i) + k3x(i)*h
  k3y(i) = yp(i)
  dumy(i) = y(i) + k3y(i)*h
  k3z(i) = zp(i)
  dumz(i) = z(i) + k3z(i)*h
 enddo
   call derivative(t+h,dumx,xp,dumy,yp,dumz,zp,n)
 do i = 1,n,1
  k4x(i) = xp(i)
  x(i) = x(i) + h/6.0d0*(k1x(i) + 2.0d0*k2x(i) + 2.0d0*k3x(i) + k4x(i)) 
  k4y(i) = yp(i)
  y(i) = y(i) + h/6.0d0*(k1y(i) + 2.0d0*k2y(i) + 2.0d0*k3y(i) + k4y(i))
  k4z(i) = zp(i)
  z(i) = z(i) + h/6.0d0*(k1z(i) + 2.0d0*k2z(i) + 2.0d0*k3z(i) + k4z(i))
 enddo
 
 return
end
