program pendulum 
implicit none
integer n,i
real*8 x,y(2),yp(2),z(2),zp(2),xst,xmax,h

open(unit=10,file='Time_Evolution.dat',status='unknown')
open(unit=20,file='Phase.dat',status='unknown')

n = 2
y(1) = 1.0d0 
y(2) = 0.0d0 
z(1) = 0.0d0
z(2) = 0.0d0

xst = 0.0d0
xmax = 10.0d0
h = 0.0010d0

do x=xst,xmax,h

  call derivative(x,y,yp,z,zp,n)
  call rk4(x,h,y,yp,z,zp,n)
  
  !write(10,*) x+h,(y(i),i=1,n),(z(i),i=1,n)
   write(10,*) x+h,y(1),z(1)
   write(20,*) y(1)*dcos(x+h), y(1)*dsin(x+h),z(1)*dcos(x+h), z(1)*dsin(x+h)
   
enddo

end program pendulum

!================ DIFFERENTIAL EQUATION =============================================

subroutine derivative(x,y,yp,z,zp,n)
implicit none
real*8 y(2),yp(2),z(2),zp(2),x,wy,wz
integer n,i

  wy = 1.0d0
  wz = 1.0d0*wy
  
  yp(1) = + y(2)
  yp(2) = - wy*y(1) + wy*z(1) !- y(1) ! SHM
  
  zp(1) = + z(2)
  zp(2) = - wz*z(1) + wz*y(1)
  
  return
end 

!==================== RK 4 SUBROUTINE ===============================================

subroutine rk4(x,h,y,yp,z,zp,n)
implicit none
real*8 x,h,y(2),yp(2),z(2),zp(2),k1y(2),k2y(2),k3y(2),k4y(2),dumy(2),k1z(2),k2z(2),k3z(2),k4z(2),dumz(2)
integer n,i

 do i = 1,n,1
  k1y(i) = yp(i)
  dumy(i) = y(i) + k1y(i)*h/2.0d0
  k1z(i) = zp(i)
  dumz(i) = z(i) + k1z(i)*h/2.0d0
 enddo
   call derivative(x+h/2.0d0,dumy,yp,dumz,zp,n)
 do i = 1,n,1
  k2y(i) = yp(i)
  dumy(i) = y(i) + k2y(i)*h/2.0d0
  k2z(i) = zp(i)
  dumz(i) = z(i) + k2z(i)*h/2.0d0
 enddo
   call derivative(x+h/2.0d0,dumy,yp,dumz,zp,n)
 do i = 1,n,1
  k3y(i) = yp(i)
  dumy(i) = y(i) + k3y(i)*h
  k3z(i) = zp(i)
  dumz(i) = z(i) + k3z(i)*h
 enddo
   call derivative(x+h,dumy,yp,dumz,zp,n)
 do i = 1,n,1
  k4y(i) = yp(i)
  y(i) = y(i) + h/6.0d0*(k1y(i) + 2.0d0*k2y(i) + 2.0d0*k3y(i) + k4y(i))
  k4z(i) = zp(i)
  z(i) = z(i) + h/6.0d0*(k1z(i) + 2.0d0*k2z(i) + 2.0d0*k3z(i) + k4z(i))
 enddo
 
 return
end
