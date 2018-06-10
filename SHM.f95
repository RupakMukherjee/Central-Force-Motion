program SHM 
implicit none
integer n,i
real*8 x,y(2),yp(2),xst,xmax,h,Error

open(unit=10,file='Time_Evolution.dat',status='unknown')
open(unit=20,file='Phase.dat',status='unknown')

n = 2
y(1) = 0.0d0 ! u
y(2) = 1.0d0 ! y

xst = 0.0d0
xmax = 1000.0d0
h = 0.0010d0

do x=xst,xmax,h
  
  call derivative(x,y,yp,n)
  call rk4(x,h,y,yp,n)
  
  Error = dsin(x+h) - y(1)
  
  !if (mod(x,10.0d0) < 1.0d-5) then
  write(10,*) x+h,(y(i),i=1,n),Error
  write(20,*) y(1)*dcos(x+h), y(1)*dsin(x+h)
  !endif
  
  call flush (10)
  call flush (20)
enddo

end program SHM

!================ DIFFERENTIAL EQUATION =============================================

subroutine derivative(x,y,yp,n)
real*8 y(2),yp(2),x
integer n,i

  yp(1) =   y(2)
  yp(2) = - y(1) ! SHM
  
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
