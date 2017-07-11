program t4

implicit none
real :: width,height,k,q,dym,dxm,dx,dy
integer :: nvx,nvy,npx,npy,n,cc1,cc2,cc3,cc4
real, dimension (100000,10) :: msh,an
real, dimension (1000000) :: t

call reader (width,height,nvx,nvy,k,q,cc1,cc2,cc3,cc4)
call mesh (width, height, nvx, nvy, npx, npy, n, dym, dxm, msh, dx, dy,an)
call calc (npx,dy,msh,n,k,width,height,an,q,t,dxm)
call gmsh (msh,t,width,npx,n)

end program
!--------------------------------------------------------------------------
subroutine reader (width,height,nvx,nvy,k,q,cc1,cc2,cc3,cc4)
real :: width, height, k, q
integer :: nvx, nvy,cc1,cc2,cc3,cc4

open(unit=1, file='reader.dat')
read(1,*)
read(1,*)
read(1,*)width
read(1,*)height
read(1,*)nvx
read(1,*)nvy
read(1,*)k
read(1,*)q
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)cc1
read(1,*)cc2
read(1,*)cc3
read(1,*)cc4


end subroutine
!--------------------------------------------------------------------------
subroutine mesh (width, height, nvx, nvy, npx, npy, n, dym, dxm, msh, dx, dy,an)
real :: width, height, dx, dy, dxm, dym
integer :: nvx, nvy, npx, npy, n
real, dimension (100000,10) :: msh, an

dx=width/nvx
dy=height/nvy
dxm=dx/2
dym=dy/2
npx=nvx+2
npy=nvy+2

n=0
msh=0.

do j=1,npy
        if (j==1) then
           y=0.
           else if ( j==2 .or. j==npy) then
                  y=y+dym
                 else if (j>2 .and. j<npy) then
                       y=y+dy
        end if

    do i=1,npx
    n=n+1
        if (i==1) then
            x=0.
            else if (i==2 .or. i==npx) then
                x=x+dxm
                else if (i>2 .and. i<npx) then
                     x=x+dx
        end if

    msh(n,1)=x
    msh(n,2)=y

    end do
    end do


    do i=1,n
      if (msh(i,1)==0) then
        an(i,1)=0.
        an(i,2)=0.
        an(i,3)=0.
        an(i,4)=0.
      else if (msh(i,2)==0) then
        an(i,1)=0.
        an(i,2)=0.
        an(i,3)=0.
        an(i,4)=0.
      else if (msh(i,1)>=width-1.e-6) then
        an(i,1)=0.
        an(i,2)=0.
        an(i,3)=0.
        an(i,4)=0.
      else if (msh(i,2)>=height-1.e-6) then
        an(i,1)=0.
        an(i,2)=0.
        an(i,3)=0.
        an(i,4)=0.
      else if ((msh(i,1)>=dxm-1.e-6).and.(msh(i,1)<=dxm+1.e-6)) then
        an(i,1)=0.
        an(i,2)=k*(dy /(msh(i+1,1)-msh(i,1)))
        an(i,3)=k*(dx /(msh(i+npx,2)-msh(i,2)))
        an(i,4)=k*(dx /(msh(i,2)-msh(i-npx,2)))
      else if ((msh(i,1)>=((width-dxm)-1.e-6)).and.(msh(i,1)<=((width-dxm)+1.e-6)))  then
        an(i,1)=k*(dy /(msh(i,1)-msh(i-1,1)))
        an(i,2)=0.
        an(i,3)=k*(dx /(msh(i+npx,2)-msh(i,2)))
        an(i,4)=k*(dx /(msh(i,2)-msh(i-npx,2)))
      else
        an(i,1)=k*(dy /(msh(i,1)-msh(i-1,1)))
        an(i,2)=k*(dy /(msh(i+1,1)-msh(i,1)))
        an(i,3)=k*(dx /(msh(i+npx,2)-msh(i,2)))
        an(i,4)=k*(dx /(msh(i,2)-msh(i-npx,2)))
      endif
      an(i,5)=an(i,1)+an(i,2)+an(i,3)+an(i,4)
    end do
end subroutine
!--------------------------------------------------------------------------
subroutine calc (npx,dy,msh,n,k,width,height,an,q,t,dxm)
implicit none
real :: dy,dxm, k, q, width, height, t(1000000), t0(1000000), erro(1000000), erro1, erro0
real, dimension (100000,10) :: msh, an
integer i,n,npx
open (unit=2, file='temperatura.dat')

do i=1,n
  if (msh(i,2)==0) then
    t(i)= 400.
  else if ((msh(i,2)<=height+1.e-6).and.(msh(i,2)>=height-1.e-6)) then
    t(i)= 100.
  else
    t(i)=0
  endif
end do

do i=1,n
t0(i)=1
erro(i)=100
end do

erro1=100

do while(erro1>=1.E-6)
erro0=0
  do i=1,n
	if (msh(i,2)==0) then
	erro(i)=0.
	 else if (msh(i,2)>=height-1.e-6) then
	 erro(i)=0.
	  else if (msh(i,1)==0) then
	  erro(i)=0.
	   else if (msh(i,1)>=width-1.e-6) then
	   erro(i)=0.
		else if ((msh(i,1)>=dxm-1.e-6).and.(msh(i,1)<=dxm+1.e-6)) then
		t0(i)=t(i)
	        t(i)= (an(i,1)*t(i-1) + an(i,2)*t(i+1) + an(i,3)*t(i+npx) + an(i,4)*t(i-npx) + (q*dy))/ an(i,5)
		erro(i)=(abs((t(i)-t0(i))/t(i)))*100
		else
		t0(i)=t(i)
	        t(i)= (an(i,1)*t(i-1) + an(i,2)*t(i+1) + an(i,3)*t(i+npx) + an(i,4)*t(i-npx))/ an(i,5)
		erro(i)=(abs((t(i)-t0(i))/t(i)))*100
	end if

	if (erro0>=erro(i)) then
		erro1=erro0
		else
		erro1=erro(i)
	end if
	erro0=erro1

	if (msh(i,1)==0) then

	t(i)=t(i+1)+(q*dy/k)

	else if (msh(i,1)>=width-1.e-6) then
 	t(i)=t(i-1)

	end if
 end do

end do


do i=1,n
print*, t(i),i
end do

end subroutine calc

!--------------------------------------------------------------------------
subroutine gmsh (msh,t,width,npx,n)
implicit none
integer :: i,npx,n
real :: width, t(1000000)
real, dimension (100000,10) :: msh

open (unit=3 , file ="./gmsh.pos",status="unknown")

!========================================================================================
!POS-PROCESSAMENTO - GMSH
write(3,*)'View "Temperature" {'

do i=1,(n-npx)
	if ((msh(i,1)<=width+1.d-6).and.(msh(i,1)>=width-1.d-6)) then
		else
		write(3,*)'SQ(',msh(i+npx,1),',',msh(i+npx,2),',0,',msh(i,1),',',msh(i,2),',0,'
		write(3,*)msh(i+1,1),',',msh(i+1,2),',0,',msh(i+npx+1,1),',',msh(i+npx+1,2),',0)'
		write(3,*)'{',t(i+npx),',',t(i),',',t(i+1),',',t(i+npx+1),'};'
	end if
end do
write(3,*)'};'

!===================================================================================

return
end
