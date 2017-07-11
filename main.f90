program t4

implicit none
real :: a,b,k,q,dym,dxm,dx,dy
integer :: nva,nvb,npx,npy,n,cc1,cc2,cc3,cc4,i
real, dimension (100000,10) :: msh,an
real, dimension (1000000) :: T

call leitor (a,b,nva,nvb,k,q,cc1,cc2,cc3,cc4)
call malha (a, b, nva, nvb, npx, npy, n, dym, dxm, msh, dx, dy)
call distribuicao2d (npx,npy,dx,dy,msh,n,k,a,b,dxm,dym,q,an)
call calculo (npx,npy,dx,dy,msh,n,k,a,b,an,q,T,dxm)
call imprime (a,b,nva,nvb,k,q,cc1,cc2,cc3,cc4)
call prints (msh,T,a,npx,n)

end program
!--------------------------------------------------------------------------
subroutine leitor (a,b,nva,nvb,k,q,cc1,cc2,cc3,cc4)
real :: a, b, k, q
integer :: nva, nvb,cc1,cc2,cc3,cc4

open(unit=1, file='leitor.dat')
read(1,*)
read(1,*)
read(1,1)a
read(1,1)b
read(1,2)nva
read(1,2)nvb
read(1,1)k
read(1,1)q
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,2)cc1
read(1,2)cc2
read(1,2)cc3
read(1,2)cc4

1 format (32x,F10.4)
2 format (32x, I5)


end subroutine

!--------------------------------------------------------------------------
subroutine imprime (a,b,nva,nvb,k,q,cc1,cc2,cc3,cc4)
real :: a, b, k, q
integer :: nva, nvb,cc1,cc2,cc3,cc4

!print*, a,b,nva,nvb,k,q,cc1,cc2,cc3,cc4

end subroutine
!--------------------------------------------------------------------------
subroutine malha (a, b, nva, nvb, npx, npy, n, dym, dxm, msh, dx, dy)
real :: a, b, dx, dy, dxm, dym
integer :: nvx, nvy, npx, npy, n
real, dimension (100000,10) :: msh

dx=a/nva
dy=b/nvb
dxm=dx/2
dym=dy/2
npx=nva+2
npy=nvb+2

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

end subroutine
!--------------------------------------------------------------------------
subroutine distribuicao2d (npx,npy,dx,dy,msh,n,k,a,b,dxm,dym,q,an)
implicit none
real  dx, dy,dxm,dym, k, a, b, q
real, dimension (100000,10) :: msh, an
integer i,n,x,y,npx,npy

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

	else if (msh(i,1)>=a-1.e-6) then
	an(i,1)=0.
	an(i,2)=0.
        an(i,3)=0.
        an(i,4)=0.

	else if (msh(i,2)>=b-1.e-6) then
	an(i,1)=0.
	an(i,2)=0.
        an(i,3)=0.
        an(i,4)=0.

!parede oeste
	else if ((msh(i,1)>=dxm-1.e-6).and.(msh(i,1)<=dxm+1.e-6)) then
	an(i,1)=0.
	an(i,2)=k*(dy /(msh(i+1,1)-msh(i,1)))
        an(i,3)=k*(dx /(msh(i+npx,2)-msh(i,2)))
        an(i,4)=k*(dx /(msh(i,2)-msh(i-npx,2)))


!parede leste
	else if ((msh(i,1)>=((a-dxm)-1.e-6)).and.(msh(i,1)<=((a-dxm)+1.e-6)))  then
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



        an(i,5)=an(i,1)+an(i,2)+an(i,3)+an(i,4)			!Ap somatório


end do

end subroutine distribuicao2d

!--------------------------------------------------------------------------
subroutine calculo (npx,npy,dx,dy,msh,n,k,a,b,an,q,T,dxm)
implicit none
real  dx, dy,dxm,dym, k, q, a, b, T(1000000), To(1000000), erro(1000000), erro1, erro0
real, dimension (100000,10) :: msh, an
integer i,n,x,y,npx,npy
open (unit=2, file='temperatura.dat')

do i=1,n
	if (msh(i,2)==0) then
	T(i)= 400.
	else if ((msh(i,2)<=b+1.e-6).and.(msh(i,2)>=b-1.e-6)) then
	T(i)= 100.
	else
	T(i)=0
	endif
end do

do i=1,n
To(i)=1
erro(i)=100
end do

erro1=100

do while(erro1>=1.E-6)
erro0=0
  do i=1,n
	if (msh(i,2)==0) then
	erro(i)=0.
	 else if (msh(i,2)>=b-1.e-6) then
	 erro(i)=0.
	  else if (msh(i,1)==0) then
	  erro(i)=0.
	   else if (msh(i,1)>=a-1.e-6) then
	   erro(i)=0.
		else if ((msh(i,1)>=dxm-1.e-6).and.(msh(i,1)<=dxm+1.e-6)) then
		To(i)=T(i)
	        T(i)= (an(i,1)*T(i-1) + an(i,2)*T(i+1) + an(i,3)*T(i+npx) + an(i,4)*T(i-npx) + (q*dy))/ an(i,5)
		erro(i)=(abs((T(i)-To(i))/T(i)))*100
		else
		To(i)=T(i)
	        T(i)= (an(i,1)*T(i-1) + an(i,2)*T(i+1) + an(i,3)*T(i+npx) + an(i,4)*T(i-npx))/ an(i,5)
		erro(i)=(abs((T(i)-To(i))/T(i)))*100
	end if

	if (erro0>=erro(i)) then
		erro1=erro0
		else
		erro1=erro(i)
	end if
	erro0=erro1

	if (msh(i,1)==0) then

	T(i)=T(i+1)+(q*dy/k)

	else if (msh(i,1)>=a-1.e-6) then
 	T(i)=T(i-1)

	end if
 end do

end do


do i=1,n
print*, T(i),i
end do

end subroutine calculo

!--------------------------------------------------------------------------
subroutine prints (msh,T,a,npx,n)
implicit none
integer :: i,npx,n
real :: a, T(1000000)
real, dimension (100000,10) :: msh

open (unit=3 , file ="./gmsh.pos",status="unknown")


!==========================================================================================
!FORMATOS
1 FORMAT (A1,E12.5,A1,E12.5,A1,E12.5,A1,E12.5,A2)
2 FORMAT (A3,E12.5,A1,E12.5,A3,E12.5,A1,E12.5,A3)
3 FORMAT (E12.5,A1,E12.5,A3,E12.5,A1,E12.5,A3)
4 FORMAT (A12,A10,A14,A14,A10)

!========================================================================================
!POS-PROCESSAMENTO - GMSH
write(3,*)'View "Temperature" {'

do i=1,(n-npx)
	if ((msh(i,1)<=a+1.d-6).and.(msh(i,1)>=a-1.d-6)) then
		else
		write(3,2)'SQ(',msh(i+npx,1),',',msh(i+npx,2),',0,',msh(i,1),',',msh(i,2),',0,'
		write(3,3)msh(i+1,1),',',msh(i+1,2),',0,',msh(i+npx+1,1),',',msh(i+npx+1,2),',0)'
		write(3,1)'{',T(i+npx),',',T(i),',',T(i+1),',',T(i+npx+1),'};'
	end if
end do
write(3,*)'};'

!===================================================================================

return
end
