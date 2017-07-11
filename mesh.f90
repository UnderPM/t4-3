subroutine mesh (width, height, nvx, nvy, npx, npy, n, dym, dxm, msh, dx, dy,an,k)
    implicit none
    real :: width, height, dx, dy, dxm, dym,k,x,y
    integer :: nvx, nvy, npx, npy, n,i,j
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
            an(i,1)=0
            an(i,2)=0
            an(i,3)=0
            an(i,4)=0
        else if (msh(i,2)==0) then
            an(i,1)=0
            an(i,2)=0
            an(i,3)=0
            an(i,4)=0
        else if (msh(i,1)>=width-1.e-6) then
            an(i,1)=0
            an(i,2)=0
            an(i,3)=0
            an(i,4)=0
        else if (msh(i,2)>=height-1.e-6) then
            an(i,1)=0
            an(i,2)=0
            an(i,3)=0
            an(i,4)=0
        else if ((msh(i,1)>=dxm-1.e-6).and.(msh(i,1)<=dxm+1.e-6)) then
            an(i,1)=0
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
