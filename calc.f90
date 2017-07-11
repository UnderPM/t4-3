subroutine calc (npx,dy,msh,n,k,width,height,an,q,t,dxm)
    implicit none
    real :: dy,dxm, k, q, width, height, t(1000000), t0(1000000), error(1000000), err_1, err_0
    real, dimension (100000,10) :: msh, an
    integer i,n,npx

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
        error(i)=100
    end do

    err_1=100

    do while(err_1>=1.E-6)
        err_0=0
        do i=1,n
            if (msh(i,2)==0) then
                error(i)=0.
            else if (msh(i,2)>=height-1.e-6) then
                error(i)=0.
            else if (msh(i,1)==0) then
                error(i)=0.
            else if (msh(i,1)>=width-1.e-6) then
                error(i)=0.
            else if ((msh(i,1)>=dxm-1.e-6).and.(msh(i,1)<=dxm+1.e-6)) then
                t0(i)=t(i)
                t(i)= (an(i,1)*t(i-1) + an(i,2)*t(i+1) + an(i,3)*t(i+npx) + an(i,4)*t(i-npx) + (q*dy))/ an(i,5)
                error(i)=(abs((t(i)-t0(i))/t(i)))*100
            else
                t0(i)=t(i)
                t(i)= (an(i,1)*t(i-1) + an(i,2)*t(i+1) + an(i,3)*t(i+npx) + an(i,4)*t(i-npx))/ an(i,5)
                error(i)=(abs((t(i)-t0(i))/t(i)))*100
            end if

            if (err_0>=error(i)) then
                err_1=err_0
            else
                err_1=error(i)
            end if

            err_0=err_1

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

end subroutine
