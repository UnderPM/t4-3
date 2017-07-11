subroutine gmsh (msh,t,width,npx,n)
    implicit none
    integer :: i,npx,n
    real :: width, t(1000000)
    real, dimension (100000,10) :: msh

    open (unit=3 , file ="./gmsh.pos")

    write(3,*)'View "Temperature" {'

    do i=1,(n-npx)
        if ((msh(i,1)<=width+1.d-6).and.(msh(i,1)>=width-1.d-6)) then
            !continue
        else
            write(3,*)'SQ(',msh(i+npx,1),',',msh(i+npx,2),',0,',msh(i,1),',',msh(i,2),',0,'
            write(3,*)msh(i+1,1),',',msh(i+1,2),',0,',msh(i+npx+1,1),',',msh(i+npx+1,2),',0)'
            write(3,*)'{',t(i+npx),',',t(i),',',t(i+1),',',t(i+npx+1),'};'
        end if
    end do
    write(3,*)'};'
    !return
end subroutine
