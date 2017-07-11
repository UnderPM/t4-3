program t4

    implicit none
    real :: width,height,k,q,dym,dxm,dx,dy
    integer :: nvx,nvy,npx,npy,n,bc_e,bc_w,bc_n,bc_s
    real, dimension (100000,10) :: msh,an
    real, dimension (1000000) :: t

    call reader (width,height,nvx,nvy,k,q,bc_e,bc_w,bc_n,bc_s)
    call mesh (width, height, nvx, nvy, npx, npy, n, dym, dxm, msh, dx, dy,an,k)
    call calc (npx,dy,msh,n,k,width,height,an,q,t,dxm)
    call gmsh (msh,t,width,npx,n)

end program
