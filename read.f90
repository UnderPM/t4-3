subroutine reader (width,height,nvx,nvy,k,q,bc_e,bc_w,bc_n,bc_s)
    implicit none
    real :: width, height, k, q
    integer :: nvx, nvy,bc_e,bc_w,bc_n,bc_s

    open(unit=1, file='dados.dat')
    read(1,*)width
    read(1,*)height
    read(1,*)nvx
    read(1,*)nvy
    read(1,*)k
    read(1,*)q
    read(1,*)bc_e
    read(1,*)bc_w
    read(1,*)bc_n
    read(1,*)bc_s

end subroutine

