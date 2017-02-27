module zwxy_kern_mod
!
implicit none
!
type, extends(kernel_type) :: zwxy_type
  type(arg), dimension(5) :: meta_args = (/ &
    arg(WRITE, 3D, CU), &
    arg(WRITE, 3D, CV), &
    arg(READ,  3D, CT, STENCIL(NE)), &
    arg(READ,  3D, CU), &
    arg(READ,  3D, CV) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => zwxy_kern
end type zwxy_type
!
contains
  !
  subroutine zwxy_kern(zwx,zwy,mydomain,umask,vmask,ji,jj,jk)
    !
    real*8,  intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real*8,  intent(in)  :: mydomain(:,:,:)
    real*8,  intent(in)  :: umask(:,:,:), vmask(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
    zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
    !
  end subroutine zwxy_kern
  !
end module zwxy_kern_mod
