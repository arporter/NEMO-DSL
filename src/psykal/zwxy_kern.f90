module zwxy_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zwxy_type
!  type(arg), dimension(5) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CU), &
!    arg(WRITE, DIMS(3), CV), &
!    arg(READ,  DIMS(3), CT, STENCIL(NE)), &
!    arg(READ,  DIMS(3), CU), &
!    arg(READ,  DIMS(3), CV) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zwxy_kern
!end type zwxy_type
!
contains
  !
  subroutine zwxy_kern(zwx,zwy,mydomain,umask,vmask,ji,jj,jk)
    !
    real(wp),  intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real(wp),  intent(in)  :: mydomain(:,:,:)
    real(wp),  intent(in)  :: umask(:,:,:), vmask(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
    zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
    !
  end subroutine zwxy_kern
  !
end module zwxy_kern_mod
