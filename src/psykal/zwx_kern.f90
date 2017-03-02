module zwx_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zwx_type
!  type(arg), dimension(3) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CU), &
!    arg(READ,  TMASK_3D), &
!    arg(READ,  DIMS(3), CT, STENCIL(D)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zwx_kern
!end type zwx_type
!
contains
  !
  subroutine zwx_kern(zwx,tmask,mydomain,ji,jj,jk)
    !
    real(wp), intent(out) :: zwx(:,:,:)
    real(wp), intent(in)  :: tmask(:,:,:), mydomain(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zwx(ji,jj,jk) = tmask(ji,jj,jk) * ( mydomain(ji,jj,jk-1) - mydomain(ji,jj,jk) )
    !
  end subroutine zwx_kern
  !
end module zwx_kern_mod
