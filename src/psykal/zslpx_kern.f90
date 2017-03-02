module zslpx_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zslpx_type
!  type(arg), dimension(2) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CU), &
!    arg(READ,  DIMS(3), CU, STENCIL(U)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zslpx_kern
!end type zslpx_type
!        
contains
  !
  subroutine zslpx_kern(zslpx,zwx,ji,jj,jk)
    !
    real(wp), intent(out) :: zslpx(:,:,:)
    real(wp), intent(in)  :: zwx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
         &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )

    !
  end subroutine zslpx_kern
  !
end module zslpx_kern_mod
