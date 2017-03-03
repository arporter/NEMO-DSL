module zslpxy_update_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zslpxy_update_type
!  type(arg), dimension(4) :: meta_args = (/ &
!    arg(READWRITE, DIMS(3), CU), &
!    arg(READWRITE, DIMS(3), CV), &
!    arg(READ,      DIMS(3), CU, STENCIL(W)), &
!    arg(READ,      DIMS(3), CV, STENCIL(S)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zslpxy_update_kern
!end type zslpxy_update_type
!
contains
  !
  subroutine zslpxy_update_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
    !
    real(wp), intent(inout) :: zslpx(:,:,:), zslpy(:,:,:)
    real(wp), intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji  ,jj,jk) ),   &
         &                                                2.d0*ABS( zwx  (ji-1,jj,jk) ),   &
         &                                                2.d0*ABS( zwx  (ji  ,jj,jk) ) )
    zslpy(ji,jj,jk) = SIGN( 1.d0, zslpy(ji,jj,jk) ) * MIN(    ABS( zslpy(ji,jj  ,jk) ),   &
         &                                                2.d0*ABS( zwy  (ji,jj-1,jk) ),   &
         &                                                2.d0*ABS( zwy  (ji,jj  ,jk) ) )
    !
  end subroutine zslpxy_update_kern
  !
end module zslpxy_update_kern_mod
