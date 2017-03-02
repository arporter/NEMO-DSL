module zslpxy_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zslpxy_type
!  type(arg), dimension(4) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CU), &
!    arg(WRITE, DIMS(3), CV), &
!    arg(READ,  DIMS(3), CU, STENCIL(W)), &
!    arg(READ,  DIMS(3), CV, STENCIL(S)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zslpxy_kern
!end type zslpxy_type
!
contains
  !
  subroutine zslpxy_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
    !
    real(wp), intent(out) :: zslpx(:,:,:), zslpy(:,:,:)
    real(wp), intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji-1,jj  ,jk) )   &
         &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji-1,jj  ,jk) ) )
    zslpy(ji,jj,jk) =                    ( zwy(ji,jj,jk) + zwy(ji  ,jj-1,jk) )   &
         &            * ( 0.25d0 + SIGN( 0.25d0, zwy(ji,jj,jk) * zwy(ji  ,jj-1,jk) ) )
    !
  end subroutine zslpxy_kern
  !
end module zslpxy_kern_mod
