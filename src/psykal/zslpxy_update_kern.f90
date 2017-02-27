module zslpxy_update_kern_mod
!
implicit none
!
type, extends(kernel_type) :: zslpxy_update_type
  type(arg), dimension(4) :: meta_args = (/ &
    arg(READWRITE, 3D, CU), &
    arg(READWRITE, 3D, CV), &
    arg(READ,      3D, CU, STENCIL(W)), &
    arg(READ,      3D, CV, STENCIL(S)) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => zslpxy_update_kern
end type zslpxy_update_type
!
contains
  !
  subroutine zslpxy_update_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
    !
    real*8, intent(inout) :: zslpx(:,:,:), zslpy(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:), zwy(:,:,:)
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
