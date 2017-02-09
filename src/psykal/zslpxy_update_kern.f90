module zslpxy_update_kern_mod
!
implicit none
!
contains
  !
  subroutine zslpxy_update_kern(zslpx,zslpy,zwx,zwy,jk,jj,ji)
    !
    real*8, intent(inout) :: zslpx(:,:,:), zslpy(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: jk,jj,ji
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
