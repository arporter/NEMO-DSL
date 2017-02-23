module zslpx_update_kern_mod
!
implicit none
!
contains
  !
  subroutine zslpx_update_kern(zslpx,zwx,ji,jj,jk)
    !
    real*8, intent(inout) :: zslpx(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN( ABS( zslpx(ji,jj,jk  ) ), &
         &                                               2.d0*ABS( zwx  (ji,jj,jk+1) ),   &
         &                                               2.d0*ABS( zwx  (ji,jj,jk  ) )  )
    !
  end subroutine zslpx_update_kern
  !
end module zslpx_update_kern_mod
