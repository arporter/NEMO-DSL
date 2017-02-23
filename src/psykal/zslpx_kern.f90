module zslpx_kern_mod
!
implicit none
!
contains
  !
  subroutine zslpx_kern(zslpx,zwx,ji,jj,jk)
    !
    real*8, intent(out) :: zslpx(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
         &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )

    !
  end subroutine zslpx_kern
  !
end module zslpx_kern_mod
