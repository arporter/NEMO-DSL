module zwx_kern_mod
!
implicit none
!
contains
  !
  subroutine zwx_kern(zwx,tmask,mydomain,jk,jj,ji)
    !
    real*8, intent(out) :: zwx(:,:,:)
    real*8, intent(in)  :: tmask(:,:,:), mydomain(:,:,:)
    integer, intent(in) :: jk,jj,ji
    !
    zwx(ji,jj,jk) = tmask(ji,jj,jk) * ( mydomain(ji,jj,jk-1) - mydomain(ji,jj,jk) )
    !
  end subroutine zwx_kern
  !
end module zwx_kern_mod
