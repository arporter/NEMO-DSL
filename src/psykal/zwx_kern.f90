module zwx_kern_mod
!
implicit none
!
contains
  !
  subroutine zwx_kern(zwx,tmask,mydomain,ji,jj,jk)
    !
    real*8, intent(out) :: zwx(:,:,:)
    real*8, intent(in)  :: tmask(:,:,:), mydomain(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zwx(ji,jj,jk) = tmask(ji,jj,jk) * ( mydomain(ji,jj,jk-1) - mydomain(ji,jj,jk) )
    !
  end subroutine zwx_kern
  !
end module zwx_kern_mod
