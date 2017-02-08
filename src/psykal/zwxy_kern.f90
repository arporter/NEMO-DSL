module zwxy_kern_mod
!
implicit none
!
contains
  !
  subroutine zwxy_kern(zwx,zwy,mydomain,umask,vmask,jk,jj,ji)
    !
    real*8,  intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real*8,  intent(in)  :: mydomain(:,:,:)
    real*8,  intent(in)  :: umask(:,:,:), vmask(:,:,:)
    integer, intent(in) :: jk,jj,ji
    !
    zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
    zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
    !
  end subroutine zwxy_kern
  !
end module zwxy_kern_mod
