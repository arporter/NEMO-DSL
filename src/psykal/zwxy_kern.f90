module zwxy_kern_mod
!
implicit none
!
contains
  !
  subroutine zwxy_kern(zwx,zwy,mydomain,umask,vmask,jpk,jk,jj,ji)
    !
    real*8,  intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real*8,  intent(in)  :: mydomain(:,:,:)
    real*8,  intent(in)  :: umask(:,:,:), vmask(:,:,:)
    integer, intent(in) :: jk,jj,ji,jpk
    !
    if (jk.eq.jpk) then
       zwx(ji,jj,jk) = 0.e0
       zwy(ji,jj,jk) = 0.e0
    else
       zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
       zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
    end if
    !
  end subroutine zwxy_kern
  !
end module zwxy_kern_mod
