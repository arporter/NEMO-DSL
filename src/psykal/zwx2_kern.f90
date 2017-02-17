module zwx2_kern_mod
!
implicit none
!
contains
  !
  subroutine zwx2_kern(zwx,pwn,mydomain,zind,zslpx,jk,jj,ji)
    !
    real*8, intent(out) :: zwx(:,:,:)
    real*8, intent(in) :: pwn(:,:,:), mydomain(:,:,:), zind(:,:,:), zslpx(:,:,:)
    integer, intent(in) :: jk,jj,ji
    ! local variables
    real*8 :: z0w, zalpha, zw, zzwx, zzwy, zdt, zbtr
    !
    zdt=1
    zbtr=1.
    z0w = SIGN( 0.5d0, pwn(ji,jj,jk+1) )
    zalpha = 0.5d0 + z0w
    zw  = z0w - 0.5d0 * pwn(ji,jj,jk+1) * zdt * zbtr
    
    zzwx = mydomain(ji,jj,jk+1) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk+1))
    zzwy = mydomain(ji,jj,jk  ) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk  ))
    
    zwx(ji,jj,jk+1) = pwn(ji,jj,jk+1) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
    !
  end subroutine zwx2_kern
  !
end module zwx2_kern_mod
