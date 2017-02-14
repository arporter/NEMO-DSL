module zwxy2_kern_mod
!
implicit none
!
contains
  !
  subroutine zwxy2_kern(zwx,zwy,zdt,pun,pvn,mydomain,zind,zslpx,zslpy,jk,jj,ji)
    !
    real*8, intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real*8, intent(in)  :: pun(:,:,:), pvn(:,:,:), mydomain(:,:,:), zind(:,:,:), zslpx(:,:,:), zslpy(:,:,:)
    real*8, intent(in)  :: zdt
    integer, intent(in) :: jk,jj,ji
    ! local variables
    real*8 :: z0u, zalpha, zu, zzwx, zzwy, z0v, zv
    !
    z0u = SIGN( 0.5d0, pun(ji,jj,jk) )
    zalpha = 0.5d0 - z0u
    zu  = z0u - 0.5d0 * pun(ji,jj,jk) * zdt

    zzwx = mydomain(ji+1,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji+1,jj,jk))
    zzwy = mydomain(ji  ,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji  ,jj,jk))

    zwx(ji,jj,jk) = pun(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
    
    z0v = SIGN( 0.5d0, pvn(ji,jj,jk) )
    zalpha = 0.5d0 - z0v
    zv  = z0v - 0.5d0 * pvn(ji,jj,jk) * zdt

    zzwx = mydomain(ji,jj+1,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj+1,jk))
    zzwy = mydomain(ji,jj  ,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj  ,jk))

    zwy(ji,jj,jk) = pvn(ji,jj,jk) * ( zalpha * zzwx + (1.d0-zalpha) * zzwy )
    !
  end subroutine zwxy2_kern
!
end module zwxy2_kern_mod
