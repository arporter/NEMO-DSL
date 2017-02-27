module zwxy2_kern_mod
!
implicit none
!
type, extends(kernel_type) :: zwxy2_type
  type(arg), dimension(8) :: meta_args = (/ &
    arg(WRITE, 3D, CU), &
    arg(WRITE, 3D, CV), &
    arg(READ,  3D, CU), &
    arg(READ,  3D, CV), &
    arg(READ,  3D, CT, STENCIL(NE)), &
    arg(READ,  3D, CT), &
    arg(READ,  3D, CU, STENCIL(E)), &
    arg(READ,  3D, CV, STENCIL(N)) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => zwxy2_kern
end type zwxy2_type
!
contains
  !
  subroutine zwxy2_kern(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy,ji,jj,jk)
    !
    real*8, intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real*8, intent(in)  :: pun(:,:,:), pvn(:,:,:), mydomain(:,:,:), zind(:,:,:), zslpx(:,:,:), zslpy(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real*8 :: z0u, zalpha, zu, zzwx, zzwy, z0v, zv, zdt
    !
    zdt = 1.d0
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
