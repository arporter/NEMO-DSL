module zwxy2_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zwxy2_type
!  type(arg), dimension(8) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CU), &
!    arg(WRITE, DIMS(3), CV), &
!    arg(READ,  DIMS(3), CU), &
!    arg(READ,  DIMS(3), CV), &
!    arg(READ,  DIMS(3), CT, STENCIL(NE)), &
!    arg(READ,  DIMS(3), CT), &
!    arg(READ,  DIMS(3), CU, STENCIL(E)), &
!    arg(READ,  DIMS(3), CV, STENCIL(N)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zwxy2_kern
!end type zwxy2_type
!
contains
  !
  subroutine zwxy2_kern(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy,ji,jj,jk)
    !
    real(wp), intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real(wp), intent(in)  :: pun(:,:,:), pvn(:,:,:), mydomain(:,:,:), zind(:,:,:), zslpx(:,:,:), zslpy(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real(wp) :: z0u, zalpha, zu, zzwx, zzwy, z0v, zv, zdt
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
