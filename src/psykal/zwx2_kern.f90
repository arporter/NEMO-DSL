module zwx2_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zslpx_type
!  type(arg), dimension(5) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CU), &
!    arg(READ,  DIMS(3), CT), &
!    arg(READ,  DIMS(3), CT, STENCIL(D)), &
!    arg(READ,  DIMS(3), CT, STENCIL(D)), &
!    arg(READ,  DIMS(3), CU, STENCIL(D)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zslpx_kern
!end type zslpx_type
!        
contains
  !
  subroutine zwx2_kern(zwx,pwn,mydomain,zind,zslpx,ji,jj,jk)
    !
    real(wp), intent(out) :: zwx(:,:,:)
    real(wp), intent(in) :: pwn(:,:,:), mydomain(:,:,:), zind(:,:,:), zslpx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real(wp) :: z0w, zalpha, zw, zzwx, zzwy, zdt, zbtr
    !
    zdt = 1.d0
    zbtr = 1.d0
    z0w = SIGN( 0.5d0, pwn(ji,jj,jk) )
    zalpha = 0.5d0 + z0w
    zw  = z0w - 0.5d0 * pwn(ji,jj,jk) * zdt * zbtr
    
    zzwx = mydomain(ji,jj,jk  ) + zind(ji,jj,jk-1) * (zw * zslpx(ji,jj,jk  ))
    zzwy = mydomain(ji,jj,jk-1) + zind(ji,jj,jk-1) * (zw * zslpx(ji,jj,jk-1))
    
    zwx(ji,jj,jk) = pwn(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
    !
  end subroutine zwx2_kern
  !
end module zwx2_kern_mod
