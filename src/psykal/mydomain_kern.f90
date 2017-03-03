module mydomain_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: mydomain_type
!  type(arg), dimension(2) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CT), &
!    arg(READ,  DIMS(3), CU, STENCIL(U)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => mydomain_kern
!end type mydomain_type
!        
contains
  !
  subroutine mydomain_kern(mydomain,zwx,ji,jj,jk)
    !
    real(wp), intent(out) :: mydomain(:,:,:)
    real(wp), intent(in)  :: zwx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real(wp) :: ztra, zbtr
    !
    zbtr = 1.d0
    ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
    mydomain(ji,jj,jk) = ztra
    !
  end subroutine mydomain_kern
  !
end module mydomain_kern_mod
