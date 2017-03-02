module mydomain_update_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: mydomain_update_type
!  type(arg), dimension(3) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CT), &
!    arg(READ,  DIMS(3), CU, STENCIL(W)), &
!    arg(READ,  DIMS(3), CV, STENCIL(S)) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => mydomain_update_kern
!end type mydomain_update_type
!
contains
  !
  subroutine mydomain_update_kern(mydomain,zwx,zwy,ji,jj,jk)
    !
    real(wp),  intent(inout) :: mydomain(:,:,:)
    real(wp),  intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real(wp) :: zbtr,ztra
    !
    zbtr = 1.d0
    ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
         &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
    mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
    !
  end subroutine mydomain_update_kern
  !
end module mydomain_update_kern_mod
