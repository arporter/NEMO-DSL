module mydomain_update_kern_mod
!
implicit none
!
type, extends(kernel_type) :: mydomain_update_type
  type(arg), dimension(3) :: meta_args = (/ &
    arg(WRITE, 3D, CT), &
    arg(READ,  3D, CU, STENCIL(W)), &
    arg(READ,  3D, CV, STENCIL(S)) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => mydomain_update_kern
end type mydomain_update_type
!
contains
  !
  subroutine mydomain_update_kern(mydomain,zwx,zwy,ji,jj,jk)
    !
    real*8,  intent(inout) :: mydomain(:,:,:)
    real*8,  intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real*8 :: zbtr,ztra
    !
    zbtr = 1.d0
    ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
         &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
    mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
    !
  end subroutine mydomain_update_kern
  !
end module mydomain_update_kern_mod
