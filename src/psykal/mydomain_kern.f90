module mydomain_kern_mod
!
implicit none
!
type, extends(kernel_type) :: mydomain_type
  type(arg), dimension(2) :: meta_args = (/ &
    arg(WRITE, 3D, CT), &
    arg(READ,  3D, CU, STENCIL(U)) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => mydomain_kern
end type mydomain_type
!        
contains
  !
  subroutine mydomain_kern(mydomain,zwx,ji,jj,jk)
    !
    real*8, intent(out) :: mydomain(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real*8 :: ztra, zbtr
    !
    zbtr = 1.d0
    ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
    mydomain(ji,jj,jk) = ztra
    !
  end subroutine mydomain_kern
  !
end module mydomain_kern_mod
