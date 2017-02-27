module zslpx_update_kern_mod
!
implicit none
!
type, extends(kernel_type) :: zslpx_update_type
  type(arg), dimension(2) :: meta_args = (/ &
    arg(READWRITE, 3D, CU), &
    arg(READ,      3D, CU, STENCIL(U)) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => zslpx_update_kern
end type zslpx_update_type
!        
contains
  !
  subroutine zslpx_update_kern(zslpx,zwx,ji,jj,jk)
    !
    real*8, intent(inout) :: zslpx(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN( ABS( zslpx(ji,jj,jk  ) ), &
         &                                               2.d0*ABS( zwx  (ji,jj,jk+1) ),   &
         &                                               2.d0*ABS( zwx  (ji,jj,jk  ) )  )
    !
  end subroutine zslpx_update_kern
  !
end module zslpx_update_kern_mod
