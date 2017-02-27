module zslpx_kern_mod
!
implicit none
!
type, extends(kernel_type) :: zslpx_type
  type(arg), dimension(2) :: meta_args = (/ &
    arg(WRITE, 3D, CU), &
    arg(READ,  3D, CU, STENCIL(U)) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => zslpx_kern
end type zslpx_type
!        
contains
  !
  subroutine zslpx_kern(zslpx,zwx,ji,jj,jk)
    !
    real*8, intent(out) :: zslpx(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:)
    integer, intent(in) :: ji,jj,jk
    !
    zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
         &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )

    !
  end subroutine zslpx_kern
  !
end module zslpx_kern_mod
