module zind_kern_mod
!
implicit none
!
type, extends(kernel_type) :: zind_type
  type(arg), dimension(7) :: meta_args = (/ &
    arg(WRITE, 3D, CT), &
    arg(READ,  3D, CT), &
    arg(READ,  2D, CT), &
    arg(READ,  2D, CT), &
    arg(READ,  1D    ), &
    arg(READ,  2D, CT), &
    arg(READ,  TMASK_3D) /)
  integer :: ITERATES_OVER = 3D
  integer :: INDEX_OFFSET = OFFSET_NE
contains
  procedure, nopass :: code => zind_kern
end type zind_type
!        
contains
  !
  subroutine zind_kern(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask,ji,jj,jk)
    !
    real*8,  intent(out) :: zind(:,:,:)
    real*8,  intent(in)  :: tsn(:,:,:)
    real*8,  intent(in)  :: tmask(:,:,:)
    real*8,  intent(in)  :: ztfreez(:,:), rnfmsk(:,:), rnfmsk_z(:), upsmsk(:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real*8 :: zice
    !
    IF( tsn(ji,jj,jk) <= ztfreez(ji,jj) + 0.1d0 ) THEN   ;   zice = 1.d0
    ELSE                                                 ;   zice = 0.d0
    ENDIF
    zind(ji,jj,jk) = MAX (   &
        rnfmsk(ji,jj) * rnfmsk_z(jk),      & 
        upsmsk(ji,jj)               ,      &
        zice                               &
        &                  ) * tmask(ji,jj,jk)
    zind(ji,jj,jk) = 1 - zind(ji,jj,jk)
    !
  end subroutine zind_kern
  !
end module zind_kern_mod
