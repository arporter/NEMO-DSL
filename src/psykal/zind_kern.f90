module zind_kern_mod
!
use kind_params_mod, only : wp
!
implicit none
!
!type, extends(kernel_type) :: zind_type
!  type(arg), dimension(7) :: meta_args = (/ &
!    arg(WRITE, DIMS(3), CT), &
!    arg(READ,  DIMS(3), CT), &
!    arg(READ,  DIMS(2), CT), &
!    arg(READ,  DIMS(2), CT), &
!    arg(READ,  DIMS(1)    ), &
!    arg(READ,  DIMS(2), CT), &
!    arg(READ,  GRID_MASK_T_3D) /)
!  integer :: ITERATES_OVER = DIMS(3)
!  integer :: INDEX_OFFSET = OFFSET_NE
!contains
!  procedure, nopass :: code => zind_kern
!end type zind_type
!        
contains
  !
  subroutine zind_kern(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask,ji,jj,jk)
    !
    real(wp),  intent(out) :: zind(:,:,:)
    real(wp),  intent(in)  :: tsn(:,:,:)
    real(wp),  intent(in)  :: tmask(:,:,:)
    real(wp),  intent(in)  :: ztfreez(:,:), rnfmsk(:,:), rnfmsk_z(:), upsmsk(:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real(wp) :: zice
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
