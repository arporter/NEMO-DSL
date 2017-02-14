module zind_kern_mod
!
implicit none
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
