module psy_mod
!
implicit none
!
contains
  !
  subroutine zind_psy(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask,jpk,jpj,jpi)
    !
    use zind_kern_mod, only : zind_kern
    !
    real*8,  intent(out) :: zind(:,:,:)
    real*8,  intent(in)  :: tsn(:,:,:)
    real*8,  intent(in)  :: tmask(:,:,:)
    real*8,  intent(in)  :: ztfreez(:,:), rnfmsk(:,:), rnfmsk_z(:), upsmsk(:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,ji,jj
    !
    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zind_kern(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask,ji,jj,jk)
        END DO
      END DO
    END DO
    !
  end subroutine zind_psy
  !
  subroutine zwxy_psy(zwx,zwy,mydomain,umask,vmask,jpk,jpj,jpi)
    !
    use zwxy_kern_mod, only : zwxy_kern
    !
    real*8,  intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real*8,  intent(in)  :: mydomain(:,:,:)
    real*8,  intent(in)  :: umask(:,:,:), vmask(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji
    !
    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zwxy_kern(zwx,zwy,mydomain,umask,vmask,jpk,jk,jj,ji)
        END DO
      END DO
    END DO
    !
  end subroutine zwxy_psy
  !
end module psy_mod
