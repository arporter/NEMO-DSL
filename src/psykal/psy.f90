module psy_mod
!
implicit none
!
contains
  !
  subroutine zero_layer(field,jpj,jpi)
    !
    real*8, intent(out) :: field(:,:)
    integer, intent(in) :: jpj,jpi
    ! local variables
    integer :: jj,ji
    !
    DO jj=1,jpj
      DO ji=1,jpi
        field(ji,jj) = 0.e0
      END DO
    END DO
    !
  end subroutine zero_layer
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
    DO jk = 1, jpk-1
      DO jj = 1, jpj-1
        DO ji = 1, jpi-1
          call zwxy_kern(zwx,zwy,mydomain,umask,vmask,jk,jj,ji)
        END DO
      END DO
    END DO
    !
  end subroutine zwxy_psy
  !
  subroutine zslpxy_psy(zslpx,zslpy,zwx,zwy,jpk,jpj,jpi)
    !
    use zslpxy_kern_mod, only : zslpxy_kern
    !
    real*8, intent(out) :: zslpx(:,:,:), zslpy(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj
          DO ji = 2, jpi
             call zslpxy_kern(zslpx,zslpy,zwx,zwy,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine zslpxy_psy
  !
end module psy_mod
