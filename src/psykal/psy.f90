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
  subroutine multiply_layer(field_out,field_in1,field_in2,jpj,jpi)
    !
    real*8, intent(out) :: field_out(:,:)
    real*8, intent(in) :: field_in1(:,:), field_in2(:,:)
    integer, intent(in) :: jpj,jpi
    ! local variables
    integer :: jj,ji
    !
    DO jj=1,jpj
      DO ji=1,jpi
        field_out(ji,jj) = field_in1(ji,jj) * field_in2(ji,jj)
      END DO
    END DO
    !
  end subroutine multiply_layer
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
  subroutine zslpxy_update_psy(zslpx,zslpy,zwx,zwy,jpk,jpj,jpi)
    !
    use zslpxy_update_kern_mod, only : zslpxy_update_kern
    !
    real*8, intent(inout) :: zslpx(:,:,:), zslpy(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji   
    !
    DO jk = 1, jpk-1    
       DO jj = 2, jpj
          DO ji = 2, jpi
             call zslpxy_update_kern(zslpx,zslpy,zwx,zwy,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine
  !
  subroutine zwxy2_psy(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy,jpk,jpj,jpi)
    !
    use zwxy2_kern_mod, only : zwxy2_kern
    !
    real*8, intent(out) :: zwx(:,:,:), zwy(:,:,:)
    real*8, intent(in)  :: pun(:,:,:), pvn(:,:,:), mydomain(:,:,:), zind(:,:,:), zslpx(:,:,:), zslpy(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji   
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1
          DO ji = 2, jpi-1
             call zwxy2_kern(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine zwxy2_psy
  !
  subroutine mydomain_update_psy(mydomain,zwx,zwy,jpk,jpj,jpi)
    !
    use mydomain_update_kern_mod, only : mydomain_update_kern
    !
    real*8, intent(inout) :: mydomain(:,:,:)
    real*8, intent(in)    :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in)   :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1     
          DO ji = 2, jpi-1
             call mydomain_update_kern(mydomain,zwx,zwy,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine mydomain_update_psy
  !
  subroutine zwx_psy(zwx,tmask,mydomain,jpk,jpj,jpi)
    !
    use zwx_kern_mod, only : zwx_kern
    !
    real*8, intent(out) :: zwx(:,:,:)
    real*8, intent(in) :: tmask(:,:,:), mydomain(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji
    !
    DO jk = 2, jpk-1
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zwx_kern(zwx,tmask,mydomain,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine zwx_psy
  !
  subroutine zslpx_psy(zslpx,zwx,jpk,jpj,jpi)
    !
    use zslpx_kern_mod, only : zslpx_kern
    !
    real*8, intent(out) :: zslpx(:,:,:)
    real*8, intent(in) :: zwx(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji
    !
    DO jk = 2, jpk-1
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zslpx_kern(zslpx,zwx,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine zslpx_psy
  !
  subroutine zslpx_update_psy(zslpx,zwx,jpk,jpj,jpi)
    !
    use zslpx_update_kern_mod, only : zslpx_update_kern
    !
    real*8, intent(inout) :: zslpx(:,:,:)
    real*8, intent(in) :: zwx(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji
    !
    DO jk = 2, jpk-1     
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zslpx_update_kern(zslpx,zwx,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine zslpx_update_psy
  !
  subroutine zwx2_psy(zwx,pwn,mydomain,zind,zslpx,jpk,jpj,jpi)
    !
    use zwx2_kern_mod, only : zwx2_kern
    !
    real*8, intent(inout) :: zwx(:,:,:)
    real*8, intent(in) :: pwn(:,:,:), mydomain(:,:,:), zind(:,:,:), zslpx(:,:,:)
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    integer :: jk,jj,ji
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1
          DO ji = 2, jpi-1
             call zwx2_kern(zwx,pwn,mydomain,zind,zslpx,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine zwx2_psy
  !
  subroutine mydomain_psy(mydomain,zbtr,zwx,jpk,jpj,jpi)
    !
    use mydomain_kern_mod, only : mydomain_kern
    !
    real*8, intent(out) :: mydomain(:,:,:)
    real*8, intent(in) :: zwx(:,:,:)
    real*8, intent(in) :: zbtr
    integer, intent(in) :: jpk,jpj,jpi
    ! local variables
    real*8 :: ztra
    integer :: jk,jj,ji
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1     
          DO ji = 2, jpi-1
             call mydomain_kern(mydomain,zbtr,zwx,jk,jj,ji)
          END DO
       END DO
    END DO
    !
  end subroutine mydomain_psy
  !
end module psy_mod
