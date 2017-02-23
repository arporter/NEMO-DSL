module psy_mod
!
implicit none
!
private
public :: set_bounds
public :: tracer_advection
!
integer :: jpi, jpj, jpk
!
contains
  !
  subroutine set_bounds(jpi_in,jpj_in,jpk_in)
    !
    integer, intent(in) :: jpi_in,jpj_in,jpk_in
    !
    jpi=jpi_in
    jpj=jpj_in
    jpk=jpk_in
    !
  end subroutine set_bounds
  !
  subroutine tracer_advection(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk, &
       tmask,zwx,zwy,mydomain,umask,vmask,zslpx,zslpy,pun,pvn,pwn)
    !
    use zind_kern_mod, only : zind_kern
    use zwxy_kern_mod, only : zwxy_kern
    use zslpxy_kern_mod, only : zslpxy_kern
    use zslpxy_update_kern_mod, only : zslpxy_update_kern
    !
    ! RF all intent out's here are actually local to this call
    real*8, intent(out) :: zind(jpi,jpj,jpk), zwx(jpi,jpj,jpk), &
         zwy(jpi,jpj,jpk), zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    real*8, intent(in)  :: tsn(jpi,jpj,jpk), ztfreez(jpi,jpj), &
         rnfmsk(jpi,jpj), rnfmsk_z(jpk), upsmsk(jpi,jpj), tmask(jpi,jpj,jpk), &
         umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk), pun(jpi,jpj,jpk), &
         pvn(jpi,jpj,jpk), pwn(jpi,jpj,jpk)
    real*8, intent(inout) :: mydomain(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    !call zind_psy(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask)
    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zind_kern(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask,ji,jj,jk)
        END DO
      END DO
    END DO
    !call zero_bottom_layer(zwx)
    DO jj=1,jpj
      DO ji=1,jpi
        zwx(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !call zero_bottom_layer(zwy)
    DO jj=1,jpj
      DO ji=1,jpi
        zwy(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !call zwxy_psy(zwx,zwy,mydomain,umask,vmask)
    DO jk = 1, jpk-1
      DO jj = 1, jpj-1
        DO ji = 1, jpi-1
          call zwxy_kern(zwx,zwy,mydomain,umask,vmask,ji,jj,jk)
        END DO
      END DO
    END DO
    !call zero_bottom_layer(zslpx)
    DO jj=1,jpj
      DO ji=1,jpi
        zslpx(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !call zero_bottom_layer(zslpy)
    DO jj=1,jpj
      DO ji=1,jpi
        zslpy(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !call zslpxy_psy(zslpx,zslpy,zwx,zwy)
    DO jk = 1, jpk-1
      DO jj = 2, jpj
        DO ji = 2, jpi
          call zslpxy_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
        END DO
      END DO
    END DO
    !call zslpxy_update_psy(zslpx,zslpy,zwx,zwy)
    DO jk = 1, jpk-1
       DO jj = 2, jpj
          DO ji = 2, jpi
             call zslpxy_update_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
          END DO
       END DO
    END DO
    call zwxy2_psy(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy)
    call mydomain_update_psy(mydomain,zwx,zwy)
    call zero_top_layer(zwx)
    call zero_bottom_layer(zwx)
    call zwx_psy(zwx,tmask,mydomain)
    call zero_top_layer(zslpx)
    call zslpx_psy(zslpx,zwx)
    call zslpx_update_psy(zslpx,zwx)
    call multiply_top_layer(zwx,pwn,mydomain)
    call zwx2_psy(zwx,pwn,mydomain,zind,zslpx)
    call mydomain_psy(mydomain,zwx)
    !
  end subroutine tracer_advection
  
  subroutine zero_top_layer(field)
    !
    real*8, intent(inout) :: field(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj
    !
    DO jj=1,jpj
      DO ji=1,jpi
        field(ji,jj,1) = 0.d0
      END DO
    END DO
    !
  end subroutine zero_top_layer
  !
  subroutine zero_bottom_layer(field)
    !
    real*8, intent(inout) :: field(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj
    !
    DO jj=1,jpj
      DO ji=1,jpi
        field(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
  end subroutine zero_bottom_layer
  !
  subroutine multiply_top_layer(field_out,field_in1,field_in2)
    !
    real*8, intent(out) :: field_out(jpi,jpj,jpk)
    real*8, intent(in) :: field_in1(jpi,jpj,jpk), field_in2(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj
    !
    DO jj=1,jpj
      DO ji=1,jpi
        field_out(ji,jj,1) = field_in1(ji,jj,1) * field_in2(ji,jj,1)
      END DO
    END DO
    !
  end subroutine multiply_top_layer
  !
  subroutine zind_psy(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask)
    !
    use zind_kern_mod, only : zind_kern
    !
    real*8,  intent(out) :: zind(jpi,jpj,jpk)
    real*8,  intent(in)  :: tsn(jpi,jpj,jpk)
    real*8,  intent(in)  :: tmask(jpi,jpj,jpk)
    real*8,  intent(in)  :: ztfreez(jpi,jpj), rnfmsk(jpi,jpj), rnfmsk_z(jpk), upsmsk(jpi,jpj)
    ! local variables
    integer :: ji,jj,jk
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
  subroutine zwxy_psy(zwx,zwy,mydomain,umask,vmask)
    !
    use zwxy_kern_mod, only : zwxy_kern
    !
    real*8,  intent(out) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    real*8,  intent(in)  :: mydomain(jpi,jpj,jpk)
    real*8,  intent(in)  :: umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 1, jpk-1
      DO jj = 1, jpj-1
        DO ji = 1, jpi-1
          call zwxy_kern(zwx,zwy,mydomain,umask,vmask,ji,jj,jk)
        END DO
      END DO
    END DO
    !
  end subroutine zwxy_psy
  !
  subroutine zslpxy_psy(zslpx,zslpy,zwx,zwy)
    !
    use zslpxy_kern_mod, only : zslpxy_kern
    !
    real*8, intent(out) :: zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj
          DO ji = 2, jpi
             call zslpxy_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine zslpxy_psy
  !
  subroutine zslpxy_update_psy(zslpx,zslpy,zwx,zwy)
    !
    use zslpxy_update_kern_mod, only : zslpxy_update_kern
    !
    real*8, intent(inout) :: zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj
          DO ji = 2, jpi
             call zslpxy_update_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine
  !
  subroutine zwxy2_psy(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy)
    !
    use zwxy2_kern_mod, only : zwxy2_kern
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    real*8, intent(in)  :: pun(jpi,jpj,jpk), pvn(jpi,jpj,jpk), &
                           mydomain(jpi,jpj,jpk), zind(jpi,jpj,jpk), &
                           zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1
          DO ji = 2, jpi-1
             call zwxy2_kern(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine zwxy2_psy
  !
  subroutine mydomain_update_psy(mydomain,zwx,zwy)
    !
    use mydomain_update_kern_mod, only : mydomain_update_kern
    !
    real*8, intent(inout) :: mydomain(jpi,jpj,jpk)
    real*8, intent(in)    :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1
          DO ji = 2, jpi-1
             call mydomain_update_kern(mydomain,zwx,zwy,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine mydomain_update_psy
  !
  subroutine zwx_psy(zwx,tmask,mydomain)
    !
    use zwx_kern_mod, only : zwx_kern
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk)
    real*8, intent(in) :: tmask(jpi,jpj,jpk), mydomain(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 2, jpk-1
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zwx_kern(zwx,tmask,mydomain,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine zwx_psy
  !
  subroutine zslpx_psy(zslpx,zwx)
    !
    use zslpx_kern_mod, only : zslpx_kern
    !
    real*8, intent(out) :: zslpx(jpi,jpj,jpk)
    real*8, intent(in) :: zwx(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 2, jpk-1
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zslpx_kern(zslpx,zwx,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine zslpx_psy
  !
  subroutine zslpx_update_psy(zslpx,zwx)
    !
    use zslpx_update_kern_mod, only : zslpx_update_kern
    !
    real*8, intent(inout) :: zslpx(jpi,jpj,jpk)
    real*8, intent(in) :: zwx(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 2, jpk-1     
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zslpx_update_kern(zslpx,zwx,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine zslpx_update_psy
  !
  subroutine zwx2_psy(zwx,pwn,mydomain,zind,zslpx)
    !
    use zwx2_kern_mod, only : zwx2_kern
    !
    real*8, intent(inout) :: zwx(jpi,jpj,jpk)
    real*8, intent(in) :: pwn(jpi,jpj,jpk), mydomain(jpi,jpj,jpk), &
                          zind(jpi,jpj,jpk), zslpx(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 2, jpk
       DO jj = 2, jpj-1
          DO ji = 2, jpi-1
             call zwx2_kern(zwx,pwn,mydomain,zind,zslpx,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine zwx2_psy
  !
  subroutine mydomain_psy(mydomain,zwx)
    !
    use mydomain_kern_mod, only : mydomain_kern
    !
    real*8, intent(out) :: mydomain(jpi,jpj,jpk)
    real*8, intent(in) :: zwx(jpi,jpj,jpk)
    ! local variables
    integer :: ji,jj,jk
    !
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1     
          DO ji = 2, jpi-1
             call mydomain_kern(mydomain,zwx,ji,jj,jk)
          END DO
       END DO
    END DO
    !
  end subroutine mydomain_psy
  !
end module psy_mod
