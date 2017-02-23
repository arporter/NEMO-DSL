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
    use zwxy2_kern_mod, only : zwxy2_kern
    use mydomain_update_kern_mod, only : mydomain_update_kern
    use zwx_kern_mod, only : zwx_kern
    use zslpx_kern_mod, only : zslpx_kern
    use zslpx_update_kern_mod, only : zslpx_update_kern
    use zwx2_kern_mod, only : zwx2_kern
    use mydomain_kern_mod, only : mydomain_kern
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
    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zind_kern(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwy(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 1, jpj-1
        DO ji = 1, jpi-1
          call zwxy_kern(zwx,zwy,mydomain,umask,vmask,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zslpx(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zslpy(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj
        DO ji = 2, jpi
          call zslpxy_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj
        DO ji = 2, jpi
          call zslpxy_update_kern(zslpx,zslpy,zwx,zwy,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj-1
        DO ji = 2, jpi-1
          call zwxy2_kern(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj-1
        DO ji = 2, jpi-1
          call mydomain_update_kern(mydomain,zwx,zwy,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx(ji,jj,1) = 0.d0
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jk = 2, jpk-1
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zwx_kern(zwx,tmask,mydomain,ji,jj,jk)
          END DO
       END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zslpx(ji,jj,1) = 0.d0
      END DO
    END DO
    !
    DO jk = 2, jpk-1
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zslpx_kern(zslpx,zwx,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 2, jpk-1     
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zslpx_update_kern(zslpx,zwx,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx(ji,jj,1) = pwn(ji,jj,1) * mydomain(ji,jj,1)
      END DO
    END DO
    !
    DO jk = 2, jpk
      DO jj = 2, jpj-1
        DO ji = 2, jpi-1
          call zwx2_kern(zwx,pwn,mydomain,zind,zslpx,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj-1     
        DO ji = 2, jpi-1
          call mydomain_kern(mydomain,zwx,ji,jj,jk)
        END DO
      END DO
    END DO
    !
  end subroutine tracer_advection
  !
end module psy_mod
