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
                              zwx,zwy,mydomain,zslpx,zslpy,pun,pvn,pwn)
    ! Infrastructure modules
    use field_mod
    use grid_mod
    ! Kernel modules
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
    implicit none
    !
    ! RF all intent out's here are actually local to this call
    type(r3d_field), intent(inout) :: zind, zwx, zwy, zslpx, zslpy
    type(r3d_field), intent(in)    :: tsn, pun, pvn, pwn
    type(r2d_field), intent(in)    :: rnfmsk, upsmsk, ztfreez
    type(r1d_field), intent(in)    :: rnfmsk_z
    type(r3d_field), intent(inout) :: mydomain
    ! local variables
    integer :: ji,jj,jk
    type(grid_type), pointer :: grid
    !
    ! Get a pointer to the grid from a read-only argument. We then use
    ! this to access grid-properties (such as masks) required by the
    ! kernels
    grid => tsn%grid

    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zind_kern(zind%data, tsn%data, ztfreez%data, rnfmsk%data, &
                         rnfmsk_z%data, upsmsk%data, grid%tmask, ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx%data(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwy%data(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 1, jpj-1
        DO ji = 1, jpi-1
          call zwxy_kern(zwx%data,zwy%data,mydomain%data, &
                         grid%umask,grid%vmask,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zslpx%data(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zslpy%data(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj
        DO ji = 2, jpi
          call zslpxy_kern(zslpx%data,zslpy%data,zwx%data,zwy%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj
        DO ji = 2, jpi
          call zslpxy_update_kern(zslpx%data,zslpy%data,zwx%data,zwy%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj-1
        DO ji = 2, jpi-1
          call zwxy2_kern(zwx%data,zwy%data,pun%data,pvn%data,mydomain%data, &
                          zind%data,zslpx%data,zslpy%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj-1
        DO ji = 2, jpi-1
          call mydomain_update_kern(mydomain%data,zwx%data,zwy%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx%data(ji,jj,1) = 0.d0
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx%data(ji,jj,jpk) = 0.d0
      END DO
    END DO
    !
    DO jk = 2, jpk-1
       DO jj = 1, jpj
          DO ji = 1, jpi
             call zwx_kern(zwx%data,grid%tmask,mydomain%data,ji,jj,jk)
          END DO
       END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zslpx%data(ji,jj,1) = 0.d0
      END DO
    END DO
    !
    DO jk = 2, jpk-1
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zslpx_kern(zslpx%data,zwx%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 2, jpk-1     
      DO jj = 1, jpj
        DO ji = 1, jpi
          call zslpx_update_kern(zslpx%data,zwx%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jj=1,jpj
      DO ji=1,jpi
        zwx%data(ji,jj,1) = pwn%data(ji,jj,1) * mydomain%data(ji,jj,1)
      END DO
    END DO
    !
    DO jk = 2, jpk
      DO jj = 2, jpj-1
        DO ji = 2, jpi-1
          call zwx2_kern(zwx%data,pwn%data,mydomain%data,zind%data, &
                         zslpx%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
    DO jk = 1, jpk-1
      DO jj = 2, jpj-1     
        DO ji = 2, jpi-1
          call mydomain_kern(mydomain%data,zwx%data,ji,jj,jk)
        END DO
      END DO
    END DO
    !
  end subroutine tracer_advection
  !
end module psy_mod
