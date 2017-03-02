module tra_adv_kern
  !
  implicit none
  private
  public :: init_3d_arrays, init_2d_arrays, init_1d_arrays, set_bounds
  !
  integer :: jpi,jpj,jpk
  integer :: ji,jj,jk
  !
contains
  !
  subroutine set_bounds(jpi_in,jpj_in,jpk_in)
    !
    integer, intent(in) :: jpi_in, jpj_in, jpk_in
    !
    jpi=jpi_in
    jpj=jpj_in
    jpk=jpk_in
    !
  end subroutine set_bounds
  !
  subroutine init_3d_arrays(umask, mydomain, pun, pvn, pwn, vmask, tsn, tmask)
    !
    real*8, intent(out), dimension(jpi,jpj,jpk) :: umask, mydomain, pun, pvn, pwn, vmask, tsn, tmask
    !
    real*8 :: r
    !
    r = jpi*jpj*jpk
    DO jk = 1, jpk
       DO jj = 1, jpj
          DO ji = 1, jpi
             umask(ji,jj,jk) = ji*jj*jk/r
             mydomain(ji,jj,jk) =ji*jj*jk/r
             pun(ji,jj,jk) =ji*jj*jk/r
             pvn(ji,jj,jk) =ji*jj*jk/r
             pwn(ji,jj,jk) =ji*jj*jk/r
             vmask(ji,jj,jk)= ji*jj*jk/r
             tsn(ji,jj,jk)= ji*jj*jk/r
             tmask(ji,jj,jk)= ji*jj*jk/r
          END DO
       END DO
    END DO
    !
  end subroutine init_3d_arrays
  !
  subroutine init_2d_arrays(ztfreez, upsmsk, rnfmsk)
    !
    real*8, intent(out), dimension(jpi,jpj) :: ztfreez, upsmsk, rnfmsk
    !
    real*8 :: r
    !
    r = jpi*jpj
    DO jj=1, jpj
       DO ji=1, jpi
          ztfreez(ji,jj) = ji*jj/r
          upsmsk(ji,jj) = ji*jj/r
          rnfmsk(ji,jj) = ji*jj/r
       END DO
    END DO
    !
  end subroutine init_2d_arrays
  !
  subroutine init_1d_arrays(rnfmsk_z)
    !
    real*8, intent(out) :: rnfmsk_z(jpk)
    !
    DO jk=1, jpk
       rnfmsk_z(jk)=jk/jpk
    END DO
    !
  end subroutine init_1d_arrays
  !
end module tra_adv_kern
