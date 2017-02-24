   !!=====================================================================================
   !! ***  traadv kernel extracted from the NEMO software (http://www.nemo-ocean.eu ) ***
   !! ***          governed by the CeCILL licence (http://www.cecill.info)            ***
   !!                                                   
   !! ***                             IS-ENES2 - CMCC/STFC                            ***
   !!=====================================================================================
PROGRAM tra_adv
   USE dl_timer, only: timer_init, timer_register, timer_start, timer_stop, timer_report
   !> Modules from the dl-esm-inf library
   USE kind_params_mod, only: wp
   USE grid_mod
   USE field_mod
   !PSyclone ... we would declare the links to the Kernel metadata here
   USE psy_mod, only : tracer_advection
   implicit none
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tmask_input
   !> 'Now' ocean temperature/salinity (in NEMO this is a 4D array with the fourth
   !! dimension holding the tracer index)
   type(r3d_field) :: tsn 
   !REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tsn 
   !> Three ocean velocity components
   type(r3d_field) :: pun, pvn, pwn
   !REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: pun, pvn, pwn
   type(r3d_field) :: mydomain, zslpx, zslpy, zwx, zwy, zind
   !REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: mydomain, zslpx, zslpy, zwx, zwy, umask, vmask, tmask, zind
   type(r2d_field) :: ztfreez, rnfmsk, upsmsk
   !REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:) :: ztfreez, rnfmsk, upsmsk
   type(r1d_field) :: rnfmsk_z
   !REAL*8, ALLOCATABLE, SAVE, DIMENSION(:) :: rnfmsk_z
   REAL*8                                        :: zice, zu, z0u, zzwx, zv, z0v, zzwy, ztra, zalpha
   REAL*8                                        :: zw, z0w
   INTEGER                                       :: jpi, jpj, jpk, jt
   INTEGER*8                                     :: it
   CHARACTER(len=10)                             :: env
   !> Timer indexes, one for initialisation, one for the 'time-stepping'
   INTEGER :: init_timer, step_timer
   !> The grid on which our fields are defined
   type(grid_type), target :: model_grid
   !> The (constant) grid resolution
   real(wp) :: dx=1.0d0, dy=1.0d0, dz=1.0d0

   CALL get_environment_variable("JPI", env)
   READ ( env, '(i10)' ) jpi
   CALL get_environment_variable("JPJ", env)
   READ ( env, '(i10)' ) jpj
   CALL get_environment_variable("JPK", env)
   READ ( env, '(i10)' ) jpk
   CALL get_environment_variable("IT", env)
   READ ( env, '(i10)' ) it

   ! Set-up our timers

   CALL timer_init()
   CALL timer_register(init_timer, label='Initialisation')
   CALL timer_register(step_timer, label='Time-stepping', num_repeats=it)

   ! Initialisation

   call timer_start(init_timer)

   ! Create the model grid. We use a NE offset (i.e. the U, V and F
   ! points immediately to the North and East of a T point all have the
   ! same i,j index).  This is the same offset scheme as used by NEMO.
   model_grid = grid_type(ARAKAWA_C, &
                          !  BC_PERIODIC, BC_NON_PERIODIC ??
                          (/BC_EXTERNAL,BC_EXTERNAL,BC_NONE/), &
                          OFFSET_NE)

   ! Generate our 'input' T mask
   allocate(tmask_input(jpi,jpj,jpk))
   tmask_input(:,:,:) = 1.0d0 ! all inner cells

   ! Having specified the T points mask, we can set up mesh parameters
   call grid_init(model_grid, jpi, jpj, jpk, dx, dy, dz, tmask_input)

   !> Example tracer field (e.g. temperature or salinity)
   tsn = r3d_field(model_grid, T_POINTS)
   mydomain = r3d_field(model_grid, T_POINTS)
   !> \todo Check that I've got the grid-point types right for these
   zwx = r3d_field(model_grid, U_POINTS)
   zwy = r3d_field(model_grid, V_POINTS)
   zslpx = r3d_field(model_grid, U_POINTS)
   zslpy = r3d_field(model_grid, V_POINTS)
   !> Ocean velocity components
   pun = r3d_field(model_grid, U_POINTS)
   pvn = r3d_field(model_grid, V_POINTS)
   !> \todo Where do vertical components live on the grid?
   pwn = r3d_field(model_grid, T_POINTS)
   !> \todo This is a workspace array; I don't yet know to what it corresponds
   zind = r3d_field(model_grid, T_POINTS)
   !> River mouth runoff mask (horiz)
   rnfmsk = r2d_field(model_grid, T_POINTS)
   !> River mouth runoff mask (vert)
   rnfmsk_z = r1d_field(model_grid)
   !> \todo Work out what these three fields represent
   ztfreez = r2d_field(model_grid, T_POINTS)
   upsmsk = r2d_field(model_grid, T_POINTS)
   tsn = r3d_field(model_grid, T_POINTS)

   call init_fields()

   call timer_stop(init_timer)

!***********************
!* Start of the symphony
!***********************
   call timer_start(step_timer)

   DO jt = 1, it

      call tracer_advection(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk, &
                            zwx,zwy,mydomain,zslpx,zslpy,pun,pvn,pwn)

      !PSyclone ... this is what the algorithm code would look like
      !call invoke (
      !  zind_kern(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask),
      !  zero_bottom_layer(zwx),
      !  zero_bottom_layer(zwy),
      !  zwxy_kern(zwx,zwy,mydomain,umask,vmask),
      !  zero_bottom_layer(zslpx),
      !  zero_bottom_layer(zslpy),
      !  zslpxy_kern(zslpx,zslpy,zwx,zwy),
      !  zslpxy_update_kern(zslpx,zslpy,zwx,zwy),
      !  zwxy2_kern(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy),
      !  mydomain_update_kern(mydomain,zwx,zwy),
      !  zero_top_layer(zwx),
      !  zero_bottom_layer(zwx),
      !  zwx_kern(zwx,tmask,mydomain),
      !  zero_top_layer(zslpx),
      !  zslpx_kern(zslpx,zwx),
      !  zslpx_update_kern(zslpx,zwx),
      !  multiply_top_layer(zwx,pwn,mydomain),
      !  zwx2_kern(zwx,pwn,mydomain,zind,zslpx),
      !  mydomain_kern(mydomain,zwx),
      !  name="tracer_advection")

  END DO

  call timer_stop(step_timer)

  call write_field(mydomain)

  CALL timer_report()

CONTAINS

  subroutine init_fields()
    REAL(wp) :: r
    integer :: ji, jj, jk

    ! Array initialization
    r = jpi*jpj*jpk

    ! the following three lines can be uncommented to randomize the array
    ! initialization
    !call random_seed()
    !call random_number(r)
    !r = r*jpi*jpj*jpk

   DO jk = 1, jpk
      DO jj = 1, jpj
          DO ji = 1, jpi
              mydomain%data(ji,jj,jk) =ji*jj*jk/r
              pun%data(ji,jj,jk) =ji*jj*jk/r
              pvn%data(ji,jj,jk) =ji*jj*jk/r
              pwn%data(ji,jj,jk) =ji*jj*jk/r
              tsn%data(ji,jj,jk)= ji*jj*jk/r
          END DO
      END DO
   END DO

   r = jpi*jpj
   DO jj=1, jpj
      DO ji=1, jpi
         ztfreez%data(ji,jj) = ji*jj/r
         upsmsk%data(ji,jj) = ji*jj/r
         rnfmsk%data(ji,jj) = ji*jj/r
      END DO
   END DO

   DO jk=1, jpk
      rnfmsk_z%data(jk)=jk/jpk
   END DO

  end subroutine init_fields

  subroutine write_field(mydomain)
    type(r3d_field), intent(in) :: mydomain
    integer :: ji, jj, jk

    OPEN(unit = 24, file = 'output.dat', form='formatted')
  
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1
          DO ji = 2, jpi-1
             write(24,*) mydomain%data(ji,jj,jk)
          END DO
       END DO
    END DO

    CLOSE(24)
  end subroutine write_field

END PROGRAM tra_adv
