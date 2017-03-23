module tra_adv_kern
  !
  implicit none
  !
  private
  public :: init_3d_arrays, init_2d_arrays, init_1d_arrays, set_bounds, zind_compute, &
       zero_bottom_layer, zwxy_compute, zslpxy_compute, zslpxy_update_compute, &
       zwxy2_compute, mydomain_update_compute, zero_top_layer, zwx_compute, &
       zslpx_compute, zslpx_update_compute, multiply_top_layer, zwx2_compute, &
       mydomain_compute
  !
  integer :: jpi,jpj,jpk
  integer :: ji,jj,jk
  !
  ! For coarse-grained OpenMP tiling
  TYPE :: tile_type
     INTEGER :: kstart, kstop
     INTEGER :: jstart, jstop
     ! For loop limits that are 2 instead of 1 and jp{j,k}-1 instead
     ! of jp{j,k}
     INTEGER :: kstartp1, kstopm1
     INTEGER :: jstartp1, jstopm1
  END TYPE tile_type

  INTEGER, SAVE                                    :: ntiles
  INTEGER, SAVE                                    :: ntilez, ntiley
  ! Whether or not to automatically compute dimensions of tiling grid
  LOGICAL, SAVE                                    :: auto_tile
  ! Extent of this 1D array will be ntiles which, by default will
  ! be set to the no. of threads
  TYPE(tile_type), ALLOCATABLE, SAVE, DIMENSION(:) :: tile
  INTEGER, SAVE                                    :: max_tile_width
  INTEGER, SAVE                                    :: max_tile_height

contains
  
  !-------------------------------------------------------------------

  subroutine set_bounds(jpi_in,jpj_in,jpk_in)
    !
    integer, intent(in) :: jpi_in, jpj_in, jpk_in
    !
    jpi=jpi_in
    jpj=jpj_in
    jpk=jpk_in
    !
    call openmp_grid_init()
    
  end subroutine set_bounds
  
  !-------------------------------------------------------------------
  
  SUBROUTINE openmp_grid_init()
     !$ use omp_lib
     INTEGER :: idy, idz, ival, jval ! For tile extent calculation
     INTEGER :: ierr, nwidth
     INTEGER :: ji,jj, ith
     INTEGER :: nthreads       ! No. of OpenMP threads being used
     INTEGER :: kover, kunder, idztmp
     INTEGER :: jover, junder, idytmp
     ! For doing stats on tile sizes
     INTEGER :: nvects, nvects_sum, nvects_min, nvects_max 
     LOGICAL, PARAMETER :: print_tiles = .TRUE.

     ! Set-up regular grid of tiles

     ! Dimensions of the grid of tiles. Done here to save coding
     ! to parse command line for the moment.
     ntilez = 2
     ntiley = 2
     ntiles = ntilez * ntiley

     nthreads = 1
!$   nthreads = omp_get_max_threads()
     WRITE (*,"(/'Have ',I3,' OpenMP threads available.')") nthreads

     ! If we've not manually specified a grid of tiles then use the no. of
     ! threads
     IF(ntiles == 1 .AND. auto_tile)ntiles = nthreads

     ALLOCATE(tile(ntiles), Stat=ierr)
     IF(ierr /= 0 )THEN
        STOP 'Harness: ERROR: failed to allocate tiling structures'
     END IF

     IF(auto_tile)THEN

        ! We make the tiles as square as possible
        ntilez = INT( SQRT(REAL(ntiles)) )
        DO WHILE(MOD(ntiles,ntilez) /= 0)
           ntilez = ntilez - 1
        END DO
        ntiley = ntiles/ntilez

     END IF ! automatic determination of tiling grid

     WRITE (*,"('OpenMP thread tiling using grid of ',I3,'x',I3)") ntiley,ntilez

     ! NEMO arrays are (i,j,k) so that innermost loops are over i. We want
     ! to leave such loops for SIMD vectorisation. Therefore only the j and k
     ! dimensions are available for tiling. Hence a 'tile' will, physically
     ! speaking, consist of a bundle of rows, contiguous in i.
     !
     ! i             j                     -----\                         
     !  \--     ----/                 ----/      ------\                  
     !     \---/                 -------\               -----\            
     !       |             -----/x x x x ------\              ------\     
     !       |        ------\ x x x x x x x x x -----\               ---- 
     !       k   ----/       ------\ x x x x x x x x x------\    ----/  | 
     !        --/\                  -----\x x x x x x x x x -----       | 
     !        |   ------\                 ------\x x x-----/x x|        | 
     !        |          -----\                  ----/ x x x x |        | 
     !        |                ------\      ----/    |x x x x x|        | 
     !        |                       -----/         | x x x x |        | 
     !        |                           |          |x x x x x|       ---
     !        |                           |          | x x x x |  ----/ | 
     !        |---\                       |          |x x x x----/      | 
     !        |x x ------\                |          | -----/  |        | 
     !        | x x x x x -----\          |       ----/        |        | 
     !        |x x x x x x x x x------\   |  ----/x x|         |        | 
     !        | x x x x x x x x x x x x-----/x x x x |         |        | 
     !        |x x x x x x x x x x x x x x| x x x x x|         |       ---
     !        | x x x x x x x x x x x x x |x x x x x |         |  ----/   
     !         ---\x x x x x x x x x x x x| x x x x x|       ----/        
     !             ------\x x x x x x x x |x x x x x | -----/             
     !                    -----\ x x x x x| x x x ----/                   
     !                          ------\ x |x ----/                        
     !                                 -----/                             
     ! Tiles at left and right of domain (or top and bottom) only have single
     ! overlap. Every other tile has two overlaps. So: 
     ! jpj = (ntiley-2)*(idy-2) + 2*(idy-1)
     !     =  ntiley*idy - 2*ntiley + 4 - 2
     !=> jpj - 2 + 2*ntiley = ntiley * idy
     !=> idy = (jpj - 2)/ntiley + 2
     ! Rearranging this gives the following expressions for idz and idy:
     idz = (jpk - 2)/ntilez + 2
     idy = (jpj - 2)/ntiley + 2

     ! Integer arithmetic means that ntiley tiles of width idy might
     ! actually span a width greater or less than jpj. If so, we try and
     ! reduce the width of each row by just one until we've accounted
     ! for the <jover> extra rows.
     nwidth = (ntiley-2)*(idy-2) + 2*(idy-1)
     IF(nwidth > jpj)THEN
        jover  = nwidth - jpj
        junder = 0
     ELSE IF(nwidth < jpj)THEN
        jover  = 0
        junder = jpj - nwidth
     ELSE
        jover  = 0
        junder = 0
     END IF
     ! Ditto for z dimension
     nwidth = (ntilez-2)*(idz-2) + 2*(idz-1)
     IF(nwidth > jpk)THEN
        kover  = nwidth - jpk
        kunder = 0
     ELSE IF(nwidth < jpk)THEN
        kover  = 0
        kunder = jpk - nwidth
     ELSE
        kover  = 0
        kunder = 0
     END IF

     ! For AVX-256 instructions, I think we want MOD(idx,4) == 0
     !idx = idx + (4 - MOD(idx,4))

     WRITE(*,"('Tile width = ',I4,', tile height = ',I4)") idy, idz
     WRITE(*,"('jover = ',I3,', junder = ',I3)") jover, junder
     WRITE(*,"('kover = ',I3,', kunder = ',I3)") kover, kunder

     ith = 1
     jval = 1

     nvects_max = 0
     nvects_min = 1000000
     nvects_sum = 0
     max_tile_width  = 0
     max_tile_height = 0

     IF(print_tiles)WRITE(*,"(/'Tile dimensions:')")

     DO jj = 1, ntiley, 1

        ! If necessary, correct the height of this tile row
        IF(jover > 0)THEN
           idytmp = idy - 1
           jover = jover - 1
        ELSE IF(junder > 0)THEN
           idytmp = idy + 1
           junder = junder - 1
        ELSE
           idytmp = idy
        END IF

        ival = 1

        DO ji = 1, ntilez, 1
         
           ! If necessary, correct the depth of this tile column
           IF(kover > 0)THEN
              idztmp = idz - 1
              kover = kover - 1
           ELSE IF(kunder > 0)THEN
              idztmp = idz + 1
              kunder = kunder - 1
           ELSE
              idztmp = idz
           END IF

           tile(ith)%kstart = ival
           tile(ith)%kstartp1 = tile(ith)%kstart + 1

           IF(ji == ntilez)THEN
              tile(ith)%kstop = jpk
              tile(ith)%kstopm1 = jpk - 1
           ELSE
              tile(ith)%kstop = MIN(ival + idztmp - 1, jpk)
              tile(ith)%kstopm1 = tile(ith)%kstop - 1
           END IF

           tile(ith)%jstart = jval
           tile(ith)%jstartp1 = tile(ith)%jstart + 1

           IF(jj == ntiley)THEN
              tile(ith)%jstop = jpj
              tile(ith)%jstopm1 = jpj - 1
           ELSE
              tile(ith)%jstop = MIN(jval + idytmp - 1, jpj)
              tile(ith)%jstopm1 = tile(ith)%jstop - 1
           END IF

           IF(print_tiles)THEN
              WRITE(*,"('tile[',I4,'](',I4,':',I4,')(',I4,':',I4,'), &
                   & interior:(',I4,':',I4,')(',I4,':',I4,') ')") &
                 ith,                                   &
                 tile(ith)%kstart, tile(ith)%kstop,     &
                 tile(ith)%jstart, tile(ith)%jstop,     &
                 tile(ith)%kstartp1, tile(ith)%kstopm1, &
                 tile(ith)%jstartp1, tile(ith)%jstopm1
           END IF

           ! Collect some data on the distribution of tile sizes for
           ! loadbalance info
           nvects = (tile(ith)%kstopm1 - tile(ith)%kstartp1 + 1) * &
                    (tile(ith)%jstopm1 - tile(ith)%jstartp1 + 1)
           nvects_sum = nvects_sum + nvects
           nvects_min = MIN(nvects_min, nvects)
           nvects_max = MAX(nvects_max, nvects)

           ! For use when allocating tile-'private' work arrays
           max_tile_width  = MAX(max_tile_width, &
                                 (tile(ith)%kstop - tile(ith)%kstart + 1) )
           max_tile_height = MAX(max_tile_height, &
                                 (tile(ith)%jstop - tile(ith)%jstart + 1) )

           ival = ival + idztmp - 2
           ith = ith + 1
        END DO
        jval = jval + idytmp - 2
     END DO

     ! Print tile-size statistics
     WRITE(*,"(/'Mean tile size = ',F6.1,' cols = ',F7.1,' KB')") &
                                   REAL(nvects_sum)/REAL(ntiles), &
                                   REAL(8*jpi*nvects_sum)/REAL(ntiles*1024)
     WRITE(*,"('Min,max tile size = ',I4,',',I4)") nvects_min,nvects_max
     WRITE(*,"('Tile load imbalance (%) =',F5.1)") &
                              100.0*(nvects_max-nvects_min)/REAL(nvects_min)
     WRITE (*,"('Max tile dims are ',I4,'x',I4)") max_tile_width, &
                                                  max_tile_height

  end subroutine openmp_grid_init
  
  !-------------------------------------------------------------------

  subroutine init_3d_arrays(umask, mydomain, pun, pvn, pwn, vmask, tsn, tmask)
    !
    real*8, intent(out), dimension(jpi,jpj,jpk) :: umask, mydomain, pun, &
         pvn, pwn, vmask, tsn, tmask
    !
    real*8 :: r
    integer :: tid
    !
    r = jpi*jpj*jpk
    !$omp parallel do default(none), private(ji, jj, jk, tid), &
    !$omp shared(jpi, ntiles, tile, r, umask, mydomain, pun, pvn, pwn, vmask, &
    !$omp tsn, tmask)
    DO tid = 1, ntiles, 1
       DO jk = tile(tid)%kstart, tile(tid)%kstop ! 1, jpk
          DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
             ! omp simd
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
    END DO
    !$omp end parallel do
    !
  end subroutine init_3d_arrays
  !
  subroutine init_2d_arrays(ztfreez, upsmsk, rnfmsk)
    !
    real*8, intent(out), dimension(jpi,jpj) :: ztfreez, upsmsk, rnfmsk
    !
    real*8  :: r
    integer :: tid
    !
    r = jpi*jpj
    !$omp parallel do default(none), private(ji, jj, jk), &
    !$omp shared(ntiles, tile, jpi, r, ztfreez, upsmsk, rnfmsk)
    DO tid = 1, ntiles, 1
       DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
          DO ji=1, jpi
             ztfreez(ji,jj) = ji*jj/r
             upsmsk(ji,jj) = ji*jj/r
             rnfmsk(ji,jj) = ji*jj/r
          END DO
       END DO
    END DO
    !$omp end parallel do
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
  subroutine zind_compute(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask)
    !
    real*8, intent(out) :: zind(jpi,jpj,jpk)
    real*8, intent(in)  :: tsn(jpi,jpj,jpk), tmask(jpi,jpj,jpk)
    real*8, intent(in)  :: ztfreez(jpi,jpj), rnfmsk(jpi,jpj), upsmsk(jpi,jpj)
    real*8, intent(in)  :: rnfmsk_z(jpk)
    !
    real*8 :: zice
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk, zice), &
    !$omp shared(ntiles, tile, jpi, tmask, rnfmsk, rnfmsk_z, upsmsk, tsn, ztfreez, zind)
    DO tid = 1, ntiles, 1
       DO jk = tile(tid)%kstart, tile(tid)%kstop ! 1, jpk
          DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
             DO ji = 1, jpi
                IF( tsn(ji,jj,jk) <= ztfreez(ji,jj) + 0.1d0 ) THEN   ;   zice = 1.d0
                ELSE                                                 ;   zice = 0.d0
                ENDIF
                zind(ji,jj,jk) = MAX (   &
                     rnfmsk(ji,jj) * rnfmsk_z(jk),      & 
                     upsmsk(ji,jj)               ,      &
                     zice                               &
                     &                  ) * tmask(ji,jj,jk)
                zind(ji,jj,jk) = 1 - zind(ji,jj,jk)
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zind_compute
  !
  subroutine zero_bottom_layer(field)
    !
    real*8, intent(out) :: field(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj), &
    !$omp shared(ntiles, tile, jpi, jpk, field)
    DO tid = 1, ntiles, 1
       DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
          DO ji = 1, jpi
             field(ji,jj,jpk) = 0.0d0
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zero_bottom_layer
  !
  subroutine zwxy_compute(zwx,zwy,umask,vmask,mydomain)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    real*8, intent(in)  :: umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk), &
                           mydomain(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk), &
    !$omp shared(jpi, ntiles, tile, umask, vmask, mydomain, zwx, zwy)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstart, tile(tid)%kstopm1 ! 1, jpk-1
          DO jj = tile(tid)%jstart, tile(tid)%jstopm1 ! 1, jpj-1
             DO ji = 1, jpi-1
                zwx(ji,jj,jk) = umask(ji,jj,jk) * &
                           ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
                zwy(ji,jj,jk) = vmask(ji,jj,jk) * &
                           ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zwxy_compute
  !
  subroutine zslpxy_compute(zslpx,zslpy,zwx,zwy)
    !
    real*8, intent(out) :: zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    real*8, intent(in) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk), &
    !$omp shared(jpi, ntiles, tile, zwx, zwy, zslpx, zslpy)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstart, tile(tid)%kstopm1 ! 1, jpk-1
          DO jj = tile(tid)%jstartp1, tile(tid)%jstop ! 2, jpj
             DO ji = 2, jpi 
                zslpx(ji,jj,jk) = ( zwx(ji,jj,jk) + zwx(ji-1,jj  ,jk) )   &
                  &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji-1,jj  ,jk) ) )
                zslpy(ji,jj,jk) = ( zwy(ji,jj,jk) + zwy(ji  ,jj-1,jk) )   &
                     &            * ( 0.25d0 + SIGN( 0.25d0, zwy(ji,jj,jk) * zwy(ji  ,jj-1,jk) ) )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zslpxy_compute
  !
  subroutine zslpxy_update_compute(zslpx,zslpy,zwx,zwy)
    !
    real*8, intent(inout) :: zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    real*8, intent(in) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk), &
    !$omp shared(jpi, ntiles, tile, zwx, zwy, zslpx, zslpy)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstart, tile(tid)%kstopm1 ! 1, jpk-1    
          DO jj = tile(tid)%jstartp1, tile(tid)%jstop ! 2, jpj
             DO ji = 2, jpi
                zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji  ,jj,jk) ),   &
                  &                           2.d0*ABS( zwx  (ji-1,jj,jk) ),   &
                  &                           2.d0*ABS( zwx  (ji  ,jj,jk) ) )
                zslpy(ji,jj,jk) = SIGN( 1.d0, zslpy(ji,jj,jk) ) * MIN(    ABS( zslpy(ji,jj  ,jk) ),   &
                  &                           2.d0*ABS( zwy  (ji,jj-1,jk) ),   &
                  &                           2.d0*ABS( zwy  (ji,jj  ,jk) ) )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zslpxy_update_compute
  !
  subroutine zwxy2_compute(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    real*8, intent(in), dimension(jpi,jpj,jpk)  :: pun, pvn, mydomain, zind, zslpx, zslpy
    !
    real*8 :: zdt, z0u, zalpha, zu, zzwx, zzwy, z0v, zv
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk, zdt, z0u, z0v, &
    !$omp zu, zv, zalpha, zzwx, zzwy), shared(jpi, ntiles, tile, zslpx, zslpy, &
    !$omp zind, mydomain, pun, pvn, zwx, zwy)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstart, tile(tid)%kstopm1 ! 1, jpk-1
          zdt  = 1
          DO jj = tile(tid)%jstartp1, tile(tid)%jstopm1 ! 2, jpj-1
             DO ji = 2, jpi-1
                z0u = SIGN( 0.5d0, pun(ji,jj,jk) )
                zalpha = 0.5d0 - z0u
                zu  = z0u - 0.5d0 * pun(ji,jj,jk) * zdt
             
                zzwx = mydomain(ji+1,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji+1,jj,jk))
                zzwy = mydomain(ji  ,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji  ,jj,jk))
             
                zwx(ji,jj,jk) = pun(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
             
                z0v = SIGN( 0.5d0, pvn(ji,jj,jk) )
                zalpha = 0.5d0 - z0v
                zv  = z0v - 0.5d0 * pvn(ji,jj,jk) * zdt
             
                zzwx = mydomain(ji,jj+1,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj+1,jk))
                zzwy = mydomain(ji,jj  ,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj  ,jk))
             
                zwy(ji,jj,jk) = pvn(ji,jj,jk) * ( zalpha * zzwx + (1.d0-zalpha) * zzwy )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zwxy2_compute
  !
  subroutine mydomain_update_compute(mydomain,zwx,zwy)
    !
    real*8, intent(inout) :: mydomain(jpi,jpj,jpk)
    real*8, intent(in)    :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    !
    real*8 :: zbtr, ztra
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk, zbtr, ztra), &
    !$omp shared(jpi, ntiles, tile, zwx, zwy, mydomain)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstart, tile(tid)%kstopm1 ! 1, jpk-1
          DO jj = tile(tid)%jstartp1, tile(tid)%jstopm1 ! 2, jpj-1     
             DO ji = 2, jpi-1
                zbtr = 1.
                ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                  &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
                mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine mydomain_update_compute
  !
  subroutine zero_top_layer(field)
    !
    real*8, intent(out) :: field(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj), &
    !$omp shared(ntiles, tile, jpi, field)
    DO tid = 1, ntiles
       DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
          DO ji = 1, jpi
             field(ji,jj,1) = 0.0d0
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zero_top_layer
  !
  subroutine zwx_compute(zwx,tmask,mydomain)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk)
    real*8, intent(in)  :: tmask(jpi,jpj,jpk), mydomain(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk), &
    !$omp shared(ntiles, tile, jpi, zwx, mydomain, tmask)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstartp1, tile(tid)%kstopm1 ! 2, jpk-1
          DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
             DO ji = 1, jpi
                zwx(ji,jj,jk) = tmask(ji,jj,jk) * ( mydomain(ji,jj,jk-1) - &
                                                    mydomain(ji,jj,jk) )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zwx_compute
  !
  subroutine zslpx_compute(zslpx,zwx)
    !
    real*8, intent(out) :: zslpx(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk), &
    !$omp shared(jpi, ntiles, tile, zwx, zslpx)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstartp1, tile(tid)%kstopm1 ! 2, jpk-1    
          DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
             DO ji = 1, jpi
                zslpx(ji,jj,jk) = ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
                     & * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zslpx_compute
  !
  subroutine zslpx_update_compute(zslpx,zwx)
    !
    real*8, intent(inout) :: zslpx(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj, jk), &
    !$omp shared(jpi, ntiles, tile, zwx, zslpx)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstartp1, tile(tid)%kstopm1 ! 2, jpk-1     
          DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
             DO ji = 1, jpi
                zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN( ABS( zslpx(ji,jj,jk  ) ), &
                  &                                               2.d0*ABS( zwx  (ji,jj,jk+1) ),   &
                  &                                               2.d0*ABS( zwx  (ji,jj,jk  ) )  )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zslpx_update_compute
  !
  subroutine multiply_top_layer(field1,field2,field3)
    !
    real*8, intent(out) :: field1(jpi,jpj,jpk)
    real*8, intent(in)  :: field2(jpi,jpj,jpk), field3(jpi,jpj,jpk)
    integer :: tid
    !
    !$omp parallel do default(none), private(tid, ji, jj), &
    !$omp shared(ntiles, tile, jpi, field1, field2, field3)
    DO tid = 1, ntiles
       DO jj = tile(tid)%jstart, tile(tid)%jstop ! 1, jpj
          DO ji = 1, jpi
             field1(ji,jj,1) = field2(ji,jj,1) * field3(ji,jj,1)
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine multiply_top_layer
  !
  subroutine zwx2_compute(zwx,pwn,mydomain,zind,zslpx)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk)
    real*8, intent(in), dimension(jpi,jpj,jpk)  :: pwn, mydomain, zind, zslpx
    !
    real*8 :: zdt, zbtr, z0w, zalpha, zw, zzwx, zzwy
    integer :: tid
    !
    zdt  = 1
    zbtr = 1.
    !$omp parallel do default(none), &
    !$omp private(tid, ji, jj, jk, z0w, zalpha, zw, zzwx, zzwy), &
    !$omp shared(zdt, zbtr, jpi, ntiles, tile, zwx, zslpx, zind, mydomain, pwn)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstartp1, tile(tid)%kstop ! 2, jpk
          DO jj = tile(tid)%jstartp1, tile(tid)%jstopm1 ! 2, jpj-1     
             DO ji = 2, jpi-1
                z0w = SIGN( 0.5d0, pwn(ji,jj,jk) )
                zalpha = 0.5d0 + z0w
                zw  = z0w - 0.5d0 * pwn(ji,jj,jk) * zdt * zbtr
                
                zzwx = mydomain(ji,jj,jk  ) + zind(ji,jj,jk-1) * (zw * zslpx(ji,jj,jk  ))
                zzwy = mydomain(ji,jj,jk-1) + zind(ji,jj,jk-1) * (zw * zslpx(ji,jj,jk-1))
                zwx(ji,jj,jk) = pwn(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine zwx2_compute
  !
  subroutine mydomain_compute(mydomain,zwx)
    !
    real*8, intent(out) :: mydomain(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk)
    !
    real*8 :: zbtr, ztra
    integer :: tid
    !
    zbtr = 1.
    !$omp parallel do default(none), private(tid, ji, jj, jk, ztra), &
    !$omp shared(jpi, ntiles, tile, zbtr, zwx, mydomain)
    DO tid = 1, ntiles
       DO jk = tile(tid)%kstart, tile(tid)%kstopm1 ! 1, jpk-1
          DO jj = tile(tid)%jstartp1, tile(tid)%jstopm1 ! 2, jpj-1     
             DO ji = 2, jpi-1
                ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
                mydomain(ji,jj,jk) = ztra
             END DO
          END DO
       END DO
    END DO
    !$omp end parallel do
    !
  end subroutine mydomain_compute
  !
end module tra_adv_kern
