   !!=====================================================================================
   !! ***  traadv kernel extracted from the NEMO software (http://www.nemo-ocean.eu ) ***
   !! ***          governed by the CeCILL licence (http://www.cecill.info)            ***
   !!                                                   
   !! ***                             IS-ENES2 - CMCC/STFC                            ***
   !!=====================================================================================
PROGRAM tra_adv
   USE dl_timer, only: timer_init, timer_register, timer_start, timer_stop, timer_report
   USE tra_adv_kern, only: init_3d_arrays, init_2d_arrays, init_1d_arrays, set_bounds, &
        zind_compute, zero_bottom_layer, zwxy_compute, zslpxy_compute, &
        zslpxy_update_compute, zwxy2_compute, mydomain_update_compute, zero_top_layer, &
        zwx_compute, zslpx_compute, zslpx_update_compute, multiply_top_layer, &
        zwx2_compute, mydomain_compute
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: t3sn, t3ns, t3ew, t3we
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tsn 
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: pun, pvn, pwn
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mydomain, zslpx, zslpy, zwx, zwy, umask, vmask, tmask, zind
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: ztfreez, rnfmsk, upsmsk
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:)       :: rnfmsk_z
   REAL*8                                        :: zice, zu, z0u, zzwx, zv, z0v, zzwy, ztra, zbtr, zdt, zalpha
   REAL*8                                        :: r
   REAL*8                                        :: zw, z0w
   INTEGER                                       :: jpi, jpj, jpk, ji, jj, jk, jt
   INTEGER*8                                     :: it
   CHARACTER(len=10)                             :: env
   !> Timer indexes, one for initialisation, one for the 'time-stepping'
   INTEGER :: init_timer, step_timer

   CALL get_environment_variable("JPI", env)
   READ ( env, '(i10)' ) jpi
   CALL get_environment_variable("JPJ", env)
   READ ( env, '(i10)' ) jpj
   CALL get_environment_variable("JPK", env)
   READ ( env, '(i10)' ) jpk
   CALL get_environment_variable("IT", env)
   READ ( env, '(i10)' ) it

   call set_bounds(jpi, jpj, jpk)

   ! Set-up our timers

   CALL timer_init()
   CALL timer_register(init_timer, label='Initialisation')
   CALL timer_register(step_timer, label='Time-stepping', num_repeats=it)

   ! Initialisation

   call timer_start(init_timer)

   ALLOCATE( mydomain (jpi,jpj,jpk))
   ALLOCATE( zwx (jpi,jpj,jpk))
   ALLOCATE( zwy (jpi,jpj,jpk))
   ALLOCATE( zslpx (jpi,jpj,jpk))
   ALLOCATE( zslpy (jpi,jpj,jpk))
   ALLOCATE( pun (jpi,jpj,jpk))
   ALLOCATE( pvn (jpi,jpj,jpk))
   ALLOCATE( pwn (jpi,jpj,jpk))
   ALLOCATE( umask (jpi,jpj,jpk))
   ALLOCATE( vmask (jpi,jpj,jpk))
   ALLOCATE( tmask (jpi,jpj,jpk))
   ALLOCATE( zind (jpi,jpj,jpk))
   ALLOCATE( ztfreez (jpi,jpj))
   ALLOCATE( rnfmsk (jpi,jpj))
   ALLOCATE( upsmsk (jpi,jpj))
   ALLOCATE( rnfmsk_z (jpk))
   ALLOCATE( tsn(jpi,jpj,jpk))

! arrays initialization

   call init_3d_arrays(umask, mydomain, pun, pvn, pwn, vmask, tsn, tmask)
   call init_2d_arrays(ztfreez, upsmsk, rnfmsk)
   call init_1d_arrays(rnfmsk_z)
   call timer_stop(init_timer)

!***********************
!* Start of the symphony
!***********************
   call timer_start(step_timer)

   DO jt = 1, it
      
      call zind_compute(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask)
      call zero_bottom_layer(zwx)
      call zero_bottom_layer(zwy)
      call zwxy_compute(zwx,zwy,umask,vmask,mydomain)
      call zero_bottom_layer(zslpx)
      call zero_bottom_layer(zslpy)
      call zslpxy_compute(zslpx,zslpy,zwx,zwy)
      call zslpxy_update_compute(zslpx,zslpy,zwx,zwy)
      call zwxy2_compute(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy)
      call mydomain_update_compute(mydomain,zwx,zwy)
      call zero_top_layer(zwx)
      call zero_bottom_layer(zwx)
      call zwx_compute(zwx,tmask,mydomain)
      call zero_top_layer(zslpx)
      call zslpx_compute(zslpx,zwx)
      call zslpx_update_compute(zslpx,zwx)
      call multiply_top_layer(zwx,pwn,mydomain)
      call zwx2_compute(zwx,pwn,mydomain,zind,zslpx)
      call mydomain_compute(mydomain,zwx)

   END DO

   call timer_stop(step_timer)
  
   OPEN(unit = 4, file = 'output.dat', form='formatted')
  
   DO jk = 1, jpk-1
      DO jj = 2, jpj-1
         DO ji = 2, jpi-1
            write(4,*) mydomain(ji,jj,jk)
         END DO
      END DO
   END DO
  
   CLOSE(4)

   DEALLOCATE( mydomain )
   DEALLOCATE( zwx )
   DEALLOCATE( zwy )
   DEALLOCATE( zslpx )
   DEALLOCATE( zslpy )
   DEALLOCATE( pun )
   DEALLOCATE( pvn )
   DEALLOCATE( pwn )
   DEALLOCATE( umask)
   DEALLOCATE( vmask)
   DEALLOCATE( tmask)
   DEALLOCATE( zind )
   DEALLOCATE( ztfreez )
   DEALLOCATE( rnfmsk)
   DEALLOCATE( upsmsk)
   DEALLOCATE( rnfmsk_z)
   DEALLOCATE( tsn)

   CALL timer_report()

END PROGRAM tra_adv
