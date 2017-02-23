   !!=====================================================================================
   !! ***  traadv kernel extracted from the NEMO software (http://www.nemo-ocean.eu ) ***
   !! ***          governed by the CeCILL licence (http://www.cecill.info)            ***
   !!                                                   
   !! ***                             IS-ENES2 - CMCC/STCF                            ***
   !!=====================================================================================
PROGRAM tra_adv
   USE dl_timer, only: timer_init, timer_register, timer_start, timer_stop, timer_report
   USE psy_mod, only : zind_psy, zwxy_psy, zslpxy_psy, zslpxy_update_psy, zwxy2_psy, &
                       mydomain_update_psy, zwx_psy, zslpx_psy, zslpx_update_psy, &
                       zwx2_psy, mydomain_psy
   USE psy_mod, only : zero_top_layer, zero_bottom_layer, multiply_top_layer
   USE psy_mod, only : set_bounds
   implicit none
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: t3sn, t3ns, t3ew, t3we
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tsn 
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: pun, pvn, pwn
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mydomain, zslpx, zslpy, zwx, zwy, umask, vmask, tmask, zind
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: ztfreez, rnfmsk, upsmsk
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:)       :: rnfmsk_z
   REAL*8                                        :: zice, zu, z0u, zzwx, zv, z0v, zzwy, ztra, zalpha
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

   r = jpi*jpj*jpk

   ! the following three lines can be uncommented to randomize arrays initialization
   !call random_seed()
   !call random_number(r)
   !r = r*jpi*jpj*jpk

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

   r = jpi*jpj
   DO jj=1, jpj
      DO ji=1, jpi
         ztfreez(ji,jj) = ji*jj/r
         upsmsk(ji,jj) = ji*jj/r
         rnfmsk(ji,jj) = ji*jj/r
      END DO
   END DO

   DO jk=1, jpk
      rnfmsk_z(jk)=jk/jpk
   END DO

   call timer_stop(init_timer)

   ! temporary way to provide dimension information to PSy layer
   call set_bounds(jpi,jpj,jpk)

!***********************
!* Start of the synphony
!***********************
   call timer_start(step_timer)

   DO jt = 1, it

      call zind_psy(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask)
      call zero_bottom_layer(zwx)
      call zero_bottom_layer(zwy)
      call zwxy_psy(zwx,zwy,mydomain,umask,vmask)
      call zero_bottom_layer(zslpx)
      call zero_bottom_layer(zslpy)
      call zslpxy_psy(zslpx,zslpy,zwx,zwy)
      call zslpxy_update_psy(zslpx,zslpy,zwx,zwy)
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
