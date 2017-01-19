   !!======================================================================
   !!                       ***  traadv NEMO kernel   ***
   !!                             IS-ENES2 - CMCC 
   !!======================================================================
PROGRAM tra_adv
#include "mpif.h"

   REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: t3sn, t3ns, t3ew, t3we
   REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tsn 
   REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: pun, pvn, pwn
   REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mydomain, zslpx, zslpy, zwx, zwy, umask, vmask, tmask, zind
   REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: ztfreez, rnfmsk, upsmsk
   REAL, ALLOCATABLE, SAVE, DIMENSION(:)       :: rnfmsk_z
   REAL                                        :: zice, zu, z0u, zzwx, zv, z0v, zzwy, ztra, zbtr, zdt, zalpha
   INTEGER                                     :: nbondi, nbondj, nono, noso, noea, nowe
   INTEGER                                     :: jpi, jpj, jpk
   INTEGER                                     :: myrank, p, ierr, realcomm
   INTEGER                                     :: px, py, cpt
   CHARACTER(len=10)                           :: env
   INTEGER                                     :: jpreci, jprecj, nreci, nrecj, nlci, nlcj
   INTEGER                                     :: i, hpmid, color, a, j, ji, jj, jk, jt

   CALL get_environment_variable("PROC_X", env)
   READ ( env, '(i10)' ) px
   CALL get_environment_variable("PROC_Y", env)
   READ ( env, '(i10)' ) py
   CALL get_environment_variable("JPI", env)
   READ ( env, '(i10)' ) jpi
   CALL get_environment_variable("JPJ", env)
   READ ( env, '(i10)' ) jpj
   CALL get_environment_variable("JPK", env)
   READ ( env, '(i10)' ) jpk

   CALL get_environment_variable("CPU_PER_TASK", env)
   READ ( env, '(i10)' ) cpt 

   CALL mpi_init( ierr )

   CALL mpi_comm_rank( MPI_COMM_WORLD, myrank, ierr )
   CALL mpi_comm_size( MPI_COMM_WORLD, p, ierr )
      WRITE (*,*) "size=",p," PROC_X=",px, "PROC_Y", py
      color = MOD(myrank,cpt)
   CALL mpi_comm_split(MPI_COMM_WORLD, color, myrank, realcomm, ierr)
   IF (color == 0) THEN
   CALL mpi_comm_rank( realcomm, myrank, ierr )
   CALL mpi_comm_size( realcomm, p, ierr )
   
   IF (myrank == 0) THEN
      WRITE (*,*) "(PROC_X, PROC_Y)(",px,",",py,")"
      WRITE (*,*) "(jpi, jpj, jpk)(",jpi,",",jpj,",",jpk,")"
   END IF

   IF (p .ne. px*py) THEN
      WRITE (*,*) "Error! the total number of processes does not match with the values of PROC_X and PROC_Y env variables"
      WRITE (*,*) "size=",p," PROC_X=",px, "PROC_Y", py
      CALL EXIT()
   ENDIF

   hpmid = myrank + 10
   nowe = myrank - px
   noea = myrank + px
   noso = myrank - 1
   nono = myrank + 1
   IF (py == 1) THEN
      nbondi = 2
   ELSE IF (myrank < px) THEN
      nbondi = -1
   ELSE IF (myrank >= px*(py-1)) THEN
      nbondi = 1
   ELSE
      nbondi = 0
   END IF
      
   IF (px == 1) THEN
      nbondj = 2
   ELSE IF (MOD(myrank, px) == px-1) THEN
      nbondj = 1
   ELSE IF (MOD(myrank, px) == 0) THEN
      nbondj = -1
   ELSE
      nbondj = 0
   END IF
     

   jpreci=1
   jprecj=1
   nlci = jpi
   nlcj = jpj
   nreci = 2*jpreci
   nrecj = 2*jprecj
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
   ALLOCATE( t3ew(jpj,jpreci,jpk,2))
   ALLOCATE( t3we(jpj,jpreci,jpk,2))
   ALLOCATE( t3ns(jpi,jprecj,jpk,2))
   ALLOCATE( t3sn(jpi,jprecj,jpk,2))



! array initialization
   DO jk = 1, jpk
      DO jj = 1, jpj
          DO ji = 1, jpi
              umask(ji,jj,jk) = ji*jj*jk
              mydomain(ji,jj,jk) =ji*jj*jk
              pun(ji,jj,jk) =ji*jj*jk
              pvn(ji,jj,jk) =ji*jj*jk
              pwn(ji,jj,jk) =ji*jj*jk
              vmask(ji,jj,jk)= ji*jj*jk
              tsn(ji,jj,jk)= ji*jj*jk
              tmask(ji,jj,jk)= ji*jj*jk
          END DO
      END DO
   END DO

   DO jj=1, jpj
      DO ji=1, jpi
         ztfreez(ji,jj)=ji*jj
         upsmsk(ji,jj)=ji*jj
         rnfmsk(ji,jj) = ji*jj
      END DO
   END DO

   DO jk=1, jpk
      rnfmsk_z(jk)=jk
   END DO


!***********************
!* Start of the synphony
!***********************
   DO jk = 1, jpk
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( tsn(ji,jj,jk) <= ztfreez(ji,jj) + 0.1 ) THEN   ;   zice = 1.e0
            ELSE                                               ;   zice = 0.e0
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

   zwx(:,:,jpk) = 0.e0   ;   zwy(:,:,jpk) = 0.e0

   DO jk = 1, jpk-1
      DO jj = 1, jpj-1
         DO ji = 1, jpi-1
             zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
             zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
         END DO
      END DO
   END DO
   CALL mpp_lnk_3d(zwx)
   CALL mpp_lnk_3d(zwy)

   zslpx(:,:,jpk) = 0.e0   ;   zslpy(:,:,jpk) = 0.e0 

   DO jk = 1, jpk-1
     DO jj = 2, jpj
        DO ji = 2, jpi 
           zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji-1,jj  ,jk) )   &
           &            * ( 0.25 + SIGN( 0.25, zwx(ji,jj,jk) * zwx(ji-1,jj  ,jk) ) )
           zslpy(ji,jj,jk) =                    ( zwy(ji,jj,jk) + zwy(ji  ,jj-1,jk) )   &
           &            * ( 0.25 + SIGN( 0.25, zwy(ji,jj,jk) * zwy(ji  ,jj-1,jk) ) )
        END DO
     END DO
  END DO

  DO jk = 1, jpk-1    
     DO jj = 2, jpj
        DO ji = 2, jpi
           zslpx(ji,jj,jk) = SIGN( 1., zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji  ,jj,jk) ),   &
           &                                                 2.*ABS( zwx  (ji-1,jj,jk) ),   &
           &                                                 2.*ABS( zwx  (ji  ,jj,jk) ) )
           zslpy(ji,jj,jk) = SIGN( 1., zslpy(ji,jj,jk) ) * MIN(    ABS( zslpy(ji,jj  ,jk) ),   &
           &                                                 2.*ABS( zwy  (ji,jj-1,jk) ),   &
           &                                                 2.*ABS( zwy  (ji,jj  ,jk) ) )
        END DO
     END DO
  END DO 

  DO jk = 1, jpk-1
     zdt  = 1
     DO jj = 2, jpj-1
        DO ji = 2, jpi-1
            z0u = SIGN( 0.5, pun(ji,jj,jk) )
            zalpha = 0.5 - z0u
            zu  = z0u - 0.5 * pun(ji,jj,jk) * zdt

            zzwx = mydomain(ji+1,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji+1,jj,jk))
            zzwy = mydomain(ji  ,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji  ,jj,jk))

            zwx(ji,jj,jk) = pun(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
            
            z0v = SIGN( 0.5, pvn(ji,jj,jk) )
            zalpha = 0.5 - z0v
            zv  = z0v - 0.5 * pvn(ji,jj,jk) * zdt

            zzwx = mydomain(ji,jj+1,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj+1,jk))
            zzwy = mydomain(ji,jj  ,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj  ,jk))

            zwy(ji,jj,jk) = pvn(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
         END DO
      END DO
  END DO

  CALL mpp_lnk_3d(zwx)
  CALL mpp_lnk_3d(zwy)

  DO jk = 1, jpk-1
     DO jj = 2, jpj-1     
        DO ji = 2, jpi-1
           zbtr = 1.
           ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
           &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
           mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
        END DO
     END DO
  END DO

  zwx (:,:, 1 ) = 0.e0    ;    zwx (:,:,jpk) = 0.e0

  DO jk = 2, jpk-1   
     zwx(:,:,jk) = tmask(:,:,jk) * ( mydomain(:,:,jk-1) - mydomain(:,:,jk) )
  END DO

  zslpx(:,:,1) = 0.e0

  DO jk = 2, jpk-1    
     DO jj = 1, jpj
        DO ji = 1, jpi
           zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
           &            * ( 0.25 + SIGN( 0.25, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )
        END DO
     END DO
  END DO

  DO jk = 2, jpk-1     
     DO jj = 1, jpj
        DO ji = 1, jpi
           zslpx(ji,jj,jk) = SIGN( 1., zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji,jj,jk  ) ),   &
           &                                                 2.*ABS( zwx  (ji,jj,jk+1) ),   &
           &                                                 2.*ABS( zwx  (ji,jj,jk  ) )  )
        END DO
     END DO
  END DO

  zwx(:,:, 1 ) = pwn(:,:,1) * mydomain(:,:,1)

  zdt  = 1
  zbtr = 1.
  DO jk = 1, jpk-1
     DO jj = 2, jpj-1     
        DO ji = 2, jpi-1
           z0w = SIGN( 0.5, pwn(ji,jj,jk+1) )
           zalpha = 0.5 + z0w
           zw  = z0w - 0.5 * pwn(ji,jj,jk+1) * zdt * zbtr

           zzwx = mydomain(ji,jj,jk+1) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk+1))
           zzwy = mydomain(ji,jj,jk  ) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk  ))

           zwx(ji,jj,jk+1) = pwn(ji,jj,jk+1) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
        END DO
     END DO
  END DO

  zbtr = 1.
  DO jk = 1, jpk-1
     DO jj = 2, jpj-1     
        DO ji = 2, jpi-1
           ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
           mydomain(ji,jj,jk) =  mydomain(ji,jj,jk) + ztra
        END DO
     END DO
  END DO

  END IF

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
  DEALLOCATE( t3ew)
  DEALLOCATE( t3we)
  DEALLOCATE( t3ns)
  DEALLOCATE( t3sn)

   CALL mpi_finalize(ierr)

CONTAINS 

   SUBROUTINE mpp_lnk_3d( ptab )
      REAL, DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptab
      !!
      INTEGER  ::   ji, jj, jk, jl
      INTEGER  ::   imigr, iihom, ijhom
      INTEGER  ::   ml_req1, ml_req2, ml_err
      INTEGER, DIMENSION(mpi_status_size) ::   ml_stat
      INTEGER :: iflag
      INTEGER :: istatus(mpi_status_size)
      DOUBLE PRECISION :: t1, t2

      SELECT CASE ( nbondi )
      CASE ( -1, 0, 1 )     
         iihom = nlci-nreci
         DO jl = 1, jpreci
            t3ew(:,jl,:,1) = ptab(jpreci+jl,:,:)
            t3we(:,jl,:,1) = ptab(iihom +jl,:,:)
         END DO
      END SELECT  
      imigr = jpreci * jpj * jpk
      SELECT CASE ( nbondi ) 
      CASE ( -1 )
         CALL mpi_isend( t3we(1,1,1,1), imigr, mpi_double_precision, noea , 0, realcomm, ml_req1, iflag )
         CALL mpi_recv( t3ew(1,1,1,2), imigr, mpi_double_precision, noea, 0, realcomm, istatus, iflag )

         CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mpi_isend( t3ew(1,1,1,1), imigr, mpi_double_precision, nowe , 0, realcomm, ml_req1, iflag )
         CALL mpi_isend( t3we(1,1,1,1), imigr, mpi_double_precision, noea , 0, realcomm, ml_req2, iflag )
         
         CALL mpi_recv( t3ew(1,1,1,2), imigr, mpi_double_precision, noea, 0, realcomm, istatus, iflag )
         CALL mpi_recv( t3we(1,1,1,2), imigr, mpi_double_precision, nowe, 0, realcomm, istatus, iflag )

         CALL mpi_wait(ml_req1, ml_stat, ml_err)
         CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 )
         CALL mpi_isend( t3ew(1,1,1,1), imigr, mpi_double_precision, nowe , 0, realcomm, ml_req1, iflag )
         CALL mpi_recv( t3we(1,1,1,2), imigr, mpi_double_precision, nowe, 0, realcomm, istatus, iflag )
         CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT

      iihom = nlci-jpreci
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, jpreci
            ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
         END DO
      CASE ( 0 ) 
         DO jl = 1, jpreci
            ptab(jl      ,:,:) = t3we(:,jl,:,2)
            ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
         END DO
      CASE ( 1 )
         DO jl = 1, jpreci
            ptab(jl      ,:,:) = t3we(:,jl,:,2)
         END DO
      END SELECT

      IF( nbondj /= 2 ) THEN
         ijhom = nlcj-nrecj
         DO jl = 1, jprecj
            t3sn(:,jl,:,1) = ptab(:,ijhom +jl,:)
            t3ns(:,jl,:,1) = ptab(:,jprecj+jl,:)
         END DO
      ENDIF
      imigr = jprecj * jpi * jpk
      SELECT CASE ( nbondj )     
      CASE ( -1 )
         CALL mpi_isend( t3sn(1,1,1,1), imigr, mpi_double_precision, nono , 0, realcomm, ml_req1, iflag )
         CALL mpi_recv( t3ns(1,1,1,2), imigr, mpi_double_precision, nono, 0, realcomm, istatus, iflag )

         CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mpi_isend( t3ns(1,1,1,1), imigr, mpi_double_precision, noso , 0, realcomm, ml_req1, iflag )
         CALL mpi_isend( t3sn(1,1,1,1), imigr, mpi_double_precision, nono , 0, realcomm, ml_req2, iflag )
         CALL mpi_recv( t3ns(1,1,1,2), imigr, mpi_double_precision, nono, 0, realcomm, istatus, iflag )
         CALL mpi_recv( t3sn(1,1,1,2), imigr, mpi_double_precision, noso, 0, realcomm, istatus, iflag )

         CALL mpi_wait(ml_req1, ml_stat, ml_err)
         CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 ) 
         CALL mpi_isend( t3ns(1,1,1,1), imigr, mpi_double_precision, noso , 0, realcomm, ml_req1, iflag )
         CALL mpi_recv( t3sn(1,1,1,2), imigr, mpi_double_precision, noso, 0, realcomm, istatus, iflag )

         CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      ijhom = nlcj-jprecj
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, jprecj
            ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
         END DO
      CASE ( 0 ) 
         DO jl = 1, jprecj
            ptab(:,jl      ,:) = t3sn(:,jl,:,2)
            ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
         END DO
      CASE ( 1 )
         DO jl = 1, jprecj
            ptab(:,jl,:) = t3sn(:,jl,:,2)
         END DO
      END SELECT

   END SUBROUTINE mpp_lnk_3d

END PROGRAM tra_adv
