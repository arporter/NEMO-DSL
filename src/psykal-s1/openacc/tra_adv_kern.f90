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
  subroutine zind_compute(zind,tsn,ztfreez,rnfmsk,rnfmsk_z,upsmsk,tmask)
    !
    real*8, intent(out) :: zind(jpi,jpj,jpk)
    real*8, intent(in)  :: tsn(jpi,jpj,jpk), tmask(jpi,jpj,jpk)
    real*8, intent(in)  :: ztfreez(jpi,jpj), rnfmsk(jpi,jpj), upsmsk(jpi,jpj)
    real*8, intent(in)  :: rnfmsk_z(jpk)
    !
    real*8 :: zice
    
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 1, jpk
       !$acc loop private(jj)
       DO jj = 1, jpj
          !$acc loop private(ji)
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
    !$acc end kernels
    !
  end subroutine zind_compute
  !
  subroutine zero_bottom_layer(field)
    !
    real*8, intent(out) :: field(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jj)
    DO jj = 1, jpj
       !$acc loop private(ji)
       DO ji = 1, jpi
          field(ji,jj,jpk) = 0.0d0
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zero_bottom_layer
  !
  subroutine zwxy_compute(zwx,zwy,umask,vmask,mydomain)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    real*8, intent(in) :: umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk), mydomain(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 1, jpk-1
       !$acc loop private(jj)
       DO jj = 1, jpj-1
          !$acc loop private(ji)
          DO ji = 1, jpi-1
             zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
             zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zwxy_compute
  !
  subroutine zslpxy_compute(zslpx,zslpy,zwx,zwy)
    !
    real*8, intent(out) :: zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    real*8, intent(in) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 1, jpk-1
       !$acc loop private(jj)
       DO jj = 2, jpj
          !$acc loop private(ji)
          DO ji = 2, jpi 
             zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji-1,jj  ,jk) )   &
                  &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji-1,jj  ,jk) ) )
             zslpy(ji,jj,jk) =                    ( zwy(ji,jj,jk) + zwy(ji  ,jj-1,jk) )   &
                  &            * ( 0.25d0 + SIGN( 0.25d0, zwy(ji,jj,jk) * zwy(ji  ,jj-1,jk) ) )
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zslpxy_compute
  !
  subroutine zslpxy_update_compute(zslpx,zslpy,zwx,zwy)
    !
    real*8, intent(inout) :: zslpx(jpi,jpj,jpk), zslpy(jpi,jpj,jpk)
    real*8, intent(in) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 1, jpk-1    
       !$acc loop private(jj)
       DO jj = 2, jpj
          !$acc loop private(ji)
          DO ji = 2, jpi
             zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji  ,jj,jk) ),   &
                  &                                                2.d0*ABS( zwx  (ji-1,jj,jk) ),   &
                  &                                                2.d0*ABS( zwx  (ji  ,jj,jk) ) )
             zslpy(ji,jj,jk) = SIGN( 1.d0, zslpy(ji,jj,jk) ) * MIN(    ABS( zslpy(ji,jj  ,jk) ),   &
                  &                                                2.d0*ABS( zwy  (ji,jj-1,jk) ),   &
                  &                                                2.d0*ABS( zwy  (ji,jj  ,jk) ) )
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zslpxy_update_compute
  !
  subroutine zwxy2_compute(zwx,zwy,pun,pvn,mydomain,zind,zslpx,zslpy)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    real*8, intent(in), dimension(jpi,jpj,jpk)  :: pun, pvn, mydomain, zind, zslpx, zslpy
    !
    real*8 :: zdt, z0u, zalpha, zu, zzwx, zzwy, z0v, zv
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 1, jpk-1
       zdt  = 1
       !$acc loop private(jj)
       DO jj = 2, jpj-1
          !$acc loop private(ji)
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
    !$acc end kernels
    !
  end subroutine zwxy2_compute
  !
  subroutine mydomain_update_compute(mydomain,zwx,zwy)
    !
    real*8, intent(inout) :: mydomain(jpi,jpj,jpk)
    real*8, intent(in)    :: zwx(jpi,jpj,jpk), zwy(jpi,jpj,jpk)
    !
    real*8 :: zbtr, ztra
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 1, jpk-1
       !$acc loop private(jj)
       DO jj = 2, jpj-1     
          !$acc loop private(ji)
          DO ji = 2, jpi-1
             zbtr = 1.
             ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                  &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
             mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine mydomain_update_compute
  !
  subroutine zero_top_layer(field)
    !
    real*8, intent(out) :: field(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jj)
    DO jj = 1, jpj
       !$acc loop private(ji)
       DO ji = 1, jpi
          field(ji,jj,1) = 0.0d0
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zero_top_layer
  !
  subroutine zwx_compute(zwx,tmask,mydomain)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk)
    real*8, intent(in)  :: tmask(jpi,jpj,jpk), mydomain(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 2, jpk-1
       !$acc loop private(jj)
       DO jj = 1, jpj
          !$acc loop private(ji)
          DO ji = 1, jpi
             zwx(ji,jj,jk) = tmask(ji,jj,jk) * ( mydomain(ji,jj,jk-1) - mydomain(ji,jj,jk) )
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zwx_compute
  !
  subroutine zslpx_compute(zslpx,zwx)
    !
    real*8, intent(out) :: zslpx(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 2, jpk-1    
       !$acc loop private(jj)
       DO jj = 1, jpj
          !$acc loop private(ji)
          DO ji = 1, jpi
             zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
                  &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zslpx_compute
  !
  subroutine zslpx_update_compute(zslpx,zwx)
    !
    real*8, intent(inout) :: zslpx(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 2, jpk-1     
       !$acc loop private(jj)
       DO jj = 1, jpj
          !$acc loop private(ji)
          DO ji = 1, jpi
             zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN( ABS( zslpx(ji,jj,jk  ) ), &
                  &                                               2.d0*ABS( zwx  (ji,jj,jk+1) ),   &
                  &                                               2.d0*ABS( zwx  (ji,jj,jk  ) )  )
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine zslpx_update_compute
  !
  subroutine multiply_top_layer(field1,field2,field3)
    !
    real*8, intent(out) :: field1(jpi,jpj,jpk)
    real*8, intent(in)  :: field2(jpi,jpj,jpk), field3(jpi,jpj,jpk)
    !
    !$acc kernels
    !$acc loop private(jj)
    DO jj = 1, jpj
       !$acc loop private(ji)
       DO ji = 1, jpi
          field1(ji,jj,1) = field2(ji,jj,1) * field3(ji,jj,1)
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine multiply_top_layer
  !
  subroutine zwx2_compute(zwx,pwn,mydomain,zind,zslpx)
    !
    real*8, intent(out) :: zwx(jpi,jpj,jpk)
    real*8, intent(in), dimension(jpi,jpj,jpk)  :: pwn, mydomain, zind, zslpx
    !
    real*8 :: zdt, zbtr, z0w, zalpha, zw, zzwx, zzwy
    !
    zdt  = 1
    zbtr = 1.
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 2, jpk
       !$acc loop private(jj)
       DO jj = 2, jpj-1     
          !$acc loop private(ji)
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
    !$acc end kernels
    !
  end subroutine zwx2_compute
  !
  subroutine mydomain_compute(mydomain,zwx)
    !
    real*8, intent(out) :: mydomain(jpi,jpj,jpk)
    real*8, intent(in)  :: zwx(jpi,jpj,jpk)
    !
    real*8 :: zbtr, ztra
    !
    zbtr = 1.
    !$acc kernels
    !$acc loop private(jk)
    DO jk = 1, jpk-1
       !$acc loop private(jj)
       DO jj = 2, jpj-1     
          !$acc loop private(ji)
          DO ji = 2, jpi-1
             ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
             mydomain(ji,jj,jk) = ztra
          END DO
       END DO
    END DO
    !$acc end kernels
    !
  end subroutine mydomain_compute
  !
end module tra_adv_kern
