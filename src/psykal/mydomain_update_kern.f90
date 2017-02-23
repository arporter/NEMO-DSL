module mydomain_update_kern_mod
!
implicit none
!
contains
  !
  subroutine mydomain_update_kern(mydomain,zwx,zwy,ji,jj,jk)
    !
    real*8,  intent(inout) :: mydomain(:,:,:)
    real*8,  intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: ji,jj,jk
    ! local variables
    real*8 :: zbtr,ztra
    !
    zbtr = 1.d0
    ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
         &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
    mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
    !
  end subroutine mydomain_update_kern
  !
end module mydomain_update_kern_mod
