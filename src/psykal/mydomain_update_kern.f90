module mydomain_update_kern_mod
!
implicit none
!
contains
  !
  subroutine mydomain_update_kern(mydomain,zwx,zwy,jk,jj,ji)
    !
    real*8,  intent(inout) :: mydomain(:,:,:)
    real*8,  intent(in)  :: zwx(:,:,:), zwy(:,:,:)
    integer, intent(in) :: jk,jj,ji
    ! local variables
    real*8 :: zbtr,ztra
    !
    zbtr = 1.
    ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
         &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
    mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
    !
  end subroutine mydomain_update_kern
  !
end module mydomain_update_kern_mod
