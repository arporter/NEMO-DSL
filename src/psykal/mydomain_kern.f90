module mydomain_kern_mod
!
implicit none
!
contains
  !
  subroutine mydomain_kern(mydomain,zbtr,zwx,jk,jj,ji)
    !
    real*8, intent(out) :: mydomain(:,:,:)
    real*8, intent(in)  :: zwx(:,:,:)
    real*8, intent(in) :: zbtr
    integer, intent(in) :: jk,jj,ji
    ! local variables
    real*8 :: ztra
    !
    ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
    mydomain(ji,jj,jk) = ztra
    !
  end subroutine mydomain_kern
  !
end module mydomain_kern_mod
