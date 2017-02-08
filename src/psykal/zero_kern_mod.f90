module zero_kern_mod
!
implicit none
!
contains
  !
  subroutine zero_kern_2d(field,jj,ji)
    !
    real*8, intent(out) :: field
    !
    field(jj,ji) = 0.0
    !
  end subroutine zero_kern_2d
!
end module zero_kern_mod
