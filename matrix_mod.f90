! Module to define the derived type matrix
!
! E. Mucciolo & C. Chamon (Aug/01/2012)
!
 module matrix_mod
!
   use double_precision_mod
!
   type matrix
      real (kind=double), allocatable, dimension(:,:,:) :: M
   end type matrix
!
 end module matrix_mod
