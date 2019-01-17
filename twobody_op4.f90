! Subroutine to implement a two-body operation on consecutive spins on
! a chain.
!
! no truncation: trmode = 'f'
! truncation at D > Dcap: trmode = 't' 
!
! E. Mucciolo & C. Chamon (February 22, 2018)
!
 subroutine twobody_op4 (A, B, MA, R, NB, Atil, Btil, Dmax, Dcap, Rtil, &
      Rseq, operat2, alpha, eps, trmode, ifail)
!
 use double_precision_mod
!
 implicit none
!
! external variables
 integer, intent(in) :: MA, R, NB, Dmax, Dcap, Rtil
 integer, intent(out) :: Rseq
 real (kind=double), intent(in), dimension(MA,R,0:1) :: A
 real (kind=double), intent(in), dimension(R,NB,0:1) :: B
 real (kind=double), intent(out), dimension(MA,Rtil,0:1) :: Atil
 real (kind=double), intent(out), dimension(Rtil,NB,0:1) :: Btil
 real (kind=double), intent(in) :: eps, alpha
 integer, intent(out) :: ifail
 character (len=1) :: trmode
!
! external subroutine containing the two-body operator definition
 external :: operat2
!
! internal variables
 integer :: Raux, Rnew, Msize, Nsize
 integer :: INFO
 integer, allocatable, dimension(:) :: seq
 real (kind=double), allocatable, dimension(:,:) :: matrix
 real (kind=double), allocatable, dimension(:) :: S
 real (kind=double), allocatable, dimension(:,:) :: U, V
 integer :: ifail_trun
 integer :: k
!
! determine enlarged matrix dimensions
 Msize = 2*MA
 Nsize = 2*NB
!
! allocate space for auxiliary matrix
 allocate (matrix(Msize,Nsize))
!
! apply operator and obtain resulting matrix
 call operat2 (MA, R, NB, A, B, alpha, matrix)
!
! allocate matrices and arrays
 allocate (U(Msize,Msize),V(Nsize,Nsize))
 allocate (S(Rtil))
!
! do the singular value decomposition
 call singvaldec (matrix, Msize, Nsize, S, U, V, INFO)
!
! check decomposition
 if (INFO /= 0) then
    ifail = 10
 else
    ifail = 0
 end if
!
! deallocate auxiliary matrix
 deallocate (matrix)
!
! find the new matrix rank (new bond dimension)
 Raux = 1
 do k = 2, Rtil
    if (S(k)/S(1) > eps) then
       Raux = Raux + 1
    end if
 end do
!
! reset negligible singular values to zero
 do k = Raux+1, Rtil
    S(k) = 0.d0
 end do
!
! check limit to bond dimension and pick the number of active singular values
 if (Raux > Dmax) then
    Rnew = Dmax
    do k = Dmax+1, Raux
       S(k) = 0.d0
    end do
 else
    Rnew = Raux
 end if
!
! no truncation
 if (trmode == 'f') then
    Rseq = Rnew
    allocate (seq(Rseq))
    do k = 1, Rseq
       seq(k) = k
    end do
!
! truncation of lowest singular values
 else if (trmode == 't') then
    Rseq = min(Rnew,Dcap)
    allocate (seq(Rseq))
    do k = 1, Rseq
       seq(k) = k
    end do
!
 else
    ifail = 15
    Rseq = 1
    allocate (seq(Rseq))
    seq(1) = 1
!
 end if
!
! recompose the matrix product using selected singular values
 call factorize_seq (S, U, V, MA, Rtil, Rseq, NB, seq, Atil, Btil, ifail_trun)
 if (ifail_trun /= 0) then
    ifail = 20 + ifail_trun
 end if
!
! deallocate matrices and arrays
 deallocate (S)
 deallocate (U,V)
 deallocate (seq)
!
! final check
 if (ifail /= 0) then
    write(54,'("ifail(2body) = ",i2," MA,R,NB =",i4,1x,i4,1x,i4,1x,i4)') &
         ifail, MA, Rnew, Rseq, NB
 end if
!
 end subroutine twobody_op4
