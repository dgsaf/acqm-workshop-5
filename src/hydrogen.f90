!>
module hydrogen

  implicit none

contains

  ! hydrogen_wf
  !
  ! Hydrogen wavefunctions and energies for given l, alpha, n_b, n_r, r_grid
  subroutine hydrogen_wf (l, alpha, n_b, n_r, r_grid, wf, energies, status)
    integer , intent(in) :: l
    double precision , intent(in) :: alpha
    integer , intent(in) :: n_b
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(out) :: wf(n_r, n_b)
    double precision , intent(out) :: energies(n_b)
    integer , intent(out) :: status
    double precision :: basis(n_r, n_b)
    double precision :: B(n_b, n_b), K(n_b, n_b), V(n_b, n_b), H(n_b, n_b)
    double precision :: eigenvectors(n_b, n_b)
    integer :: ii, jj, kk

    ! check if arguments are valid
    status = 0

    if ((l < 0) .or. (n_b < 1) .or. (alpha < 0.0d0) .or. (n_r < 1)) then
      status = 1
      return
    end if

    ! setup basis
    call radial_basis(l, alpha, n_b, n_r, r_grid, basis, status)
    if (status /= 0) return

    ! setup matrices
    call overlap_matrix(l, n_b, B, status)
    if (status /= 0) return

    call kinetic_matrix(l, alpha, n_b, K, status)
    if (status /= 0) return

    call coulomb_matrix(l, alpha, n_b, V, status)
    if (status /= 0) return

    H(:, :) = K(:, :) + V(:, :)

    ! diagonalise
    energies(:) = 0d0
    eigenvectors(:, :) = 0d0

    call rsg(n_b, n_b, H, B, energies, 1, eigenvectors, status)
    if (status /= 0) return

    ! calculate wavefunctions
    wf(:, :) = 0d0

    do jj = 1, n_b
      do ii = 1, n_r
        do kk = 1, n_b
          wf(ii, jj) = wf(ii, jj)  + (eigenvectors(kk, jj) * basis(ii, kk))
        end do
      end do
    end do

  end subroutine hydrogen_wf

  ! radial_basis
  !
  ! For given l, alpha, and r_grid, yields the functions varphi_{k, l}(r) for
  ! k = 1, ..., n_b, on the radial values specified in the grid.
  subroutine radial_basis (l, alpha, n_b, n_r, r_grid, basis, status)
    integer , intent(in) :: l
    double precision , intent(in) :: alpha
    integer , intent(in) :: n_b
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(out) :: basis(n_r, n_b)
    integer , intent(out) :: status
    double precision :: norm(n_b)
    double precision :: alpha_grid(n_r)
    integer :: kk

    ! check if arguments are valid
    status = 0

    if ((l < 0) .or. (n_b < 1) .or. (n_r < 1)) then
      status = 1
      return
    end if

    ! recurrence relation for basis normalisation constants
    norm(1) = sqrt(alpha / dble((l + 1d0) * gamma(dble((2d0 * l) + 2d0))))

    if (n_b >= 2) then
      do kk = 2, n_b
        norm(kk) = norm(kk-1) * sqrt(dble((kk - 1d0) * (kk - 1d0 + l)) / &
            dble((kk + l) * (kk + (2d0 * l))))
      end do
    end if

    ! in-lined array since r_grid(:) on its own is never used
    alpha_grid(:) = alpha * r_grid(:)

    ! recurrence relation for basis functions
    basis(:, 1) = ((2.0d0 * alpha_grid(:)) ** (l + 1)) * &
        exp(-alpha_grid(:))

    if (n_b >= 2) then
      basis(:, 2) = 2.0d0 * (dble(l + 1d0) - alpha_grid(:)) * basis(:, 1)
    end if

    if (n_b >= 3) then
      do kk = 3, n_b
        basis(:, kk) = &
            ((2.0d0 * (dble(kk - 1d0 + l) - alpha_grid(:)) * basis(:, kk-1)) &
            - dble(kk + (2d0 * l) - 1d0) * basis(:, kk-2)) / dble(kk - 1)
      end do
    end if

    ! scaling basis functions by normalisation constants
    do kk = 1, n_b
      basis(:, kk) = basis(:, kk) * norm(kk)
    end do

  end subroutine radial_basis

  ! overlap_matrix
  !
  ! Overlap matrix elements for given l.
  subroutine overlap_matrix(l, n_b, B, status)
    integer , intent(in) :: l
    integer , intent(in) :: n_b
    double precision , intent(out) :: B(n_b, n_b)
    integer , intent(out) :: status
    integer :: kk

    ! check if arguments are valid
    status = 0

    if ((l < 0) .or. (n_b < 1)) then
      status = 1
      return
    end if

    ! initialise overlap matrix to zero
    B(:, :) = 0.0d0

    ! determine tri-diagonal overlap matrix elements
    do kk = 1, n_b-1
      B(kk, kk) = 1.0d0

      B(kk, kk+1) = - 0.5d0 * sqrt(1 - &
          (dble(l * (l + 1)) / dble((kk + l) * (kk + l + 1))))

      B(kk+1, kk) = B(kk, kk+1)
    end do

    ! last term (not covered by loop)
    B(n_b, n_b) = 1.0d0

  end subroutine overlap_matrix

  ! kinetic_matrix
  !
  ! Kinetic matrix elements for given l, alpha.
  subroutine kinetic_matrix(l, alpha, n_b, K, status)
    integer , intent(in) :: l
    integer , intent(in) :: n_b
    double precision , intent(in) :: alpha
    double precision , intent(out) :: K(n_b, n_b)
    integer , intent(out) :: status
    integer :: kk

    ! check if arguments are valid
    status = 0

    if ((l < 0) .or. (n_b < 1)) then
      status = 1
      return
    end if

    ! initialise kinetic matrix to zero
    K(:, :) = 0.0d0

    ! determine tri-diagonal kinetic matrix elements
    do kk = 1, n_b-1
      K(kk, kk) = 0.5d0 * (alpha ** 2)

      K(kk, kk+1) = (alpha ** 2) * 0.25d0 * sqrt(1 - &
          (dble(l * (l + 1)) / dble((kk + l) * (kk + l + 1))))

      K(kk+1, kk) = K(kk, kk+1)
    end do

    ! last term (not covered by loop)
    K(n_b, n_b) = 0.5d0 * (alpha ** 2)

  end subroutine kinetic_matrix

  ! potential_matrix
  !
  ! Hydrogen potential matrix elements for given l, alpha.
  subroutine potential_matrix(l, alpha, n_b, V, status)
    integer , intent(in) :: l
    integer , intent(in) :: n_b
    double precision , intent(in) :: alpha
    double precision , intent(out) :: V(n_b, n_b)
    integer , intent(out) :: status
    integer :: kk

    ! check if arguments are valid
    status = 0

    if ((l < 0) .or. (n_b < 1)) then
      status = 1
      return
    end if

    ! initialise coulomb matrix to zero
    V(:, :) = 0.0d0

    ! determine diagonal coulomb matrix elements
    do kk = 1, n_b
      V(kk, kk) = - alpha / dble(kk + l)
    end do

  end subroutine potential_matrix

end module hydrogen
