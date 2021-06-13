!>
module vmatrix

  implicit none

contains

  ! direct_matrix
  !
  ! Potential direct matrix elements.
  ! Note:
  ! - delta_i_f indicates if i==f
  subroutine direct_matrix (n_r, r_grid, r_weights, delta_i_f, wf_i, wf_f, &
      n_k, k_grid, V_direct, status)
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: r_weights(n_r)
    logical , intent(in) :: delta_i_f
    double precision , intent(in) :: wf_i(n_r)
    double precision , intent(in) :: wf_f(n_r)
    integer , intent(in) :: n_k
    double precision , intent(in) :: k_grid(n_k)
    double precision , intent(out) :: V_direct(n_k, n_k)
    integer , intent(out) :: status
    ! local variables
    double precision :: f(n_r), g(n_r)
    double precision :: k_i(n_r), k_f(n_r)
    integer :: ii, ff

    ! check if arguments are valid
    if ((n_r < 1) .or. (n_k < 1)) then
      status = 1
      return
    end if

    ! initialise V_direct
    V_direct(:, :) = 0.0d0

    ! determine direct matrix elements numerically
    do ii = 1, n_k
      do ff = 1, n_k
        call direct_matrix_element(n_r, r_grid, r_weights, delta_i_f, &
            wf_i, wf_f, k_grid(ii), k_grid(ff), V_direct(ff, ii), status)
      end do
    end do

  end subroutine direct_matrix

  ! direct_matrix
  !
  ! Potential direct matrix elements for given ki, kf.
  ! Note:
  ! - delta_i_f indicates if i==f
  subroutine direct_matrix_element (n_r, r_grid, r_weights, &
      delta_i_f, wf_i, wf_f, ki, kf, V_direct, status)
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: r_weights(n_r)
    logical , intent(in) :: delta_i_f
    double precision , intent(in) :: wf_i(n_r)
    double precision , intent(in) :: wf_f(n_r)
    double precision , intent(in) :: ki
    double precision , intent(in) :: kf
    double precision , intent(out) :: V_direct
    integer , intent(out) :: status
    ! local variables
    ! double precision :: k_i(n_r), k_f(n_r)
    double precision :: f(n_r), g(n_r)
    integer :: ii, ff

    ! check if arguments are valid
    if (n_r < 1) then
      status = 1
      return
    end if

    ! determine direct matrix elements numerically
    f(:) = sin(kf * r_grid(:)) * sin(ki * r_grid(:))
    g(:) = wf_f(:) * wf_i(:)

    V_direct = two_e_integral(n_r, r_grid, r_weights, f, g)

    ! if i == f, then subtract (1/r) term
    if (delta_i_f) then
      V_direct = V_direct - integrate(n_r, r_weights, f(:) / r_grid(:))
    end if

  end subroutine direct_matrix_element

  ! exchange_matrix
  !
  ! Potential exchange matrix elements.
  subroutine exchange_matrix (n_r, r_grid, r_weights, energy, wf_i, wf_f, &
      n_k, k_grid, V_exchange, status)
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: r_weights(n_r)
    double precision , intent(in) :: energy
    double precision , intent(in) :: wf_i(n_r)
    double precision , intent(in) :: wf_f(n_r)
    integer , intent(in) :: n_k
    double precision , intent(in) :: k_grid(n_k)
    double precision , intent(out) :: V_exchange(n_k, n_k)
    integer , intent(out) :: status
    ! local variables
    double precision :: v(n_r)
    double precision :: a, b, c, d
    double precision :: k_i(n_r), k_f(n_r)
    double precision :: ki, kf
    integer :: ii, ff

    ! check if arguments are valid
    if ((n_r < 1) .or. (n_k < 1)) then
      status = 1
      return
    end if

    ! initialise V_exchange
    V_exchange(:, :) = 0.0d0

    ! determine direct matrix elements numerically
    do ii = 1, n_k
      do ff = 1, n_k
        call exchange_matrix_element(n_r, r_grid, r_weights, energy, &
            wf_i, wf_f, k_grid(ii), k_grid(ff), V_exchange(ff, ii), status)
      end do
    end do

  end subroutine exchange_matrix

  ! exchange_matrix_element
  !
  ! Potential exchange matrix element for given ki, kf
  subroutine exchange_matrix_element (n_r, r_grid, r_weights, &
      energy, wf_i, wf_f, ki, kf, V_exchange, status)
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: r_weights(n_r)
    double precision , intent(in) :: energy
    double precision , intent(in) :: wf_i(n_r)
    double precision , intent(in) :: wf_f(n_r)
    double precision , intent(in) :: ki, kf
    double precision , intent(out) :: V_exchange
    integer , intent(out) :: status
    ! local variables
    double precision :: v(n_r)
    double precision :: a, b, c, d
    ! double precision :: k_i(n_r), k_f(n_r)
    double precision :: f(n_r), g(n_r)
    integer :: ii, ff

    ! check if arguments are valid
    if (n_r < 1) then
      status = 1
      return
    end if

    ! cache one-electron potential
    v(:) = - 1.0d0 / r_grid(:)

    f(:) = sin(kf * r_grid(:)) * wf_i(:)
    g(:) = sin(ki * r_grid(:)) * wf_f(:)

    V_exchange = &
        (energy - (0.5d0 * ((kf ** 2) + (ki ** 2)))) * intg(f) * intg(g) &
        - intg(v * f) * intg(g) &
        - intg(f) * intg(v * g) &
        - two_e_integral(n_r, r_grid, r_weights, f, g)

    ! a = (energy - (0.5d0 * ((kf ** 2) + (ki ** 2)))) * intg(f) * intg(g)
    ! b = - intg(v * f) * intg(g)
    ! c = - intg(f) * intg(v * g)
    ! d = - two_e_integral(n_r, r_grid, r_weights, f, g)

    ! V_exchange = a + b + c + d

  contains

    ! in-lining integrate() for compactness
    double precision function intg (f)
      double precision , intent(in) :: f(n_r)
      intg = integrate(n_r, r_weights, f(:))
    end function intg

  end subroutine exchange_matrix_element

  ! two_e_integral
  !
  ! Numerically evaluates l=0 term of two-electron multi-pole integrals
  ! > \int_{0}^{\infty} f(r_1)
  ! >   (   (1/r_1) * \int_{0}^{r_1}                g(r_2) dr_2
  ! >     +           \int_{r_1}^{\infty} (1/r_2) * g(r_2) dr_2
  ! >   )
  function two_e_integral (n_r, r_grid, r_weights, f, g) result (integral)
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: r_weights(n_r)
    double precision , intent(in) :: f(n_r)
    double precision , intent(in) :: g(n_r)
    double precision :: integral
    ! local variables
    double precision :: a(n_r)
    double precision :: b(n_r)
    integer :: ii

    ! check if arguments are valid
    if (n_r < 1) then
      return
    end if

    ! populate arrays
    a(1) = r_weights(1) * g(1)
    do ii = 2, n_r
      a(ii) = a(ii-1) + (r_weights(ii) * g(ii))
    end do

    b(n_r) = (r_weights(n_r) * g(n_r)) / r_grid(n_r)
    do ii = n_r-1, 1, -1
      b(ii) = b(ii+1) + ((r_weights(ii) * g(ii)) / r_grid(ii))
    end do

    ! evaluate the integral
    integral = integrate(n_r, r_weights, &
        f(:) * ((a(:) / r_grid(:)) + b(:)))

  end function two_e_integral

  ! integrate
  !
  ! Numerically evaluates
  ! > \int_{0}^{\infty} f(r) dr = \sum_{i = 1}^{n_r} w_i + f_i
  ! where
  ! > f_i = f(r_i)
  function integrate (n_r, r_weights, f) result (integral)
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_weights(n_r)
    double precision , intent(in) :: f(n_r)
    double precision :: integral
    ! local variables
    integer :: ii

    ! check if arguments are valid
    if (n_r < 1) then
      return
    end if

    ! evaluate the integral
    integral = 0.0d0
    do ii = 1, n_r
      integral = integral + (r_weights(ii) * f(ii))
    end do

  end function integrate

end module vmatrix
