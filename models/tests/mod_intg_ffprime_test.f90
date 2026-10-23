!> James Carpenter, Oct 2026
!> Test numerical integration of FF' to F
!> We start here with F, construct FF', then use
!> JOREK routines to reconstruct F. Note in tests we only have psi_n = 0,..,1
!> so x-point handling is untested here

module mod_intg_ffprime_test
    use fruit

    !> Same phys_module vars as integrate_ffprime routine
    use phys_module, only: xpoint, F0, FF_0, FF_1, FF_coef, num_Fprofile, freeboundary_equil
    use phys_module, only: num_ffprime, num_ffprime_len, num_ffprime_x, num_ffprime_y0
    use phys_module, only: intg_ffprime_y0, intg_ffprime_y1, intg_ffprime_y2, intg_ffprime_y3
    use constants,   only: TWOPI
    implicit none

    private
    public :: run_fruit_intg_ffprime

    !> Vars to change for testing -----------------------------
    integer, parameter :: len_f   = 500
    real*8,  parameter :: f_0     = -0.8d0
    real*8,  parameter :: psi_ax  = 0.03d0
    real*8,  parameter :: psi_bnd = -1.42d0
    real*8,  parameter :: tol     = 0.001d0

    !> internal vars
    real*8  :: x(len_f), f_c(len_f), f_g(len_f), f_e(len_f), f_s(len_f)
    real*8  :: psi_arr(len_f)
    logical :: write_prof = .true.

contains
subroutine run_fruit_intg_ffprime()
    implicit none
    write(*, '(A)') " ... setting-up: ffprime integration tests"
    call setup
    write(*, '(A)') " ... running: ffprime integration tests"
    call run_test_case(test_intg_ffprime_cubic, "test_intg_ffprime_cubic")
    call run_test_case(test_intg_ffprime_guass, "test_intg_ffprime_guass")
    call run_test_case(test_intg_ffprime_exp,   "test_intg_ffprime_exp")
    call run_test_case(test_intg_ffprime_sin,   "test_intg_ffprime_sin")
    write(*,'(A)') " ... tearing-down: ffprime integration tests"
    call teardown
end subroutine run_fruit_intg_ffprime

subroutine test_intg_ffprime_cubic()
    implicit none
    real*8 :: f_calc(len_f)
    call test_intg_ffprime(f_c, f_calc, psi_arr)
    if (write_prof) call write_profile("intg_ffprime_cubic_test.dat", psi_arr, f_c, num_ffprime_y0, f_calc)
    call assert_equals(f_c, f_calc, len_f, tol, &
        "ERROR: ffprime integration failed on cubic profile")
end subroutine test_intg_ffprime_cubic

subroutine test_intg_ffprime_guass()
    implicit none
    real*8 :: f_calc(len_f)
    call test_intg_ffprime(f_g, f_calc, psi_arr)
    if (write_prof) call write_profile("intg_ffprime_gauss_test.dat", psi_arr, f_g, num_ffprime_y0, f_calc)
    call assert_equals(f_g, f_calc, len_f, tol, &
        "ERROR: ffprime integration failed on guass profile")
end subroutine test_intg_ffprime_guass

subroutine test_intg_ffprime_exp()
    implicit none
    real*8 :: f_calc(len_f)
    call test_intg_ffprime(f_e, f_calc, psi_arr)
    if (write_prof) call write_profile("intg_ffprime_exp_test.dat", psi_arr, f_e, num_ffprime_y0, f_calc)
    call assert_equals(f_e, f_calc, len_f, tol, &
        "ERROR: ffprime integration failed on exp profile")
end subroutine test_intg_ffprime_exp

subroutine test_intg_ffprime_sin()
    implicit none
    real*8 :: f_calc(len_f)
    call test_intg_ffprime(f_s, f_calc, psi_arr)
    if (write_prof) call write_profile("intg_ffprime_sin_test.dat", psi_arr, f_s, num_ffprime_y0, f_calc)
    call assert_equals(f_s, f_calc, len_f, tol, &
        "ERROR: ffprime integration failed on sin profile")
end subroutine test_intg_ffprime_sin

subroutine setup()
    implicit none

    real*8  :: ran(11)
    real*8  :: c(3), g(3), e(2), s(3)
    real*8  :: fc0, fg0, fe0, fs0
    integer :: i

    !> Makes 4 profiles : cubic / gaussian / exp / sinusoidal
    !>        N. coeffs :   3   /     3    /  2  /     3 
    !> We make 10 arbitrary coefficients, then generate profiles
    call random_number(ran)

    !> cubic coeffs, all in range [-10,10]
    c(1:3) = 20*ran(1:3) - 10

    !> gaussian coeffs
    g(1) = ran(4)*20 - 10
    g(3) = 1.d0 + ran(5)*9
    g(2) = 0.5*g(3)*(2*ran(6)-1)   ! choose position of guassian st non-zero on 0,1 region

    !> exp coeffs
    e(1) = 20.d0*ran(7)
    e(2) = 1.d0 + 19.d0*ran(8)

    !> sinusoidal coeffs
    s(1:2) = 20*ran(9:10)
    s(3)   = TWOPI*ran(11)

    !> generate functions
    do i = 1, len_f
        x(i)   = dble(i-1)/dble(len_f - 1)
        f_c(i) = c(1)*x(i)**3 + c(2)*x(i)**2 + c(3)*x(i)
        f_g(i) = g(1)*exp(-(x(i)-g(2))**2/(2*g(3)**2))
        f_e(i) = e(1)*exp(x(i)/e(2))
        f_s(i) = s(1)*sin(s(2)*x(i) + s(3))
    enddo

    !> find values on bnd
    fc0 = f_c(len_f); fg0 = f_g(len_f); fe0 = f_e(len_f); fs0 = f_s(len_f)

    !> scale profiles to be F0 on boundary
    f_c = f_c - fc0 + f_0
    f_g = f_g - fg0 + f_0
    f_e = f_e - fe0 + f_0
    f_s = f_s - fs0 + f_0
end subroutine setup

subroutine teardown()
    implicit none
    x = 0.d0; f_c = 0.d0; f_g = 0.d0; f_e = 0.d0; f_s = 0.d0
end subroutine teardown

subroutine safe_alloc(arr, n)
    implicit none
    real*8,  allocatable, intent(inout) :: arr(:)
    integer,              intent(in)    :: n

    if (allocated(arr)) then
        if (size(arr) == n) return
        deallocate(arr)
    endif
    allocate(arr(n))
endsubroutine safe_alloc

subroutine test_intg_ffprime(f, f_calc, psi_arr)

    use mod_F_profile
    use profiles, only: derivProf
    implicit none

    real*8,  intent(in)  :: f(:)
    real*8,  intent(out) :: f_calc(:)
    real*8,  intent(out) :: psi_arr(:)

    real*8, allocatable :: f_prime(:)
    real*8              :: dp, F_prof
    integer             :: n, i

    !> allocatable copies
    real*8, allocatable :: psi_loc(:), f_loc(:)

    !> fluff vars, but needed for F_profile routine
    real*8 :: dF_dpsi, dF_dz, dF_dpsi2, dF_dx2, dF_dpsi_dz, FFprime_prof
    real*8 :: dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz

    !> get length of input f prof
    n = size(f)

    !> Set all phys_module vars
    xpoint      = .false.; F0 = f(n); FF_0 = 0.d0; FF_1 = 0.d0; FF_coef = 0.d0
    num_ffprime = .true. ; num_ffprime_len = n
    num_Fprofile = .false.; freeboundary_equil = .false.

    !> allocate arrs
    allocate(psi_loc(n), f_loc(n), f_prime(n))
    call safe_alloc(num_ffprime_x,   n)
    call safe_alloc(num_ffprime_y0,  n)
    call safe_alloc(intg_ffprime_y0, n)
    call safe_alloc(intg_ffprime_y1, n)
    call safe_alloc(intg_ffprime_y2, n)
    call safe_alloc(intg_ffprime_y3, n)
    intg_ffprime_y0 = 0.d0; intg_ffprime_y1 = 0.d0; intg_ffprime_y2 = 0.d0; intg_ffprime_y3 = 0.d0
    psi_arr         = 0.d0; f_prime         = 0.d0; f_calc          = 0.d0
    num_ffprime_x   = 0.d0; num_ffprime_y0  = 0.d0

    !> dpsi / dpsi_n = dp (same convention as F_profile: psi_bnd - psi_ax)
    dp = psi_bnd - psi_ax

    !> set up psi and psi_n arrays
    do i = 1, n
        num_ffprime_x(i) = dble(i-1)/dble(n - 1)
        psi_arr(i)       = dp*num_ffprime_x(i) + psi_ax
    enddo

    !> copy psi and f
    psi_loc = psi_arr; f_loc = f

    !> take deriv of F and find FF'
    call derivProf(psi_loc, f_loc, n, f_prime)
    do i = 1, n
        num_ffprime_y0(i) = f(i)*f_prime(i)
    enddo

    ! write(*,*) "DEBUG : num_ffprime_x  : ", num_ffprime_x
    ! write(*,*) "DEBUG : num_ffprime_y0 : ", num_ffprime_y0

    !> integrate ffprime prof in psi_n
    call integrate_ffprime(.false.)

    ! write(*,*) "DEBUG : intg_ffprime_y0 : ", intg_ffprime_y0

    !> calc derivs of integrated ffprime prof
    call derivProf(num_ffprime_x, intg_ffprime_y0, n, intg_ffprime_y1)
    call derivProf(num_ffprime_x, intg_ffprime_y1, n, intg_ffprime_y2)
    call derivProf(num_ffprime_x, intg_ffprime_y2, n, intg_ffprime_y3)

    !> Find F at all values of psi
    do i = 1, n
        call F_profile(.false., 1, 0.d0, [0.d0,0.d0], psi_arr(i), psi_ax, psi_bnd,    &
                       F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dx2, dF_dpsi_dz, &
                       FFprime_prof, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz)
        f_calc(i) = F_prof
    enddo
    deallocate(psi_loc, f_loc, f_prime)
end subroutine test_intg_ffprime

subroutine write_profile(fname, psi, f, ffp, f_calc)
  character(len=*), intent(in) :: fname
  real*8, intent(in) :: psi(:), f(:), ffp(:), f_calc(:)
  integer :: iu, i
  open(newunit=iu, file=fname, status="replace", action="write")
  write(iu,"(A)") "#  psi  f  ffp  f_calc  err"
  do i = 1, size(f)
    write(iu,"(5ES16.7)") psi(i), f(i), ffp(i), f_calc(i), abs(f_calc(i) - f(i))
  enddo
  close(iu)
end subroutine

end module mod_intg_ffprime_test
    