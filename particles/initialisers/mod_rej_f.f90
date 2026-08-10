!> Module defining abstract interfaces and implemtation for 
!> real space rejection functions used during particle initialisation.
!>
!> A rejection funciton (rej_f) maps field values at a sample point to
!> an acceptance probability in [0,1]. The spatial_pdf derived type bundles 
!> a procedure pointer with the variable list it expects, so both pieces of
!> information always travel together
!>
!> ## Variable index convention (vars(:) entries)
!>
!>  Positive index -> JOREK variable number (model-dependent)
!>  -1             -> R   (major radius)
!>  -2             -> Z   (height)
!>  -3             -> Phi (toroidal angle)
!>
!> The P(:) array passed to a rej_f function is ordered to match vars(:),
!> so a function must be paired with the vars array it was written for. 
!> Use spatial_pdf_from_name() to get a correctly-paired spatial_pdf.

module mod_rej_f
    use equil_info       ! (R_axis,Z_axis), etc..
    use mod_model_settings

    implicit none
    private

    public :: spatial_pdf
    public :: rej_f
    public :: itpa_tae_pdf, RZ_pdf, analytical_pdf, current_pdf
    public :: spatial_pdf_from_name

    ! =============================================================================
    ! Interface
    ! =============================================================================
    
    !> Primary interface for spatial rejection functions used in particle init
    !>
    !> @param n     Number of field values in P and gradP
    !> @param P     Field values at the sample point, ordered as vars(:)
    !> @param gradP Gradients (3,n); may be unused
    !> @return      Acceptance probability in [0,1]
    abstract interface
        function rej_f(n, P, gradP)
            implicit none
            integer,                intent(in) :: n
            real*8, dimension(n),   intent(in) :: P
            real*8, dimension(3,n), intent(in) :: gradP
            real*4                             :: rej_f
        end function rej_f
    end interface

    ! =============================================================================
    ! spatial_pdf -> Type bundling rejection f and vars it expects
    ! =============================================================================

    !> Bundle type pairing a spatial rejection function with
    !> the JOREK / (R,Z,Phi) variable list it expects
    !>
    !> Always construct via spatial_pdf_from_name() or by setting
    !> both components together - the procedure f is written to 
    !> expect exactly size(vars) values in P(:), ordered to match vars(:).
    !>
    !> Example (custom profile)
    !>
    !>      type(spatial_pdf) :: pdf
    !>      pdf%f    => my_rej_function
    !>      pdf%vars = [var_psi, -1] ! psi and R
    type :: spatial_pdf
        procedure(rej_f), nopass, pointer :: f => null()
        integer,          allocatable     :: vars(:)
    end type spatial_pdf
    

    contains
    ! =============================================================================
    ! Constructor
    ! =============================================================================

    !> Returns a correctly-paired spatial_pdf type
    !>
    !> Recongnised names:
    !>    'itpa_tae'        - reproduce EP dist in ITPA TAE benchmark
    !>    'RZ'              - weight by 1/(1+r_minor)^2
    !>    'analytical'      - weight by (1-(r/a)^2)^nu
    !>    'current'         - weight by normalised j_tor
    !>    'none'            - no rejection (f => null, vars unallocated)
    function spatial_pdf_from_name(name) result(pdf)
        character(len=*), intent(in) :: name
        type(spatial_pdf)            :: pdf

        select case(trim(name))

        case('itpa_tae')
            pdf%f    => itpa_tae_pdf
            pdf%vars =  [var_psi]
        
        case('RZ')
            pdf%f    => RZ_pdf
            pdf%vars =  [-2, -1]

        case('analytical')
            pdf%f    => analytical_pdf
            pdf%vars =  [-2, -1]

        case('current')
            !> NOTE : var_zj = 0 in fullMHD (j_tor is not a stored variable)
            !> This rej_f is only valid for models where var_zj > 0 (e.g. model600)
            if (var_zj == 0) then
                write(*,*) "ERROR (mod_rej_f): 'current' spatial PDF requires model"
                write(*,*) "  with var_zj > 0. j_tor is not a stored variable in this model"
                stop 1
            endif
            pdf%f    => current_pdf
            pdf%vars =  [-1, var_zj]

        case('none')
            !> Leave f => null() and vars unallocated
            !> caller should handle this case by skipping use of rej_f
            !> no rej_f should be passed to initialiser, sampling will be uniform in space

        case default
            write(*,*) "ERROR (mod_rej_f): Unknown spatial PDF name '", trim(name), "'"
            write(*,*) "      Valid names: 'itpa_tae', 'RZ', 'analytical', "
            write(*,*) "                   'current', 'none'"
            stop 1
        end select
    end function spatial_pdf_from_name


    ! =============================================================================
    ! Rejection functions
    ! =============================================================================


    !> rejection function to produce EP spatial, distribution from ITPA TAE benchmark
    !> A. Könies et al 2018 Nucl. Fusion 58 126027, https://doi.org/10.1088/1741-4326/aae4e6
    !>
    !> Expected vars : [var_psi]
    !>                 P(1) = psi
    pure function itpa_tae_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f
        
        !> internal vars
        real*8              :: s, psi_norm, coeff(0:3)

        coeff(0)=0.49123
        coeff(1)=0.298228
        coeff(2)=0.198739
        coeff(3)=0.521298

        psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)
        s        = 0.957 * psi_norm + 0.043 * psi_norm**2

        f = real(coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2)))), 4)
    end function itpa_tae_pdf

    !> Weight as 1/(1 + r_minor)^2 - broad profile peaked on axis
    !>
    !> Expected vars : [-2, -1]
    !>                 P(1) = Z
    !>                 P(2) = R
    pure function RZ_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f

        real*8 :: minor_r

        minor_r = sqrt((P(1) - ES%Z_axis)**2 + (P(2) - ES%R_axis)**2)
        f = real(1.d0 / (1.d0 + minor_r)**2, 4)
    end function RZ_pdf
    
    !> Analytically prescribed profile (1 - (r_minor/a)^2)^nu, nu=2
    !>
    !> Expected vars : [-2, -1]
    !>                 P(1) = Z
    !>                 P(2) = R
    pure function analytical_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f

        real*8, parameter :: nu = 2.d0
        real*8            :: minor_r

        minor_r = sqrt((P(1) - ES%Z_axis)**2 + (P(2) - ES%R_axis)**2)
        f = real((1.d0 - (minor_r/ES%LCFS_a)**2)**nu, 4)
        f = max(f, 0.0e0)  !> clamp outside LCFS
    end function analytical_pdf

    !> Weight proportional to normalised toroidal current density j_tor
    !>
    !> Expected vars : [-1, var_zj]
    !>                 P(1) = R
    !>                 P(2) = j_tor
    !>
    !> NOTE : Won't work for fMHD models, see note in constructor above
    !> NOTE : jzmin and jzmax are currently hardcoded, needs better handling
    pure function current_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f

        !> TODO: make jzmin.jzmax namelist params
        real*8, parameter :: jzmax = 3.0d0 / 10.0d0
        real*8, parameter :: jzmin = 1.239d-4 / 11.0d0

        f = real((P(2)/P(1) - jzmin) / (jzmax - jzmin), 4)
        f = max(f, 0.0e0)
    end function current_pdf
end module mod_rej_f
