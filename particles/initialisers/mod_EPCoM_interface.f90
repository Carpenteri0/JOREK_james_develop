!>
!> Adapted from Guillaume Broachard's interface in XTOR-K
!>
!> Reads netcdf file named 'EPCoM_CoM_pdf.nc' constructs spline and
!> samples particles based on f(P_phi, muB_0/E, E, sigma) distribution
!> contained within netcdf file
!>
!> James Carpenter
!> zhpw96@durham.ac.uk
!> Aug, 2026

module mod_EPCoM_interface
    !$ use omp_lib
    use netcdf
    use mpi
    use mod_particle_sim
    use mod_particle_types
    use initialisers_base
    use mod_pusher_tools
    use mod_rng
    use mod_pcg32_rng
    use mod_random_seed
    use data_structure
    use mod_interp
    use constants, only: TWOPI, EL_CHG, ATOMIC_MASS_UNIT
    use equil_info

    implicit none
    
    !> EPCoM params
    character(len=*), parameter                   :: EPCoM_file_CoM = 'EPCoM_CoM_pdf.nc'
    integer                                       :: N_P, N_L, N_E
    real*8, dimension(:),             allocatable :: P_ar, L_ar, E_ar
    real*8, dimension(:,:,:,:,:),     allocatable :: F_CoM
    real*8, dimension(:,:,:,:),       allocatable :: Jac_CoM
    real*8, dimension(:),             allocatable :: max_F_CoM
    real*8                                        :: mEP, qEP
    real*8                                        :: mass, atomic_no
    real*8                                        :: scaling_pdf

    contains

    !> Master subroutine for EPCoM initialisation
    !> Reads file, generates real/phase-space coordinates
    !> from uniform sampling. Compares against actual pdf (with spline)
    !> accepts if probability is sufficient. Recreates given
    !> distribution in CoM space. 
    subroutine initialise_from_EPCoM(sim, group_num)
        implicit none

        !> i/o vars
        class(particle_sim), intent(inout) :: sim
        integer,             intent(in)    :: group_num
        
        !> internal vars
        !class(pcg32_rng)                          :: rng_base
        class(particle_kinetic_leapfrog), pointer :: particles(:)
        !type(pcg32_rng), allocatable              :: rng
        class(type_rng), dimension(:), allocatable :: rng
        integer :: n_particles
        real*8  :: RBox(2), ZBox(2), EBox(2), LBox(2)
        integer :: i_part, i_elm
        logical :: accepted, is_outside
        real*8  :: ran(7)
        real*8  :: R, Z, phi, Energy, Pitch, P_phi, gyro_angle, Lambda
        real*8  :: v(3), v_squared, v_par, v_perp, v_phi, sig
        real*8  :: s, t, psi, dummy_R, dummy_Z, E(3), B(3), U
        real*8  :: B_norm(3), e1(3), e2(3)
        real*8  :: E_axis(3), B_axis(3), psi_axis, U_axis, B_axis_mag
        real*8  :: Jac_Ep, Jac_RZ
        real*8  :: F_loc, H_loc, randvarf
        real*8  :: max_H
        integer :: ifail, ierr
        integer :: tid
        real*8, allocatable :: stored_x(:,:), stored_v(:,:), stored_st(:,:)
        integer, allocatable:: stored_ielm(:)

        if (sim%my_id .eq. 0) then
            write(*,*) ""
            write(*,*) "=============================="
            write(*,*) "       EPCoM interface        "
            write(*,*) "=============================="
            write(*,*) ""
        endif

        mass      = sim%groups(group_num)%mass * ATOMIC_MASS_UNIT
        atomic_no = sim%groups(group_num)%Z

        !> Calc B_norm on axis for normalised magnetic moment later
        call sim%fields%calc_EBpsiU(sim%time, ES%i_elm_axis, [ES%s_axis, ES%t_axis], 0.d0, E_axis, B_axis, psi_axis, U_axis)
        B_axis_mag = norm2(B_axis)

        !allocate(pcg32_rng :: rng)
        !call rng%initialize(7, random_seed(), 1, 1, ifail)
        call setup_shared_rngs(7, random_seed(), pcg32_rng(), rng)

        if (sim%my_id .eq. 0) then
            write(*,*) "mass   : ", mass
            write(*,*) "Z      : ", atomic_no
        endif
        
        n_particles = size(sim%groups(group_num)%particles(:))
        allocate(stored_x(3, n_particles))
        allocate(stored_v(3, n_particles))
        allocate(stored_st(2, n_particles))
        allocate(stored_ielm(n_particles))

        !> Read EPCoM file and load arrays
        call read_EPCoM_file(sim%my_id)

        !> Boxes for rejection sampling
        call domain_bounding_box(sim%fields%node_list, sim%fields%element_list, RBox(1), RBox(2), ZBox(1), ZBox(2))
        EBox(1) = minval(E_ar); Ebox(2) = maxval(E_ar)


        !> calculate maximum of F*Jacobian to scale pdf in rej sampling
        !> Jac_Ep = 2pi*sqrt(2E/m) - max at EBox(2)
        Jac_RZ = TWOPI
        max_H  = maxval(max_F_CoM) * (TWOPI*sqrt(2.0d0*EBox(2)*EL_CHG / mass)) * Jac_RZ
        if (sim%my_id .eq. 0) then
            write(*,*) "RBox = ", RBox
            write(*,*) "ZBox = ", ZBox
            write(*,*) "EBox = ", EBox
            write(*,*) "LBox = ", minval(L_ar), maxval(L_ar)
            write(*,*) "PBox = ", minval(P_ar), maxval(P_ar)
            write(*,*) "max_H= ", max_H
        endif

       !$omp parallel do default(none) &
       !$omp shared(sim, mass, atomic_no, E_ar, L_ar, P_ar, N_E, N_L, N_P, psi_axis,&
       !$omp        F_CoM, Jac_RZ, max_H, RBox, ZBox, EBox, group_num, n_particles, &
       !$omp        stored_x, stored_v, stored_st, stored_ielm, rng, B_axis_mag)    &
       !$omp private(i_part, R, Z, phi, Energy, Pitch, P_phi, gyro_angle, Lambda,   &
       !$omp         ran, s, t, psi, U, E, B, B_norm, e1, e2, v, v_squared, v_par,  &
       !$omp         v_perp, v_phi, sig, Jac_Ep, F_loc, H_loc, randvarf, ifail,     &
       !$omp         ierr, i_elm, dummy_R, dummy_Z, accepted, is_outside, tid)

        !> omp loop over total number of marker particles to find
        do i_part = 1, n_particles

            !$ tid = omp_get_thread_num() + 1 !> thread-id for rng

            !> loop untill found acceptable particle (rej sample loop)
            accepted = .false.
            do while (.not. accepted)

        ! ---- 1. Draw (unif) random number and convert to phase space point ----
                call rng(tid)%next(ran)  !> 7 uniform values in [0,1]
                R           = RBox(1) + (RBox(2) - RBox(1)) * ran(1) 
                Z           = ZBox(1) + (ZBox(2) - ZBox(1)) * ran(2)
                phi         = TWOPI * ran(3) 
                Energy      = EBox(1) + (EBox(2) - EBox(1)) * ran(4)
                Pitch       = 2.0d0 * (ran(5) - 0.5d0) 
                gyro_angle  = TWOPI * ran(6)
                randvarf    = max_H * ran(7)

        ! ---- 2. Check (R,Z) is inside domain & find local coords ----
                call find_RZ(sim%fields%node_list, sim%fields%element_list, &
                             R, Z, dummy_R, dummy_Z, i_elm, s, t, ifail)
                if (i_elm .le. 0) cycle  !> start inner loop again if not in domain

        ! ---- 3. Find Psi and B-field at this (i_elm, s, t, phi) ----
                call sim%fields%calc_EBpsiU(sim%time, i_elm, [s, t], phi, E, B, psi, U)

        ! ---- 4. Compute velocity components, P_phi, and normalised magnetic moment ----
                B_norm    = B / sqrt(dot_product(B, B))
                v_squared = 2.0d0 * (Energy*EL_CHG) / mass
                v_par     = Pitch * sqrt(v_squared)
                v_perp    = sqrt(max(0.d0,v_squared - v_par**2))
                sig       = sign(1.0d0, Pitch)

                call get_orthonormals(B_norm, e1, e2)
                v      = v_par*B_norm + v_perp*(cos(gyro_angle)*e1 + sin(gyro_angle)*e2)
                v_phi  = v(3)
                P_phi  = mass*R*v_phi/(atomic_no*EL_CHG) + (psi-psi_axis) !> Normalised by q=Z*e
                Lambda = (v_perp**2 * B_axis_mag) / (v_squared*norm2(B))

        ! ---- 5. Evalutate pdf(P_phi, L, E, sgn(v_||)) with EPCoM spline ----
                call EPCoM_spline(P_phi, Lambda, Energy, sig, F_loc)

        ! ---- 6. Jacobians (phase-space volume element) ----
                Jac_Ep = TWOPI * sqrt(v_squared)
                H_loc  = F_loc * Jac_Ep * Jac_RZ

        ! ---- 7. Rejection test ----
                !> Unphysical spline overshoot
                if (H_loc .lt. 0.d0) cycle
                !> reject if prob of being at this rand position is too high
                if (randvarf .le. H_loc) accepted = .true.

            enddo !> rejection loop

        ! ---- 8. Store accepted particle
            stored_x(:,  i_part) = [R, Z, phi]
            stored_v(:,  i_part) = v
            stored_st(:, i_part) = [s, t]
            stored_ielm( i_part) = i_elm
        
        enddo !> particle loop
        !$omp end parallel do

        ! ---- 9. Pass stored values to particle object
        select type(p => sim%groups(group_num)%particles)
        type is (particle_kinetic_leapfrog)
            do i_part = 1, n_particles
                p(i_part)%x      = stored_x(:,  i_part)
                p(i_part)%v      = stored_v(:,  i_part)
                p(i_part)%st     = stored_st(:, i_part)
                p(i_part)%i_elm  = stored_ielm( i_part)
                p(i_part)%weight = 1.0
            enddo
        class default
            if (sim%my_id .eq. 0) then
                write(*,*) "ERROR : initialise_from_EPCoM only works for particles"
                write(*,*) "        of type particle_kinetic_leapfrog"
                write(*,*) "        check type of group_num : ", group_num
            endif
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end select

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        deallocate(stored_x, stored_v, stored_st, stored_ielm, rng)

    end subroutine initialise_from_EPCoM

    !> Read file on root proc and communicate to all other procs
    subroutine read_EPCoM_file(my_id)
        implicit none

        !> i/o vars
        integer, intent(in) :: my_id
        
        !> netcdf ids
        integer                               :: ncid, varid, icount
        integer                               :: state, status
        integer, dimension(nf90_max_var_dims) :: dimids

        !> mpi error code
        integer :: ierr

        !> loop indecies
        integer :: i, j, k, l

        !> read file on root proc
        if (my_id == 0) then
            write(*,*) "Reading EPCoM file : ", trim(EPCoM_file_CoM)

            !> open netcdf file
            status = nf90_open(EPCoM_file_CoM, nf90_nowrite, ncid)
            if (status /= NF90_NOERR) then
                write(*,*) "ERROR: Could not open EPCoM file: ", trim(EPCoM_file_CoM)
                write(*,*) "Exiting...."
                stop
            end if

            write(*,*) "Successfully openened EPCoM file : ", trim(EPCoM_file_CoM)

            !> get CoM array dimensions
            status = nf90_inq_varid(ncid, 'CoM_pdf', varid)
            status = nf90_inquire_variable(ncid, varid, dimids=dimids)
            status = nf90_inquire_dimension(ncid, dimids(1), len = N_P)
            status = nf90_inquire_dimension(ncid, dimids(2), len = N_L)
            status = nf90_inquire_dimension(ncid, dimids(3), len = N_E)

            write(*,*) "Pzeta dimension  N_P = ", N_P
            write(*,*) "Lambda dimension N_L = ", N_L
            write(*,*) "Energy dimension N_E = ", N_E
        endif

        !> Broadcast CoM dimensions
        icount = 1
        call MPI_BCAST(N_P, icount, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(N_L, icount, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(N_E, icount, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        !> Allocate CoM arrays
        allocate(P_ar(N_P), L_ar(N_L), E_ar(N_E), STAT=state)
        allocate(F_CoM(N_P,N_L,N_E,2,27), Jac_CoM(N_P,N_L,N_E,2), STAT=state)
        allocate(max_F_CoM(N_E), STAT=state)

        max_F_CoM = -1.0
        scaling_pdf = 1.0d0

        if (my_id == 0) then

            !> Load CoM pdf (units : eV^-1.s.m^-4)
            status = nf90_get_var(ncid, varid, F_CoM)
            call check_status(status, 'F_CoM')
            F_CoM  = scaling_pdf*F_CoM                !> Apply rescaling factor

            !> Load the Pphi array (units : T.m^2)
            status = nf90_inq_varid(ncid, 'P_ar', varid)
            status = nf90_get_var(ncid, varid, P_ar)
            call check_status(status, 'P_ar')

            !> Load the lambda = mu Baxis/E array (units : 1)
            status = nf90_inq_varid(ncid, 'L_ar', varid)
            status = nf90_get_var(ncid, varid, L_ar)
            call check_status(status, 'L_ar')

            !> Load the energy array (units : eV)
            status = nf90_inq_varid(ncid, 'E_ar', varid)
            status = nf90_get_var(ncid, varid, E_ar)
            call check_status(status, 'E_ar')

            !> Load the CoM jacobian matrix (units : m^2.s^-1.T^-1)
            status = nf90_inq_varid(ncid,'Jac_CoM', varid)
            status = nf90_get_var(ncid,varid,Jac_CoM)
            call check_status(status, 'Jac_CoM')

            ! Load the EP charge number (units : 1)
            status = nf90_inq_varid(ncid,'Z', varid)
            status = nf90_get_var(ncid,varid,qEP)
            call check_status(status, 'qEP')

            ! Load the EP charge number (units : 1)
            status = nf90_inq_varid(ncid,'mA', varid)
            status = nf90_get_var(ncid,varid,mEP)
            call check_status(status, 'mEP')

            !> TODO: Check mass/charge specified in nml is the same as EPCoM

            !> Close netCDF file
            status = nf90_close(ncid)
            if (status /= NF90_NOERR) then
                write(*,*) "ERROR : Could not close netcdf file", nf90_strerror(status)
            endif

            write(*,*) "Successfuly loaded EPCoM file"

            do k = 1, N_E
                do j = 1, N_L
                    do i = 1, N_P
                        do l = 1, 2
                            if (Jac_CoM(i,j,k,l)>1e-3 .and. F_CoM(i,j,k,l,1)>max_F_CoM(k)) then
                                max_F_CoM(k) = F_CoM(i,j,k,l,1)
                            endif
                        enddo
                    enddo
                enddo
            enddo

        endif
        deallocate(Jac_CoM)

        !> Broadcast CoM arrays
        icount=N_P*N_L*N_E*2*27
        call MPI_BCAST(F_CoM,     icount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(max_F_CoM, N_E,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(P_ar,      N_P,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(L_ar,      N_L,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(E_ar,      N_E,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end subroutine read_EPCoM_file

    !> Evaluate if variable was read from netcdf correctly
    subroutine check_status(status, var)
        implicit none
        
        integer,          intent(in) :: status
        character(len=*), intent(in) :: var

        if (status /= NF90_NOERR) then
            write(*,*) "ERROR: Could not read variable: ", trim(var)
            write(*,*) "Exiting...."
            stop
        else
            write(*,*) "Successfully read variable : ", trim(var)
        end if

    end subroutine check_status

    !> Evaluate spline at given CoM position
    subroutine EPCoM_spline(P_phi, Lambda, Energy, sig, F_loc)
        implicit none

        !> i/o vars
        real*8, intent(in)  :: P_phi, Lambda, Energy, sig
        real*8, intent(out) :: F_loc

        !> internal vars
        integer                :: ind_P, ind_L, ind_E, l, i
        real*8                 :: dP, dL, dE, dx, dy, dz
        real*8, dimension(27)  :: dF_CoM 

        !> cell widths
        dx = P_ar(2) - P_ar(1)
        dy = L_ar(2) - L_ar(1)
        dz = E_ar(2) - E_ar(1)

        !> index of nearest cell corner
        ind_P = max(1, min(N_P-1, ceiling((P_phi  - P_ar(1))/dx)))
        ind_L = max(1, min(N_L-1, ceiling((Lambda - L_ar(1))/dy)))
        ind_E = max(1, min(N_E-1, ceiling((Energy - E_ar(1))/dz)))

        !> distance to nearest cell corner
        dP = P_phi  - P_ar(ind_P)
        dL = lambda - L_ar(ind_L)
        dE = Energy - E_ar(ind_E)

        !> Spline coeffs
        dF_CoM(1)     = 1.d0
        dF_CoM(2)     = dP
        dF_CoM(3)     = dP**2
        dF_CoM(4:6)   = dF_CoM(1:3)*dL
        dF_CoM(7:9)   = dF_CoM(1:3)*dL**2
        dF_CoM(10:18) = dF_CoM(1:9)*dE
        dF_CoM(19:27) = dF_CoM(1:9)*dE**2

        !> Choose sigma branch
        if (sig .gt. 0) then
            l = 1
        else
            l = 2
        endif

        !> Evaluate spline at provided position
        F_loc = 0.0d0
        do i = 1, 27
            F_loc = F_loc + F_CoM(ind_P, ind_L, ind_E, l, i)*dF_CoM(i)
        enddo
    end subroutine EPCoM_spline

end module mod_EPCoM_interface