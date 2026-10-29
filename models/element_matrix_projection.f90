subroutine element_matrix_projection(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
! itype = 0: projection of a profile onto variable ivar_out
! itype = 5: weak-form projection of the current density component ivar_in (1=J_R, 2=J_Z, 3=J_phi) of
!            harmonic i_harm (full MHD only), i.e. int v J R dA with the curl moved onto the test function v,
!            so that only first derivatives of A are needed (used by jorek2vtk)
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use mod_basisfunctions, only: basisfunctions
use equil_info, only: ES
use phys_module, only: xpoint, xcase, F0, n_jropes, R_jropes, Z_jropes, w_jropes, rho_jropes, T_jropes, normalized_velocity_profile, forceSDN, mode
use mod_F_profile
use mod_interp
use nodes_elements
use mod_sources
use constants

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: factor(n_gauss,n_gauss)
real*8     :: eq_g(n_gauss,n_gauss),   eq_s(n_gauss,n_gauss),   eq_t(n_gauss,n_gauss)
real*8     :: eq2_g(n_gauss,n_gauss),  eq2_s(n_gauss,n_gauss),  eq2_t(n_gauss,n_gauss)
real*8     :: ELM(n_vertex_max*n_degrees,n_vertex_max*n_degrees), RHS(n_vertex_max*n_degrees)

real*8     :: xjac, wst
real*8     :: v, psi, rhs_ij
integer    :: ms, mt, i, j, k, l, index_ij, index_kl, itype, ivar_in, ivar_out, i_harm
real*8     :: F_prof        ,dF_dpsi      ,dF_dz      , dF_dpsi2      ,dF_dz2       ,dF_dpsi_dz
real*8     :: zFFprime      ,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi2,dFFprime_dz2 ,dFFprime_dpsi_dz
real*8     :: zn,  dn_dpsi,  dn_dpsi2,  dn_dz,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi2_dz,  dn_dpsi_dz2
real*8     :: zT,  dT_dpsi,  dT_dpsi2,  dT_dz,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi2_dz,  dT_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: zrn, drn_dpsi, drn_dpsi2, drn_dz, drn_dz2, drn_dpsi_dz, drn_dpsi3, drn_dpsi2_dz, drn_dpsi_dz2
real*8     :: zV,  dV_dpsi,  dV_dpsi2,  dV_dz,  dV_dz2,  dV_dpsi_dz,  dV_dpsi3,  dV_dpsi2_dz,  dV_dpsi_dz2
real*8     :: var_RHS
real*8     :: R_out,Z_out,s_out,t_out
integer    :: ielm_out,ifail
real*8     :: dd1,dd2,dd3,dd4,dd5
integer    :: nj
real*8     :: rr,ww
real*8     :: R,Z
real*8     :: G(n_vertex_max,n_degrees), G_s(n_vertex_max,n_degrees), G_t(n_vertex_max,n_degrees)
real*8     :: Bfield(3), Bfield_p(2), R_s, R_t, Z_s, Z_t, v_R, v_Z, nrm_R, nrm_Z, s_side, t_side, dphi_fact
integer    :: i_partner, iside

ELM=0.d0
RHS=0.d0

! --- itype=5: the toroidal derivative of harmonic i_harm is given by its cos/sin partner harmonic
i_partner = i_harm
dphi_fact = 0.d0
if (i_harm .gt. 1) then
  if (mod(i_harm,2) .eq. 0) then
    i_partner = i_harm + 1 ; dphi_fact = + float(mode(i_harm))
  else
    i_partner = i_harm - 1 ; dphi_fact = - float(mode(i_harm))
  endif
endif

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g(:,:)   = 0.d0; x_s(:,:)   = 0.d0; x_t(:,:)   = 0.d0;
y_g(:,:)   = 0.d0; y_s(:,:)   = 0.d0; y_t(:,:)   = 0.d0;
eq_g(:,:)  = 0.d0; eq_s(:,:)  = 0.d0; eq_t(:,:)  = 0.d0;
eq2_g(:,:) = 0.d0; eq2_s(:,:) = 0.d0; eq2_t(:,:) = 0.d0;

do i=1,n_vertex_max
 do j=1,n_degrees
   do ms=1, n_gauss
     do mt=1, n_gauss

       x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
       y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

       x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
       x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
       y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
       y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

       eq_g(ms,mt)  = eq_g(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * element%size(i,j) * H(i,j,ms,mt)
       eq_s(ms,mt)  = eq_s(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * element%size(i,j) * H_s(i,j,ms,mt)
       eq_t(ms,mt)  = eq_t(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * element%size(i,j) * H_t(i,j,ms,mt)

       if (ivar_out .eq. 710) then
#ifdef fullmhd
         eq2_g(ms,mt)  = eq2_g(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H(i,j,ms,mt)
         eq2_s(ms,mt)  = eq2_s(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_s(i,j,ms,mt)
         eq2_t(ms,mt)  = eq2_t(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_t(i,j,ms,mt)
#endif
       else
         eq2_g(ms,mt)  = eq2_g(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H(i,j,ms,mt)
         eq2_s(ms,mt)  = eq2_s(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H_s(i,j,ms,mt)
         eq2_t(ms,mt)  = eq2_t(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H_t(i,j,ms,mt)
       endif

     enddo
   enddo
 enddo
enddo


factor =  x_g                                   ! Poisson

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

  do mt=1, n_gauss
 
    wst = wgauss(ms)*wgauss(mt)

    R = x_g(ms,mt) 
    Z = y_g(ms,mt) 
    xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

    if (itype .eq. 5) call field_at_point(H(:,:,ms,mt), H_s(:,:,ms,mt), H_t(:,:,ms,mt), R, R_s, R_t, Z_s, Z_t, Bfield, Bfield_p)
    
    if (ivar_out .eq. 710) then
      ! --- note: no need to use psi_axis_init, psi_bnd_init etc., since this routine should only be called at t=0
      call F_profile      (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           F_prof        ,dF_dpsi      ,dF_dz      , dF_dpsi2      ,dF_dz2       ,dF_dpsi_dz , &
                           zFFprime      ,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi2,dFFprime_dz2 ,dFFprime_dpsi_dz)
      var_RHS = F_prof
    endif
    if ( (ivar_out .eq. var_rho) .or. (ivar_out .eq. var_zj) ) then
      call density        (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
      if (ivar_out .eq. var_rho) then
        var_RHS = zn
        ! --- Define density blobs for model002 and model003
        if ((jorek_model .eq. 002) .or. (jorek_model .eq. 003) .or. (jorek_model .eq. 004)) then
          do nj=1,n_jropes
            rr = sqrt((R-R_jropes(nj))**2 + (Z-Z_jropes(nj))**2)
            ww = w_jropes(nj)
            if (rr .lt. 1.00*ww) then
              var_RHS = var_RHS + rho_jropes(nj) * (1.0 - (rr/ww)**2 )**2
            endif
          enddo
        endif
      endif
    endif
    if ( (ivar_out .eq. var_T) .or. ( (ivar_out .eq. var_zj) .and. (.not. with_TiTe)) ) then
      call temperature    (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
      if (ivar_out .eq. var_T) then
        var_RHS = zT
        ! --- Define temperature blobs for model003
        if (jorek_model .eq. 003) then
          do nj=1,n_jropes
            rr = sqrt((R-R_jropes(nj))**2 + (Z-Z_jropes(nj))**2)
            ww = w_jropes(nj)
            if (rr .lt. 1.00*ww) then
              var_RHS = var_RHS + T_jropes(nj) * (1.0 - (rr/ww)**2 )**2
            endif
          enddo
        endif
      endif
    endif
    if ( (ivar_out .eq. var_Ti) .or. ( (ivar_out .eq. var_zj) .and. with_TiTe) ) then
      call temperature_i  (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
      if (ivar_out .eq. var_Ti) var_RHS = zTi
    endif
    if ( (ivar_out .eq. var_Te) .or. ( (ivar_out .eq. var_zj) .and. with_TiTe) ) then
      call temperature_e  (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
      if (ivar_out .eq. var_Te) var_RHS = zTe
    endif
    if (ivar_out .eq. var_zj) then
      call FFprime        (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)
    endif
#ifdef WITH_Vpar
    if (ivar_out .eq. var_Vpar) then
      call velocity       (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zV,dV_dpsi,dV_dz,dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
      if (normalized_velocity_profile) then
        var_RHS = zV
      else
        var_RHS = 2.d0 * PI / F0 * R**2 * zV
      endif
    endif
#endif
    if (ivar_out .eq. var_rhon) then
      call neutral_density(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zrn,drn_dpsi,drn_dz,drn_dpsi2,drn_dz2,drn_dpsi_dz,drn_dpsi3,drn_dpsi_dz2, drn_dpsi2_dz)
      var_RHS = zrn
    endif
    
    ! --- Special case for psi, we interpolate from the old grid copied into grid_xpoint_data module
    if (ivar_out .eq. var_psi) then

      if (forceSDN) then
        call find_RZ(node_list,element_list,x_g(ms,mt),-y_g(ms,mt),R_out,Z_out,ielm_out,s_out,t_out,ifail)
        if (ifail .ne. 0) then ! if we can't find the point on the previous grid, use whatever was set before (should be the eqdsk value)
          var_RHS = eq2_g(ms,mt)
        else
          call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,psi,dd1,dd2,dd3,dd4,dd5)
          var_RHS = (psi+eq2_g(ms,mt))/2.0
        end if

      else
        call find_RZ(node_list,element_list,x_g(ms,mt),y_g(ms,mt),R_out,Z_out,ielm_out,s_out,t_out,ifail)
        if (ifail .ne. 0) then ! if we can't find the point on the previous grid, use whatever was set before (should be the eqdsk value)
          var_RHS = eq2_g(ms,mt)
        else
          call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,psi,dd1,dd2,dd3,dd4,dd5)
          var_RHS = psi
        end if

      end if

    end if

    ! --- Special case for current
    if (ivar_out .eq. var_zj) then
      if (with_TiTe) then
        zT = zTi + zTe
        dT_dpsi = dTi_dpsi + dTe_dpsi
      endif
      var_RHS = zFFprime - x_g(ms,mt)**2 * (dn_dpsi * zT + zn * dT_dpsi)
    endif
 
    do i=1,n_vertex_max
 
      do j=1,n_degrees
 
        index_ij = (i-1)*n_degrees + j
 
        v   = h(i,j,ms,mt)  * element%size(i,j)
 
        rhs_ij = + var_RHS
 
        if (itype .eq. 5) then  ! weak form of int v J R dA: the curl is moved onto the test function v
          v_R = (   y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
          v_Z = ( - x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
          select case (ivar_in)
          case (1)  ! J_R   = dB_phi/dZ - dB_Z/dphi / R
            rhs_ij = - R * v_Z * Bfield(3) - v * Bfield_p(2)
          case (2)  ! J_Z   = ( dB_R/dphi - d(R B_phi)/dR ) / R
            rhs_ij = + R * v_R * Bfield(3) + v * Bfield_p(1)
          case (3)  ! J_phi = dB_Z/dR - dB_R/dZ
            rhs_ij = + R * v_Z * Bfield(1) - (R * v_R + v) * Bfield(2)
          end select
          RHS(index_ij) = RHS(index_ij) - rhs_ij * xjac * wst
        else
          RHS(index_ij) = RHS(index_ij) - v * rhs_ij       * factor(ms,mt) * xjac * wst
        endif
        RHS(index_ij) = RHS(index_ij) + v * eq2_g(ms,mt) * factor(ms,mt) * xjac * wst    ! solve for perturbation only
 
        do k=1,n_vertex_max
 
          do l=1,n_degrees
 
            psi   = h(k,l,ms,mt)  * element%size(k,l)
 
            index_kl = (k-1)*n_degrees + l
 
            ELM(index_ij,index_kl) =  ELM(index_ij,index_kl) - psi * v * factor(ms,mt) * xjac * wst
 
          enddo
        enddo
 
      enddo
    enddo
 
  enddo
enddo

!--------------------------------------------------- itype=5: boundary term of the integration by parts
if (itype .eq. 5) then
  do iside = 1, 4
    ! --- sides on the domain boundary have two boundary nodes (as in boundary_from_grid)
    if ( (nodes(iside)%boundary .le. 0) .or. (nodes(mod(iside,4)+1)%boundary .le. 0) ) cycle
    do ms = 1, n_gauss
      select case (iside)  ! (s,t) along the side, counter-clockwise (as elm_coords in mod_edge_elements)
      case (1) ; s_side = xgauss(ms)        ; t_side = 0.d0
      case (2) ; s_side = 1.d0              ; t_side = xgauss(ms)
      case (3) ; s_side = 1.d0 - xgauss(ms) ; t_side = 1.d0
      case (4) ; s_side = 0.d0              ; t_side = 1.d0 - xgauss(ms)
      end select
      call basisfunctions(s_side, t_side, G, G_s, G_t)
      call field_at_point(G, G_s, G_t, R, R_s, R_t, Z_s, Z_t, Bfield, Bfield_p)
      ! --- outward normal times line element; its sign(xjac) cancels with that of the volume terms above
      if (mod(iside,2) .eq. 1) then
        nrm_R = + Z_s ; nrm_Z = - R_s
      else
        nrm_R = + Z_t ; nrm_Z = - R_t
      endif
      if (iside .ge. 3) then
        nrm_R = - nrm_R ; nrm_Z = - nrm_Z
      endif
      select case (ivar_in)
      case (1) ; rhs_ij = + Bfield(3) * nrm_Z
      case (2) ; rhs_ij = - Bfield(3) * nrm_R
      case (3) ; rhs_ij = + Bfield(2) * nrm_R - Bfield(1) * nrm_Z
      end select
      do i=1,n_vertex_max
        do j=1,n_degrees
          index_ij = (i-1)*n_degrees + j
          RHS(index_ij) = RHS(index_ij) - G(i,j) * element%size(i,j) * R * rhs_ij * wgauss(ms)
        enddo
      enddo
    enddo
  enddo
endif

return

contains

!> Magnetic field (B_R, B_Z, B_phi) of harmonic i_harm and its toroidal derivatives (dB_R/dphi, dB_Z/dphi) at one point
!> of the element, from the basis functions G, G_s, G_t evaluated there. Same definitions as in the full MHD models:
!> B_R = (A3_Z - AZ_p)/R,  B_Z = (AR_p - A3_R)/R,  B_phi = AZ_R - AR_Z + F/R
subroutine field_at_point(G, G_s, G_t, R, R_s, R_t, Z_s, Z_t, B, B_p)
  real*8, intent(in)  :: G(n_vertex_max,n_degrees), G_s(n_vertex_max,n_degrees), G_t(n_vertex_max,n_degrees)
  real*8, intent(out) :: R, R_s, R_t, Z_s, Z_t, B(3), B_p(2)
  real*8  :: A(3,2,3), A_R(3,2), A_Z(3,2), F, sz, xj, n2
  integer :: iv, k, m, ivars(3), iharms(2)

  ! --- A(variable: A3/AR/AZ, harmonic: i_harm/partner, value/d_s/d_t)
  ivars  = (/ var_A3, var_AR, var_AZ /)
  iharms = (/ i_harm, i_partner /)
  R = 0.d0; R_s = 0.d0; R_t = 0.d0; Z_s = 0.d0; Z_t = 0.d0; A = 0.d0; F = 0.d0
  do iv = 1, n_vertex_max
    do k = 1, n_degrees
      sz  = element%size(iv,k)
      R   = R   + nodes(iv)%x(1,k,1) * sz * G(iv,k)
      R_s = R_s + nodes(iv)%x(1,k,1) * sz * G_s(iv,k)
      R_t = R_t + nodes(iv)%x(1,k,1) * sz * G_t(iv,k)
      Z_s = Z_s + nodes(iv)%x(1,k,2) * sz * G_s(iv,k)
      Z_t = Z_t + nodes(iv)%x(1,k,2) * sz * G_t(iv,k)
      do m = 1, 2
        A(:,m,1) = A(:,m,1) + nodes(iv)%values(iharms(m),k,ivars) * sz * G(iv,k)
        A(:,m,2) = A(:,m,2) + nodes(iv)%values(iharms(m),k,ivars) * sz * G_s(iv,k)
        A(:,m,3) = A(:,m,3) + nodes(iv)%values(iharms(m),k,ivars) * sz * G_t(iv,k)
      enddo
#ifdef fullmhd
      if (i_harm .eq. 1) F = F + nodes(iv)%Fprof_eq(k) * sz * G(iv,k)
#endif
    enddo
  enddo
  xj  = R_s * Z_t - R_t * Z_s
  A_R = (   Z_t * A(:,:,2) - Z_s * A(:,:,3) ) / xj
  A_Z = ( - R_t * A(:,:,2) + R_s * A(:,:,3) ) / xj
  n2  = float(mode(i_harm))**2

  B(1)   = ( A_Z(1,1) - dphi_fact * A(3,2,1) ) / R
  B(2)   = ( dphi_fact * A(2,2,1) - A_R(1,1) ) / R
  B(3)   = A_R(3,1) - A_Z(2,1) + F / R
  B_p(1) = ( dphi_fact * A_Z(1,2) + n2 * A(3,1,1) ) / R
  B_p(2) = ( - n2 * A(2,1,1) - dphi_fact * A_R(1,2) ) / R

end subroutine field_at_point

end
