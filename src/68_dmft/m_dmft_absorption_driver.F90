!!****m* ABINIT/m_dmft_absorption_driver
!! NAME
!!  m_dmft_absorption_driver
!!
!! FUNCTION
!!  High-level driver for DFT+DMFT absorption spectra with spin-flip excitations.
!!  Orchestrates the workflow: two-particle measurement -> vertex extraction ->
!!  lattice BSE -> optical response -> spectral attribution.
!!
!!  This module implements the Matsubara-axis response (dmft_resp_mode=1) as a
!!  first release. Real-frequency absorption (dmft_resp_mode=2) requires a
!!  real-axis impurity backend that is not yet implemented.
!!
!!  Spectral attribution decomposes the response into spin-conserving and
!!  spin-flip (S+S-, S-S+) channels with orbital resolution, enabling
!!  identification of the physical origin of absorption features.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2026 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_dmft_absorption_driver

 use defs_basis
 use m_abicore
 use m_errors

 use m_dtset, only : dataset_type
 use m_paw_dmft, only : paw_dmft_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type
 use m_dmft_two_particle, only : chi_loc_type, init_chi_loc, destroy_chi_loc, &
                                & compute_chi0_imp, write_chi_loc
 use m_dmft_vertex, only : vertex_irr_type, init_vertex_irr, destroy_vertex_irr, &
                          & extract_vertex_irr
 use m_dmft_lattice_bse, only : lattice_bse_type, init_lattice_bse, destroy_lattice_bse, &
                               & compute_chi0_lattice, solve_lattice_bse
 use m_dmft_optic_kernel, only : optic_kernel_type, init_optic_kernel, destroy_optic_kernel, &
                                & compute_bubble_conductivity, write_optic_kernel
 use m_dmft_spectral_attribution, only : spectral_attribution_type, &
   & init_spectral_attribution, destroy_spectral_attribution, &
   & compute_spectral_attribution, write_spectral_attribution

 implicit none

 private

 public :: dmft_absorption_run

!!***

CONTAINS

!!****f* m_dmft_absorption_driver/dmft_absorption_run
!! NAME
!!  dmft_absorption_run
!!
!! FUNCTION
!!  Main entry point for the DMFT absorption spectrum calculation.
!!  Called after the DFT+DMFT self-consistent cycle has converged,
!!  while the converged Green function is still available.
!!
!! INPUTS
!!  dtset=dataset structure containing input variables
!!  paw_dmft=DFT+DMFT data structure with converged one-particle quantities
!!  cryst_struc=crystal structure
!!  green_imp=converged impurity Green function (from DMFT loop)
!!
!! NOTES
!!  This routine checks preconditions (dmft_resp_mode, dmft_solv, nspinor)
!!  and dispatches to the appropriate calculation stages.
!!
!!  The converged Green function is passed directly from dmft_solve before
!!  it is destroyed, since paw_dmft_type does not store the Green function.
!!
!! SOURCE

subroutine dmft_absorption_run(dtset, paw_dmft, cryst_struc, green_imp)

 type(dataset_type), intent(in) :: dtset
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(crystal_t), intent(in) :: cryst_struc
 type(green_type), intent(in) :: green_imp

!Local variables
 character(len=500) :: msg
 type(chi_loc_type) :: chi_loc
 type(chi_loc_type) :: chi0_loc
 type(vertex_irr_type) :: vertex_irr
 type(lattice_bse_type) :: latt_bse
 type(optic_kernel_type) :: optic_kern
 type(spectral_attribution_type) :: attrib_chi0, attrib_chi
 integer :: nboson, niw_vertex, norb_corr, resp_mode
 integer :: ndim_orb
 real(dp) :: beta

! *********************************************************************

 resp_mode = dtset%dmft_resp_mode
 if (resp_mode == 0) return

 write(msg,'(a)') ' === DFT+DMFT Absorption Spectrum Driver ==='
 call wrtout(std_out, msg)

 ! --- Precondition checks (belt-and-suspenders; m_chkinp should have caught these) ---
 if (dtset%usedmft /= 1) then
   ABI_ERROR('dmft_resp_mode>0 requires usedmft=1')
 end if

 if (dtset%dmft_resp_spinflip == 1) then
   if (dtset%dmft_solv /= 7) then
     ABI_ERROR('dmft_resp_spinflip=1 requires dmft_solv=7 (rotationally invariant Slater)')
   end if
 end if

 if (resp_mode == 2) then
   if (dtset%dmft_resp_realaxis_backend /= 1) then
     ABI_ERROR('dmft_resp_mode=2 requires dmft_resp_realaxis_backend=1 (not yet implemented)')
   end if
   write(msg,'(3a)') &
   'dmft_resp_mode=2 (real-frequency absorption) is not yet implemented.',ch10,&
   'Use dmft_resp_mode=1 for Matsubara-axis response.'
   ABI_ERROR(msg)
 end if

 nboson = dtset%dmft_resp_nboson
 niw_vertex = dtset%dmft_resp_niw_vertex
 beta = one / paw_dmft%temp

 ! --- Determine correlated orbital count ---
 ndim_orb = 2 * paw_dmft%maxlpawu + 1
 norb_corr = paw_dmft%nspinor * ndim_orb

 write(msg,'(a,i4,a,i4,a,i4,a,es14.6)') &
 ' Response parameters: nboson=', nboson, ' niw_vertex=', niw_vertex, &
 ' norb_corr=', norb_corr, ' beta=', beta
 call wrtout(std_out, msg)

 ! --- Stage 1: Compute local impurity bubble chi0 ---
 write(msg,'(a)') ' Stage 1: Computing local impurity bubble chi0_imp'
 call wrtout(std_out, msg)

 call init_chi_loc(chi0_loc, norb_corr, niw_vertex, nboson)
 call compute_chi0_imp(chi0_loc, green_imp, paw_dmft, norb_corr, niw_vertex, nboson)

 ! Write chi0 diagnostics
 call write_chi_loc(chi0_loc, 'DMFT_chi0_imp.dat')

 ! --- Stage 1b: Spectral attribution of chi0 (bubble) ---
 if (dtset%dmft_resp_spinflip == 1 .and. paw_dmft%nspinor == 2) then
   write(msg,'(a)') ' Stage 1b: Spectral attribution of impurity bubble chi0'
   call wrtout(std_out, msg)

   call init_spectral_attribution(attrib_chi0, nboson, ndim_orb, paw_dmft%nspinor)
   call compute_spectral_attribution(attrib_chi0, chi0_loc, paw_dmft%nspinor)
   call write_spectral_attribution(attrib_chi0, 'DMFT_attrib_chi0_imp.dat', beta)
   call destroy_spectral_attribution(attrib_chi0)
 end if

 ! --- Stage 2: Two-particle measurement from impurity solver ---
 write(msg,'(3a)') &
 ' Stage 2: Local two-particle correlation function measurement.',ch10,&
 ' WARNING: TRIQS two-particle measurement interface not yet connected.'
 call wrtout(std_out, msg)

 call init_chi_loc(chi_loc, norb_corr, niw_vertex, nboson)
 ! chi_loc%chi_mat would be filled by the TRIQS interface extension.
 ! For now, chi_loc remains zero (placeholder for future TRIQS integration).

 write(msg,'(3a)') &
 ' WARNING: chi_loc is zero (TRIQS two-particle interface not connected).',ch10,&
 ' Vertex extraction and BSE results will be trivial until this is implemented.'
 call wrtout(std_out, msg)

 ! --- Stage 3: Extract irreducible vertex ---
 write(msg,'(a)') ' Stage 3: Extracting local irreducible vertex Gamma_imp'
 call wrtout(std_out, msg)

 call init_vertex_irr(vertex_irr, norb_corr, niw_vertex, nboson)
 call extract_vertex_irr(vertex_irr, chi_loc, chi0_loc, norb_corr, niw_vertex, nboson)

 ! --- Stage 4: Lattice BSE ---
 write(msg,'(a)') ' Stage 4: Solving lattice Bethe-Salpeter equation'
 call wrtout(std_out, msg)

 call init_lattice_bse(latt_bse, norb_corr, niw_vertex, nboson, paw_dmft%nkpt)
 call compute_chi0_lattice(latt_bse, paw_dmft, norb_corr, niw_vertex, nboson)
 call solve_lattice_bse(latt_bse, vertex_irr, norb_corr, niw_vertex, nboson)

 ! --- Stage 5: Optical kernel (bubble conductivity on Matsubara axis) ---
 if (resp_mode >= 1) then
   write(msg,'(a)') ' Stage 5: Computing Matsubara optical conductivity'
   call wrtout(std_out, msg)

   call init_optic_kernel(optic_kern, nboson, 3)
   call compute_bubble_conductivity(optic_kern, paw_dmft, nboson)

   ! Write optical kernel output
   call write_optic_kernel(optic_kern, 'DMFT_optic_kernel.dat', beta)

   write(msg,'(a)') ' Matsubara-axis response computation complete.'
   call wrtout(std_out, msg)

   call destroy_optic_kernel(optic_kern)
 end if

 ! --- Cleanup ---
 call destroy_lattice_bse(latt_bse)
 call destroy_vertex_irr(vertex_irr)
 call destroy_chi_loc(chi_loc)
 call destroy_chi_loc(chi0_loc)

 write(msg,'(a)') ' === DFT+DMFT Absorption Spectrum Driver Complete ==='
 call wrtout(std_out, msg)

end subroutine dmft_absorption_run

!!***

END MODULE m_dmft_absorption_driver
