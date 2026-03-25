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
   & compute_spectral_attribution, write_spectral_attribution, &
   & write_attribution_comparison
 use m_dmft_spinor_proj, only : spinor_proj_type, init_spinor_proj, destroy_spinor_proj, &
   & populate_from_chipsi, check_spinor_completeness

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
 character(len=fnlen) :: fname
 type(chi_loc_type) :: chi_loc
 type(chi_loc_type) :: chi0_loc
 type(chi_loc_type) :: chi0_atom
 type(vertex_irr_type) :: vertex_irr
 type(lattice_bse_type) :: latt_bse
 type(optic_kernel_type) :: optic_kern
 type(spinor_proj_type) :: sproj
 type(spectral_attribution_type) :: attrib_tmp
 type(spectral_attribution_type) :: attrib_imp
 type(spectral_attribution_type) :: attrib_latt
 type(spectral_attribution_type) :: attrib_bse
 logical :: do_spinflip_attrib
 integer :: nboson, niw_vertex, norb_corr, resp_mode
 integer :: ndim_orb, iatom, ncorr_atoms, iatom_latt_count
 integer :: ncorr_atoms_latt
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

 ! Count correlated atoms (total and those usable for lattice bubble)
 ncorr_atoms = 0
 ncorr_atoms_latt = 0
 do iatom = 1, paw_dmft%natom
   if (paw_dmft%lpawu(iatom) >= 0) then
     ncorr_atoms = ncorr_atoms + 1
     if (paw_dmft%lpawu(iatom) == paw_dmft%maxlpawu) then
       ncorr_atoms_latt = ncorr_atoms_latt + 1
     end if
   end if
 end do

 ! Determine if spin-flip attribution should be computed
 do_spinflip_attrib = (dtset%dmft_resp_spinflip == 1 .and. paw_dmft%nspinor == 2)

 write(msg,'(a,i4,a,i4,a,i4,a,es14.6,a,i4,a,i4)') &
 ' Response parameters: nboson=', nboson, ' niw_vertex=', niw_vertex, &
 ' norb_corr=', norb_corr, ' beta=', beta, &
 ' ncorr_atoms=', ncorr_atoms, ' ncorr_atoms_latt=', ncorr_atoms_latt
 call wrtout(std_out, msg)

 if (ncorr_atoms_latt < ncorr_atoms) then
   write(msg,'(a,i4,a,i4,a)') &
   ' WARNING: ', ncorr_atoms - ncorr_atoms_latt, ' of ', ncorr_atoms, &
   ' correlated atoms have lpawu /= maxlpawu and will be excluded from lattice bubble.'
   call wrtout(std_out, msg)
 end if

 ! =====================================================================
 ! Stage 1: Compute local impurity bubble chi0 for each correlated atom
 ! =====================================================================
 write(msg,'(a)') ' Stage 1: Computing local impurity bubble chi0_imp (all correlated atoms)'
 call wrtout(std_out, msg)

 ! Initialize the total (summed over atoms) chi0_imp
 call init_chi_loc(chi0_loc, norb_corr, niw_vertex, nboson)

 do iatom = 1, paw_dmft%natom
   if (paw_dmft%lpawu(iatom) < 0) cycle

   write(msg,'(a,i4,a,i2)') &
   '   Computing chi0_imp for atom ', iatom, ' lpawu=', paw_dmft%lpawu(iatom)
   call wrtout(std_out, msg)

   ! Compute atom-resolved chi0
   call init_chi_loc(chi0_atom, norb_corr, niw_vertex, nboson)
   call compute_chi0_imp(chi0_atom, green_imp, paw_dmft, norb_corr, niw_vertex, nboson, &
     & iatom_index=iatom)

   ! Write atom-resolved chi0 diagnostics
   write(fname,'(a,i3.3,a)') 'DMFT_chi0_imp_atom', iatom, '.dat'
   call write_chi_loc(chi0_atom, trim(fname))

   ! Stage 1b: Atom-resolved spectral attribution
   if (do_spinflip_attrib) then
     write(msg,'(a,i4)') '   Stage 1b: Spectral attribution for atom ', iatom
     call wrtout(std_out, msg)

     call init_spectral_attribution(attrib_tmp, nboson, ndim_orb, paw_dmft%nspinor)
     call compute_spectral_attribution(attrib_tmp, chi0_atom, paw_dmft%nspinor)
     write(fname,'(a,i3.3,a)') 'DMFT_attrib_chi0_imp_atom', iatom, '.dat'
     call write_spectral_attribution(attrib_tmp, trim(fname), beta)
     call destroy_spectral_attribution(attrib_tmp)
   end if

   ! Accumulate into total chi0
   chi0_loc%chi_mat(:,:,:) = chi0_loc%chi_mat(:,:,:) + chi0_atom%chi_mat(:,:,:)

   call destroy_chi_loc(chi0_atom)
 end do

 ! Write total chi0 diagnostics
 call write_chi_loc(chi0_loc, 'DMFT_chi0_imp_total.dat')

 ! Total spectral attribution (sum over atoms) — kept alive for cross-level comparison
 if (do_spinflip_attrib) then
   write(msg,'(a)') '   Stage 1b: Spectral attribution of total chi0_imp'
   call wrtout(std_out, msg)

   call init_spectral_attribution(attrib_imp, nboson, ndim_orb, paw_dmft%nspinor)
   call compute_spectral_attribution(attrib_imp, chi0_loc, paw_dmft%nspinor)
   call write_spectral_attribution(attrib_imp, 'DMFT_attrib_chi0_imp_total.dat', beta)
   ! attrib_imp is kept alive for Stage 6 cross-level comparison
 end if

 ! =====================================================================
 ! Stage 2: Two-particle measurement from impurity solver
 ! =====================================================================
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

 ! =====================================================================
 ! Stage 3: Extract irreducible vertex
 ! =====================================================================
 write(msg,'(a)') ' Stage 3: Extracting local irreducible vertex Gamma_imp'
 call wrtout(std_out, msg)

 call init_vertex_irr(vertex_irr, norb_corr, niw_vertex, nboson)
 call extract_vertex_irr(vertex_irr, chi_loc, chi0_loc, norb_corr, niw_vertex, nboson)

 ! =====================================================================
 ! Stage 4: Lattice BSE (multi-atom lattice bubble accumulation)
 ! =====================================================================
 write(msg,'(a)') ' Stage 4: Solving lattice Bethe-Salpeter equation'
 call wrtout(std_out, msg)

 call init_lattice_bse(latt_bse, norb_corr, niw_vertex, nboson, paw_dmft%nkpt)
 call init_spinor_proj(sproj, paw_dmft)

 ! --- Multi-atom lattice bubble: loop over all correlated atoms ---
 ! Each atom contributes through its own spinor projectors P_atom(k).
 ! The total lattice bubble is the sum of contributions from all atoms:
 !   chi0^latt_total = sum_atom chi0^latt_atom
 ! where each atom's contribution uses G^loc_atom(k,iw) = P_atom G^KS P_atom^dag.
 !
 ! Atoms with lpawu /= maxlpawu are skipped because the spinor projector
 ! dimension (nspinor*(2*lpawu+1)) would not match norb_corr. This is a
 ! structural constraint, not a heuristic: the composite index space is
 ! defined by maxlpawu, and atoms with different lpawu have incompatible
 ! orbital spaces that cannot be directly accumulated.

 iatom_latt_count = 0
 do iatom = 1, paw_dmft%natom
   if (paw_dmft%lpawu(iatom) < 0) cycle
   if (paw_dmft%lpawu(iatom) /= paw_dmft%maxlpawu) then
     write(msg,'(a,i4,a,i2,a,i2,a)') &
     '   Skipping atom ', iatom, ' for lattice bubble: lpawu=', paw_dmft%lpawu(iatom), &
     ' /= maxlpawu=', paw_dmft%maxlpawu, ' (incompatible orbital dimension)'
     call wrtout(std_out, msg)
     cycle
   end if

   iatom_latt_count = iatom_latt_count + 1

   write(msg,'(a,i4,a,i2,a,i4,a,i4)') &
   '   Computing lattice bubble for atom ', iatom, &
   ' lpawu=', paw_dmft%lpawu(iatom), &
   ' (', iatom_latt_count, ' of ', ncorr_atoms_latt, ')'
   call wrtout(std_out, msg)

   ! Populate spinor projectors for this atom
   call populate_from_chipsi(sproj, paw_dmft, iatom)
   call check_spinor_completeness(sproj, tol4)

   ! Compute lattice bubble contribution from this atom
   ! First atom: ladd=.false. (zeros chi0_latt then fills)
   ! Subsequent atoms: ladd=.true. (accumulates into existing chi0_latt)
   call compute_chi0_lattice(lbse=latt_bse, green_imp=green_imp, paw_dmft=paw_dmft, &
     & sproj=sproj, norb_corr=norb_corr, niw_vertex=niw_vertex, nboson=nboson, &
     & ladd=(iatom_latt_count > 1))

 end do

 if (iatom_latt_count == 0) then
   write(msg,'(a)') &
   ' WARNING: No atoms contributed to the lattice bubble (no correlated atoms with lpawu=maxlpawu).'
   call wrtout(std_out, msg)
 else
   write(msg,'(a,i4,a)') &
   ' Lattice bubble accumulated from ', iatom_latt_count, ' correlated atom(s).'
   call wrtout(std_out, msg)
 end if

 ! Write lattice bubble diagnostics (reuse chi_loc write format)
 call write_chi_mat_as_chi_loc(latt_bse%chi0_latt, norb_corr, niw_vertex, nboson, &
   & 'DMFT_chi0_lattice.dat')

 ! Stage 4b: Spectral attribution of total lattice bubble — kept alive for comparison
 if (do_spinflip_attrib) then
   write(msg,'(a)') '   Stage 4b: Spectral attribution of lattice bubble chi0_latt'
   call wrtout(std_out, msg)

   call init_spectral_attribution(attrib_latt, nboson, ndim_orb, paw_dmft%nspinor)
   call init_chi_loc(chi0_atom, norb_corr, niw_vertex, nboson)
   chi0_atom%chi_mat(:,:,:) = latt_bse%chi0_latt(:,:,:)
   call compute_spectral_attribution(attrib_latt, chi0_atom, paw_dmft%nspinor)
   call destroy_chi_loc(chi0_atom)
   call write_spectral_attribution(attrib_latt, 'DMFT_attrib_chi0_lattice.dat', beta)
   ! attrib_latt is kept alive for Stage 6 cross-level comparison
 end if

 ! --- Solve BSE ---
 call solve_lattice_bse(latt_bse, vertex_irr, norb_corr, niw_vertex, nboson)

 ! Stage 4c: Spectral attribution of BSE-corrected chi — kept alive for comparison
 if (do_spinflip_attrib) then
   write(msg,'(a)') '   Stage 4c: Spectral attribution of BSE chi_full'
   call wrtout(std_out, msg)

   call init_spectral_attribution(attrib_bse, nboson, ndim_orb, paw_dmft%nspinor)
   call init_chi_loc(chi0_atom, norb_corr, niw_vertex, nboson)
   chi0_atom%chi_mat(:,:,:) = latt_bse%chi_full(:,:,:)
   call compute_spectral_attribution(attrib_bse, chi0_atom, paw_dmft%nspinor)
   call destroy_chi_loc(chi0_atom)
   call write_spectral_attribution(attrib_bse, 'DMFT_attrib_chi_full.dat', beta)
   ! attrib_bse is kept alive for Stage 6 cross-level comparison
 end if

 ! =====================================================================
 ! Stage 5: Optical kernel (bubble conductivity on Matsubara axis)
 ! =====================================================================
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

 ! =====================================================================
 ! Stage 6: Cross-level attribution comparison
 ! =====================================================================
 ! Compare attributions across three levels: chi0_imp → chi0_latt → chi_full
 ! to quantify the effects of k-point dispersion and vertex corrections.
 if (do_spinflip_attrib) then
   write(msg,'(a)') ' Stage 6: Cross-level attribution comparison'
   call wrtout(std_out, msg)

   call write_attribution_comparison(attrib_imp, attrib_latt, attrib_bse, &
     & 'DMFT_attrib_comparison.dat', beta)

   write(msg,'(a)') ' Cross-level comparison written to DMFT_attrib_comparison.dat'
   call wrtout(std_out, msg)
 end if

 ! =====================================================================
 ! Cleanup
 ! =====================================================================
 if (do_spinflip_attrib) then
   call destroy_spectral_attribution(attrib_imp)
   call destroy_spectral_attribution(attrib_latt)
   call destroy_spectral_attribution(attrib_bse)
 end if
 call destroy_spinor_proj(sproj)
 call destroy_lattice_bse(latt_bse)
 call destroy_vertex_irr(vertex_irr)
 call destroy_chi_loc(chi_loc)
 call destroy_chi_loc(chi0_loc)

 write(msg,'(a)') ' === DFT+DMFT Absorption Spectrum Driver Complete ==='
 call wrtout(std_out, msg)

end subroutine dmft_absorption_run

!!***

!!****f* m_dmft_absorption_driver/write_chi_mat_as_chi_loc
!! NAME
!!  write_chi_mat_as_chi_loc
!!
!! FUNCTION
!!  Write a chi_mat array to file using the chi_loc_type write format.
!!  Creates a temporary chi_loc_type, copies data, writes, and destroys.
!!
!! SOURCE

subroutine write_chi_mat_as_chi_loc(chi_mat, norb_corr, niw_vertex, nboson, fname)

 complex(dp), intent(in) :: chi_mat(:,:,:)
 integer, intent(in) :: norb_corr, niw_vertex, nboson
 character(len=*), intent(in) :: fname

!Local variables
 type(chi_loc_type) :: chi_tmp

! *********************************************************************

 call init_chi_loc(chi_tmp, norb_corr, niw_vertex, nboson)
 chi_tmp%chi_mat(:,:,:) = chi_mat(:,:,:)
 call write_chi_loc(chi_tmp, trim(fname))
 call destroy_chi_loc(chi_tmp)

end subroutine write_chi_mat_as_chi_loc

!!***

!!****f* m_dmft_absorption_driver/attribute_chi_mat
!! NAME
!!  attribute_chi_mat
!!
!! FUNCTION
!!  Perform spectral attribution on a chi_mat array and write to file.
!!  Creates temporary chi_loc_type and spectral_attribution_type objects.
!!
!! SOURCE

subroutine attribute_chi_mat(chi_mat, norb_corr, niw_vertex, nboson, &
  & ndim_orb, nspinor, beta, fname)

 complex(dp), intent(in) :: chi_mat(:,:,:)
 integer, intent(in) :: norb_corr, niw_vertex, nboson, ndim_orb, nspinor
 real(dp), intent(in) :: beta
 character(len=*), intent(in) :: fname

!Local variables
 type(chi_loc_type) :: chi_tmp
 type(spectral_attribution_type) :: attrib

! *********************************************************************

 call init_spectral_attribution(attrib, nboson, ndim_orb, nspinor)
 call init_chi_loc(chi_tmp, norb_corr, niw_vertex, nboson)
 chi_tmp%chi_mat(:,:,:) = chi_mat(:,:,:)
 call compute_spectral_attribution(attrib, chi_tmp, nspinor)
 call write_spectral_attribution(attrib, trim(fname), beta)
 call destroy_spectral_attribution(attrib)
 call destroy_chi_loc(chi_tmp)

end subroutine attribute_chi_mat

!!***

END MODULE m_dmft_absorption_driver
