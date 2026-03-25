!!****m* ABINIT/m_dmft_spectral_attribution
!! NAME
!!  m_dmft_spectral_attribution
!!
!! FUNCTION
!!  Spectral attribution (decomposition) of two-particle response functions
!!  into physically meaningful channels for DFT+DMFT absorption spectra.
!!
!!  Implements decomposition of the local susceptibility and its bubble
!!  by spin channel (spin-conserving vs spin-flip) and by orbital pair,
!!  enabling identification of the physical origin of absorption peaks.
!!
!!  Spin channel decomposition:
!!    chi_total = chi_{spin-conserving} + chi_{spin-flip}
!!
!!  Spin-flip channels (S+S- and S-S+):
!!    chi^{+-}(iOm) = sum_{m,m'} chi_{(m,up)(m,down),(m',down)(m',up)}(iOm)
!!    chi^{-+}(iOm) = sum_{m,m'} chi_{(m,down)(m,up),(m',up)(m',down)}(iOm)
!!
!!  Orbital-resolved spin-flip:
!!    chi^{+-}_{mm'}(iOm) = chi_{(m,up)(m,down),(m',down)(m',up)}(iOm)
!!
!!  Composite spinor-orbital index convention (for nspinor=2):
!!    alpha = 1..ndim_orb               : spin up,   m = alpha
!!    alpha = ndim_orb+1..2*ndim_orb    : spin down,  m = alpha - ndim_orb
!!  where ndim_orb = 2*lpawu + 1.
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

MODULE m_dmft_spectral_attribution

 use defs_basis
 use m_abicore
 use m_errors

 use m_dmft_two_particle, only : chi_loc_type

 implicit none

 private

 public :: spectral_attribution_type
 public :: init_spectral_attribution
 public :: destroy_spectral_attribution
 public :: compute_spectral_attribution
 public :: write_spectral_attribution
 public :: write_attribution_comparison
 public :: write_attribution_summary
 public :: write_frequency_profile

!!***

!!****t* m_dmft_spectral_attribution/spectral_attribution_type
!! NAME
!!  spectral_attribution_type
!!
!! FUNCTION
!!  Stores the decomposition of a two-particle correlation function
!!  into spin and orbital channels.
!!
!!  All traces are summed over fermionic Matsubara frequencies.
!!
!! SOURCE

 type, public :: spectral_attribution_type

   integer :: nboson = 0
   integer :: ndim_orb = 0       ! Number of orbital indices (2*lpawu+1)
   integer :: norb_corr = 0      ! Full spinor-orbital dimension (ndim_orb * nspinor)
   integer :: nspinor = 0

   ! Channel-resolved susceptibility traces as function of bosonic frequency
   complex(dp), allocatable :: chi_total(:)
   ! chi_total(nboson) : total trace = sum_{n,a,b} chi_{ab,ba}(n,n;iOm)

   complex(dp), allocatable :: chi_spin_conserving(:)
   ! chi_spin_conserving(nboson) : spin-conserving part (sigma_a = sigma_b)

   complex(dp), allocatable :: chi_spin_flip_pm(:)
   ! chi_spin_flip_pm(nboson) : S+S- channel

   complex(dp), allocatable :: chi_spin_flip_mp(:)
   ! chi_spin_flip_mp(nboson) : S-S+ channel

   ! Orbital-resolved S+S- susceptibility
   complex(dp), allocatable :: chi_pm_orbital(:,:,:)
   ! chi_pm_orbital(nboson, ndim_orb, ndim_orb) :
   ! chi^{+-}_{mm'}(iOm) = sum_n chi_{(m,up)(m,down),(m',down)(m',up)}(n,n;iOm)

   ! Orbital-resolved S-S+ susceptibility
   complex(dp), allocatable :: chi_mp_orbital(:,:,:)
   ! chi_mp_orbital(nboson, ndim_orb, ndim_orb) :
   ! chi^{-+}_{mm'}(iOm) = sum_n chi_{(m,down)(m,up),(m',up)(m',down)}(n,n;iOm)

 end type spectral_attribution_type

!!***

CONTAINS

!!****f* m_dmft_spectral_attribution/init_spectral_attribution
!! NAME
!!  init_spectral_attribution
!!
!! FUNCTION
!!  Initialize the spectral attribution structure.
!!
!! INPUTS
!!  nboson = number of bosonic Matsubara frequencies
!!  ndim_orb = number of orbital indices (2*lpawu+1)
!!  nspinor = number of spinor components (1 or 2)
!!
!! SOURCE

subroutine init_spectral_attribution(attrib, nboson, ndim_orb, nspinor)

 type(spectral_attribution_type), intent(inout) :: attrib
 integer, intent(in) :: nboson, ndim_orb, nspinor

!Local variables
 character(len=500) :: msg

! *********************************************************************

 attrib%nboson = nboson
 attrib%ndim_orb = ndim_orb
 attrib%nspinor = nspinor
 attrib%norb_corr = ndim_orb * nspinor

 write(msg,'(a,i4,a,i4,a,i2)') &
 ' init_spectral_attribution: nboson=', nboson, &
 ' ndim_orb=', ndim_orb, ' nspinor=', nspinor
 call wrtout(std_out, msg)

 ABI_MALLOC(attrib%chi_total, (nboson))
 ABI_MALLOC(attrib%chi_spin_conserving, (nboson))
 ABI_MALLOC(attrib%chi_spin_flip_pm, (nboson))
 ABI_MALLOC(attrib%chi_spin_flip_mp, (nboson))
 ABI_MALLOC(attrib%chi_pm_orbital, (nboson, ndim_orb, ndim_orb))
 ABI_MALLOC(attrib%chi_mp_orbital, (nboson, ndim_orb, ndim_orb))

 attrib%chi_total = czero
 attrib%chi_spin_conserving = czero
 attrib%chi_spin_flip_pm = czero
 attrib%chi_spin_flip_mp = czero
 attrib%chi_pm_orbital = czero
 attrib%chi_mp_orbital = czero

end subroutine init_spectral_attribution

!!***

!!****f* m_dmft_spectral_attribution/destroy_spectral_attribution
!! NAME
!!  destroy_spectral_attribution
!!
!! FUNCTION
!!  Deallocate the spectral attribution structure.
!!
!! SOURCE

subroutine destroy_spectral_attribution(attrib)

 type(spectral_attribution_type), intent(inout) :: attrib

! *********************************************************************

 if (allocated(attrib%chi_total)) then
   ABI_FREE(attrib%chi_total)
 end if
 if (allocated(attrib%chi_spin_conserving)) then
   ABI_FREE(attrib%chi_spin_conserving)
 end if
 if (allocated(attrib%chi_spin_flip_pm)) then
   ABI_FREE(attrib%chi_spin_flip_pm)
 end if
 if (allocated(attrib%chi_spin_flip_mp)) then
   ABI_FREE(attrib%chi_spin_flip_mp)
 end if
 if (allocated(attrib%chi_pm_orbital)) then
   ABI_FREE(attrib%chi_pm_orbital)
 end if
 if (allocated(attrib%chi_mp_orbital)) then
   ABI_FREE(attrib%chi_mp_orbital)
 end if

 attrib%nboson = 0
 attrib%ndim_orb = 0
 attrib%norb_corr = 0
 attrib%nspinor = 0

end subroutine destroy_spectral_attribution

!!***

!!****f* m_dmft_spectral_attribution/compute_spectral_attribution
!! NAME
!!  compute_spectral_attribution
!!
!! FUNCTION
!!  Compute the spectral attribution (spin/orbital decomposition)
!!  from a two-particle correlation function chi_loc_type.
!!
!!  For each bosonic frequency iOm, the susceptibility trace is decomposed:
!!
!!  Total:
!!    chi_total(iOm) = sum_{n} sum_{a,b} chi_{ab,ba}(n,n;iOm)
!!
!!  Spin-conserving (sigma_a = sigma_b):
!!    chi_sc(iOm) = sum_{n} sum_{a,b with same spin} chi_{ab,ba}(n,n;iOm)
!!
!!  S+S- channel (alpha=(m,up), beta=(m,down), gamma=(m',down), delta=(m',up)):
!!    chi^{+-}(iOm) = sum_{n} sum_{m,m'} chi_{(m,up)(m,down),(m',down)(m',up)}(n,n;iOm)
!!
!!  S-S+ channel (alpha=(m,down), beta=(m,up), gamma=(m',up), delta=(m',down)):
!!    chi^{-+}(iOm) = sum_{n} sum_{m,m'} chi_{(m,down)(m,up),(m',up)(m',down)}(n,n;iOm)
!!
!! INPUTS
!!  chi = two-particle correlation function
!!  nspinor = number of spinor components (must be 2 for spin-flip decomposition)
!!
!! SIDE EFFECTS
!!  attrib = on output, filled with the decomposed susceptibility traces
!!
!! NOTES
!!  For nspinor=2, the spinor-orbital index alpha is packed as:
!!    alpha = 1..ndim_orb          -> spin up,   orbital m = alpha
!!    alpha = ndim_orb+1..2*ndim_orb -> spin down, orbital m = alpha - ndim_orb
!!  where ndim_orb = 2*lpawu+1.
!!
!!  The trace chi_{ab,ba} means:
!!    I = (n-1)*norb^2 + (a-1)*norb + b
!!    J = (n-1)*norb^2 + (b-1)*norb + a
!!
!! SOURCE

subroutine compute_spectral_attribution(attrib, chi, nspinor)

 type(spectral_attribution_type), intent(inout) :: attrib
 type(chi_loc_type), intent(in) :: chi
 integer, intent(in) :: nspinor

!Local variables
 integer :: iom, iw, ialpha, ibeta, idx_i, idx_j
 integer :: norb_corr, norb_sq, ndim_orb, niw_vertex, nboson
 integer :: ispin_a, ispin_b, im_a, im_b
 character(len=500) :: msg

! *********************************************************************

 norb_corr = chi%norb_corr
 niw_vertex = chi%niw_vertex
 nboson = chi%nboson
 norb_sq = norb_corr * norb_corr

 if (nspinor /= 2) then
   write(msg,'(a,i2)') &
   'compute_spectral_attribution: spin-flip decomposition requires nspinor=2, got ', nspinor
   ABI_ERROR(msg)
 end if

 ndim_orb = norb_corr / nspinor

 write(msg,'(a,i4,a,i4)') &
 ' compute_spectral_attribution: norb_corr=', norb_corr, ' ndim_orb=', ndim_orb
 call wrtout(std_out, msg)

 attrib%chi_total = czero
 attrib%chi_spin_conserving = czero
 attrib%chi_spin_flip_pm = czero
 attrib%chi_spin_flip_mp = czero
 attrib%chi_pm_orbital = czero
 attrib%chi_mp_orbital = czero

 ! --- Compute decomposition ---
 ! For each bosonic frequency, sum over fermionic frequencies and orbital indices.
 ! The trace chi_{ab,ba} uses:
 !   I = (n-1)*norb^2 + (a-1)*norb + b
 !   J = (n-1)*norb^2 + (b-1)*norb + a

 do iom = 1, nboson
   do iw = 1, niw_vertex
     do ialpha = 1, norb_corr
       do ibeta = 1, norb_corr
         idx_i = (iw - 1) * norb_sq + (ialpha - 1) * norb_corr + ibeta
         idx_j = (iw - 1) * norb_sq + (ibeta - 1) * norb_corr + ialpha

         ! Total trace
         attrib%chi_total(iom) = attrib%chi_total(iom) + chi%chi_mat(iom, idx_i, idx_j)

         ! Determine spin indices from composite index
         ! alpha = 1..ndim_orb -> spin up (ispin=1), alpha = ndim_orb+1..2*ndim_orb -> spin down (ispin=2)
         if (ialpha <= ndim_orb) then
           ispin_a = 1
           im_a = ialpha
         else
           ispin_a = 2
           im_a = ialpha - ndim_orb
         end if

         if (ibeta <= ndim_orb) then
           ispin_b = 1
           im_b = ibeta
         else
           ispin_b = 2
           im_b = ibeta - ndim_orb
         end if

         ! Classify by spin channel
         if (ispin_a == ispin_b) then
           ! Spin-conserving: both indices have same spin
           attrib%chi_spin_conserving(iom) = attrib%chi_spin_conserving(iom) + &
             chi%chi_mat(iom, idx_i, idx_j)

         else if (ispin_a == 1 .and. ispin_b == 2) then
           ! S+S- channel: alpha=(m,up), beta=(m',down)
           ! The trace chi_{ab,ba} with a=(m,up), b=(m',down) gives the S+S- component.
           ! In the trace, we contract: c†_{m,up} c_{m',down} * c†_{m',down} c_{m,up}
           ! which is S+_{mm'} * S-_{m'm} at the operator level.
           attrib%chi_spin_flip_pm(iom) = attrib%chi_spin_flip_pm(iom) + &
             chi%chi_mat(iom, idx_i, idx_j)
           attrib%chi_pm_orbital(iom, im_a, im_b) = &
             attrib%chi_pm_orbital(iom, im_a, im_b) + chi%chi_mat(iom, idx_i, idx_j)

         else if (ispin_a == 2 .and. ispin_b == 1) then
           ! S-S+ channel: alpha=(m,down), beta=(m',up)
           attrib%chi_spin_flip_mp(iom) = attrib%chi_spin_flip_mp(iom) + &
             chi%chi_mat(iom, idx_i, idx_j)
           attrib%chi_mp_orbital(iom, im_a, im_b) = &
             attrib%chi_mp_orbital(iom, im_a, im_b) + chi%chi_mat(iom, idx_i, idx_j)
         end if

       end do
     end do
   end do
 end do

 write(msg,'(a)') ' compute_spectral_attribution: Spin/orbital decomposition complete'
 call wrtout(std_out, msg)

end subroutine compute_spectral_attribution

!!***

!!****f* m_dmft_spectral_attribution/write_spectral_attribution
!! NAME
!!  write_spectral_attribution
!!
!! FUNCTION
!!  Write the spectral attribution to file in human-readable format.
!!
!! INPUTS
!!  attrib = spectral attribution data
!!  fname = output file name
!!  beta = inverse temperature (for bosonic frequency values)
!!
!! SOURCE

subroutine write_spectral_attribution(attrib, fname, beta)

 type(spectral_attribution_type), intent(in) :: attrib
 character(len=*), intent(in) :: fname
 real(dp), intent(in) :: beta

!Local variables
 integer :: unt, iom, im, imp, ios
 real(dp) :: omega_boson
 character(len=500) :: msg

! *********************************************************************

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT Spectral Attribution (Spin/Orbital Decomposition)'
 write(unt,'(a,i6)') '# nboson = ', attrib%nboson
 write(unt,'(a,i4)') '# ndim_orb = ', attrib%ndim_orb
 write(unt,'(a,i2)') '# nspinor = ', attrib%nspinor
 write(unt,'(a,es14.6)') '# beta = ', beta
 write(unt,'(a)')

 ! --- Channel-resolved susceptibility traces ---
 write(unt,'(a)') '# === Channel-resolved susceptibility traces ==='
 write(unt,'(a)') '# iOm  Omega_boson  Re(total)  Im(total)  ' // &
   'Re(spin_conserv)  Im(spin_conserv)  Re(S+S-)  Im(S+S-)  Re(S-S+)  Im(S-S+)'

 do iom = 1, attrib%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   write(unt,'(i6,es14.6,8es18.8)') iom, omega_boson, &
     real(attrib%chi_total(iom)), aimag(attrib%chi_total(iom)), &
     real(attrib%chi_spin_conserving(iom)), aimag(attrib%chi_spin_conserving(iom)), &
     real(attrib%chi_spin_flip_pm(iom)), aimag(attrib%chi_spin_flip_pm(iom)), &
     real(attrib%chi_spin_flip_mp(iom)), aimag(attrib%chi_spin_flip_mp(iom))
 end do

 write(unt,'(a)')

 ! --- Orbital-resolved S+S- susceptibility ---
 write(unt,'(a)') '# === Orbital-resolved S+S- susceptibility ==='
 write(unt,'(a)') '# iOm  Omega_boson  m  m_prime  Re(chi_pm)  Im(chi_pm)'

 do iom = 1, attrib%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   do im = 1, attrib%ndim_orb
     do imp = 1, attrib%ndim_orb
       write(unt,'(i6,es14.6,2i4,2es18.8)') iom, omega_boson, im, imp, &
         real(attrib%chi_pm_orbital(iom, im, imp)), &
         aimag(attrib%chi_pm_orbital(iom, im, imp))
     end do
   end do
 end do

 write(unt,'(a)')

 ! --- Orbital-resolved S-S+ susceptibility ---
 write(unt,'(a)') '# === Orbital-resolved S-S+ susceptibility ==='
 write(unt,'(a)') '# iOm  Omega_boson  m  m_prime  Re(chi_mp)  Im(chi_mp)'

 do iom = 1, attrib%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   do im = 1, attrib%ndim_orb
     do imp = 1, attrib%ndim_orb
       write(unt,'(i6,es14.6,2i4,2es18.8)') iom, omega_boson, im, imp, &
         real(attrib%chi_mp_orbital(iom, im, imp)), &
         aimag(attrib%chi_mp_orbital(iom, im, imp))
     end do
   end do
 end do

 close(unt)

 write(msg,'(3a)') ' write_spectral_attribution: Written to ', trim(fname)
 call wrtout(std_out, msg)

end subroutine write_spectral_attribution

!!***

!!****f* m_dmft_spectral_attribution/write_attribution_comparison
!! NAME
!!  write_attribution_comparison
!!
!! FUNCTION
!!  Write a cross-level comparison of spectral attributions to a single file.
!!  Shows how the susceptibility decomposes at three levels:
!!    1. Impurity bubble chi0_imp (local, per-site)
!!    2. Lattice bubble chi0_latt (k-summed, no vertex correction)
!!    3. BSE chi_full (k-summed, with vertex correction from Gamma_imp)
!!
!!  The comparison enables direct identification of:
!!    - Effect of k-point dispersion: chi0_latt - chi0_imp
!!    - Effect of vertex corrections: chi_full - chi0_latt
!!
!!  Output includes total susceptibility and all spin channels at each level.
!!
!! INPUTS
!!  attrib_imp = attribution of total chi0_imp
!!  attrib_latt = attribution of chi0_latt
!!  attrib_bse = attribution of chi_full (BSE solution)
!!  fname = output file name
!!  beta = inverse temperature
!!
!! SOURCE

subroutine write_attribution_comparison(attrib_imp, attrib_latt, attrib_bse, fname, beta)

 type(spectral_attribution_type), intent(in) :: attrib_imp
 type(spectral_attribution_type), intent(in) :: attrib_latt
 type(spectral_attribution_type), intent(in) :: attrib_bse
 character(len=*), intent(in) :: fname
 real(dp), intent(in) :: beta

!Local variables
 integer :: unt, iom, im, imp, ios
 real(dp) :: omega_boson
 complex(dp) :: delta_disp, delta_vtx
 character(len=500) :: msg

! *********************************************************************

 ! Validate that all attributions have the same dimensions
 if (attrib_imp%nboson /= attrib_latt%nboson .or. &
     attrib_imp%nboson /= attrib_bse%nboson) then
   ABI_ERROR('write_attribution_comparison: nboson mismatch between attribution levels')
 end if
 if (attrib_imp%ndim_orb /= attrib_latt%ndim_orb .or. &
     attrib_imp%ndim_orb /= attrib_bse%ndim_orb) then
   ABI_ERROR('write_attribution_comparison: ndim_orb mismatch between attribution levels')
 end if

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT Cross-Level Attribution Comparison'
 write(unt,'(a)') '# Compares spectral decomposition at three levels:'
 write(unt,'(a)') '#   IMP  = impurity bubble chi0_imp (local, per-site sum)'
 write(unt,'(a)') '#   LATT = lattice bubble chi0_latt (k-summed, no vertex)'
 write(unt,'(a)') '#   BSE  = full chi from BSE (k-summed, with vertex)'
 write(unt,'(a)') '# Delta_disp = LATT - IMP  (effect of k-point dispersion)'
 write(unt,'(a)') '# Delta_vtx  = BSE - LATT  (effect of vertex corrections)'
 write(unt,'(a,i6)') '# nboson = ', attrib_imp%nboson
 write(unt,'(a,i4)') '# ndim_orb = ', attrib_imp%ndim_orb
 write(unt,'(a,es14.6)') '# beta = ', beta
 write(unt,'(a)')

 ! --- Section 1: Total susceptibility comparison ---
 write(unt,'(a)') '# === Section 1: Total susceptibility trace ==='
 write(unt,'(a)') '# iOm  Omega_boson  Re(IMP)  Im(IMP)  Re(LATT)  Im(LATT)  ' // &
   'Re(BSE)  Im(BSE)  Re(Delta_disp)  Im(Delta_disp)  Re(Delta_vtx)  Im(Delta_vtx)'

 do iom = 1, attrib_imp%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   delta_disp = attrib_latt%chi_total(iom) - attrib_imp%chi_total(iom)
   delta_vtx = attrib_bse%chi_total(iom) - attrib_latt%chi_total(iom)
   write(unt,'(i6,es14.6,10es18.8)') iom, omega_boson, &
     real(attrib_imp%chi_total(iom)), aimag(attrib_imp%chi_total(iom)), &
     real(attrib_latt%chi_total(iom)), aimag(attrib_latt%chi_total(iom)), &
     real(attrib_bse%chi_total(iom)), aimag(attrib_bse%chi_total(iom)), &
     real(delta_disp), aimag(delta_disp), &
     real(delta_vtx), aimag(delta_vtx)
 end do

 write(unt,'(a)')

 ! --- Section 2: Spin-conserving channel comparison ---
 write(unt,'(a)') '# === Section 2: Spin-conserving channel ==='
 write(unt,'(a)') '# iOm  Omega_boson  Re(IMP)  Im(IMP)  Re(LATT)  Im(LATT)  ' // &
   'Re(BSE)  Im(BSE)  Re(Delta_disp)  Im(Delta_disp)  Re(Delta_vtx)  Im(Delta_vtx)'

 do iom = 1, attrib_imp%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   delta_disp = attrib_latt%chi_spin_conserving(iom) - attrib_imp%chi_spin_conserving(iom)
   delta_vtx = attrib_bse%chi_spin_conserving(iom) - attrib_latt%chi_spin_conserving(iom)
   write(unt,'(i6,es14.6,10es18.8)') iom, omega_boson, &
     real(attrib_imp%chi_spin_conserving(iom)), aimag(attrib_imp%chi_spin_conserving(iom)), &
     real(attrib_latt%chi_spin_conserving(iom)), aimag(attrib_latt%chi_spin_conserving(iom)), &
     real(attrib_bse%chi_spin_conserving(iom)), aimag(attrib_bse%chi_spin_conserving(iom)), &
     real(delta_disp), aimag(delta_disp), &
     real(delta_vtx), aimag(delta_vtx)
 end do

 write(unt,'(a)')

 ! --- Section 3: S+S- channel comparison ---
 write(unt,'(a)') '# === Section 3: S+S- (spin-flip) channel ==='
 write(unt,'(a)') '# iOm  Omega_boson  Re(IMP)  Im(IMP)  Re(LATT)  Im(LATT)  ' // &
   'Re(BSE)  Im(BSE)  Re(Delta_disp)  Im(Delta_disp)  Re(Delta_vtx)  Im(Delta_vtx)'

 do iom = 1, attrib_imp%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   delta_disp = attrib_latt%chi_spin_flip_pm(iom) - attrib_imp%chi_spin_flip_pm(iom)
   delta_vtx = attrib_bse%chi_spin_flip_pm(iom) - attrib_latt%chi_spin_flip_pm(iom)
   write(unt,'(i6,es14.6,10es18.8)') iom, omega_boson, &
     real(attrib_imp%chi_spin_flip_pm(iom)), aimag(attrib_imp%chi_spin_flip_pm(iom)), &
     real(attrib_latt%chi_spin_flip_pm(iom)), aimag(attrib_latt%chi_spin_flip_pm(iom)), &
     real(attrib_bse%chi_spin_flip_pm(iom)), aimag(attrib_bse%chi_spin_flip_pm(iom)), &
     real(delta_disp), aimag(delta_disp), &
     real(delta_vtx), aimag(delta_vtx)
 end do

 write(unt,'(a)')

 ! --- Section 4: S-S+ channel comparison ---
 write(unt,'(a)') '# === Section 4: S-S+ (reverse spin-flip) channel ==='
 write(unt,'(a)') '# iOm  Omega_boson  Re(IMP)  Im(IMP)  Re(LATT)  Im(LATT)  ' // &
   'Re(BSE)  Im(BSE)  Re(Delta_disp)  Im(Delta_disp)  Re(Delta_vtx)  Im(Delta_vtx)'

 do iom = 1, attrib_imp%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   delta_disp = attrib_latt%chi_spin_flip_mp(iom) - attrib_imp%chi_spin_flip_mp(iom)
   delta_vtx = attrib_bse%chi_spin_flip_mp(iom) - attrib_latt%chi_spin_flip_mp(iom)
   write(unt,'(i6,es14.6,10es18.8)') iom, omega_boson, &
     real(attrib_imp%chi_spin_flip_mp(iom)), aimag(attrib_imp%chi_spin_flip_mp(iom)), &
     real(attrib_latt%chi_spin_flip_mp(iom)), aimag(attrib_latt%chi_spin_flip_mp(iom)), &
     real(attrib_bse%chi_spin_flip_mp(iom)), aimag(attrib_bse%chi_spin_flip_mp(iom)), &
     real(delta_disp), aimag(delta_disp), &
     real(delta_vtx), aimag(delta_vtx)
 end do

 write(unt,'(a)')

 ! --- Section 5: Orbital-resolved S+S- comparison (dominant contributions) ---
 write(unt,'(a)') '# === Section 5: Orbital-resolved S+S- comparison ==='
 write(unt,'(a)') '# iOm  Omega_boson  m  m_prime  Re(IMP)  Im(IMP)  Re(LATT)  Im(LATT)  Re(BSE)  Im(BSE)'

 do iom = 1, attrib_imp%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   do im = 1, attrib_imp%ndim_orb
     do imp = 1, attrib_imp%ndim_orb
       write(unt,'(i6,es14.6,2i4,6es18.8)') iom, omega_boson, im, imp, &
         real(attrib_imp%chi_pm_orbital(iom, im, imp)), &
         aimag(attrib_imp%chi_pm_orbital(iom, im, imp)), &
         real(attrib_latt%chi_pm_orbital(iom, im, imp)), &
         aimag(attrib_latt%chi_pm_orbital(iom, im, imp)), &
         real(attrib_bse%chi_pm_orbital(iom, im, imp)), &
         aimag(attrib_bse%chi_pm_orbital(iom, im, imp))
     end do
   end do
 end do

 close(unt)

 write(msg,'(3a)') ' write_attribution_comparison: Written to ', trim(fname)
 call wrtout(std_out, msg)

end subroutine write_attribution_comparison

!!***

!!****f* m_dmft_spectral_attribution/write_attribution_summary
!! NAME
!!  write_attribution_summary
!!
!! FUNCTION
!!  Write a physical interpretation summary of spectral attribution data.
!!
!!  Extracts directly physical quantities from the Matsubara-axis
!!  attribution data:
!!
!!  1. Static susceptibilities chi(iOm=0) for each spin channel.
!!     These are directly measurable (e.g., by neutron scattering for
!!     the transverse channel).
!!
!!  2. Channel fractions: relative contribution of each spin channel
!!     to the total static susceptibility.
!!
!!  3. Matsubara convergence diagnostics: ratio of the highest-frequency
!!     value to the static value. If this ratio is not small, the
!!     Matsubara frequency grid is insufficient and results are unreliable.
!!
!!  4. Dominant orbital pairs in the S+S- channel at the static limit,
!!     sorted by magnitude. This directly identifies which d-orbital
!!     transitions dominate the spin-flip susceptibility.
!!
!! INPUTS
!!  attrib = spectral attribution data (must be computed beforehand)
!!  fname = output file name
!!  beta = inverse temperature
!!
!! NOTES
!!  No heuristic processing is performed. All quantities are exact
!!  within the approximation level used to compute the attribution.
!!  The convergence diagnostic is an assessment, not a correction.
!!
!! SOURCE

subroutine write_attribution_summary(attrib, fname, beta)

 type(spectral_attribution_type), intent(in) :: attrib
 character(len=*), intent(in) :: fname
 real(dp), intent(in) :: beta

!Local variables
 integer :: unt, ios, nboson, ndim_orb, npairs
 real(dp) :: omega_last
 real(dp) :: abs_total_0, abs_pm_0, abs_mp_0, abs_sc_0
 real(dp) :: abs_total_last, abs_pm_last, abs_mp_last, abs_sc_last
 real(dp) :: ratio_conv
 complex(dp) :: chi_static_total, chi_static_pm, chi_static_mp, chi_static_sc
 character(len=500) :: msg

! *********************************************************************

 nboson = attrib%nboson
 ndim_orb = attrib%ndim_orb

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT Attribution Summary — Physical Quantity Extraction'
 write(unt,'(a)') '# All quantities are exact on the Matsubara axis (no heuristics applied).'
 write(unt,'(a,es14.6)') '# beta (inverse temperature) = ', beta
 write(unt,'(a,es14.6)') '# T (temperature in Ha) = ', one/beta
 write(unt,'(a,i6)') '# nboson = ', nboson
 write(unt,'(a,i4)') '# ndim_orb = ', ndim_orb
 write(unt,'(a)')

 ! --- Section 1: Static susceptibilities chi(iOm=0) ---
 ! iOm=0 corresponds to iom=1 (the first bosonic frequency)
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 1: Static susceptibilities chi(iOm=0)'
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# These are the zero-frequency (static) limits of each channel.'
 write(unt,'(a)') '# chi_total = chi_sc + chi_pm + chi_mp must hold exactly.'
 write(unt,'(a)')

 chi_static_total = attrib%chi_total(1)
 chi_static_sc = attrib%chi_spin_conserving(1)
 chi_static_pm = attrib%chi_spin_flip_pm(1)
 chi_static_mp = attrib%chi_spin_flip_mp(1)

 write(unt,'(a,2es22.12)') '# chi_total(0)           Re,Im = ', &
   real(chi_static_total), aimag(chi_static_total)
 write(unt,'(a,2es22.12)') '# chi_spin_conserving(0) Re,Im = ', &
   real(chi_static_sc), aimag(chi_static_sc)
 write(unt,'(a,2es22.12)') '# chi_S+S-(0)            Re,Im = ', &
   real(chi_static_pm), aimag(chi_static_pm)
 write(unt,'(a,2es22.12)') '# chi_S-S+(0)            Re,Im = ', &
   real(chi_static_mp), aimag(chi_static_mp)

 ! Verify sum rule: total = sc + pm + mp
 write(unt,'(a,2es22.12)') '# sum_check (sc+pm+mp)   Re,Im = ', &
   real(chi_static_sc + chi_static_pm + chi_static_mp), &
   aimag(chi_static_sc + chi_static_pm + chi_static_mp)
 write(unt,'(a,2es22.12)') '# difference (total-sum) Re,Im = ', &
   real(chi_static_total - chi_static_sc - chi_static_pm - chi_static_mp), &
   aimag(chi_static_total - chi_static_sc - chi_static_pm - chi_static_mp)
 write(unt,'(a)')

 ! --- Section 2: Channel fractions at static limit ---
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 2: Channel fractions at static limit'
 write(unt,'(a)') '# ================================================================='

 abs_total_0 = abs(chi_static_total)
 abs_sc_0 = abs(chi_static_sc)
 abs_pm_0 = abs(chi_static_pm)
 abs_mp_0 = abs(chi_static_mp)

 if (abs_total_0 > tol14) then
   write(unt,'(a,es14.6)') '# |chi_sc(0)| / |chi_total(0)|  = ', abs_sc_0/abs_total_0
   write(unt,'(a,es14.6)') '# |chi_pm(0)| / |chi_total(0)|  = ', abs_pm_0/abs_total_0
   write(unt,'(a,es14.6)') '# |chi_mp(0)| / |chi_total(0)|  = ', abs_mp_0/abs_total_0
   write(unt,'(a)') '#'
   write(unt,'(a)') '# Interpretation:'
   write(unt,'(a)') '#   Large chi_pm fraction -> spin-flip excitations are significant'
   write(unt,'(a)') '#   chi_pm ~ chi_mp -> time-reversal symmetry approximately holds'
 else
   write(unt,'(a)') '# chi_total(0) is zero or negligible. Channel fractions undefined.'
   write(unt,'(a)') '# This may indicate: (1) system has no susceptibility at this level,'
   write(unt,'(a)') '# or (2) cancellation between channels.'
 end if
 write(unt,'(a)')

 ! --- Section 3: Matsubara convergence diagnostics ---
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 3: Matsubara convergence diagnostics'
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Compares the magnitude at the highest bosonic frequency to the static'
 write(unt,'(a)') '# value. A large ratio indicates insufficient number of bosonic frequencies.'

 if (nboson > 1) then
   omega_last = two_pi * dble(nboson - 1) / beta

   abs_total_last = abs(attrib%chi_total(nboson))
   abs_pm_last = abs(attrib%chi_spin_flip_pm(nboson))
   abs_mp_last = abs(attrib%chi_spin_flip_mp(nboson))
   abs_sc_last = abs(attrib%chi_spin_conserving(nboson))

   write(unt,'(a,es14.6)') '# Highest bosonic Matsubara frequency = ', omega_last
   write(unt,'(a)')
   write(unt,'(a)') '# Channel            |chi(0)|         |chi(max)|       ratio'

   if (abs_total_0 > tol14) then
     ratio_conv = abs_total_last / abs_total_0
     write(unt,'(a,3es16.6)') '# total           ', abs_total_0, abs_total_last, ratio_conv
   else
     write(unt,'(a,3es16.6)') '# total           ', abs_total_0, abs_total_last, zero
   end if

   if (abs_sc_0 > tol14) then
     write(unt,'(a,3es16.6)') '# spin-conserving ', abs_sc_0, abs_sc_last, abs_sc_last/abs_sc_0
   else
     write(unt,'(a,3es16.6)') '# spin-conserving ', abs_sc_0, abs_sc_last, zero
   end if

   if (abs_pm_0 > tol14) then
     write(unt,'(a,3es16.6)') '# S+S-            ', abs_pm_0, abs_pm_last, abs_pm_last/abs_pm_0
   else
     write(unt,'(a,3es16.6)') '# S+S-            ', abs_pm_0, abs_pm_last, zero
   end if

   if (abs_mp_0 > tol14) then
     write(unt,'(a,3es16.6)') '# S-S+            ', abs_mp_0, abs_mp_last, abs_mp_last/abs_mp_0
   else
     write(unt,'(a,3es16.6)') '# S-S+            ', abs_mp_0, abs_mp_last, zero
   end if

   write(unt,'(a)')
   if (abs_total_0 > tol14) then
     ratio_conv = abs_total_last / abs_total_0
     if (ratio_conv > 0.01_dp) then
       write(unt,'(a)') '# WARNING: total ratio > 0.01. Matsubara summation may NOT be converged.'
       write(unt,'(a)') '#          Consider increasing dmft_resp_nboson.'
     else
       write(unt,'(a)') '# Matsubara convergence appears satisfactory (total ratio < 0.01).'
     end if
   end if
 else
   write(unt,'(a)') '# Only 1 bosonic frequency; convergence cannot be assessed.'
   write(unt,'(a)') '# At minimum 2 bosonic frequencies are needed for diagnostics.'
 end if
 write(unt,'(a)')

 ! --- Section 4: Dominant orbital pairs in S+S- channel at static limit ---
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 4: Dominant orbital pairs in S+S- at iOm=0'
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Orbital pairs (m, m_prime) sorted by |chi^{+-}_{mm_prime}(0)|.'
 write(unt,'(a)') '# Large values identify the d-orbital transitions that dominate'
 write(unt,'(a)') '# the spin-flip susceptibility.'
 write(unt,'(a)') '#'
 write(unt,'(a)') '# Rank   m   m_prime   |chi_pm|        Re(chi_pm)         Im(chi_pm)'

 npairs = ndim_orb * ndim_orb
 call write_orbital_ranking(unt, attrib%chi_pm_orbital(1,:,:), ndim_orb, npairs)

 write(unt,'(a)')

 ! --- Section 5: Dominant orbital pairs in S-S+ channel at static limit ---
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 5: Dominant orbital pairs in S-S+ at iOm=0'
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Rank   m   m_prime   |chi_mp|        Re(chi_mp)         Im(chi_mp)'

 call write_orbital_ranking(unt, attrib%chi_mp_orbital(1,:,:), ndim_orb, npairs)

 close(unt)

 write(msg,'(3a)') ' write_attribution_summary: Written to ', trim(fname)
 call wrtout(std_out, msg)

end subroutine write_attribution_summary

!!***

!!****f* m_dmft_spectral_attribution/write_frequency_profile
!! NAME
!!  write_frequency_profile
!!
!! FUNCTION
!!  Write a frequency-by-frequency attribution profile identifying
!!  the dominant spin channel and orbital pair at each bosonic frequency.
!!
!!  For each bosonic Matsubara frequency iOm_m, this output shows:
!!  - Which spin channel (conserving, S+S-, S-S+) carries the largest weight
!!  - Which orbital pair (m, m') has the largest S+S- and S-S+ magnitude
!!  - The ratio of spin-flip to total susceptibility as function of frequency
!!
!!  This enables identifying energy-scale-dependent attribution:
!!  e.g., whether spin-flip processes dominate at low vs high frequencies,
!!  indicating different excitation energy scales.
!!
!! INPUTS
!!  attrib = spectral attribution data (must be computed beforehand)
!!  fname = output file name
!!  beta = inverse temperature
!!
!! NOTES
!!  No heuristic processing. All quantities are computed from the
!!  Matsubara-axis data without analytic continuation.
!!
!! SOURCE

subroutine write_frequency_profile(attrib, fname, beta)

 type(spectral_attribution_type), intent(in) :: attrib
 character(len=*), intent(in) :: fname
 real(dp), intent(in) :: beta

!Local variables
 integer :: unt, ios, iom, im, imp, ndim_orb, nboson
 integer :: max_pm_m, max_pm_mp, max_mp_m, max_mp_mp
 real(dp) :: omega_boson, abs_total, abs_sc, abs_pm, abs_mp
 real(dp) :: frac_sc, frac_pm, frac_mp
 real(dp) :: max_pm_val, max_mp_val, cur_val
 integer :: nchannel_changes
 character(len=20) :: dominant_channel, current_channel
 character(len=500) :: msg

! *********************************************************************

 nboson = attrib%nboson
 ndim_orb = attrib%ndim_orb

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT Frequency-Dependent Attribution Profile'
 write(unt,'(a)') '# For each bosonic frequency, identifies the dominant spin channel'
 write(unt,'(a)') '# and orbital pair contributing to the susceptibility.'
 write(unt,'(a)') '# No analytic continuation or heuristic processing applied.'
 write(unt,'(a,es14.6)') '# beta = ', beta
 write(unt,'(a,i6)') '# nboson = ', nboson
 write(unt,'(a,i4)') '# ndim_orb = ', ndim_orb
 write(unt,'(a)')

 ! --- Section 1: Channel fractions vs frequency ---
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 1: Channel fractions vs bosonic frequency'
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Shows how the relative weight of each spin channel changes with'
 write(unt,'(a)') '# bosonic frequency. A frequency-dependent fraction indicates that'
 write(unt,'(a)') '# different excitation energy scales are governed by different channels.'
 write(unt,'(a)') '#'
 write(unt,'(a)') '# iOm  Omega_boson  |total|  |sc|  |pm|  |mp|  ' // &
   'frac_sc  frac_pm  frac_mp  dominant_channel'

 do iom = 1, nboson
   omega_boson = two_pi * dble(iom - 1) / beta

   abs_total = abs(attrib%chi_total(iom))
   abs_sc = abs(attrib%chi_spin_conserving(iom))
   abs_pm = abs(attrib%chi_spin_flip_pm(iom))
   abs_mp = abs(attrib%chi_spin_flip_mp(iom))

   if (abs_total > tol14) then
     frac_sc = abs_sc / abs_total
     frac_pm = abs_pm / abs_total
     frac_mp = abs_mp / abs_total
   else
     frac_sc = zero
     frac_pm = zero
     frac_mp = zero
   end if

   ! Determine dominant channel
   if (abs_sc >= abs_pm .and. abs_sc >= abs_mp) then
     dominant_channel = 'spin-conserving'
   else if (abs_pm >= abs_sc .and. abs_pm >= abs_mp) then
     dominant_channel = 'S+S-'
   else
     dominant_channel = 'S-S+'
   end if

   write(unt,'(i6,es14.6,4es14.4,3f8.4,2x,a)') &
     iom, omega_boson, abs_total, abs_sc, abs_pm, abs_mp, &
     frac_sc, frac_pm, frac_mp, trim(dominant_channel)
 end do

 write(unt,'(a)')

 ! --- Section 2: Dominant orbital pair at each frequency ---
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 2: Dominant orbital pair at each bosonic frequency'
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# For each frequency, the orbital pair (m, m_prime) with the largest'
 write(unt,'(a)') '# |chi^{+-}_{mm_prime}| and |chi^{-+}_{mm_prime}| is identified.'
 write(unt,'(a)') '# This reveals whether the orbital character of spin-flip excitations'
 write(unt,'(a)') '# changes with excitation energy scale.'
 write(unt,'(a)') '#'
 write(unt,'(a)') '# iOm  Omega_boson  pm_m  pm_mp  |chi_pm|  mp_m  mp_mp  |chi_mp|'

 do iom = 1, nboson
   omega_boson = two_pi * dble(iom - 1) / beta

   ! Find dominant S+S- orbital pair at this frequency
   max_pm_val = zero
   max_pm_m = 1
   max_pm_mp = 1
   do im = 1, ndim_orb
     do imp = 1, ndim_orb
       cur_val = abs(attrib%chi_pm_orbital(iom, im, imp))
       if (cur_val > max_pm_val) then
         max_pm_val = cur_val
         max_pm_m = im
         max_pm_mp = imp
       end if
     end do
   end do

   ! Find dominant S-S+ orbital pair at this frequency
   max_mp_val = zero
   max_mp_m = 1
   max_mp_mp = 1
   do im = 1, ndim_orb
     do imp = 1, ndim_orb
       cur_val = abs(attrib%chi_mp_orbital(iom, im, imp))
       if (cur_val > max_mp_val) then
         max_mp_val = cur_val
         max_mp_m = im
         max_mp_mp = imp
       end if
     end do
   end do

   write(unt,'(i6,es14.6,2i4,es16.6,2i4,es16.6)') &
     iom, omega_boson, &
     max_pm_m, max_pm_mp, max_pm_val, &
     max_mp_m, max_mp_mp, max_mp_val
 end do

 write(unt,'(a)')

 ! --- Section 3: Frequency stability assessment ---
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Section 3: Frequency stability of dominant channels'
 write(unt,'(a)') '# ================================================================='
 write(unt,'(a)') '# Checks whether the dominant channel remains consistent across'
 write(unt,'(a)') '# all bosonic frequencies. A channel change with frequency indicates'
 write(unt,'(a)') '# multiple excitation mechanisms at different energy scales.'
 write(unt,'(a)') '#'

 if (nboson >= 2) then
   ! Check if dominant channel at iOm=0 is the same at other frequencies
   abs_sc = abs(attrib%chi_spin_conserving(1))
   abs_pm = abs(attrib%chi_spin_flip_pm(1))
   abs_mp = abs(attrib%chi_spin_flip_mp(1))

   if (abs_sc >= abs_pm .and. abs_sc >= abs_mp) then
     dominant_channel = 'spin-conserving'
   else if (abs_pm >= abs_sc .and. abs_pm >= abs_mp) then
     dominant_channel = 'S+S-'
   else
     dominant_channel = 'S-S+'
   end if

   write(unt,'(a,a)') '# Dominant channel at iOm=0: ', trim(dominant_channel)

   ! Count frequency points where dominant channel differs
   nchannel_changes = 0
   do iom = 2, nboson
     abs_sc = abs(attrib%chi_spin_conserving(iom))
     abs_pm = abs(attrib%chi_spin_flip_pm(iom))
     abs_mp = abs(attrib%chi_spin_flip_mp(iom))

     if (abs_sc >= abs_pm .and. abs_sc >= abs_mp) then
       current_channel = 'spin-conserving'
     else if (abs_pm >= abs_sc .and. abs_pm >= abs_mp) then
       current_channel = 'S+S-'
     else
       current_channel = 'S-S+'
     end if

     if (trim(current_channel) /= trim(dominant_channel)) then
       nchannel_changes = nchannel_changes + 1
     end if
   end do

   if (nchannel_changes == 0) then
     write(unt,'(a)') '# Dominant channel is CONSISTENT across all bosonic frequencies.'
     write(unt,'(a)') '# This indicates a single excitation mechanism dominates at all energy scales.'
   else
     write(unt,'(a,i4,a,i4,a)') '# Dominant channel CHANGES at ', nchannel_changes, ' of ', nboson - 1, &
       ' frequency points.'
     write(unt,'(a)') '# This indicates MULTIPLE excitation mechanisms at different energy scales.'
     write(unt,'(a)') '# Inspect Section 1 for the frequency-resolved breakdown.'
   end if
 else
   write(unt,'(a)') '# Only 1 bosonic frequency; stability cannot be assessed.'
 end if

 close(unt)

 write(msg,'(3a)') ' write_frequency_profile: Written to ', trim(fname)
 call wrtout(std_out, msg)

end subroutine write_frequency_profile

!!***

! Private helper: write orbital pair ranking sorted by magnitude
subroutine write_orbital_ranking(unt, chi_orbital, ndim_orb, npairs)

 integer, intent(in) :: unt, ndim_orb, npairs
 complex(dp), intent(in) :: chi_orbital(ndim_orb, ndim_orb)

!Local variables
 integer :: im, imp, irank, ii, jj, tmp_m, tmp_mp
 real(dp) :: tmp_val
 integer, allocatable :: rank_m(:), rank_mp(:)
 real(dp), allocatable :: rank_val(:)

! *********************************************************************

 ABI_MALLOC(rank_m, (npairs))
 ABI_MALLOC(rank_mp, (npairs))
 ABI_MALLOC(rank_val, (npairs))

 irank = 0
 do im = 1, ndim_orb
   do imp = 1, ndim_orb
     irank = irank + 1
     rank_m(irank) = im
     rank_mp(irank) = imp
     rank_val(irank) = abs(chi_orbital(im, imp))
   end do
 end do

 ! Sort by magnitude in descending order (insertion sort; npairs is small)
 do ii = 1, npairs - 1
   do jj = ii + 1, npairs
     if (rank_val(jj) > rank_val(ii)) then
       tmp_val = rank_val(ii); rank_val(ii) = rank_val(jj); rank_val(jj) = tmp_val
       tmp_m = rank_m(ii); rank_m(ii) = rank_m(jj); rank_m(jj) = tmp_m
       tmp_mp = rank_mp(ii); rank_mp(ii) = rank_mp(jj); rank_mp(jj) = tmp_mp
     end if
   end do
 end do

 do irank = 1, npairs
   write(unt,'(i6,2i6,es16.6,2es20.10)') irank, rank_m(irank), rank_mp(irank), &
     rank_val(irank), &
     real(chi_orbital(rank_m(irank), rank_mp(irank))), &
     aimag(chi_orbital(rank_m(irank), rank_mp(irank)))
 end do

 ABI_FREE(rank_m)
 ABI_FREE(rank_mp)
 ABI_FREE(rank_val)

end subroutine write_orbital_ranking

!!***

END MODULE m_dmft_spectral_attribution
