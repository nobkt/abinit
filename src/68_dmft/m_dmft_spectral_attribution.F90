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

END MODULE m_dmft_spectral_attribution
