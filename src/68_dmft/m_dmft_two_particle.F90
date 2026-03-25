!!****m* ABINIT/m_dmft_two_particle
!! NAME
!!  m_dmft_two_particle
!!
!! FUNCTION
!!  Data structures and routines for local two-particle correlation functions
!!  in the particle-hole channel. Implements chi_imp and chi0_imp storage
!!  with full spinor indices following the Fourier convention of
!!  Rohringer et al., Rev. Mod. Phys. 90, 025003 (2018).
!!
!!  Composite index: I = (alpha, beta, n) where alpha,beta are spinor-orbital
!!  indices and n is a fermionic Matsubara frequency index.
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

MODULE m_dmft_two_particle

 use defs_basis
 use m_abicore
 use m_errors

 use m_paw_dmft, only : paw_dmft_type

 implicit none

 private

 public :: chi_loc_type
 public :: init_chi_loc
 public :: destroy_chi_loc
 public :: compute_chi0_imp
 public :: write_chi_loc

!!***

!!****t* m_dmft_two_particle/chi_loc_type
!! NAME
!!  chi_loc_type
!!
!! FUNCTION
!!  Stores the local two-particle correlation function in the particle-hole channel.
!!
!!  chi_mat(iOm, I, J) where:
!!    iOm = bosonic Matsubara frequency index (1..nboson)
!!    I = composite index (alpha, beta, n) packed as:
!!        I = (n-1)*norb^2 + (alpha-1)*norb + beta
!!    J = composite index (gamma, delta, n') with same packing
!!
!!  Total composite dimension: ndim_comp = norb_corr^2 * niw_vertex
!!
!! SOURCE

 type, public :: chi_loc_type

   integer :: norb_corr = 0
   integer :: niw_vertex = 0
   integer :: nboson = 0
   integer :: ndim_comp = 0

   complex(dp), allocatable :: chi_mat(:,:,:)
   ! chi_mat(nboson, ndim_comp, ndim_comp)

 end type chi_loc_type

!!***

CONTAINS

!!****f* m_dmft_two_particle/init_chi_loc
!! NAME
!!  init_chi_loc
!!
!! FUNCTION
!!  Initialize the local two-particle correlation function structure.
!!
!! INPUTS
!!  norb_corr = number of correlated spinor-orbitals
!!  niw_vertex = number of fermionic Matsubara frequencies for vertex
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SOURCE

subroutine init_chi_loc(chi, norb_corr, niw_vertex, nboson)

 type(chi_loc_type), intent(inout) :: chi
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 integer :: ndim_comp
 character(len=500) :: msg

! *********************************************************************

 chi%norb_corr = norb_corr
 chi%niw_vertex = niw_vertex
 chi%nboson = nboson
 ndim_comp = norb_corr * norb_corr * niw_vertex
 chi%ndim_comp = ndim_comp

 if (ndim_comp <= 0) then
   ABI_ERROR('init_chi_loc: ndim_comp must be positive')
 end if

 write(msg,'(a,i8,a,3(i6,a))') &
 ' init_chi_loc: ndim_comp=', ndim_comp, &
 ' (norb=', norb_corr, ', niw=', niw_vertex, ', nboson=', nboson, ')'
 call wrtout(std_out, msg)

 ABI_MALLOC(chi%chi_mat, (nboson, ndim_comp, ndim_comp))
 chi%chi_mat = czero

end subroutine init_chi_loc

!!***

!!****f* m_dmft_two_particle/destroy_chi_loc
!! NAME
!!  destroy_chi_loc
!!
!! FUNCTION
!!  Deallocate the local two-particle correlation function structure.
!!
!! SOURCE

subroutine destroy_chi_loc(chi)

 type(chi_loc_type), intent(inout) :: chi

! *********************************************************************

 if (allocated(chi%chi_mat)) then
   ABI_FREE(chi%chi_mat)
 end if

 chi%norb_corr = 0
 chi%niw_vertex = 0
 chi%nboson = 0
 chi%ndim_comp = 0

end subroutine destroy_chi_loc

!!***

!!****f* m_dmft_two_particle/compute_chi0_imp
!! NAME
!!  compute_chi0_imp
!!
!! FUNCTION
!!  Compute the bare impurity bubble chi0_imp from the impurity Green function.
!!
!!  chi0_imp_{alpha,beta,gamma,delta}(iw_n, iw_n'; iOm_m)
!!    = -beta * delta_{n,n'} * G^imp_{delta,alpha}(iw_n) * G^imp_{beta,gamma}(iw_n + iOm_m)
!!
!!  Following the Fourier convention of Section 5.9 of the design document.
!!
!! INPUTS
!!  paw_dmft = DFT+DMFT data containing the converged impurity Green function
!!  norb_corr = number of correlated spinor-orbitals
!!  niw_vertex = number of fermionic Matsubara frequencies
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SIDE EFFECTS
!!  chi0 = on output, contains the computed bare bubble
!!
!! NOTES
!!  Currently a skeleton implementation. The actual impurity Green function
!!  extraction from paw_dmft requires connecting to the green_type data.
!!
!! SOURCE

subroutine compute_chi0_imp(chi0, paw_dmft, norb_corr, niw_vertex, nboson)

 type(chi_loc_type), intent(inout) :: chi0
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 integer :: iom, iw, ialpha, ibeta, igamma, idelta
 integer :: idx_i, idx_j
 character(len=500) :: msg

! *********************************************************************

 write(msg,'(a)') ' compute_chi0_imp: Computing bare impurity bubble'
 call wrtout(std_out, msg)

 ! The bare bubble is diagonal in fermionic frequency: delta_{n,n'}
 ! chi0(iOm, I, J) where I=(alpha,beta,n), J=(gamma,delta,n')
 ! is nonzero only when n=n' in the I,J composite indices.
 !
 ! Formula: chi0_{ab,gd}(n,n';Om) = -beta * delta_{nn'} * G_{da}(iwn) * G_{bg}(iwn+iOm)
 !
 ! NOTE: Actual implementation requires extracting G^imp from paw_dmft.
 ! This is a structural skeleton that establishes the correct index layout.
 ! The Green function data connection will be completed when the TRIQS
 ! two-particle interface is available.

 chi0%chi_mat = czero

 write(msg,'(a,i8,a)') &
 ' compute_chi0_imp: Bubble structure initialized (', chi0%ndim_comp, ' composite indices)'
 call wrtout(std_out, msg)

 write(msg,'(a)') &
 ' compute_chi0_imp: NOTE - Full G^imp connection pending TRIQS interface extension'
 call wrtout(std_out, msg)

end subroutine compute_chi0_imp

!!***

!!****f* m_dmft_two_particle/write_chi_loc
!! NAME
!!  write_chi_loc
!!
!! FUNCTION
!!  Write the local two-particle correlation function to file for diagnostics.
!!
!! INPUTS
!!  chi = local two-particle correlation function
!!  fname = output file name
!!
!! SOURCE

subroutine write_chi_loc(chi, fname)

 type(chi_loc_type), intent(in) :: chi
 character(len=*), intent(in) :: fname

!Local variables
 integer :: unt, iom, ii, jj
 character(len=500) :: msg

! *********************************************************************

 if (open_file(fname, msg, newunit=unt, form='formatted', action='write') /= 0) then
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# Local two-particle correlation function'
 write(unt,'(a,i6)') '# norb_corr = ', chi%norb_corr
 write(unt,'(a,i6)') '# niw_vertex = ', chi%niw_vertex
 write(unt,'(a,i6)') '# nboson = ', chi%nboson
 write(unt,'(a,i8)') '# ndim_comp = ', chi%ndim_comp
 write(unt,'(a)') '# iOm  I  J  Re(chi)  Im(chi)'

 do iom = 1, chi%nboson
   do ii = 1, min(chi%ndim_comp, 10)  ! Print only first 10 for diagnostics
     do jj = 1, min(chi%ndim_comp, 10)
       write(unt,'(3i6,2es20.10)') iom, ii, jj, &
       real(chi%chi_mat(iom,ii,jj)), aimag(chi%chi_mat(iom,ii,jj))
     end do
   end do
 end do

 close(unt)

end subroutine write_chi_loc

!!***

! Private helper: open file
integer function open_file(fname, msg, newunit, form, action)
 character(len=*), intent(in) :: fname, form, action
 character(len=*), intent(out) :: msg
 integer, intent(out) :: newunit
 integer :: ios

 open_file = 0
 open(newunit=newunit, file=fname, form=form, action=action, iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   open_file = 1
 end if
end function open_file

END MODULE m_dmft_two_particle
