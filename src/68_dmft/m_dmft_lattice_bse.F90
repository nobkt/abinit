!!****m* ABINIT/m_dmft_lattice_bse
!! NAME
!!  m_dmft_lattice_bse
!!
!! FUNCTION
!!  Lattice Bethe-Salpeter equation solver for DFT+DMFT response.
!!
!!  Implements the lattice particle-hole bubble:
!!    chi0_latt(q,iw,iOm) = -(1/Nk) sum_k G(k,iw) * G(k+q, iw+iOm)
!!
!!  and the BSE:
!!    chi(q,iOm) = [ chi0(q,iOm)^{-1} - Gamma_imp(iOm) ]^{-1}
!!
!!  For optical absorption, q -> 0.
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

MODULE m_dmft_lattice_bse

 use defs_basis
 use m_abicore
 use m_errors

 use m_hide_lapack, only : xginv
 use m_paw_dmft, only : paw_dmft_type
 use m_dmft_vertex, only : vertex_irr_type

 implicit none

 private

 public :: lattice_bse_type
 public :: init_lattice_bse
 public :: destroy_lattice_bse
 public :: compute_chi0_lattice
 public :: solve_lattice_bse

!!***

!!****t* m_dmft_lattice_bse/lattice_bse_type
!! NAME
!!  lattice_bse_type
!!
!! FUNCTION
!!  Stores the lattice bubble and full susceptibility at q=0.
!!
!! SOURCE

 type, public :: lattice_bse_type

   integer :: norb_corr = 0
   integer :: niw_vertex = 0
   integer :: nboson = 0
   integer :: nkpt = 0
   integer :: ndim_comp = 0

   complex(dp), allocatable :: chi0_latt(:,:,:)
   ! chi0_latt(nboson, ndim_comp, ndim_comp) : lattice bubble at q=0

   complex(dp), allocatable :: chi_full(:,:,:)
   ! chi_full(nboson, ndim_comp, ndim_comp) : full lattice susceptibility at q=0

 end type lattice_bse_type

!!***

CONTAINS

!!****f* m_dmft_lattice_bse/init_lattice_bse
!! NAME
!!  init_lattice_bse
!!
!! FUNCTION
!!  Initialize the lattice BSE data structure.
!!
!! SOURCE

subroutine init_lattice_bse(lbse, norb_corr, niw_vertex, nboson, nkpt)

 type(lattice_bse_type), intent(inout) :: lbse
 integer, intent(in) :: norb_corr, niw_vertex, nboson, nkpt

!Local variables
 integer :: ndim_comp
 character(len=500) :: msg

! *********************************************************************

 lbse%norb_corr = norb_corr
 lbse%niw_vertex = niw_vertex
 lbse%nboson = nboson
 lbse%nkpt = nkpt
 ndim_comp = norb_corr * norb_corr * niw_vertex
 lbse%ndim_comp = ndim_comp

 if (ndim_comp <= 0) then
   ABI_ERROR('init_lattice_bse: ndim_comp must be positive')
 end if

 write(msg,'(a,i8,a,i6)') &
 ' init_lattice_bse: ndim_comp=', ndim_comp, ' nkpt=', nkpt
 call wrtout(std_out, msg)

 ABI_MALLOC(lbse%chi0_latt, (nboson, ndim_comp, ndim_comp))
 ABI_MALLOC(lbse%chi_full, (nboson, ndim_comp, ndim_comp))
 lbse%chi0_latt = czero
 lbse%chi_full = czero

end subroutine init_lattice_bse

!!***

!!****f* m_dmft_lattice_bse/destroy_lattice_bse
!! NAME
!!  destroy_lattice_bse
!!
!! FUNCTION
!!  Deallocate the lattice BSE data structure.
!!
!! SOURCE

subroutine destroy_lattice_bse(lbse)

 type(lattice_bse_type), intent(inout) :: lbse

! *********************************************************************

 if (allocated(lbse%chi0_latt)) then
   ABI_FREE(lbse%chi0_latt)
 end if
 if (allocated(lbse%chi_full)) then
   ABI_FREE(lbse%chi_full)
 end if

 lbse%norb_corr = 0
 lbse%niw_vertex = 0
 lbse%nboson = 0
 lbse%nkpt = 0
 lbse%ndim_comp = 0

end subroutine destroy_lattice_bse

!!***

!!****f* m_dmft_lattice_bse/compute_chi0_lattice
!! NAME
!!  compute_chi0_lattice
!!
!! FUNCTION
!!  Compute the lattice particle-hole bubble at q=0:
!!
!!    chi0_{alpha,beta,gamma,delta}(q=0, iw_n, iOm_m)
!!      = -(1/Nk) sum_k G_{delta,alpha}(k, iw_n) * G_{beta,gamma}(k, iw_n + iOm_m)
!!
!!  The result is stored in the composite-index matrix representation.
!!
!! INPUTS
!!  paw_dmft = DFT+DMFT data with converged lattice Green function
!!  norb_corr = number of correlated spinor-orbitals
!!  niw_vertex = number of fermionic Matsubara frequencies
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SIDE EFFECTS
!!  lbse%chi0_latt = on output, contains the lattice bubble
!!
!! NOTES
!!  Skeleton implementation. Full lattice Green function connection
!!  requires accessing G(k,iw) from paw_dmft and performing the k-sum.
!!  This will be completed when the one-particle infrastructure is connected.
!!
!! SOURCE

subroutine compute_chi0_lattice(lbse, paw_dmft, norb_corr, niw_vertex, nboson)

 type(lattice_bse_type), intent(inout) :: lbse
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 character(len=500) :: msg

! *********************************************************************

 write(msg,'(a)') ' compute_chi0_lattice: Computing lattice bubble at q=0'
 call wrtout(std_out, msg)

 ! NOTE: Full implementation requires:
 ! 1. Extracting G(k, iw_n) from paw_dmft for all k-points
 ! 2. Computing the product G(k,iw)*G(k,iw+iOm) for each (k, iw, iOm)
 ! 3. Projecting onto correlated subspace via spinor projectors
 ! 4. Summing over k-points with 1/Nk normalization
 !
 ! The lattice Green function at each k-point is:
 !   [G(k,iw)]^{-1} = (iw+mu)*S(k) - H^KS(k) - Sigma^emb(k,iw)
 !
 ! This is a structural skeleton. The complete implementation will
 ! follow when one-particle quantities are properly exported.

 lbse%chi0_latt = czero

 write(msg,'(a)') ' compute_chi0_lattice: Lattice bubble structure initialized (skeleton)'
 call wrtout(std_out, msg)

end subroutine compute_chi0_lattice

!!***

!!****f* m_dmft_lattice_bse/solve_lattice_bse
!! NAME
!!  solve_lattice_bse
!!
!! FUNCTION
!!  Solve the lattice Bethe-Salpeter equation:
!!
!!    chi(q=0, iOm) = [ chi0_latt(q=0, iOm)^{-1} - Gamma_imp(iOm) ]^{-1}
!!
!! INPUTS
!!  vertex = local irreducible vertex Gamma_imp
!!  norb_corr = number of correlated spinor-orbitals
!!  niw_vertex = number of fermionic Matsubara frequencies
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SIDE EFFECTS
!!  lbse%chi_full = on output, contains the full lattice susceptibility
!!
!! SOURCE

subroutine solve_lattice_bse(lbse, vertex, norb_corr, niw_vertex, nboson)

 type(lattice_bse_type), intent(inout) :: lbse
 type(vertex_irr_type), intent(in) :: vertex
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 integer :: iom, ndim_comp, ierr
 complex(dp), allocatable :: work_mat(:,:)
 character(len=500) :: msg

! *********************************************************************

 ndim_comp = norb_corr * norb_corr * niw_vertex

 write(msg,'(a,i4,a,i8)') &
 ' solve_lattice_bse: Solving BSE for ', nboson, ' bosonic frequencies, dim=', ndim_comp
 call wrtout(std_out, msg)

 ABI_MALLOC(work_mat, (ndim_comp, ndim_comp))

 do iom = 1, nboson

   ! Compute chi0^{-1} - Gamma
   work_mat(:,:) = lbse%chi0_latt(iom,:,:)
   call xginv(work_mat, ndim_comp, ierr)
   if (ierr /= 0) then
     write(msg,'(a,i4)') &
     ' solve_lattice_bse: chi0_latt inversion failed at bosonic freq ', iom
     ABI_WARNING(msg)
     lbse%chi_full(iom,:,:) = czero
     cycle
   end if

   work_mat(:,:) = work_mat(:,:) - vertex%gamma_mat(iom,:,:)

   ! Invert to get chi_full
   call xginv(work_mat, ndim_comp, ierr)
   if (ierr /= 0) then
     write(msg,'(a,i4)') &
     ' solve_lattice_bse: BSE matrix inversion failed at bosonic freq ', iom
     ABI_WARNING(msg)
     lbse%chi_full(iom,:,:) = czero
     cycle
   end if

   lbse%chi_full(iom,:,:) = work_mat(:,:)

 end do

 ABI_FREE(work_mat)

 write(msg,'(a)') ' solve_lattice_bse: Lattice BSE solution complete'
 call wrtout(std_out, msg)

end subroutine solve_lattice_bse

!!***

END MODULE m_dmft_lattice_bse
