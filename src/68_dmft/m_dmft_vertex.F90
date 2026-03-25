!!****m* ABINIT/m_dmft_vertex
!! NAME
!!  m_dmft_vertex
!!
!! FUNCTION
!!  Construction and storage of the local irreducible vertex in the
!!  particle-hole channel for DFT+DMFT response calculations.
!!
!!  The vertex is defined via:
!!    Gamma_imp(iOm) = [chi0_imp(iOm)]^{-1} - [chi_imp(iOm)]^{-1}
!!
!!  where the inverse is over the composite index I=(alpha,beta,n).
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

MODULE m_dmft_vertex

 use defs_basis
 use m_abicore
 use m_errors

 use m_hide_lapack, only : xginv
 use m_dmft_two_particle, only : chi_loc_type

 implicit none

 private

 public :: vertex_irr_type
 public :: init_vertex_irr
 public :: destroy_vertex_irr
 public :: extract_vertex_irr
 public :: write_vertex_irr

!!***

!!****t* m_dmft_vertex/vertex_irr_type
!! NAME
!!  vertex_irr_type
!!
!! FUNCTION
!!  Stores the local particle-hole irreducible vertex Gamma_imp(iOm_m).
!!
!!  gamma_mat(iOm, I, J) where I,J are composite indices
!!  with the same packing as chi_loc_type.
!!
!! SOURCE

 type, public :: vertex_irr_type

   integer :: norb_corr = 0
   integer :: niw_vertex = 0
   integer :: nboson = 0
   integer :: ndim_comp = 0

   complex(dp), allocatable :: gamma_mat(:,:,:)
   ! gamma_mat(nboson, ndim_comp, ndim_comp)

 end type vertex_irr_type

!!***

CONTAINS

!!****f* m_dmft_vertex/init_vertex_irr
!! NAME
!!  init_vertex_irr
!!
!! FUNCTION
!!  Initialize the irreducible vertex structure.
!!
!! SOURCE

subroutine init_vertex_irr(vertex, norb_corr, niw_vertex, nboson)

 type(vertex_irr_type), intent(inout) :: vertex
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 integer :: ndim_comp
 character(len=500) :: msg

! *********************************************************************

 vertex%norb_corr = norb_corr
 vertex%niw_vertex = niw_vertex
 vertex%nboson = nboson
 ndim_comp = norb_corr * norb_corr * niw_vertex
 vertex%ndim_comp = ndim_comp

 if (ndim_comp <= 0) then
   ABI_ERROR('init_vertex_irr: ndim_comp must be positive')
 end if

 write(msg,'(a,i8)') ' init_vertex_irr: ndim_comp=', ndim_comp
 call wrtout(std_out, msg)

 ABI_MALLOC(vertex%gamma_mat, (nboson, ndim_comp, ndim_comp))
 vertex%gamma_mat = czero

end subroutine init_vertex_irr

!!***

!!****f* m_dmft_vertex/destroy_vertex_irr
!! NAME
!!  destroy_vertex_irr
!!
!! FUNCTION
!!  Deallocate the irreducible vertex structure.
!!
!! SOURCE

subroutine destroy_vertex_irr(vertex)

 type(vertex_irr_type), intent(inout) :: vertex

! *********************************************************************

 if (allocated(vertex%gamma_mat)) then
   ABI_FREE(vertex%gamma_mat)
 end if

 vertex%norb_corr = 0
 vertex%niw_vertex = 0
 vertex%nboson = 0
 vertex%ndim_comp = 0

end subroutine destroy_vertex_irr

!!***

!!****f* m_dmft_vertex/extract_vertex_irr
!! NAME
!!  extract_vertex_irr
!!
!! FUNCTION
!!  Extract the irreducible vertex from the full and bare local
!!  two-particle correlation functions:
!!
!!    Gamma(iOm) = [chi0(iOm)]^{-1} - [chi(iOm)]^{-1}
!!
!!  For each bosonic frequency iOm, this involves two matrix inversions
!!  of ndim_comp x ndim_comp complex matrices.
!!
!! INPUTS
!!  chi_loc = full impurity two-particle correlation function
!!  chi0_loc = bare impurity bubble
!!  norb_corr = number of correlated spinor-orbitals
!!  niw_vertex = number of fermionic Matsubara frequencies
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SIDE EFFECTS
!!  vertex = on output, contains the extracted irreducible vertex
!!
!! NOTES
!!  The matrix inversion is numerically sensitive due to QMC noise.
!!  No smoothing or regularization is applied here (forbidden by design).
!!  If the inversion fails, the routine reports the error explicitly.
!!
!! SOURCE

subroutine extract_vertex_irr(vertex, chi_loc, chi0_loc, norb_corr, niw_vertex, nboson)

 type(vertex_irr_type), intent(inout) :: vertex
 type(chi_loc_type), intent(in) :: chi_loc
 type(chi_loc_type), intent(in) :: chi0_loc
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 integer :: iom, ndim_comp, ierr
 complex(dp), allocatable :: chi0_inv(:,:), chi_inv(:,:)
 character(len=500) :: msg

! *********************************************************************

 ndim_comp = norb_corr * norb_corr * niw_vertex

 write(msg,'(a,i4,a,i8)') &
 ' extract_vertex_irr: Processing ', nboson, ' bosonic frequencies, matrix dim=', ndim_comp
 call wrtout(std_out, msg)

 ABI_MALLOC(chi0_inv, (ndim_comp, ndim_comp))
 ABI_MALLOC(chi_inv, (ndim_comp, ndim_comp))

 do iom = 1, nboson

   ! Copy chi0 for this bosonic frequency and invert
   chi0_inv(:,:) = chi0_loc%chi_mat(iom,:,:)
   call xginv(chi0_inv, ndim_comp, ierr)
   if (ierr /= 0) then
     write(msg,'(a,i4,a)') &
     ' extract_vertex_irr: chi0 inversion failed at bosonic freq ', iom, &
     '. Matrix may be singular (check QMC statistics or frequency truncation).'
     ABI_WARNING(msg)
     vertex%gamma_mat(iom,:,:) = czero
     cycle
   end if

   ! Copy chi for this bosonic frequency and invert
   chi_inv(:,:) = chi_loc%chi_mat(iom,:,:)
   call xginv(chi_inv, ndim_comp, ierr)
   if (ierr /= 0) then
     write(msg,'(a,i4,a)') &
     ' extract_vertex_irr: chi inversion failed at bosonic freq ', iom, &
     '. Matrix may be singular (check QMC statistics or frequency truncation).'
     ABI_WARNING(msg)
     vertex%gamma_mat(iom,:,:) = czero
     cycle
   end if

   ! Gamma = chi0^{-1} - chi^{-1}
   vertex%gamma_mat(iom,:,:) = chi0_inv(:,:) - chi_inv(:,:)

 end do

 ABI_FREE(chi0_inv)
 ABI_FREE(chi_inv)

 write(msg,'(a)') ' extract_vertex_irr: Irreducible vertex extraction complete'
 call wrtout(std_out, msg)

end subroutine extract_vertex_irr

!!***

!!****f* m_dmft_vertex/write_vertex_irr
!! NAME
!!  write_vertex_irr
!!
!! FUNCTION
!!  Write the irreducible vertex to file for diagnostics.
!!
!! SOURCE

subroutine write_vertex_irr(vertex, fname)

 type(vertex_irr_type), intent(in) :: vertex
 character(len=*), intent(in) :: fname

!Local variables
 integer :: unt, iom, ii, jj, ios
 character(len=500) :: msg

! *********************************************************************

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# Local irreducible vertex Gamma_imp'
 write(unt,'(a,i6)') '# norb_corr = ', vertex%norb_corr
 write(unt,'(a,i6)') '# niw_vertex = ', vertex%niw_vertex
 write(unt,'(a,i6)') '# nboson = ', vertex%nboson
 write(unt,'(a,i8)') '# ndim_comp = ', vertex%ndim_comp
 write(unt,'(a)') '# iOm  I  J  Re(Gamma)  Im(Gamma)'

 do iom = 1, vertex%nboson
   do ii = 1, min(vertex%ndim_comp, 10)
     do jj = 1, min(vertex%ndim_comp, 10)
       write(unt,'(3i6,2es20.10)') iom, ii, jj, &
       real(vertex%gamma_mat(iom,ii,jj)), aimag(vertex%gamma_mat(iom,ii,jj))
     end do
   end do
 end do

 close(unt)

end subroutine write_vertex_irr

!!***

END MODULE m_dmft_vertex
