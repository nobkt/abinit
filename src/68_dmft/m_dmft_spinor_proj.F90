!!****m* ABINIT/m_dmft_spinor_proj
!! NAME
!!  m_dmft_spinor_proj
!!
!! FUNCTION
!!  Manages spinor projectors for the correlated subspace in DFT+DMFT
!!  response calculations.
!!
!!  The projector P_{alpha, a*sigma}(k) maps from the full Kohn-Sham
!!  spinor basis {a,sigma} to the correlated subspace {alpha}, where
!!  alpha = (m, sigma) is a composite orbital-spin index.
!!
!!  For spin-flip absorption calculations, the projector must preserve
!!  the full spinor structure without collapsing spin indices.
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

MODULE m_dmft_spinor_proj

 use defs_basis
 use m_abicore
 use m_errors

 use m_paw_dmft, only : paw_dmft_type

 implicit none

 private

 public :: spinor_proj_type
 public :: init_spinor_proj
 public :: destroy_spinor_proj
 public :: populate_from_chipsi
 public :: check_spinor_completeness

!!***

!!****t* m_dmft_spinor_proj/spinor_proj_type
!! NAME
!!  spinor_proj_type
!!
!! FUNCTION
!!  Stores spinor projectors for the correlated subspace.
!!
!!  proj(norb_corr, nband_ks, nkpt) where:
!!    norb_corr = (2*lpawu+1) * nspinor (composite spinor-orbital index)
!!    nband_ks = number of Kohn-Sham bands included
!!    nkpt = number of k-points
!!
!!  The projector includes both spin-up and spin-down components
!!  for each correlated orbital m, preserving off-diagonal spin mixing.
!!
!! SOURCE

 type, public :: spinor_proj_type

   integer :: norb_corr = 0
   integer :: nband_ks = 0
   integer :: nkpt = 0
   integer :: nspinor = 0
   integer :: maxlpawu = 0

   complex(dp), allocatable :: proj(:,:,:)
   ! proj(norb_corr, nband_ks, nkpt)

 end type spinor_proj_type

!!***

CONTAINS

!!****f* m_dmft_spinor_proj/init_spinor_proj
!! NAME
!!  init_spinor_proj
!!
!! FUNCTION
!!  Initialize spinor projectors from the DFT+DMFT data.
!!  Allocates the projector array and sets dimensions from paw_dmft.
!!  Actual population of projector values from chipsi is pending
!!  (requires chipsi-to-composite-index mapping not yet implemented).
!!
!! INPUTS
!!  paw_dmft = DFT+DMFT data with chipsi projections
!!
!! SIDE EFFECTS
!!  sproj = on output, contains initialized spinor projectors
!!
!! NOTES
!!  The existing chipsi in paw_dmft has dimensions:
!!    chipsi(nspinor*(2*maxlpawu+1), mbandc, nkpt, nsppol, natom)
!!  This routine reorganizes into the composite index alpha=(m,sigma).
!!
!! SOURCE

subroutine init_spinor_proj(sproj, paw_dmft)

 type(spinor_proj_type), intent(inout) :: sproj
 type(paw_dmft_type), intent(in) :: paw_dmft

!Local variables
 integer :: norb_corr, nband_ks, nkpt, nspinor, maxlpawu
 character(len=500) :: msg

! *********************************************************************

 nspinor = paw_dmft%nspinor
 maxlpawu = paw_dmft%maxlpawu
 norb_corr = nspinor * (2 * maxlpawu + 1)
 nband_ks = paw_dmft%mbandc
 nkpt = paw_dmft%nkpt

 sproj%norb_corr = norb_corr
 sproj%nband_ks = nband_ks
 sproj%nkpt = nkpt
 sproj%nspinor = nspinor
 sproj%maxlpawu = maxlpawu

 write(msg,'(a,i4,a,i4,a,i6,a,i4)') &
 ' init_spinor_proj: norb_corr=', norb_corr, ' nband_ks=', nband_ks, &
 ' nkpt=', nkpt, ' nspinor=', nspinor
 call wrtout(std_out, msg)

 ABI_MALLOC(sproj%proj, (norb_corr, nband_ks, nkpt))
 sproj%proj = czero

 ! NOTE: The actual population of projectors from chipsi requires
 ! careful index mapping between the paw_dmft chipsi layout and
 ! the composite spinor-orbital index alpha=(m,sigma).
 ! This will be connected in the next implementation phase.

 write(msg,'(a)') ' init_spinor_proj: Projector structure initialized'
 call wrtout(std_out, msg)

end subroutine init_spinor_proj

!!***

!!****f* m_dmft_spinor_proj/populate_from_chipsi
!! NAME
!!  populate_from_chipsi
!!
!! FUNCTION
!!  Populate the spinor projector array from paw_dmft%chipsi.
!!
!!  chipsi has dimensions:
!!    chipsi(nspinor*(2*maxlpawu+1), mbandc, nkpt, nsppol, natom)
!!
!!  For nspinor=2, nsppol=1 (required for spin-flip calculations),
!!  the first dimension of chipsi already contains the composite
!!  spinor-orbital index alpha = (m, sigma), with the same ordering
!!  as used by matlu%mat.
!!
!!  This routine copies chipsi(:,:,:,1,iatom) into sproj%proj(:,:,:)
!!  for a single correlated atom.
!!
!! INPUTS
!!  paw_dmft = DFT+DMFT data with chipsi projections
!!  iatom = atom index for which to populate projectors
!!
!! SIDE EFFECTS
!!  sproj = on output, proj array contains the chipsi data
!!
!! NOTES
!!  The atom must have lpawu >= 0 (correlated atom).
!!  For nspinor=2, nsppol must be 1.
!!  The ordering of the composite index in chipsi is assumed to match
!!  the matlu convention: the same index ordering that defines
!!  G^loc in green%oper(iw)%matlu(iatom)%mat.
!!
!! SOURCE

subroutine populate_from_chipsi(sproj, paw_dmft, iatom)

 type(spinor_proj_type), intent(inout) :: sproj
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: iatom

!Local variables
 integer :: norb_corr_atom, ialpha, ib, ik
 character(len=500) :: msg

! *********************************************************************

 if (iatom < 1 .or. iatom > paw_dmft%natom) then
   write(msg,'(a,i4,a,i4)') &
   'populate_from_chipsi: iatom=', iatom, ' out of range [1,', paw_dmft%natom
   ABI_ERROR(trim(msg)//']')
 end if

 if (paw_dmft%lpawu(iatom) < 0) then
   write(msg,'(a,i4,a)') &
   'populate_from_chipsi: atom ', iatom, ' is not correlated (lpawu < 0)'
   ABI_ERROR(msg)
 end if

 if (paw_dmft%nspinor == 2 .and. paw_dmft%nsppol /= 1) then
   ABI_ERROR('populate_from_chipsi: nspinor=2 requires nsppol=1')
 end if

 ! Validate dimension consistency
 norb_corr_atom = paw_dmft%nspinor * (2 * paw_dmft%lpawu(iatom) + 1)
 if (norb_corr_atom /= sproj%norb_corr) then
   write(msg,'(a,i4,a,i4)') &
   'populate_from_chipsi: dimension mismatch. atom norb_corr=', norb_corr_atom, &
   ' but sproj%norb_corr=', sproj%norb_corr
   ABI_ERROR(msg)
 end if

 if (.not. allocated(paw_dmft%chipsi)) then
   ABI_ERROR('populate_from_chipsi: paw_dmft%chipsi is not allocated')
 end if

 write(msg,'(a,i4,a,i4,a,i4)') &
 ' populate_from_chipsi: atom=', iatom, ' norb_corr=', norb_corr_atom, &
 ' nband_ks=', sproj%nband_ks
 call wrtout(std_out, msg)

 ! Copy chipsi to sproj%proj
 ! chipsi(alpha, band, kpt, isppol, iatom) -> proj(alpha, band, kpt)
 ! For nspinor=2, nsppol=1: isppol index is always 1
 do ik = 1, sproj%nkpt
   do ib = 1, sproj%nband_ks
     do ialpha = 1, sproj%norb_corr
       sproj%proj(ialpha, ib, ik) = paw_dmft%chipsi(ialpha, ib, ik, 1, iatom)
     end do
   end do
 end do

 write(msg,'(a,i4)') &
 ' populate_from_chipsi: Projectors populated for atom ', iatom
 call wrtout(std_out, msg)

end subroutine populate_from_chipsi

!!***

!!****f* m_dmft_spinor_proj/destroy_spinor_proj
!! NAME
!!  destroy_spinor_proj
!!
!! FUNCTION
!!  Deallocate spinor projector structure.
!!
!! SOURCE

subroutine destroy_spinor_proj(sproj)

 type(spinor_proj_type), intent(inout) :: sproj

! *********************************************************************

 if (allocated(sproj%proj)) then
   ABI_FREE(sproj%proj)
 end if

 sproj%norb_corr = 0
 sproj%nband_ks = 0
 sproj%nkpt = 0
 sproj%nspinor = 0
 sproj%maxlpawu = 0

end subroutine destroy_spinor_proj

!!***

!!****f* m_dmft_spinor_proj/check_spinor_completeness
!! NAME
!!  check_spinor_completeness
!!
!! FUNCTION
!!  Verify that the spinor projectors satisfy the completeness relation
!!  within the correlated subspace:
!!    sum_a P_{alpha,a}(k) * P^*_{beta,a}(k) approx delta_{alpha,beta}
!!
!!  This is a diagnostic check, not a physical approximation.
!!
!! INPUTS
!!  sproj = spinor projectors
!!  tol = tolerance for deviation from identity
!!
!! SOURCE

subroutine check_spinor_completeness(sproj, tol)

 type(spinor_proj_type), intent(in) :: sproj
 real(dp), intent(in) :: tol

!Local variables
 integer :: ik, ialpha, ibeta, ib
 real(dp) :: max_offdiag, max_diag_dev
 complex(dp) :: overlap
 character(len=500) :: msg

! *********************************************************************

 max_offdiag = zero
 max_diag_dev = zero

 do ik = 1, sproj%nkpt
   do ialpha = 1, sproj%norb_corr
     do ibeta = 1, sproj%norb_corr
       overlap = czero
       do ib = 1, sproj%nband_ks
         overlap = overlap + sproj%proj(ialpha, ib, ik) * conjg(sproj%proj(ibeta, ib, ik))
       end do
       if (ialpha == ibeta) then
         max_diag_dev = max(max_diag_dev, abs(overlap - cone))
       else
         max_offdiag = max(max_offdiag, abs(overlap))
       end if
     end do
   end do
 end do

 write(msg,'(a,es10.3,a,es10.3)') &
 ' check_spinor_completeness: max diagonal deviation=', max_diag_dev, &
 ' max off-diagonal=', max_offdiag
 call wrtout(std_out, msg)

 if (max_diag_dev > tol .or. max_offdiag > tol) then
   write(msg,'(a,es10.3,a)') &
   ' WARNING: Spinor projector completeness check failed (tol=', tol, ')'
   call wrtout(std_out, msg)
 end if

end subroutine check_spinor_completeness

!!***

END MODULE m_dmft_spinor_proj
