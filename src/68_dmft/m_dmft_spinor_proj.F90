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
