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
 use m_green, only : green_type
 use m_dmft_vertex, only : vertex_irr_type
 use m_dmft_spinor_proj, only : spinor_proj_type

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
!!      = -beta * (1/Nk) sum_k G^loc_{delta,alpha}(k, iw_n)
!!                            * G^loc_{beta,gamma}(k, iw_n + iOm_m)
!!
!!  where G^loc_{alpha,beta}(k, iw) is the lattice Green function
!!  projected onto the correlated subspace:
!!
!!    G^loc_{alpha,beta}(k, iw) = sum_{a,b} P_{alpha,a}(k) G^KS_{a,b}(k, iw)
!!                                          P*_{beta,b}(k)
!!
!!  P_{alpha,a}(k) = chipsi(alpha, a, k) is the spinor projector.
!!  G^KS_{a,b}(k, iw) = green_imp%oper(iw)%ks(a, b, ik, 1) is the
!!  lattice Green function in the KS band basis (including self-energy).
!!
!!  The result is stored in the composite-index matrix representation.
!!  The -beta factor matches the impurity bubble convention (Rohringer et al.),
!!  while 1/Nk comes from k-point averaging via wtk weights.
!!
!! INPUTS
!!  green_imp = converged Green function with KS-basis data at all k-points
!!  paw_dmft = DFT+DMFT data with system parameters
!!  sproj = spinor projectors populated from chipsi
!!  norb_corr = number of correlated spinor-orbitals
!!  niw_vertex = number of fermionic Matsubara frequencies
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SIDE EFFECTS
!!  lbse%chi0_latt = on output, contains the lattice bubble
!!
!! NOTES
!!  Requires green_imp%oper(iw)%has_operks = 1 for all needed frequencies.
!!  If KS data is not available, the routine reports an error.
!!
!!  The k-point weights paw_dmft%wtk(ik) are used for the averaging,
!!  satisfying sum_k wtk(ik) = 1.
!!
!! SOURCE

subroutine compute_chi0_lattice(lbse, green_imp, paw_dmft, sproj, &
  & norb_corr, niw_vertex, nboson)

 type(lattice_bse_type), intent(inout) :: lbse
 type(green_type), intent(in) :: green_imp
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(spinor_proj_type), intent(in) :: sproj
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 integer :: iom, iw, iw_shifted, ik
 integer :: ialpha, ibeta, igamma, idelta
 integer :: ia, ib, norb_sq, idx_i, idx_j
 integer :: mbandc
 real(dp) :: beta, wk
 complex(dp), allocatable :: gloc_n(:,:), gloc_np(:,:)
 complex(dp), allocatable :: temp_mat(:,:)
 character(len=500) :: msg

! *********************************************************************

 write(msg,'(a)') ' compute_chi0_lattice: Computing lattice bubble at q=0'
 call wrtout(std_out, msg)

 ! --- Validate KS data availability ---
 if (green_imp%oper(1)%has_operks /= 1) then
   write(msg,'(3a)') &
   'compute_chi0_lattice: Green function does not have KS basis data.',ch10,&
   'Cannot compute lattice bubble without G(k,iw). Leaving chi0_latt as zero.'
   ABI_WARNING(msg)
   lbse%chi0_latt = czero
   return
 end if

 ! --- Validate frequency bounds ---
 if (niw_vertex + nboson - 1 > green_imp%nw) then
   write(msg,'(a,i6,a,i6,a,i6)') &
   'compute_chi0_lattice: frequency bounds exceeded. niw_vertex=', niw_vertex, &
   ' nboson=', nboson, ' but green has nw=', green_imp%nw
   ABI_ERROR(msg)
 end if

 beta = one / paw_dmft%temp
 mbandc = sproj%nband_ks
 norb_sq = norb_corr * norb_corr

 write(msg,'(a,es14.6,a,i4,a,i4,a,i6)') &
 ' compute_chi0_lattice: beta=', beta, ' niw_vertex=', niw_vertex, &
 ' nboson=', nboson, ' nkpt=', lbse%nkpt
 call wrtout(std_out, msg)

 ABI_MALLOC(gloc_n, (norb_corr, norb_corr))
 ABI_MALLOC(gloc_np, (norb_corr, norb_corr))
 ABI_MALLOC(temp_mat, (norb_corr, mbandc))

 lbse%chi0_latt = czero

 ! --- Compute lattice bubble ---
 ! chi0^latt_{(ab,n),(gd,n')}(iOm) = -beta * delta_{nn'} *
 !   (1/Nk) sum_k G^loc_{da}(k,iwn) * G^loc_{bg}(k,iwn+iOm)
 !
 ! G^loc_{ab}(k,iw) = sum_{c,d} P_{ac}(k) G^KS_{cd}(k,iw) P*_{bd}(k)
 !   = [ chipsi * G^KS * chipsi^dagger ]_{ab}
 !
 ! Computed using two matrix multiplies:
 !   temp = chipsi(k) * G^KS(k,iw)    [norb_corr x mbandc]
 !   G^loc = temp * chipsi(k)^dagger   [norb_corr x norb_corr]

 do iom = 1, nboson
   do iw = 1, niw_vertex
     iw_shifted = iw + (iom - 1)

     do ik = 1, lbse%nkpt
       wk = paw_dmft%wtk(ik)

       ! --- Project G^KS(k, iw_n) to local basis ---
       ! temp = chipsi * G^KS  (norb_corr x mbandc)
       do ialpha = 1, norb_corr
         do ib = 1, mbandc
           temp_mat(ialpha, ib) = czero
           do ia = 1, mbandc
             temp_mat(ialpha, ib) = temp_mat(ialpha, ib) + &
               sproj%proj(ialpha, ia, ik) * green_imp%oper(iw)%ks(ia, ib, ik, 1)
           end do
         end do
       end do
       ! G^loc = temp * chipsi^dagger  (norb_corr x norb_corr)
       do ialpha = 1, norb_corr
         do ibeta = 1, norb_corr
           gloc_n(ialpha, ibeta) = czero
           do ib = 1, mbandc
             gloc_n(ialpha, ibeta) = gloc_n(ialpha, ibeta) + &
               temp_mat(ialpha, ib) * conjg(sproj%proj(ibeta, ib, ik))
           end do
         end do
       end do

       ! --- Project G^KS(k, iw_n + iOm_m) to local basis ---
       do ialpha = 1, norb_corr
         do ib = 1, mbandc
           temp_mat(ialpha, ib) = czero
           do ia = 1, mbandc
             temp_mat(ialpha, ib) = temp_mat(ialpha, ib) + &
               sproj%proj(ialpha, ia, ik) * green_imp%oper(iw_shifted)%ks(ia, ib, ik, 1)
           end do
         end do
       end do
       do ialpha = 1, norb_corr
         do ibeta = 1, norb_corr
           gloc_np(ialpha, ibeta) = czero
           do ib = 1, mbandc
             gloc_np(ialpha, ibeta) = gloc_np(ialpha, ibeta) + &
               temp_mat(ialpha, ib) * conjg(sproj%proj(ibeta, ib, ik))
           end do
         end do
       end do

       ! --- Accumulate into lattice bubble ---
       ! chi0^latt += -beta * wk * G^loc_{da}(k,iw) * G^loc_{bg}(k,iw+iOm)
       do ialpha = 1, norb_corr
         do ibeta = 1, norb_corr
           idx_i = (iw - 1) * norb_sq + (ialpha - 1) * norb_corr + ibeta
           do igamma = 1, norb_corr
             do idelta = 1, norb_corr
               idx_j = (iw - 1) * norb_sq + (igamma - 1) * norb_corr + idelta
               lbse%chi0_latt(iom, idx_i, idx_j) = &
                 lbse%chi0_latt(iom, idx_i, idx_j) &
                 - beta * wk * gloc_n(idelta, ialpha) * gloc_np(ibeta, igamma)
             end do
           end do
         end do
       end do

     end do  ! ik
   end do  ! iw
 end do  ! iom

 ABI_FREE(gloc_n)
 ABI_FREE(gloc_np)
 ABI_FREE(temp_mat)

 write(msg,'(a,i8,a)') &
 ' compute_chi0_lattice: Lattice bubble computed (', lbse%ndim_comp, ' composite indices)'
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
