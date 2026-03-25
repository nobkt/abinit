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
 use m_abi_linalg, only : abi_zgemm_2dd
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
 public :: kpoint_chi0_attrib_type
 public :: init_kpoint_chi0_attrib
 public :: destroy_kpoint_chi0_attrib
 public :: write_kpoint_chi0_attrib

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

!!****t* m_dmft_lattice_bse/kpoint_chi0_attrib_type
!! NAME
!!  kpoint_chi0_attrib_type
!!
!! FUNCTION
!!  Stores k-point resolved spin-channel traces of the lattice bubble.
!!
!!  For each k-point and bosonic frequency, the susceptibility trace
!!  is decomposed into spin channels:
!!    chi0_k_total(iOm, ik) = -beta * wtk * sum_{n,a,b}
!!                             G^loc_{ba}(k,iwn) G^loc_{ab}(k,iwn+iOm)
!!    chi0_k_pm, chi0_k_mp, chi0_k_sc : spin-channel components
!!
!!  This enables identification of which k-points in the Brillouin zone
!!  contribute most to spin-flip susceptibility.
!!
!! SOURCE

 type, public :: kpoint_chi0_attrib_type

   integer :: nboson = 0
   integer :: nkpt = 0
   integer :: ndim_orb = 0
   integer :: nspinor = 0

   complex(dp), allocatable :: chi0_k_total(:,:)
   ! chi0_k_total(nboson, nkpt) : total trace per k-point

   complex(dp), allocatable :: chi0_k_sc(:,:)
   ! chi0_k_sc(nboson, nkpt) : spin-conserving trace per k-point

   complex(dp), allocatable :: chi0_k_pm(:,:)
   ! chi0_k_pm(nboson, nkpt) : S+S- trace per k-point

   complex(dp), allocatable :: chi0_k_mp(:,:)
   ! chi0_k_mp(nboson, nkpt) : S-S+ trace per k-point

   complex(dp), allocatable :: chi0_k_pm_orb(:,:,:,:)
   ! chi0_k_pm_orb(nboson, nkpt, ndim_orb, ndim_orb) :
   ! Per-k orbital-resolved S+S- susceptibility.
   ! chi0_k_pm_orb(iOm, ik, m, m') =
   !   sum_n (-beta * wk) G^loc_{(m',down),(m,up)}(k,iwn) G^loc_{(m,up),(m',down)}(k,iwn+iOm)

   complex(dp), allocatable :: chi0_k_mp_orb(:,:,:,:)
   ! chi0_k_mp_orb(nboson, nkpt, ndim_orb, ndim_orb) :
   ! Per-k orbital-resolved S-S+ susceptibility.

 end type kpoint_chi0_attrib_type

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
!!  When kpt_attrib is present and nspinor=2, per-k spin-channel trace
!!  decomposition is computed on-the-fly during the k-point loop.
!!  The spin index convention is: indices 1..ndim_orb are spin-up,
!!  indices ndim_orb+1..2*ndim_orb are spin-down, where ndim_orb=norb_corr/nspinor.
!!  This matches the convention in m_dmft_spectral_attribution.
!!  If nspinor/=2, kpt_attrib is ignored (spin-channel decomposition is undefined).
!!
!! SOURCE

subroutine compute_chi0_lattice(lbse, green_imp, paw_dmft, sproj, &
  & norb_corr, niw_vertex, nboson, ladd, kpt_attrib)

 type(lattice_bse_type), intent(inout) :: lbse
 type(green_type), intent(in) :: green_imp
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(spinor_proj_type), intent(in) :: sproj
 integer, intent(in) :: norb_corr, niw_vertex, nboson
 logical, intent(in), optional :: ladd
   !! When .true., accumulate into chi0_latt instead of zeroing it first.
   !! Default is .false. (zero and compute from scratch).
 type(kpoint_chi0_attrib_type), intent(inout), optional :: kpt_attrib
   !! When present, per-k spin-channel traces are computed and accumulated.
   !! Must be initialized before calling this routine.
   !! If ladd=.true., existing values are preserved and accumulated into.

!Local variables
 integer :: iom, iw, iw_shifted, ik, iw_max
 integer :: ialpha, ibeta, igamma, idelta
 integer :: norb_sq, idx_i, idx_j
 integer :: mbandc
 integer :: ispin_a, ispin_b, ndim_orb, im_a, im_b
 logical :: ladd_local, do_kpt_attrib
 real(dp) :: beta, wk
 complex(dp) :: bubble_contrib
 complex(dp), allocatable :: gloc_n(:,:), gloc_np(:,:)
 complex(dp), allocatable :: temp_mat(:,:), proj_k(:,:)
 character(len=500) :: msg

! *********************************************************************

 ladd_local = .false.
 if (present(ladd)) ladd_local = ladd

 do_kpt_attrib = present(kpt_attrib)
 ndim_orb = norb_corr / paw_dmft%nspinor

 write(msg,'(a,l2,a,l2)') ' compute_chi0_lattice: Computing lattice bubble at q=0, accumulate=', &
   ladd_local, ' kpt_attrib=', do_kpt_attrib
 call wrtout(std_out, msg)

 ! --- Validate KS data availability across all needed frequencies ---
 iw_max = niw_vertex + nboson - 1
 if (green_imp%oper(1)%has_operks /= 1) then
   write(msg,'(3a)') &
   'compute_chi0_lattice: Green function does not have KS basis data.',ch10,&
   'Cannot compute lattice bubble without G(k,iw). Skipping this atom contribution.'
   ABI_WARNING(msg)
   if (.not. ladd_local) lbse%chi0_latt = czero
   return
 end if
 ! Verify KS data exists for the highest frequency needed
 if (green_imp%oper(iw_max)%has_operks /= 1) then
   write(msg,'(a,i6,a)') &
   'compute_chi0_lattice: KS data not available at frequency index ', iw_max, &
   '. Skipping this atom contribution.'
   ABI_WARNING(msg)
   if (.not. ladd_local) lbse%chi0_latt = czero
   return
 end if

 ! --- Validate frequency bounds ---
 if (iw_max > green_imp%nw) then
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
 ABI_MALLOC(proj_k, (norb_corr, mbandc))

 if (.not. ladd_local) then
   lbse%chi0_latt = czero
   if (do_kpt_attrib) then
     kpt_attrib%chi0_k_total = czero
     kpt_attrib%chi0_k_sc = czero
     kpt_attrib%chi0_k_pm = czero
     kpt_attrib%chi0_k_mp = czero
     kpt_attrib%chi0_k_pm_orb = czero
     kpt_attrib%chi0_k_mp_orb = czero
   end if
 end if

 ! --- Compute lattice bubble ---
 ! chi0^latt_{(ab,n),(gd,n')}(iOm) = -beta * delta_{nn'} *
 !   sum_k wk * G^loc_{da}(k,iwn) * G^loc_{bg}(k,iwn+iOm)
 !
 ! G^loc(k,iw) = P(k) * G^KS(k,iw) * P(k)^dagger
 !
 ! Computed using BLAS (abi_zgemm_2dd):
 !   temp = P(k) * G^KS(k,iw)          [norb_corr x mbandc]
 !   G^loc = temp * P(k)^dagger         [norb_corr x norb_corr]

 do iom = 1, nboson
   do iw = 1, niw_vertex
     iw_shifted = iw + (iom - 1)

     do ik = 1, lbse%nkpt
       wk = paw_dmft%wtk(ik)

       ! Extract projector for this k-point
       proj_k(:,:) = sproj%proj(:,:,ik)

       ! --- Project G^KS(k, iw_n) to local basis using BLAS ---
       ! temp = P(k) * G^KS(k,iw)  [norb_corr x mbandc]
       call abi_zgemm_2dd("N", "N", norb_corr, mbandc, mbandc, cone, &
         proj_k, norb_corr, &
         green_imp%oper(iw)%ks(1, 1, ik, 1), mbandc, &
         czero, temp_mat, norb_corr)
       ! G^loc_n = temp * P(k)^dagger  [norb_corr x norb_corr]
       call abi_zgemm_2dd("N", "C", norb_corr, norb_corr, mbandc, cone, &
         temp_mat, norb_corr, &
         proj_k, norb_corr, &
         czero, gloc_n, norb_corr)

       ! --- Project G^KS(k, iw_n + iOm_m) to local basis using BLAS ---
       call abi_zgemm_2dd("N", "N", norb_corr, mbandc, mbandc, cone, &
         proj_k, norb_corr, &
         green_imp%oper(iw_shifted)%ks(1, 1, ik, 1), mbandc, &
         czero, temp_mat, norb_corr)
       call abi_zgemm_2dd("N", "C", norb_corr, norb_corr, mbandc, cone, &
         temp_mat, norb_corr, &
         proj_k, norb_corr, &
         czero, gloc_np, norb_corr)

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


        ! --- Per-k spin-channel trace attribution ---
        ! Trace chi0_{ab,ba} = -beta * wk * G_{ba}(k,iw) * G_{ab}(k,iw+iOm)
        ! Classified by spin channel of (alpha, beta).
        if (do_kpt_attrib .and. paw_dmft%nspinor == 2) then
          do ialpha = 1, norb_corr
            do ibeta = 1, norb_corr
              bubble_contrib = -beta * wk * gloc_n(ibeta, ialpha) * gloc_np(ialpha, ibeta)

              ! Total trace
              kpt_attrib%chi0_k_total(iom, ik) = &
                kpt_attrib%chi0_k_total(iom, ik) + bubble_contrib

              ! Determine spin indices from composite spinor-orbital index.
              ! Convention: alpha=1..ndim_orb -> spin up (ispin=1)
              !             alpha=ndim_orb+1..2*ndim_orb -> spin down (ispin=2)
              ! This matches m_dmft_spectral_attribution (compute_spectral_attribution).
              if (ialpha <= ndim_orb) then
                ispin_a = 1
              else
                ispin_a = 2
              end if
              if (ibeta <= ndim_orb) then
                ispin_b = 1
              else
                ispin_b = 2
              end if

              ! Classify by spin channel
              if (ispin_a == ispin_b) then
                kpt_attrib%chi0_k_sc(iom, ik) = &
                  kpt_attrib%chi0_k_sc(iom, ik) + bubble_contrib
              else if (ispin_a == 1 .and. ispin_b == 2) then
                kpt_attrib%chi0_k_pm(iom, ik) = &
                  kpt_attrib%chi0_k_pm(iom, ik) + bubble_contrib
                ! Orbital-resolved S+S-: alpha=(m,up), beta=(m',down)
                im_a = ialpha              ! orbital index for spin-up part
                im_b = ibeta - ndim_orb    ! orbital index for spin-down part
                kpt_attrib%chi0_k_pm_orb(iom, ik, im_a, im_b) = &
                  kpt_attrib%chi0_k_pm_orb(iom, ik, im_a, im_b) + bubble_contrib
              else if (ispin_a == 2 .and. ispin_b == 1) then
                kpt_attrib%chi0_k_mp(iom, ik) = &
                  kpt_attrib%chi0_k_mp(iom, ik) + bubble_contrib
                ! Orbital-resolved S-S+: alpha=(m,down), beta=(m',up)
                im_a = ialpha - ndim_orb   ! orbital index for spin-down part
                im_b = ibeta               ! orbital index for spin-up part
                kpt_attrib%chi0_k_mp_orb(iom, ik, im_a, im_b) = &
                  kpt_attrib%chi0_k_mp_orb(iom, ik, im_a, im_b) + bubble_contrib
              end if
            end do
          end do
        end if
     end do  ! ik
   end do  ! iw
 end do  ! iom

 ABI_FREE(gloc_n)
 ABI_FREE(gloc_np)
 ABI_FREE(temp_mat)
 ABI_FREE(proj_k)

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

!!****f* m_dmft_lattice_bse/init_kpoint_chi0_attrib
!! NAME
!!  init_kpoint_chi0_attrib
!!
!! FUNCTION
!!  Initialize the k-point resolved attribution structure.
!!
!! INPUTS
!!  nboson = number of bosonic Matsubara frequencies
!!  nkpt = number of k-points
!!  ndim_orb = number of orbital indices (2*lpawu+1)
!!  nspinor = number of spinor components
!!
!! SOURCE

subroutine init_kpoint_chi0_attrib(kattr, nboson, nkpt, ndim_orb, nspinor)

 type(kpoint_chi0_attrib_type), intent(inout) :: kattr
 integer, intent(in) :: nboson, nkpt, ndim_orb, nspinor

!Local variables
 character(len=500) :: msg

! *********************************************************************

 kattr%nboson = nboson
 kattr%nkpt = nkpt
 kattr%ndim_orb = ndim_orb
 kattr%nspinor = nspinor

 write(msg,'(a,i4,a,i6,a,i4,a,i2)') &
 ' init_kpoint_chi0_attrib: nboson=', nboson, ' nkpt=', nkpt, &
 ' ndim_orb=', ndim_orb, ' nspinor=', nspinor
 call wrtout(std_out, msg)

 ABI_MALLOC(kattr%chi0_k_total, (nboson, nkpt))
 ABI_MALLOC(kattr%chi0_k_sc, (nboson, nkpt))
 ABI_MALLOC(kattr%chi0_k_pm, (nboson, nkpt))
 ABI_MALLOC(kattr%chi0_k_mp, (nboson, nkpt))
 ABI_MALLOC(kattr%chi0_k_pm_orb, (nboson, nkpt, ndim_orb, ndim_orb))
 ABI_MALLOC(kattr%chi0_k_mp_orb, (nboson, nkpt, ndim_orb, ndim_orb))

 kattr%chi0_k_total = czero
 kattr%chi0_k_sc = czero
 kattr%chi0_k_pm = czero
 kattr%chi0_k_mp = czero
 kattr%chi0_k_pm_orb = czero
 kattr%chi0_k_mp_orb = czero

end subroutine init_kpoint_chi0_attrib

!!***

!!****f* m_dmft_lattice_bse/destroy_kpoint_chi0_attrib
!! NAME
!!  destroy_kpoint_chi0_attrib
!!
!! FUNCTION
!!  Deallocate the k-point resolved attribution structure.
!!
!! SOURCE

subroutine destroy_kpoint_chi0_attrib(kattr)

 type(kpoint_chi0_attrib_type), intent(inout) :: kattr

! *********************************************************************

 if (allocated(kattr%chi0_k_total)) then
   ABI_FREE(kattr%chi0_k_total)
 end if
 if (allocated(kattr%chi0_k_sc)) then
   ABI_FREE(kattr%chi0_k_sc)
 end if
 if (allocated(kattr%chi0_k_pm)) then
   ABI_FREE(kattr%chi0_k_pm)
 end if
 if (allocated(kattr%chi0_k_mp)) then
   ABI_FREE(kattr%chi0_k_mp)
 end if
 if (allocated(kattr%chi0_k_pm_orb)) then
   ABI_FREE(kattr%chi0_k_pm_orb)
 end if
 if (allocated(kattr%chi0_k_mp_orb)) then
   ABI_FREE(kattr%chi0_k_mp_orb)
 end if

 kattr%nboson = 0
 kattr%nkpt = 0
 kattr%ndim_orb = 0
 kattr%nspinor = 0

end subroutine destroy_kpoint_chi0_attrib

!!***

!!****f* m_dmft_lattice_bse/write_kpoint_chi0_attrib
!! NAME
!!  write_kpoint_chi0_attrib
!!
!! FUNCTION
!!  Write k-point resolved lattice bubble attribution to file.
!!
!!  For each bosonic frequency, outputs the spin-channel trace
!!  contribution from each k-point. This enables identification of
!!  which regions of the Brillouin zone dominate the spin-flip
!!  susceptibility.
!!
!! INPUTS
!!  kattr = k-point resolved attribution data
!!  paw_dmft = DFT+DMFT data (for k-point coordinates and weights)
!!  fname = output file name
!!  beta = inverse temperature
!!
!! SOURCE

subroutine write_kpoint_chi0_attrib(kattr, paw_dmft, kpt_coords, fname, beta)

 type(kpoint_chi0_attrib_type), intent(in) :: kattr
 type(paw_dmft_type), intent(in) :: paw_dmft
 real(dp), intent(in) :: kpt_coords(:,:)
   !! kpt_coords(3, nkpt) : k-point coordinates in reduced coordinates
 character(len=*), intent(in) :: fname
 real(dp), intent(in) :: beta

!Local variables
 integer :: unt, iom, ik, ios, im, imp
 integer :: nentries, nrank, ientry, ii, jj, tmp_int
 integer, allocatable :: rank_ik(:), rank_im(:), rank_imp(:)
 real(dp) :: omega_boson, tmp_val
 real(dp), allocatable :: rank_val(:)
 character(len=500) :: msg

! *********************************************************************

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT k-point Resolved Lattice Bubble Attribution'
 write(unt,'(a)') '# Spin-channel trace per k-point: chi0_{ab,ba}(iOm) contributions'
 write(unt,'(a,i6)') '# nboson = ', kattr%nboson
 write(unt,'(a,i6)') '# nkpt = ', kattr%nkpt
 write(unt,'(a,i4)') '# ndim_orb = ', kattr%ndim_orb
 write(unt,'(a,es14.6)') '# beta = ', beta
 write(unt,'(a)')

 ! --- Section 1: Per-k spin channel traces ---
 write(unt,'(a)') '# === Section 1: Per-k spin channel traces ==='
 write(unt,'(a)') '# iOm  ik  wtk  kx  ky  kz  Re(total)  Im(total)  ' // &
   'Re(sc)  Im(sc)  Re(pm)  Im(pm)  Re(mp)  Im(mp)'

 do iom = 1, kattr%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   do ik = 1, kattr%nkpt
     write(unt,'(2i6,es12.4,3f10.5,8es18.8)') &
       iom, ik, paw_dmft%wtk(ik), &
       kpt_coords(1,ik), kpt_coords(2,ik), kpt_coords(3,ik), &
       real(kattr%chi0_k_total(iom, ik)), aimag(kattr%chi0_k_total(iom, ik)), &
       real(kattr%chi0_k_sc(iom, ik)), aimag(kattr%chi0_k_sc(iom, ik)), &
       real(kattr%chi0_k_pm(iom, ik)), aimag(kattr%chi0_k_pm(iom, ik)), &
       real(kattr%chi0_k_mp(iom, ik)), aimag(kattr%chi0_k_mp(iom, ik))
   end do
 end do

 write(unt,'(a)')

 ! --- Section 2: Summary — dominant k-points at static limit (iOm=0) ---
 write(unt,'(a)') '# === Section 2: Dominant k-points for S+S- at iOm=0 ==='
 write(unt,'(a)') '# ik  wtk  kx  ky  kz  |chi0_k_pm|  Re(chi0_k_pm)  Im(chi0_k_pm)'

 do ik = 1, kattr%nkpt
   write(unt,'(i6,es12.4,3f10.5,es16.6,2es18.8)') &
     ik, paw_dmft%wtk(ik), &
     kpt_coords(1,ik), kpt_coords(2,ik), kpt_coords(3,ik), &
     abs(kattr%chi0_k_pm(1, ik)), &
     real(kattr%chi0_k_pm(1, ik)), aimag(kattr%chi0_k_pm(1, ik))
 end do

 ! --- Section 3: Per-k orbital-resolved S+S- at iOm=0 ---
 if (allocated(kattr%chi0_k_pm_orb)) then
   write(unt,'(a)')
   write(unt,'(a)') '# === Section 3: Per-k orbital-resolved S+S- at iOm=0 ==='
   write(unt,'(a)') '# For each k-point, the orbital-pair decomposition of the S+S-'
   write(unt,'(a)') '# spin-flip susceptibility at the static limit.'
   write(unt,'(a)') '# ik  kx  ky  kz  m  m_prime  Re(chi0_pm_orb)  Im(chi0_pm_orb)  |chi0_pm_orb|'

   do ik = 1, kattr%nkpt
     do im = 1, kattr%ndim_orb
       do imp = 1, kattr%ndim_orb
         write(unt,'(i6,3f10.5,2i4,2es18.8,es16.6)') &
           ik, kpt_coords(1,ik), kpt_coords(2,ik), kpt_coords(3,ik), &
           im, imp, &
           real(kattr%chi0_k_pm_orb(1, ik, im, imp)), &
           aimag(kattr%chi0_k_pm_orb(1, ik, im, imp)), &
           abs(kattr%chi0_k_pm_orb(1, ik, im, imp))
       end do
     end do
   end do
 end if

 ! --- Section 4: Global ranking of dominant (k, m, m') for S+S- at iOm=0 ---
 if (allocated(kattr%chi0_k_pm_orb)) then
   write(unt,'(a)')
   write(unt,'(a)') '# === Section 4: Global ranking of dominant (k, m, m_prime) for S+S- at iOm=0 ==='
   write(unt,'(a)') '# Identifies which k-point and orbital transition contributes most'
   write(unt,'(a)') '# to the spin-flip susceptibility.'
   write(unt,'(a)') '# Rank  ik  kx  ky  kz  m  m_prime  |chi0_pm_orb|  Re  Im'

   nentries = kattr%nkpt * kattr%ndim_orb * kattr%ndim_orb
   nrank = min(nentries, 20)

   ABI_MALLOC(rank_ik, (nentries))
   ABI_MALLOC(rank_im, (nentries))
   ABI_MALLOC(rank_imp, (nentries))
   ABI_MALLOC(rank_val, (nentries))

   ientry = 0
   do ik = 1, kattr%nkpt
     do im = 1, kattr%ndim_orb
       do imp = 1, kattr%ndim_orb
         ientry = ientry + 1
         rank_ik(ientry) = ik
         rank_im(ientry) = im
         rank_imp(ientry) = imp
         rank_val(ientry) = abs(kattr%chi0_k_pm_orb(1, ik, im, imp))
       end do
     end do
   end do

   ! Sort by magnitude in descending order (partial insertion sort for top nrank)
   do ii = 1, nrank
     do jj = ii + 1, nentries
       if (rank_val(jj) > rank_val(ii)) then
         tmp_val = rank_val(ii); rank_val(ii) = rank_val(jj); rank_val(jj) = tmp_val
         tmp_int = rank_ik(ii); rank_ik(ii) = rank_ik(jj); rank_ik(jj) = tmp_int
         tmp_int = rank_im(ii); rank_im(ii) = rank_im(jj); rank_im(jj) = tmp_int
         tmp_int = rank_imp(ii); rank_imp(ii) = rank_imp(jj); rank_imp(jj) = tmp_int
       end if
     end do
   end do

   do ientry = 1, nrank
     ik = rank_ik(ientry)
     write(unt,'(i6,i6,3f10.5,2i4,es16.6,2es18.8)') &
       ientry, ik, &
       kpt_coords(1,ik), kpt_coords(2,ik), kpt_coords(3,ik), &
       rank_im(ientry), rank_imp(ientry), rank_val(ientry), &
       real(kattr%chi0_k_pm_orb(1, ik, rank_im(ientry), rank_imp(ientry))), &
       aimag(kattr%chi0_k_pm_orb(1, ik, rank_im(ientry), rank_imp(ientry)))
   end do

   ABI_FREE(rank_ik)
   ABI_FREE(rank_im)
   ABI_FREE(rank_imp)
   ABI_FREE(rank_val)
 end if

 close(unt)

 write(msg,'(3a)') ' write_kpoint_chi0_attrib: Written to ', trim(fname)
 call wrtout(std_out, msg)

end subroutine write_kpoint_chi0_attrib

!!***

END MODULE m_dmft_lattice_bse
