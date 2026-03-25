!!****m* ABINIT/m_dmft_current_vertex
!! NAME
!!  m_dmft_current_vertex
!!
!! FUNCTION
!!  Computes momentum (velocity) matrix elements for the DFT+DMFT
!!  optical response on the Matsubara axis.
!!
!!  The momentum matrix elements <psi_a(k)| -i nabla_mu |psi_b(k)> are
!!  computed for all correlated bands (a,b) at each k-point.
!!  These are used for the bubble optical conductivity:
!!    Pi_mu_nu^bubble(iOm) = -(1/(beta*Nk)) sum_{k,n}
!!        Tr[ j_mu(k) G(k,iw_n) j_nu(k) G(k,iw_n+iOm) ]
!!
!!  The computation has two contributions:
!!    1. Kinetic part: sum_G (k+G)_mu c*_a(G) c_b(G)
!!    2. PAW augmentation: sum_{ij,atom} cprj*_i nabla_ij cprj_j
!!       where nabla_ij = <phi_i|nabla|phi_j> - <tphi_i|nabla|tphi_j>
!!
!!  The PAW augmentation is included when pawtab%nabla_ij is available.
!!  If not available, only the kinetic part is computed and a warning
!!  is issued. This is documented, not a heuristic omission.
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

MODULE m_dmft_current_vertex

 use defs_basis
 use m_abicore
 use m_errors

 use m_dtset, only : dataset_type
 use m_paw_dmft, only : paw_dmft_type
 use m_pawtab, only : pawtab_type
 use m_crystal, only : crystal_t
 use m_pawcprj, only : pawcprj_type, pawcprj_get, pawcprj_alloc, pawcprj_free

 implicit none

 private

 public :: compute_psinablapsi_dmft

!!***

CONTAINS

!!****f* m_dmft_current_vertex/compute_psinablapsi_dmft
!! NAME
!!  compute_psinablapsi_dmft
!!
!! FUNCTION
!!  Compute momentum matrix elements <psi_a(k)|-i nabla_mu|psi_b(k)>
!!  for all correlated band pairs (a,b) at all k-points.
!!  The result is stored in paw_dmft%psinablapsi_dmft.
!!
!!  This function should be called from vtorho after datafordmft,
!!  where the plane-wave coefficients (cg) and PAW projector
!!  coefficients (cprj) are available.
!!
!! INPUTS
!!  cg(2,mcg) = plane-wave coefficients of wavefunctions
!!  cprj(natom,nspinor*mband_cprj*mkmem*nsppol*usecprj)
!!            = PAW projector coefficients
!!  kg(3,mpw*mkmem) = integer G-vector indices for each k-point
!!  gprimd(3,3) = reciprocal lattice vectors (Bohr^{-1})
!!  dtset = dataset structure with k-point and band info
!!  pawtab(ntypat) = PAW tabulated data (nabla_ij)
!!  cryst_struc = crystal structure (for atindx1)
!!  dimcprj(natom) = dimensions of cprj for each atom
!!  mcg = dimension of cg array
!!  mband_cprj = dimension of cprj band index
!!  my_nspinor = local number of spinor components
!!  usecprj = 1 if cprj is available
!!  mpi_comm = MPI communicator
!!  mpi_proc_distrb(nkpt,mband,nsppol) = processor distribution
!!
!! SIDE EFFECTS
!!  paw_dmft = on output, psinablapsi_dmft is filled and
!!             has_psinablapsi_dmft is set to 1
!!
!! SOURCE

subroutine compute_psinablapsi_dmft(paw_dmft, cg, cprj, kg, gprimd, dtset, &
  & pawtab, cryst_struc, dimcprj, mcg, mband_cprj, my_nspinor, usecprj, &
  & mpi_comm, mpi_proc_distrb)

 type(paw_dmft_type), intent(inout) :: paw_dmft
 real(dp), intent(in) :: cg(2,mcg)
 type(pawcprj_type), intent(in) :: cprj(paw_dmft%natom, &
   & my_nspinor*mband_cprj*paw_dmft%mkmem*paw_dmft%nsppol*usecprj)
 integer, intent(in) :: kg(3,*)
 real(dp), intent(in) :: gprimd(3,3)
 type(dataset_type), intent(in) :: dtset
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
 type(crystal_t), intent(in) :: cryst_struc
 integer, intent(in) :: dimcprj(paw_dmft%natom)
 integer, intent(in) :: mcg, mband_cprj, my_nspinor, usecprj
 integer, intent(in) :: mpi_comm
 integer, intent(in) :: mpi_proc_distrb(:,:,:)

!Local variables
 character(len=500) :: msg
 integer :: isppol, ikpt, npw_k, nband_k
 integer :: iband_a, iband_b, ia, ib, ipw, idir, ispinor
 integer :: icg, ikg, ibg
 integer :: iwf_a, iwf_b, ig_a, ig_b
 integer :: iatom, itypat, ilmn, jlmn, lmn_size
 integer :: ibsp_a, ibsp_b
 real(dp) :: cgnm1, cgnm2, cpnm1, cpnm2
 real(dp), allocatable :: kpg_k(:,:)
 real(dp) :: nabla_tmp(3)
 logical :: has_nabla, paw_augmentation_done
 type(pawcprj_type), allocatable :: cprj_k(:,:)
 integer :: mbandc, nband_k_cprj

! *********************************************************************

 write(msg,'(a)') ' compute_psinablapsi_dmft: Computing momentum matrix elements for DMFT response'
 call wrtout(std_out, msg)

 mbandc = paw_dmft%mbandc

 ! Allocate output array if needed
 if (.not. allocated(paw_dmft%psinablapsi_dmft)) then
   ABI_MALLOC(paw_dmft%psinablapsi_dmft, (2, 3, mbandc, mbandc, &
     & paw_dmft%nkpt, paw_dmft%nsppol))
 end if
 paw_dmft%psinablapsi_dmft = zero

 ! Check if PAW augmentation data (nabla_ij) is available
 has_nabla = .true.
 paw_augmentation_done = .false.
 do itypat = 1, paw_dmft%ntypat
   if (.not. allocated(pawtab(itypat)%nabla_ij)) then
     has_nabla = .false.
     exit
   end if
 end do

 if (.not. has_nabla) then
   write(msg,'(3a)') &
   ' WARNING: pawtab%nabla_ij not available.',ch10,&
   '   Only the kinetic (plane-wave) part of momentum matrix elements will be computed.'
   call wrtout(std_out, msg)
   write(msg,'(a)') '   To include PAW augmentation, set prtnabla>0 or call pawnabla_init.'
   call wrtout(std_out, msg)
 end if

 ! Main loop over spin polarizations and k-points
 icg = 0  ! plane-wave coefficient offset
 ikg = 0  ! G-vector index offset
 ibg = 0  ! cprj index offset

 do isppol = 1, paw_dmft%nsppol
   do ikpt = 1, paw_dmft%nkpt

     npw_k = paw_dmft%npwarr(ikpt)
     nband_k = dtset%nband(ikpt + (isppol-1)*paw_dmft%nkpt)

     ! --- Compute (k+G)_mu in Cartesian coordinates ---
     ! Allocated for nspinor copies so the loop ipw=1..npw*nspinor works
     ABI_MALLOC(kpg_k, (3, npw_k * my_nspinor))
     kpg_k = zero
     do ipw = 1, npw_k
       do idir = 1, 3
         kpg_k(idir, ipw) = gprimd(idir,1)*(dtset%kpt(1,ikpt)+dble(kg(1,ikg+ipw))) &
                           + gprimd(idir,2)*(dtset%kpt(2,ikpt)+dble(kg(2,ikg+ipw))) &
                           + gprimd(idir,3)*(dtset%kpt(3,ikpt)+dble(kg(3,ikg+ipw)))
       end do
     end do
     ! For nspinor=2, copy (k+G) vectors for second spinor component
     if (my_nspinor == 2) then
       kpg_k(:, npw_k+1:2*npw_k) = kpg_k(:, 1:npw_k)
     end if

     ! --- Kinetic part: sum_G (k+G)_mu * c*_a(G) * c_b(G) ---
     do iband_a = paw_dmft%dmftbandi, min(paw_dmft%dmftbandf, nband_k)
       ia = iband_a - paw_dmft%dmftbandi + 1
       iwf_a = icg + (iband_a - 1) * npw_k * my_nspinor

       do iband_b = paw_dmft%dmftbandi, min(paw_dmft%dmftbandf, nband_k)
         ib = iband_b - paw_dmft%dmftbandi + 1
         iwf_b = icg + (iband_b - 1) * npw_k * my_nspinor

         ! Sum over all plane waves (including spinor components)
         do ipw = 1, npw_k * my_nspinor
           ig_a = iwf_a + ipw
           ig_b = iwf_b + ipw

           ! Re(c*_a * c_b) = Re(a)*Re(b) + Im(a)*Im(b)
           cgnm1 = cg(1,ig_a)*cg(1,ig_b) + cg(2,ig_a)*cg(2,ig_b)
           ! Im(c*_a * c_b) = Re(a)*Im(b) - Im(a)*Re(b)
           cgnm2 = cg(1,ig_a)*cg(2,ig_b) - cg(2,ig_a)*cg(1,ig_b)

           do idir = 1, 3
             paw_dmft%psinablapsi_dmft(1, idir, ia, ib, ikpt, isppol) = &
               paw_dmft%psinablapsi_dmft(1, idir, ia, ib, ikpt, isppol) + cgnm1 * kpg_k(idir, ipw)
             paw_dmft%psinablapsi_dmft(2, idir, ia, ib, ikpt, isppol) = &
               paw_dmft%psinablapsi_dmft(2, idir, ia, ib, ikpt, isppol) + cgnm2 * kpg_k(idir, ipw)
           end do
         end do ! ipw

       end do ! iband_b
     end do ! iband_a

     ABI_FREE(kpg_k)

     ! --- PAW augmentation: sum_{ij,atom} cprj*(atom,a)_i * nabla_ij * cprj(atom,b)_j ---
     if (has_nabla .and. usecprj == 1) then

       nband_k_cprj = nband_k
       ! Extract cprj for ALL correlated bands at this k-point
       ABI_MALLOC(cprj_k, (paw_dmft%natom, my_nspinor * nband_k))
       call pawcprj_alloc(cprj_k, 0, dimcprj)
       call pawcprj_get(cryst_struc%atindx1, cprj_k, cprj, paw_dmft%natom, &
         & 1, ibg, ikpt, 0, isppol, mband_cprj, &
         & paw_dmft%mkmem, paw_dmft%natom, nband_k, nband_k_cprj, &
         & my_nspinor, paw_dmft%nsppol, paw_dmft%unpaw, &
         & mpicomm=mpi_comm, proc_distrb=mpi_proc_distrb)

       do iband_a = paw_dmft%dmftbandi, min(paw_dmft%dmftbandf, nband_k)
         ia = iband_a - paw_dmft%dmftbandi + 1

         do iband_b = paw_dmft%dmftbandi, min(paw_dmft%dmftbandf, nband_k)
           ib = iband_b - paw_dmft%dmftbandi + 1

           do iatom = 1, paw_dmft%natom
             itypat = paw_dmft%typat(iatom)
             lmn_size = pawtab(itypat)%lmn_size

             do jlmn = 1, lmn_size
               do ilmn = 1, lmn_size
                 nabla_tmp(:) = pawtab(itypat)%nabla_ij(:, ilmn, jlmn)

                 ! Sum over spinor components
                 cpnm1 = zero
                 cpnm2 = zero
                 do ispinor = 1, my_nspinor
                   ibsp_a = (iband_a - 1) * my_nspinor + ispinor
                   ibsp_b = (iband_b - 1) * my_nspinor + ispinor

                   ! Re(cprj*_a * cprj_b)
                   cpnm1 = cpnm1 + &
                     cprj_k(iatom,ibsp_a)%cp(1,ilmn)*cprj_k(iatom,ibsp_b)%cp(1,jlmn) + &
                     cprj_k(iatom,ibsp_a)%cp(2,ilmn)*cprj_k(iatom,ibsp_b)%cp(2,jlmn)
                   ! Im(cprj*_a * cprj_b)
                   cpnm2 = cpnm2 + &
                     cprj_k(iatom,ibsp_a)%cp(1,ilmn)*cprj_k(iatom,ibsp_b)%cp(2,jlmn) - &
                     cprj_k(iatom,ibsp_a)%cp(2,ilmn)*cprj_k(iatom,ibsp_b)%cp(1,jlmn)
                 end do ! ispinor

                 ! PAW nabla contribution: -i * (cprj*_a cprj_b) * nabla_ij
                 ! <psi_a| (-i nabla) |psi_b> += -i * (cpnm1 + i*cpnm2) * nabla_ij
                 !   = (cpnm2 - i*cpnm1) * nabla_ij
                 ! Real part: cpnm2 * nabla_ij
                 ! Imag part: -cpnm1 * nabla_ij
                 do idir = 1, 3
                   paw_dmft%psinablapsi_dmft(1, idir, ia, ib, ikpt, isppol) = &
                     paw_dmft%psinablapsi_dmft(1, idir, ia, ib, ikpt, isppol) + cpnm2 * nabla_tmp(idir)
                   paw_dmft%psinablapsi_dmft(2, idir, ia, ib, ikpt, isppol) = &
                     paw_dmft%psinablapsi_dmft(2, idir, ia, ib, ikpt, isppol) - cpnm1 * nabla_tmp(idir)
                 end do

               end do ! ilmn
             end do ! jlmn
           end do ! iatom

         end do ! iband_b
       end do ! iband_a

       call pawcprj_free(cprj_k)
       ABI_FREE(cprj_k)
       paw_augmentation_done = .true.
     end if ! has_nabla

     ! Advance indices
     ikg = ikg + npw_k
     icg = icg + nband_k * npw_k * my_nspinor
     ibg = ibg + nband_k * my_nspinor

   end do ! ikpt
 end do ! isppol

 paw_dmft%has_psinablapsi_dmft = 1

 write(msg,'(a,i6,a,i4,a,l1)') &
 ' compute_psinablapsi_dmft: Done. mbandc=', mbandc, &
 ' nkpt=', paw_dmft%nkpt, ' PAW_augmentation=', paw_augmentation_done
 call wrtout(std_out, msg)

end subroutine compute_psinablapsi_dmft

!!***

END MODULE m_dmft_current_vertex
