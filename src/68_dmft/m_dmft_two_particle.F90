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
 use m_green, only : green_type

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
!!  Following the Fourier convention of Section 5.9 of the design document
!!  (Rohringer et al., Rev. Mod. Phys. 90, 025003 (2018)).
!!
!!  The Green function is extracted from the converged green_type object,
!!  which stores G^imp(iw_n) in the local basis (matlu) for each frequency.
!!
!! INPUTS
!!  green_imp = converged impurity Green function (green_type)
!!  paw_dmft = DFT+DMFT data with system parameters (temperature, atom info)
!!  norb_corr = number of correlated spinor-orbitals = (2*lpawu+1)*nspinor
!!  niw_vertex = number of fermionic Matsubara frequencies for vertex
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SIDE EFFECTS
!!  chi0 = on output, contains the computed bare impurity bubble
!!
!! NOTES
!!  Requires nspinor=2 for spin-flip calculations (nsppol=1 in this case).
!!  The frequency shift iw_n + iOm_m = iw_{n+m}, so we need
!!  niw_vertex + nboson - 1 <= green_imp%nw.
!!
!! SOURCE

subroutine compute_chi0_imp(chi0, green_imp, paw_dmft, norb_corr, niw_vertex, nboson)

 type(chi_loc_type), intent(inout) :: chi0
 type(green_type), intent(in) :: green_imp
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: norb_corr, niw_vertex, nboson

!Local variables
 integer :: iom, iw, iw_shifted, ialpha, ibeta, igamma, idelta
 integer :: idx_i, idx_j, iatom, iatom_corr
 integer :: norb_sq, ndim_matlu
 real(dp) :: beta
 complex(dp), allocatable :: gmat_n(:,:), gmat_np(:,:)
 character(len=500) :: msg

! *********************************************************************

 write(msg,'(a)') ' compute_chi0_imp: Computing bare impurity bubble from Green function'
 call wrtout(std_out, msg)

 beta = one / paw_dmft%temp

 ! --- Validate frequency bounds ---
 ! Shifted frequency index: iw + (iom - 1) must not exceed green_imp%nw
 if (niw_vertex + nboson - 1 > green_imp%nw) then
   write(msg,'(a,i6,a,i6,a,i6)') &
   'compute_chi0_imp: frequency bounds exceeded. niw_vertex=', niw_vertex, &
   ' nboson=', nboson, ' but green has nw=', green_imp%nw
   ABI_ERROR(msg)
 end if

 ! --- Find the first correlated atom ---
 iatom_corr = 0
 do iatom = 1, paw_dmft%natom
   if (paw_dmft%lpawu(iatom) >= 0) then
     iatom_corr = iatom
     exit
   end if
 end do

 if (iatom_corr == 0) then
   ABI_ERROR('compute_chi0_imp: no correlated atom found (no atoms with lpawu >= 0)')
 end if

 write(msg,'(a,i4,a,i2)') &
 ' compute_chi0_imp: Using correlated atom ', iatom_corr, &
 ' with lpawu=', paw_dmft%lpawu(iatom_corr)
 call wrtout(std_out, msg)

 ! --- Validate Green function data ---
 if (green_imp%oper(1)%has_opermatlu /= 1) then
   ABI_ERROR('compute_chi0_imp: Green function does not have matlu (local basis) data')
 end if

 ! Validate orbital dimension: matlu stores (2*lpawu+1)*nspinor x (2*lpawu+1)*nspinor
 ndim_matlu = (2 * paw_dmft%lpawu(iatom_corr) + 1) * paw_dmft%nspinor
 if (ndim_matlu /= norb_corr) then
   write(msg,'(a,i4,a,i4)') &
   'compute_chi0_imp: matlu dimension mismatch. ndim_matlu=', ndim_matlu, &
   ' norb_corr=', norb_corr
   ABI_ERROR(msg)
 end if

 norb_sq = norb_corr * norb_corr

 ABI_MALLOC(gmat_n, (norb_corr, norb_corr))
 ABI_MALLOC(gmat_np, (norb_corr, norb_corr))

 chi0%chi_mat = czero

 write(msg,'(a,es14.6,a,i4,a,i4)') &
 ' compute_chi0_imp: beta=', beta, ' niw_vertex=', niw_vertex, ' nboson=', nboson
 call wrtout(std_out, msg)

 ! --- Compute chi0 ---
 ! chi0_{ab,gd}(n,n';Om) = -beta * delta_{nn'} * G_{da}(iwn) * G_{bg}(iwn+iOm)
 ! Packed index: I = (n-1)*norb^2 + (a-1)*norb + b  (1-based)
 !
 ! The bubble is block-diagonal in fermionic frequency (delta_{nn'}).
 ! For each bosonic frequency iOm_m (m = iom-1, so iom=1 gives iOm_0=0):
 !   shifted fermionic index = iw + (iom - 1)
 !   i.e., iw_n + iOm_m = iw_{n + m}

 do iom = 1, nboson
   do iw = 1, niw_vertex
     iw_shifted = iw + (iom - 1)

     ! Extract G^imp matrices at iw_n and iw_n + iOm_m
     ! For nspinor=2: nsppol=1, so isppol index is always 1
     gmat_n(:,:) = green_imp%oper(iw)%matlu(iatom_corr)%mat(:,:,1)
     gmat_np(:,:) = green_imp%oper(iw_shifted)%matlu(iatom_corr)%mat(:,:,1)

     do ialpha = 1, norb_corr
       do ibeta = 1, norb_corr
         idx_i = (iw - 1) * norb_sq + (ialpha - 1) * norb_corr + ibeta
         do igamma = 1, norb_corr
           do idelta = 1, norb_corr
             ! J uses the same frequency n as I due to delta_{nn'}
             idx_j = (iw - 1) * norb_sq + (igamma - 1) * norb_corr + idelta
             ! chi0 = -beta * G_{delta,alpha}(iwn) * G_{beta,gamma}(iwn+iOm)
             chi0%chi_mat(iom, idx_i, idx_j) = &
               -beta * gmat_n(idelta, ialpha) * gmat_np(ibeta, igamma)
           end do
         end do
       end do
     end do
   end do
 end do

 ABI_FREE(gmat_n)
 ABI_FREE(gmat_np)

 write(msg,'(a,i8,a)') &
 ' compute_chi0_imp: Bubble computed (', chi0%ndim_comp, ' composite indices)'
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
