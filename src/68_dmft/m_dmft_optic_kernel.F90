!!****m* ABINIT/m_dmft_optic_kernel
!! NAME
!!  m_dmft_optic_kernel
!!
!! FUNCTION
!!  Optical kernel for DFT+DMFT absorption spectra.
!!  Implements the current-current correlation function Pi_mu_nu on the
!!  Matsubara axis and the dressed current vertex.
!!
!!  Pi_mu_nu^bubble(iOm) = -(1/beta*Nk) sum_{k,n} Tr[ j_mu(k) G(k,iw_n) j_nu(k) G(k,iw_n+iOm) ]
!!
!!  The dressed vertex satisfies:
!!    j_tilde_nu = [ 1 - Gamma_imp * chi0_bar ]^{-1} * j_bar_nu
!!
!!  For optical absorption, tensor components mu,nu in {x,y,z} are all computed
!!  without isotropic approximation (required for anisotropic crystals like MnF2).
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

MODULE m_dmft_optic_kernel

 use defs_basis
 use m_abicore
 use m_errors

 use m_paw_dmft, only : paw_dmft_type

 implicit none

 private

 public :: optic_kernel_type
 public :: init_optic_kernel
 public :: destroy_optic_kernel
 public :: compute_bubble_conductivity
 public :: write_optic_kernel

!!***

!!****t* m_dmft_optic_kernel/optic_kernel_type
!! NAME
!!  optic_kernel_type
!!
!! FUNCTION
!!  Stores the Matsubara-axis optical kernel Pi_mu_nu(iOm_m).
!!
!!  pi_matsubara(iOm, mu, nu) : current-current correlation on Matsubara axis
!!  ndir = 3 for {x,y,z} spatial directions
!!
!! SOURCE

 type, public :: optic_kernel_type

   integer :: nboson = 0
   integer :: ndir = 0

   complex(dp), allocatable :: pi_bubble(:,:,:)
   ! pi_bubble(nboson, ndir, ndir) : bubble contribution

   complex(dp), allocatable :: pi_vertex(:,:,:)
   ! pi_vertex(nboson, ndir, ndir) : vertex correction

   complex(dp), allocatable :: pi_total(:,:,:)
   ! pi_total(nboson, ndir, ndir) : total Pi_mu_nu

 end type optic_kernel_type

!!***

CONTAINS

!!****f* m_dmft_optic_kernel/init_optic_kernel
!! NAME
!!  init_optic_kernel
!!
!! FUNCTION
!!  Initialize the optical kernel structure.
!!
!! SOURCE

subroutine init_optic_kernel(optic, nboson, ndir)

 type(optic_kernel_type), intent(inout) :: optic
 integer, intent(in) :: nboson, ndir

!Local variables
 character(len=500) :: msg

! *********************************************************************

 optic%nboson = nboson
 optic%ndir = ndir

 write(msg,'(a,i4,a,i2)') &
 ' init_optic_kernel: nboson=', nboson, ' ndir=', ndir
 call wrtout(std_out, msg)

 ABI_MALLOC(optic%pi_bubble, (nboson, ndir, ndir))
 ABI_MALLOC(optic%pi_vertex, (nboson, ndir, ndir))
 ABI_MALLOC(optic%pi_total, (nboson, ndir, ndir))
 optic%pi_bubble = czero
 optic%pi_vertex = czero
 optic%pi_total = czero

end subroutine init_optic_kernel

!!***

!!****f* m_dmft_optic_kernel/destroy_optic_kernel
!! NAME
!!  destroy_optic_kernel
!!
!! FUNCTION
!!  Deallocate the optical kernel structure.
!!
!! SOURCE

subroutine destroy_optic_kernel(optic)

 type(optic_kernel_type), intent(inout) :: optic

! *********************************************************************

 if (allocated(optic%pi_bubble)) then
   ABI_FREE(optic%pi_bubble)
 end if
 if (allocated(optic%pi_vertex)) then
   ABI_FREE(optic%pi_vertex)
 end if
 if (allocated(optic%pi_total)) then
   ABI_FREE(optic%pi_total)
 end if

 optic%nboson = 0
 optic%ndir = 0

end subroutine destroy_optic_kernel

!!***

!!****f* m_dmft_optic_kernel/compute_bubble_conductivity
!! NAME
!!  compute_bubble_conductivity
!!
!! FUNCTION
!!  Compute the bubble (vertex-uncorrected) part of the current-current
!!  correlation function on the Matsubara axis:
!!
!!  Pi_mu_nu^bubble(iOm_m) = -(1/(beta*Nk)) sum_{k,n}
!!    Tr[ j_mu(k) G(k,iw_n) j_nu(k) G(k,iw_n+iOm_m) ]
!!
!!  where the trace is over spinor-orbital indices.
!!
!! INPUTS
!!  paw_dmft = DFT+DMFT data with converged lattice Green function
!!  nboson = number of bosonic Matsubara frequencies
!!
!! SIDE EFFECTS
!!  optic%pi_bubble = on output, contains the bubble conductivity
!!
!! NOTES
!!  Skeleton implementation. Requires:
!!  1. Current vertex j_mu(k) computation (PAW velocity matrix elements)
!!  2. Lattice Green function G(k,iw) access
!!  3. Spinor trace over product j*G*j*G
!!
!!  The current vertex should include PAW nonlocal corrections:
!!    j_mu^PAW = j_mu^local + (i/hbar)[V_NL, r_mu]
!!  as implemented in src/65_paw/m_paw_optics.F90 for IPA.
!!
!! SOURCE

subroutine compute_bubble_conductivity(optic, paw_dmft, nboson)

 type(optic_kernel_type), intent(inout) :: optic
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: nboson

!Local variables
 character(len=500) :: msg

! *********************************************************************

 write(msg,'(a)') ' compute_bubble_conductivity: Computing Matsubara bubble Pi_mu_nu'
 call wrtout(std_out, msg)

 ! NOTE: Full implementation requires:
 ! 1. Computing current (velocity) matrix elements j_mu(k) for each k-point
 !    using PAW momentum matrix elements (cf. m_paw_optics.F90)
 !    with spinor components j_mu;a_up,b_down(k) preserved
 ! 2. For each k, iw_n, iOm_m:
 !    Compute Tr[ j_mu(k) * G(k,iw_n) * j_nu(k) * G(k,iw_n+iOm_m) ]
 ! 3. Sum over k and n with proper normalization
 !
 ! The current matrix elements for spinor basis must include:
 ! - Local contribution: -e * dH^KS/dk_mu (including spin-mixed components)
 ! - PAW nonlocal correction: (i/hbar)[V_NL, r_mu]
 !
 ! This is a structural skeleton establishing the interface.
 ! Actual computation will be connected after the one-particle Green
 ! function export infrastructure is in place.

 optic%pi_bubble = czero

 write(msg,'(a)') ' compute_bubble_conductivity: Bubble structure initialized (skeleton)'
 call wrtout(std_out, msg)

end subroutine compute_bubble_conductivity

!!***

!!****f* m_dmft_optic_kernel/write_optic_kernel
!! NAME
!!  write_optic_kernel
!!
!! FUNCTION
!!  Write the optical kernel to file for diagnostics.
!!
!! SOURCE

subroutine write_optic_kernel(optic, fname, beta)

 type(optic_kernel_type), intent(in) :: optic
 character(len=*), intent(in) :: fname
 real(dp), intent(in) :: beta

!Local variables
 integer :: unt, iom, mu, nu, ios
 real(dp) :: omega_boson
 character(len=500) :: msg

! *********************************************************************

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT Optical kernel Pi_mu_nu(iOm) on Matsubara axis'
 write(unt,'(a,i6)') '# nboson = ', optic%nboson
 write(unt,'(a,es14.6)') '# beta = ', beta
 write(unt,'(a)') '# iOm  Omega_boson  mu  nu  Re(Pi_bubble)  Im(Pi_bubble)  Re(Pi_total)  Im(Pi_total)'

 do iom = 1, optic%nboson
   omega_boson = two_pi * dble(iom - 1) / beta  ! iOm_m = 2*pi*m/beta for bosons (m=0,1,2,...)
   do mu = 1, optic%ndir
     do nu = 1, optic%ndir
       write(unt,'(i6,es14.6,2i3,4es20.10)') iom, omega_boson, mu, nu, &
       real(optic%pi_bubble(iom,mu,nu)), aimag(optic%pi_bubble(iom,mu,nu)), &
       real(optic%pi_total(iom,mu,nu)), aimag(optic%pi_total(iom,mu,nu))
     end do
   end do
 end do

 close(unt)

end subroutine write_optic_kernel

!!***

END MODULE m_dmft_optic_kernel
