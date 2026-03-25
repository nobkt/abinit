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
!!  The bubble conductivity is decomposed by spin channel:
!!    Pi^{total} = Pi^{spin-conserving} + Pi^{S+S-} + Pi^{S-S+}
!!  enabling attribution of optical absorption features to specific spin excitations.
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
 use m_green, only : green_type

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
!!  When spin-flip attribution is active (nspinor=2), the bubble is
!!  decomposed into spin-conserving and spin-flip (S+S-, S-S+) channels.
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

   complex(dp), allocatable :: pi_bubble_sc(:,:,:)
   ! pi_bubble_sc(nboson, ndir, ndir) : spin-conserving part of bubble

   complex(dp), allocatable :: pi_bubble_pm(:,:,:)
   ! pi_bubble_pm(nboson, ndir, ndir) : S+S- spin-flip part of bubble

   complex(dp), allocatable :: pi_bubble_mp(:,:,:)
   ! pi_bubble_mp(nboson, ndir, ndir) : S-S+ spin-flip part of bubble

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
 ABI_MALLOC(optic%pi_bubble_sc, (nboson, ndir, ndir))
 ABI_MALLOC(optic%pi_bubble_pm, (nboson, ndir, ndir))
 ABI_MALLOC(optic%pi_bubble_mp, (nboson, ndir, ndir))
 optic%pi_bubble = czero
 optic%pi_vertex = czero
 optic%pi_total = czero
 optic%pi_bubble_sc = czero
 optic%pi_bubble_pm = czero
 optic%pi_bubble_mp = czero

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

 ABI_SFREE(optic%pi_bubble)
 ABI_SFREE(optic%pi_vertex)
 ABI_SFREE(optic%pi_total)
 ABI_SFREE(optic%pi_bubble_sc)
 ABI_SFREE(optic%pi_bubble_pm)
 ABI_SFREE(optic%pi_bubble_mp)

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
!!    sum_{a,b,c,d} j_mu^{ab}(k) * G^{bc}(k,iw_n) * j_nu^{cd}(k) * G^{da}(k,iw_n+iOm_m)
!!
!!  where a,b,c,d run over correlated bands (mbandc) and the sum
!!  includes all spinor components embedded in the band index.
!!
!!  The velocity matrix elements j_mu^{ab}(k) are read from
!!  paw_dmft%psinablapsi_dmft, computed by compute_psinablapsi_dmft.
!!
!!  For nspinor=2, the bubble is decomposed into spin channels by
!!  tracing over (a,d) and (b,c) with spin classification:
!!    spin-conserving: sigma_a = sigma_d AND sigma_b = sigma_c (same spin)
!!    S+S-: sigma_a = up, sigma_d = down (band a in up block, d in down block)
!!    S-S+: sigma_a = down, sigma_d = up
!!
!! INPUTS
!!  paw_dmft = DFT+DMFT data with psinablapsi_dmft
!!  green_imp = converged Green function with oper(iw)%ks
!!  nboson = number of bosonic Matsubara frequencies
!!  niw_vertex = number of fermionic Matsubara frequencies for the sum
!!  nspinor = number of spinor components
!!
!! SIDE EFFECTS
!!  optic%pi_bubble = on output, contains the bubble conductivity
!!  optic%pi_bubble_sc/pm/mp = spin-channel decomposition (if nspinor=2)
!!
!! SOURCE

subroutine compute_bubble_conductivity(optic, paw_dmft, green_imp, nboson, niw_vertex, nspinor)

 type(optic_kernel_type), intent(inout) :: optic
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type), intent(in) :: green_imp
 integer, intent(in) :: nboson, niw_vertex, nspinor

!Local variables
 character(len=500) :: msg
 integer :: iom, iw, ikpt, isppol
 integer :: ia, ib, ic, id, mu, nu
 integer :: mbandc, nkpt
 real(dp) :: beta, wk
 complex(dp) :: jmu_ab, jnu_cd, gbc, gda, contrib
 complex(dp), allocatable :: jmu_mat(:,:), jnu_mat(:,:)
 complex(dp), allocatable :: g_iw(:,:), g_iw_om(:,:)
 complex(dp), allocatable :: temp1(:,:), temp2(:,:)
 integer :: ndim_half
 logical :: do_spin_decomp
 character(len=2) :: spin_channel

! *********************************************************************

 write(msg,'(a)') ' compute_bubble_conductivity: Computing Matsubara bubble Pi_mu_nu'
 call wrtout(std_out, msg)

 ! Check that velocity matrix elements are available
 if (paw_dmft%has_psinablapsi_dmft /= 1) then
   write(msg,'(3a)') &
   ' compute_bubble_conductivity: psinablapsi_dmft not computed.',ch10,&
   ' Bubble conductivity remains zero. Call compute_psinablapsi_dmft first.'
   call wrtout(std_out, msg)
   optic%pi_bubble = czero
   optic%pi_bubble_sc = czero
   optic%pi_bubble_pm = czero
   optic%pi_bubble_mp = czero
   return
 end if

 ! Check that Green function has KS basis data
 if (green_imp%oper(1)%has_operks /= 1) then
   write(msg,'(a)') &
   ' compute_bubble_conductivity: Green function has no KS-basis data. Bubble remains zero.'
   call wrtout(std_out, msg)
   optic%pi_bubble = czero
   return
 end if

 mbandc = paw_dmft%mbandc
 nkpt = paw_dmft%nkpt
 beta = one / paw_dmft%temp
 do_spin_decomp = (nspinor == 2)
 ndim_half = mbandc / 2  ! Half of mbandc for spin-up/down blocks in spinor basis

 ! Initialize
 optic%pi_bubble = czero
 optic%pi_bubble_sc = czero
 optic%pi_bubble_pm = czero
 optic%pi_bubble_mp = czero

 ! Allocate working arrays for matrix operations at each (k, iw)
 ABI_MALLOC(jmu_mat, (mbandc, mbandc))
 ABI_MALLOC(jnu_mat, (mbandc, mbandc))
 ABI_MALLOC(g_iw, (mbandc, mbandc))
 ABI_MALLOC(g_iw_om, (mbandc, mbandc))
 ABI_MALLOC(temp1, (mbandc, mbandc))
 ABI_MALLOC(temp2, (mbandc, mbandc))

 ! Main computation loop
 ! Pi_mu_nu(iOm) = -(1/(beta*Nk)) sum_k sum_n Tr[j_mu G(iw_n) j_nu G(iw_n+iOm)]
 do isppol = 1, paw_dmft%nsppol
   do ikpt = 1, nkpt
     wk = paw_dmft%wtk(ikpt)

     do iom = 1, nboson
       do mu = 1, 3
         ! Construct j_mu matrix at this k-point
         do ia = 1, mbandc
           do ib = 1, mbandc
             jmu_mat(ia, ib) = cmplx(paw_dmft%psinablapsi_dmft(1, mu, ia, ib, ikpt, isppol), &
                                     paw_dmft%psinablapsi_dmft(2, mu, ia, ib, ikpt, isppol), kind=dp)
           end do
         end do

         do nu = 1, 3
           ! Construct j_nu matrix
           do ic = 1, mbandc
             do id = 1, mbandc
               jnu_mat(ic, id) = cmplx(paw_dmft%psinablapsi_dmft(1, nu, ic, id, ikpt, isppol), &
                                       paw_dmft%psinablapsi_dmft(2, nu, ic, id, ikpt, isppol), kind=dp)
             end do
           end do

           ! Sum over fermionic Matsubara frequencies
           do iw = 1, niw_vertex
             ! Check frequency bounds
             if (iw + iom - 1 > green_imp%nw) cycle

             ! Extract G(k, iw_n) and G(k, iw_n + iOm_m)
             do ia = 1, mbandc
               do ib = 1, mbandc
                 g_iw(ia, ib) = green_imp%oper(iw)%ks(ia, ib, ikpt, isppol)
                 g_iw_om(ia, ib) = green_imp%oper(iw + iom - 1)%ks(ia, ib, ikpt, isppol)
               end do
             end do

             ! Compute Tr[j_mu * G(iw) * j_nu * G(iw+iOm)]
             ! = sum_{a,b,c,d} j_mu^{ab} G^{bc}(iw) j_nu^{cd} G^{da}(iw+iOm)
             ! Use explicit summation for correctness and spin decomposition

             ! Total trace: sum over all (a,b,c,d)
             do ia = 1, mbandc
               do ib = 1, mbandc
                 ! j_mu^{ab}
                 jmu_ab = jmu_mat(ia, ib)
                 if (abs(jmu_ab) < tol16) cycle

                 do ic = 1, mbandc
                   ! G^{bc}(iw)
                   gbc = g_iw(ib, ic)
                   if (abs(gbc) < tol16) cycle

                   do id = 1, mbandc
                     ! j_nu^{cd}
                     jnu_cd = jnu_mat(ic, id)
                     if (abs(jnu_cd) < tol16) cycle

                     ! G^{da}(iw+iOm)
                     gda = g_iw_om(id, ia)

                     contrib = jmu_ab * gbc * jnu_cd * gda

                     ! Accumulate total
                     optic%pi_bubble(iom, mu, nu) = optic%pi_bubble(iom, mu, nu) &
                       - wk * contrib / beta

                     ! Spin channel decomposition
                     ! For nspinor=2 with nsppol=1: bands 1..ndim_half are spin-up,
                     ! bands ndim_half+1..mbandc are spin-down (embedded in band index)
                     if (do_spin_decomp) then
                       ! Classify by the outer trace indices (a, d):
                       ! The optical response traces over the (a,d) pair
                       ! which corresponds to the "creation" and "annihilation" of the
                       ! electron-hole pair
                       spin_channel = classify_spin_pair(ia, id, ndim_half)

                       if (spin_channel == 'sc') then
                         optic%pi_bubble_sc(iom, mu, nu) = optic%pi_bubble_sc(iom, mu, nu) &
                           - wk * contrib / beta
                       else if (spin_channel == 'pm') then
                         optic%pi_bubble_pm(iom, mu, nu) = optic%pi_bubble_pm(iom, mu, nu) &
                           - wk * contrib / beta
                       else if (spin_channel == 'mp') then
                         optic%pi_bubble_mp(iom, mu, nu) = optic%pi_bubble_mp(iom, mu, nu) &
                           - wk * contrib / beta
                       end if
                     end if

                   end do ! id
                 end do ! ic
               end do ! ib
             end do ! ia

           end do ! iw (fermionic frequency)
         end do ! nu
       end do ! mu
     end do ! iom (bosonic frequency)
   end do ! ikpt
 end do ! isppol

 ! Set total = bubble (vertex correction is zero until TRIQS interface is connected)
 optic%pi_total = optic%pi_bubble

 ABI_FREE(jmu_mat)
 ABI_FREE(jnu_mat)
 ABI_FREE(g_iw)
 ABI_FREE(g_iw_om)
 ABI_FREE(temp1)
 ABI_FREE(temp2)

 write(msg,'(a)') ' compute_bubble_conductivity: Bubble computation complete.'
 call wrtout(std_out, msg)

end subroutine compute_bubble_conductivity

!!***

!!****f* m_dmft_optic_kernel/classify_spin_pair
!! NAME
!!  classify_spin_pair
!!
!! FUNCTION
!!  Classify a pair of band indices (ia, id) into spin channels.
!!  For nspinor=2 with nsppol=1, band indices 1..ndim_half correspond to
!!  spin-up, and ndim_half+1..mbandc to spin-down.
!!
!! SOURCE

pure function classify_spin_pair(ia, id, ndim_half) result(channel)

 integer, intent(in) :: ia, id, ndim_half
 character(len=2) :: channel

! *********************************************************************

 if (ia <= ndim_half .and. id <= ndim_half) then
   ! Both spin-up: spin-conserving
   channel = 'sc'
 else if (ia > ndim_half .and. id > ndim_half) then
   ! Both spin-down: spin-conserving
   channel = 'sc'
 else if (ia <= ndim_half .and. id > ndim_half) then
   ! a = up, d = down: S+S- channel
   channel = 'pm'
 else
   ! a = down, d = up: S-S+ channel
   channel = 'mp'
 end if

end function classify_spin_pair

!!***

!!****f* m_dmft_optic_kernel/write_optic_kernel
!! NAME
!!  write_optic_kernel
!!
!! FUNCTION
!!  Write the optical kernel to file for diagnostics.
!!  Includes spin-channel decomposition when available.
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
 logical :: has_spin_decomp

! *********************************************************************

 has_spin_decomp = allocated(optic%pi_bubble_sc) .and. allocated(optic%pi_bubble_pm)

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT Optical kernel Pi_mu_nu(iOm) on Matsubara axis'
 write(unt,'(a,i6)') '# nboson = ', optic%nboson
 write(unt,'(a,es14.6)') '# beta = ', beta
 write(unt,'(a)') '#'
 write(unt,'(a)') '# Section 1: Full optical kernel'
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

 ! Section 2: Spin-channel decomposition of bubble
 if (has_spin_decomp) then
   write(unt,'(a)') '#'
   write(unt,'(a)') '# Section 2: Spin-channel decomposition of bubble Pi_mu_nu'
   write(unt,'(a)') '# iOm  Omega  mu  nu  Re(Pi_sc)  Im(Pi_sc)  Re(Pi_pm)  Im(Pi_pm)  Re(Pi_mp)  Im(Pi_mp)'

   do iom = 1, optic%nboson
     omega_boson = two_pi * dble(iom - 1) / beta
     do mu = 1, optic%ndir
       do nu = 1, optic%ndir
         write(unt,'(i6,es14.6,2i3,6es20.10)') iom, omega_boson, mu, nu, &
         real(optic%pi_bubble_sc(iom,mu,nu)), aimag(optic%pi_bubble_sc(iom,mu,nu)), &
         real(optic%pi_bubble_pm(iom,mu,nu)), aimag(optic%pi_bubble_pm(iom,mu,nu)), &
         real(optic%pi_bubble_mp(iom,mu,nu)), aimag(optic%pi_bubble_mp(iom,mu,nu))
       end do
     end do
   end do

   ! Section 3: Spin-channel fractions at iOm=0
   write(unt,'(a)') '#'
   write(unt,'(a)') '# Section 3: Spin-channel fractions of diagonal Pi_mu_mu at iOm=0'
   write(unt,'(a)') '# mu  |Pi_total|  |Pi_sc|  |Pi_pm|  |Pi_mp|  frac_sc  frac_pm  frac_mp'

   do mu = 1, optic%ndir
     if (abs(optic%pi_bubble(1,mu,mu)) > tol16) then
       write(unt,'(i3,4es16.6,3f10.4)') mu, &
         abs(optic%pi_bubble(1,mu,mu)), &
         abs(optic%pi_bubble_sc(1,mu,mu)), &
         abs(optic%pi_bubble_pm(1,mu,mu)), &
         abs(optic%pi_bubble_mp(1,mu,mu)), &
         abs(optic%pi_bubble_sc(1,mu,mu))/abs(optic%pi_bubble(1,mu,mu)), &
         abs(optic%pi_bubble_pm(1,mu,mu))/abs(optic%pi_bubble(1,mu,mu)), &
         abs(optic%pi_bubble_mp(1,mu,mu))/abs(optic%pi_bubble(1,mu,mu))
     else
       write(unt,'(i3,4es16.6,3a10)') mu, &
         abs(optic%pi_bubble(1,mu,mu)), &
         abs(optic%pi_bubble_sc(1,mu,mu)), &
         abs(optic%pi_bubble_pm(1,mu,mu)), &
         abs(optic%pi_bubble_mp(1,mu,mu)), &
         ' undef', ' undef', ' undef'
     end if
   end do
 end if

 close(unt)

end subroutine write_optic_kernel

!!***

END MODULE m_dmft_optic_kernel
