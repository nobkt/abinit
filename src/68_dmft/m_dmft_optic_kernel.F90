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
 use m_dmft_spinor_proj, only : spinor_proj_type

 implicit none

 private

 public :: optic_kernel_type
 public :: init_optic_kernel
 public :: destroy_optic_kernel
 public :: compute_bubble_conductivity
 public :: write_optic_kernel
 public :: optic_orbital_attrib_type
 public :: init_optic_orbital_attrib
 public :: destroy_optic_orbital_attrib
 public :: compute_optic_orbital_attrib
 public :: write_optic_orbital_attrib

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

!!****t* m_dmft_optic_kernel/optic_orbital_attrib_type
!! NAME
!!  optic_orbital_attrib_type
!!
!! FUNCTION
!!  Stores the orbital-pair-resolved decomposition of the projected
!!  bubble optical conductivity Pi_mu_nu.
!!
!!  The velocity matrix j_mu and Green function G are projected from the
!!  KS band basis (a,b) into the correlated orbital basis (alpha,beta)
!!  using chipsi projectors:
!!    J_mu^{alpha,beta}(k) = sum_{a,b} P_{alpha,a}(k) j_mu^{ab}(k) P*_{beta,b}(k)
!!    G^loc_{alpha,beta}(k,iw) = sum_{a,b} P_{alpha,a}(k) G^{ab}(k,iw) P*_{beta,b}(k)
!!
!!  The projected bubble is:
!!    Pi_proj(iOm) = -(1/beta*Nk) sum_{k,n} Tr[J_mu G^loc(iw) J_nu G^loc(iw+iOm)]
!!
!!  The (alpha,delta) resolved contribution is then classified by spin
!!  channel and orbital index for the S+S- and S-S+ channels.
!!
!!  Note: Pi_proj covers only the intra-correlated-subspace part.
!!  Its total trace is generally smaller than the full KS-basis Pi_total.
!!  The ratio pi_projected_total / pi_bubble gives the fraction of
!!  optical response captured by correlated orbital transitions.
!!
!! SOURCE

 type, public :: optic_orbital_attrib_type

   integer :: nboson = 0
   integer :: ndir = 0
   integer :: ndim_orb = 0
   integer :: nspinor = 0

   complex(dp), allocatable :: pi_orb_pm(:,:,:,:,:)
   ! pi_orb_pm(nboson, ndir, ndir, ndim_orb, ndim_orb) :
   ! S+S- orbital pair (m,m') resolved projected bubble conductivity

   complex(dp), allocatable :: pi_orb_mp(:,:,:,:,:)
   ! pi_orb_mp(nboson, ndir, ndir, ndim_orb, ndim_orb) :
   ! S-S+ orbital pair (m,m') resolved projected bubble conductivity

   complex(dp), allocatable :: pi_projected_total(:,:,:)
   ! pi_projected_total(nboson, ndir, ndir) :
   ! Total projected bubble (sum over all alpha,delta)

   complex(dp), allocatable :: pi_projected_sc(:,:,:)
   ! pi_projected_sc(nboson, ndir, ndir) :
   ! Spin-conserving part of projected bubble

 end type optic_orbital_attrib_type

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

 ! Validate that mbandc is even for spin decomposition
 if (do_spin_decomp .and. mod(mbandc, 2) /= 0) then
   write(msg,'(a,i6,a)') &
   ' compute_bubble_conductivity: mbandc=', mbandc, ' is odd but nspinor=2.'
   ABI_ERROR(msg)
 end if

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
             ! Frequency index convention (1-based):
             !   iw indexes fermionic frequency omega_n, n = iw-1
             !   iom indexes bosonic frequency Omega_m, m = iom-1
             !   Shifted index: iw + iom - 1 gives omega_n + Omega_m
             !   Valid when iw + iom - 1 <= green_imp%nw
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

!!****f* m_dmft_optic_kernel/init_optic_orbital_attrib
!! NAME
!!  init_optic_orbital_attrib
!!
!! FUNCTION
!!  Initialize the orbital-resolved optical attribution structure.
!!
!! SOURCE

subroutine init_optic_orbital_attrib(orb_attrib, nboson, ndir, ndim_orb, nspinor)

 type(optic_orbital_attrib_type), intent(inout) :: orb_attrib
 integer, intent(in) :: nboson, ndir, ndim_orb, nspinor

!Local variables
 character(len=500) :: msg

! *********************************************************************

 orb_attrib%nboson = nboson
 orb_attrib%ndir = ndir
 orb_attrib%ndim_orb = ndim_orb
 orb_attrib%nspinor = nspinor

 write(msg,'(a,i4,a,i2,a,i4,a,i2)') &
 ' init_optic_orbital_attrib: nboson=', nboson, ' ndir=', ndir, &
 ' ndim_orb=', ndim_orb, ' nspinor=', nspinor
 call wrtout(std_out, msg)

 ABI_MALLOC(orb_attrib%pi_orb_pm, (nboson, ndir, ndir, ndim_orb, ndim_orb))
 ABI_MALLOC(orb_attrib%pi_orb_mp, (nboson, ndir, ndir, ndim_orb, ndim_orb))
 ABI_MALLOC(orb_attrib%pi_projected_total, (nboson, ndir, ndir))
 ABI_MALLOC(orb_attrib%pi_projected_sc, (nboson, ndir, ndir))

 orb_attrib%pi_orb_pm = czero
 orb_attrib%pi_orb_mp = czero
 orb_attrib%pi_projected_total = czero
 orb_attrib%pi_projected_sc = czero

end subroutine init_optic_orbital_attrib

!!***

!!****f* m_dmft_optic_kernel/destroy_optic_orbital_attrib
!! NAME
!!  destroy_optic_orbital_attrib
!!
!! FUNCTION
!!  Deallocate the orbital-resolved optical attribution structure.
!!
!! SOURCE

subroutine destroy_optic_orbital_attrib(orb_attrib)

 type(optic_orbital_attrib_type), intent(inout) :: orb_attrib

! *********************************************************************

 ABI_SFREE(orb_attrib%pi_orb_pm)
 ABI_SFREE(orb_attrib%pi_orb_mp)
 ABI_SFREE(orb_attrib%pi_projected_total)
 ABI_SFREE(orb_attrib%pi_projected_sc)

 orb_attrib%nboson = 0
 orb_attrib%ndir = 0
 orb_attrib%ndim_orb = 0
 orb_attrib%nspinor = 0

end subroutine destroy_optic_orbital_attrib

!!***

!!****f* m_dmft_optic_kernel/compute_optic_orbital_attrib
!! NAME
!!  compute_optic_orbital_attrib
!!
!! FUNCTION
!!  Compute orbital-pair resolved bubble optical conductivity by projecting
!!  velocity matrix elements and Green functions into the correlated orbital
!!  subspace using chipsi projectors.
!!
!!  The projection is:
!!    J_mu(k) = P(k) * j_mu(k) * P(k)^dag    [norb_corr x norb_corr]
!!    G^loc(k,iw) = P(k) * G^KS(k,iw) * P(k)^dag  [norb_corr x norb_corr]
!!
!!  where P_{alpha,a}(k) = chipsi projector from correlated orbital alpha
!!  to KS band a.
!!
!!  Then for each (alpha,delta) pair:
!!    contrib(alpha,delta) = -(wk/beta) * sum_n sum_{beta,gamma}
!!      J_mu(alpha,beta) * G^loc(beta,gamma;iw) * J_nu(gamma,delta) * G^loc(delta,alpha;iw+iOm)
!!
!!  This is computed as:
!!    prod(alpha,delta) = [J_mu * G^loc(iw) * J_nu](alpha,delta) * G^loc(delta,alpha;iw+iOm)
!!
!!  For the spin-flip channels:
!!    S+S-: alpha=(m,up) [1..ndim_orb], delta=(m',down) [ndim_orb+1..2*ndim_orb]
!!    S-S+: alpha=(m,down), delta=(m',up)
!!
!! INPUTS
!!  paw_dmft = DFT+DMFT data with psinablapsi_dmft and system parameters
!!  green_imp = converged Green function with oper(iw)%ks
!!  sproj = spinor projectors (from chipsi), must be populated for the atom
!!  nboson, niw_vertex = frequency counts
!!  nspinor = number of spinor components
!!  ndim_orb = number of orbital indices (2*lpawu+1)
!!
!! SIDE EFFECTS
!!  orb_attrib = on output, contains orbital-resolved attribution
!!
!! SOURCE

subroutine compute_optic_orbital_attrib(orb_attrib, paw_dmft, green_imp, sproj, &
  & nboson, niw_vertex, nspinor, ndim_orb, ladd)

 type(optic_orbital_attrib_type), intent(inout) :: orb_attrib
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type), intent(in) :: green_imp
 type(spinor_proj_type), intent(in) :: sproj
 integer, intent(in) :: nboson, niw_vertex, nspinor, ndim_orb
 logical, intent(in), optional :: ladd
   !! When .true., accumulate into existing values instead of zeroing first.

!Local variables
 character(len=500) :: msg
 integer :: iom, iw, ikpt, isppol, mu, nu
 integer :: ialpha, idelta, ispin_a, ispin_d, im_a, im_d
 integer :: mbandc, nkpt, norb_corr
 real(dp) :: beta, wk
 complex(dp) :: contrib
 logical :: ladd_local
 complex(dp), allocatable :: proj_k(:,:)
 complex(dp), allocatable :: jmu_ks(:,:), jnu_ks(:,:)
 complex(dp), allocatable :: jmu_loc(:,:), jnu_loc(:,:)
 complex(dp), allocatable :: gloc_iw(:,:), gloc_iw_om(:,:)
 complex(dp), allocatable :: temp_proj(:,:)
 complex(dp), allocatable :: prod1(:,:), prod2(:,:)

! *********************************************************************

 write(msg,'(a)') ' compute_optic_orbital_attrib: Computing orbital-resolved optical attribution'
 call wrtout(std_out, msg)

 ladd_local = .false.
 if (present(ladd)) ladd_local = ladd

 ! Check prerequisites
 if (paw_dmft%has_psinablapsi_dmft /= 1) then
   write(msg,'(3a)') &
   ' compute_optic_orbital_attrib: psinablapsi_dmft not computed.',ch10,&
   ' Orbital attribution remains unchanged.'
   call wrtout(std_out, msg)
   return
 end if

 if (green_imp%oper(1)%has_operks /= 1) then
   write(msg,'(a)') &
   ' compute_optic_orbital_attrib: Green function has no KS-basis data. Attribution remains unchanged.'
   call wrtout(std_out, msg)
   return
 end if

 if (nspinor /= 2) then
   write(msg,'(a)') &
   ' compute_optic_orbital_attrib: Orbital decomposition requires nspinor=2. Skipping.'
   call wrtout(std_out, msg)
   return
 end if

 mbandc = paw_dmft%mbandc
 nkpt = paw_dmft%nkpt
 norb_corr = nspinor * ndim_orb
 beta = one / paw_dmft%temp

 ! Validate frequency bounds
 if (niw_vertex + nboson - 1 > green_imp%nw) then
   write(msg,'(a,i6,a,i6,a,i6)') &
   ' compute_optic_orbital_attrib: frequency bounds exceeded. niw_vertex=', niw_vertex, &
   ' nboson=', nboson, ' but green has nw=', green_imp%nw
   ABI_ERROR(msg)
 end if

 ! Initialize or accumulate
 if (.not. ladd_local) then
   orb_attrib%pi_orb_pm = czero
   orb_attrib%pi_orb_mp = czero
   orb_attrib%pi_projected_total = czero
   orb_attrib%pi_projected_sc = czero
 end if

 ! Allocate working arrays
 ABI_MALLOC(proj_k, (norb_corr, mbandc))
 ABI_MALLOC(jmu_ks, (mbandc, mbandc))
 ABI_MALLOC(jnu_ks, (mbandc, mbandc))
 ABI_MALLOC(jmu_loc, (norb_corr, norb_corr))
 ABI_MALLOC(jnu_loc, (norb_corr, norb_corr))
 ABI_MALLOC(gloc_iw, (norb_corr, norb_corr))
 ABI_MALLOC(gloc_iw_om, (norb_corr, norb_corr))
 ABI_MALLOC(temp_proj, (norb_corr, mbandc))
 ABI_MALLOC(prod1, (norb_corr, norb_corr))
 ABI_MALLOC(prod2, (norb_corr, norb_corr))

 do isppol = 1, paw_dmft%nsppol
   do ikpt = 1, nkpt
     wk = paw_dmft%wtk(ikpt)

     ! Extract projector for this k-point: P(alpha, a)
     proj_k(:,:) = sproj%proj(:,:,ikpt)

     do iom = 1, nboson
       do mu = 1, 3
         ! Build j_mu in complex form and project to correlated space
         ! J_mu = P * j_mu_ks * P^dag
         do ialpha = 1, mbandc
           do idelta = 1, mbandc
             jmu_ks(ialpha, idelta) = &
               cmplx(paw_dmft%psinablapsi_dmft(1, mu, ialpha, idelta, ikpt, isppol), &
                     paw_dmft%psinablapsi_dmft(2, mu, ialpha, idelta, ikpt, isppol), kind=dp)
           end do
         end do
         call zgemm_wrapper(norb_corr, mbandc, mbandc, proj_k, jmu_ks, temp_proj)
         call zgemm_ct_wrapper(norb_corr, norb_corr, mbandc, temp_proj, proj_k, jmu_loc)

         do nu = 1, 3
           ! Build j_nu in complex form and project
           do ialpha = 1, mbandc
             do idelta = 1, mbandc
               jnu_ks(ialpha, idelta) = &
                 cmplx(paw_dmft%psinablapsi_dmft(1, nu, ialpha, idelta, ikpt, isppol), &
                       paw_dmft%psinablapsi_dmft(2, nu, ialpha, idelta, ikpt, isppol), kind=dp)
             end do
           end do
           call zgemm_wrapper(norb_corr, mbandc, mbandc, proj_k, jnu_ks, temp_proj)
           call zgemm_ct_wrapper(norb_corr, norb_corr, mbandc, temp_proj, proj_k, jnu_loc)

           ! Sum over fermionic Matsubara frequencies
           do iw = 1, niw_vertex
             if (iw + iom - 1 > green_imp%nw) cycle

             ! --- Project G^KS(k, iw) to local basis ---
             call zgemm_wrapper(norb_corr, mbandc, mbandc, proj_k, &
               green_imp%oper(iw)%ks(:,:,ikpt,isppol), temp_proj)
             call zgemm_ct_wrapper(norb_corr, norb_corr, mbandc, temp_proj, proj_k, gloc_iw)

             ! --- Project G^KS(k, iw+iOm) to local basis ---
             call zgemm_wrapper(norb_corr, mbandc, mbandc, proj_k, &
               green_imp%oper(iw+iom-1)%ks(:,:,ikpt,isppol), temp_proj)
             call zgemm_ct_wrapper(norb_corr, norb_corr, mbandc, temp_proj, proj_k, gloc_iw_om)

             ! --- Compute contrib(alpha,delta) ---
             ! prod1 = J_mu * G^loc(iw)
             call zgemm_nn_small(norb_corr, jmu_loc, gloc_iw, prod1)
             ! prod2 = prod1 * J_nu
             call zgemm_nn_small(norb_corr, prod1, jnu_loc, prod2)

             ! For each (alpha, delta):
             !   contrib = prod2(alpha,delta) * G^loc_iw_om(delta,alpha)
             do ialpha = 1, norb_corr
               do idelta = 1, norb_corr
                 contrib = -wk / beta * prod2(ialpha, idelta) * gloc_iw_om(idelta, ialpha)

                 ! Accumulate projected total
                 orb_attrib%pi_projected_total(iom, mu, nu) = &
                   orb_attrib%pi_projected_total(iom, mu, nu) + contrib

                 ! Classify spin channel
                 if (ialpha <= ndim_orb) then
                   ispin_a = 1  ! spin-up
                 else
                   ispin_a = 2  ! spin-down
                 end if
                 if (idelta <= ndim_orb) then
                   ispin_d = 1
                 else
                   ispin_d = 2
                 end if

                 if (ispin_a == ispin_d) then
                   ! Spin-conserving
                   orb_attrib%pi_projected_sc(iom, mu, nu) = &
                     orb_attrib%pi_projected_sc(iom, mu, nu) + contrib
                 else if (ispin_a == 1 .and. ispin_d == 2) then
                   ! S+S-: alpha=(m,up), delta=(m',down)
                   im_a = ialpha               ! orbital index m
                   im_d = idelta - ndim_orb    ! orbital index m'
                   orb_attrib%pi_orb_pm(iom, mu, nu, im_a, im_d) = &
                     orb_attrib%pi_orb_pm(iom, mu, nu, im_a, im_d) + contrib
                 else
                   ! S-S+: alpha=(m,down), delta=(m',up)
                   im_a = ialpha - ndim_orb    ! orbital index m
                   im_d = idelta               ! orbital index m'
                   orb_attrib%pi_orb_mp(iom, mu, nu, im_a, im_d) = &
                     orb_attrib%pi_orb_mp(iom, mu, nu, im_a, im_d) + contrib
                 end if
               end do ! idelta
             end do ! ialpha

           end do ! iw
         end do ! nu
       end do ! mu
     end do ! iom
   end do ! ikpt
 end do ! isppol

 ABI_FREE(proj_k)
 ABI_FREE(jmu_ks)
 ABI_FREE(jnu_ks)
 ABI_FREE(jmu_loc)
 ABI_FREE(jnu_loc)
 ABI_FREE(gloc_iw)
 ABI_FREE(gloc_iw_om)
 ABI_FREE(temp_proj)
 ABI_FREE(prod1)
 ABI_FREE(prod2)

 write(msg,'(a)') ' compute_optic_orbital_attrib: Orbital-resolved attribution complete.'
 call wrtout(std_out, msg)

end subroutine compute_optic_orbital_attrib

!!***

!!****f* m_dmft_optic_kernel/zgemm_wrapper
!! NAME
!!  zgemm_wrapper
!!
!! FUNCTION
!!  Simple wrapper for C = A * B (no transpose) for small matrices.
!!  Performs: C(m,n) = sum_k A(m,k) * B(k,n)
!!
!! SOURCE

subroutine zgemm_wrapper(m, n, k, A, B, C)

 integer, intent(in) :: m, n, k
 complex(dp), intent(in) :: A(m, k), B(k, n)
 complex(dp), intent(out) :: C(m, n)

!Local variables
 integer :: ii, jj, kk
 complex(dp) :: acc

! *********************************************************************

 do jj = 1, n
   do ii = 1, m
     acc = czero
     do kk = 1, k
       acc = acc + A(ii, kk) * B(kk, jj)
     end do
     C(ii, jj) = acc
   end do
 end do

end subroutine zgemm_wrapper

!!***

!!****f* m_dmft_optic_kernel/zgemm_ct_wrapper
!! NAME
!!  zgemm_ct_wrapper
!!
!! FUNCTION
!!  Wrapper for C = A * B^dag (conjugate transpose of B) for small matrices.
!!  Performs: C(m,n) = sum_k A(m,k) * conjg(B(n,k))
!!
!! SOURCE

subroutine zgemm_ct_wrapper(m, n, k, A, B, C)

 integer, intent(in) :: m, n, k
 complex(dp), intent(in) :: A(m, k), B(n, k)
 complex(dp), intent(out) :: C(m, n)

!Local variables
 integer :: ii, jj, kk
 complex(dp) :: acc

! *********************************************************************

 do jj = 1, n
   do ii = 1, m
     acc = czero
     do kk = 1, k
       acc = acc + A(ii, kk) * conjg(B(jj, kk))
     end do
     C(ii, jj) = acc
   end do
 end do

end subroutine zgemm_ct_wrapper

!!***

!!****f* m_dmft_optic_kernel/zgemm_nn_small
!! NAME
!!  zgemm_nn_small
!!
!! FUNCTION
!!  C = A * B for small square matrices of dimension n.
!!
!! SOURCE

subroutine zgemm_nn_small(n, A, B, C)

 integer, intent(in) :: n
 complex(dp), intent(in) :: A(n, n), B(n, n)
 complex(dp), intent(out) :: C(n, n)

!Local variables
 integer :: ii, jj, kk
 complex(dp) :: acc

! *********************************************************************

 do jj = 1, n
   do ii = 1, n
     acc = czero
     do kk = 1, n
       acc = acc + A(ii, kk) * B(kk, jj)
     end do
     C(ii, jj) = acc
   end do
 end do

end subroutine zgemm_nn_small

!!***

!!****f* m_dmft_optic_kernel/write_optic_orbital_attrib
!! NAME
!!  write_optic_orbital_attrib
!!
!! FUNCTION
!!  Write orbital-resolved optical conductivity attribution to file.
!!
!!  Output sections:
!!    Section 1: Coverage fraction — fraction of Pi_bubble captured by projected Pi
!!    Section 2: Orbital-pair-resolved S+S- contribution at iOm=0 for diagonal mu=nu
!!    Section 3: Orbital-pair-resolved S-S+ contribution at iOm=0 for diagonal mu=nu
!!    Section 4: Dominant (m,m') orbital pairs ranked by |Pi^{S+S-}_{mm'}| at iOm=0
!!    Section 5: Frequency dependence of dominant orbital pairs for S+S-
!!
!! INPUTS
!!  orb_attrib = orbital attribution data
!!  pi_bubble_total(nboson,ndir,ndir) = full KS-space bubble for coverage fraction
!!  fname = output file name
!!  beta = inverse temperature
!!
!! SOURCE

subroutine write_optic_orbital_attrib(orb_attrib, pi_bubble_total, fname, beta)

 type(optic_orbital_attrib_type), intent(in) :: orb_attrib
 complex(dp), intent(in) :: pi_bubble_total(:,:,:)
 character(len=*), intent(in) :: fname
 real(dp), intent(in) :: beta

!Local variables
 integer :: unt, iom, mu, im, imp, ios, ipair, npairs, irank
 real(dp) :: omega_boson, abs_proj, abs_full, coverage
 real(dp) :: abs_pm_total, abs_mp_total
 character(len=500) :: msg
 character(len=3) :: dir_label(3)
 ! For ranking
 integer :: max_rank
 real(dp), allocatable :: pair_mag(:)
 integer, allocatable :: pair_im(:), pair_imp(:)
 real(dp) :: mag_tmp
 integer :: itmp

! *********************************************************************

 dir_label = (/ ' xx', ' yy', ' zz' /)
 max_rank = min(20, orb_attrib%ndim_orb * orb_attrib%ndim_orb)
 npairs = orb_attrib%ndim_orb * orb_attrib%ndim_orb

 open(newunit=unt, file=fname, form='formatted', action='write', iostat=ios)
 if (ios /= 0) then
   write(msg,'(3a)') 'Cannot open file: ', trim(fname), ' for writing.'
   ABI_ERROR(msg)
 end if

 write(unt,'(a)') '# DFT+DMFT Orbital-Resolved Optical Conductivity Attribution'
 write(unt,'(a,i6)') '# nboson = ', orb_attrib%nboson
 write(unt,'(a,i4)') '# ndim_orb = ', orb_attrib%ndim_orb
 write(unt,'(a,i2)') '# nspinor = ', orb_attrib%nspinor
 write(unt,'(a,es14.6)') '# beta = ', beta
 write(unt,'(a)') '#'

 ! === Section 1: Coverage fraction ===
 write(unt,'(a)') '# Section 1: Coverage fraction of projected vs full bubble'
 write(unt,'(a)') '# Shows what fraction of the total optical conductivity comes from'
 write(unt,'(a)') '# correlated orbital transitions (captured by chipsi projection).'
 write(unt,'(a)') '# mu  |Pi_proj_total(0)|  |Pi_full_bubble(0)|  coverage_fraction'

 do mu = 1, orb_attrib%ndir
   abs_proj = abs(orb_attrib%pi_projected_total(1, mu, mu))
   abs_full = abs(pi_bubble_total(1, mu, mu))
   if (abs_full > tol16) then
     coverage = abs_proj / abs_full
   else
     coverage = zero
   end if
   write(unt,'(i3,a3,2es20.10,f12.6)') mu, dir_label(mu), abs_proj, abs_full, coverage
 end do

 ! === Section 2: Orbital-pair-resolved S+S- at iOm=0 ===
 write(unt,'(a)') '#'
 write(unt,'(a)') '# Section 2: Orbital-pair S+S- contribution to Pi_mu_mu at iOm=0'
 write(unt,'(a)') '# For each diagonal direction mu, shows |Pi^{S+S-}_{mm''}(0)| per orbital pair.'
 write(unt,'(a)') '# mu  m  m''  Re(Pi^{S+S-})  Im(Pi^{S+S-})  |Pi^{S+S-}|'

 do mu = 1, orb_attrib%ndir
   do im = 1, orb_attrib%ndim_orb
     do imp = 1, orb_attrib%ndim_orb
       write(unt,'(i3,a3,2i4,3es20.10)') mu, dir_label(mu), im, imp, &
         real(orb_attrib%pi_orb_pm(1, mu, mu, im, imp)), &
         aimag(orb_attrib%pi_orb_pm(1, mu, mu, im, imp)), &
         abs(orb_attrib%pi_orb_pm(1, mu, mu, im, imp))
     end do
   end do
 end do

 ! === Section 3: Orbital-pair-resolved S-S+ at iOm=0 ===
 write(unt,'(a)') '#'
 write(unt,'(a)') '# Section 3: Orbital-pair S-S+ contribution to Pi_mu_mu at iOm=0'
 write(unt,'(a)') '# mu  m  m''  Re(Pi^{S-S+})  Im(Pi^{S-S+})  |Pi^{S-S+}|'

 do mu = 1, orb_attrib%ndir
   do im = 1, orb_attrib%ndim_orb
     do imp = 1, orb_attrib%ndim_orb
       write(unt,'(i3,a3,2i4,3es20.10)') mu, dir_label(mu), im, imp, &
         real(orb_attrib%pi_orb_mp(1, mu, mu, im, imp)), &
         aimag(orb_attrib%pi_orb_mp(1, mu, mu, im, imp)), &
         abs(orb_attrib%pi_orb_mp(1, mu, mu, im, imp))
     end do
   end do
 end do

 ! === Section 4: Dominant S+S- orbital pairs ranked by magnitude ===
 write(unt,'(a)') '#'
 write(unt,'(a)') '# Section 4: Dominant S+S- orbital pairs at iOm=0 (ranked by |Pi^{S+S-}_{mm''}|)'
 write(unt,'(a)') '# Averaged over diagonal directions (xx, yy, zz).'
 write(unt,'(a)') '# rank  m  m''  |Pi^{S+S-}|_avg  fraction_of_total_S+S-'

 ABI_MALLOC(pair_mag, (npairs))
 ABI_MALLOC(pair_im, (npairs))
 ABI_MALLOC(pair_imp, (npairs))

 ! Compute direction-averaged magnitudes
 abs_pm_total = zero
 ipair = 0
 do im = 1, orb_attrib%ndim_orb
   do imp = 1, orb_attrib%ndim_orb
     ipair = ipair + 1
     pair_im(ipair) = im
     pair_imp(ipair) = imp
     pair_mag(ipair) = zero
     do mu = 1, orb_attrib%ndir
       pair_mag(ipair) = pair_mag(ipair) + abs(orb_attrib%pi_orb_pm(1, mu, mu, im, imp))
     end do
     pair_mag(ipair) = pair_mag(ipair) / dble(orb_attrib%ndir)
     abs_pm_total = abs_pm_total + pair_mag(ipair)
   end do
 end do

 ! Simple selection sort for top max_rank pairs
 do irank = 1, min(max_rank, npairs)
   do ipair = irank + 1, npairs
     if (pair_mag(ipair) > pair_mag(irank)) then
       mag_tmp = pair_mag(irank); pair_mag(irank) = pair_mag(ipair); pair_mag(ipair) = mag_tmp
       itmp = pair_im(irank); pair_im(irank) = pair_im(ipair); pair_im(ipair) = itmp
       itmp = pair_imp(irank); pair_imp(irank) = pair_imp(ipair); pair_imp(ipair) = itmp
     end if
   end do
   if (abs_pm_total > tol16) then
     write(unt,'(i4,2i4,es20.10,f12.6)') irank, pair_im(irank), pair_imp(irank), &
       pair_mag(irank), pair_mag(irank) / abs_pm_total
   else
     write(unt,'(i4,2i4,es20.10,a12)') irank, pair_im(irank), pair_imp(irank), &
       pair_mag(irank), '  undef'
   end if
 end do

 ABI_FREE(pair_mag)
 ABI_FREE(pair_im)
 ABI_FREE(pair_imp)

 ! === Section 5: Frequency dependence of spin-flip orbital pairs ===
 write(unt,'(a)') '#'
 write(unt,'(a)') '# Section 5: Frequency dependence of S+S- spin-flip by orbital pair'
 write(unt,'(a)') '# For diagonal mu=1 (xx), shows |Pi^{S+S-}_{mm''}(iOm)| vs bosonic frequency.'
 write(unt,'(a)') '# iOm  Omega_boson  m  m''  |Pi^{S+S-}|'

 do iom = 1, orb_attrib%nboson
   omega_boson = two_pi * dble(iom - 1) / beta
   do im = 1, orb_attrib%ndim_orb
     do imp = 1, orb_attrib%ndim_orb
       write(unt,'(i6,es14.6,2i4,es20.10)') iom, omega_boson, im, imp, &
         abs(orb_attrib%pi_orb_pm(iom, 1, 1, im, imp))
     end do
   end do
 end do

 close(unt)

 write(msg,'(3a)') ' write_optic_orbital_attrib: Written to ', trim(fname)
 call wrtout(std_out, msg)

end subroutine write_optic_orbital_attrib

!!***

END MODULE m_dmft_optic_kernel
