!
! Copyright 1996-2017 the Authors
!
! Licensed under the EUPL, Version 1.1 only (the "Licence");
!
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!   https://joinup.ec.europa.eu/software/page/eupl
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an
! "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
! either express or implied.
!
! See the Licence for the specific language governing permissions
! and limitations under the Licence.


! additional subroutines for KPP version
! --------------------------------------

! Author of this file:
! --------------------
!    Roland von Glasow

! Modifications :
! -------------
!

! =======================================================================
! =======================================================================

subroutine initc (box,n_bl)

! Description :
! -----------
  ! initialisation of chemistry module

! Modifications :
! -------------
  ! jjb: missing initialisations

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       cmechdir,     &
       iaertyp,      &
       iod,          &
       lpJoyce14bc,  &
       neula

  USE constants, ONLY : &
! Imported Parameters:
       Avogadro, &
       m_air

  USE file_unit, ONLY : &
       jpfunerr,   &
       jpfuneul,   &
       jpfunout,   &
       jpfunprofc, &
       jpfunsg1,   &
       jpfunsr1

  USE gas_common, ONLY : &
! Imported Parameters:
       j1, &
       j5, &
! Imported Array Variables with intent(in):
       ind_gas, ind_gas_rev, &
       gas_is_halo, &
       gas_name, &
! Imported Array Variables with intent(inout):
       s1, &
       s1_init_grd, &
       s1_init_top, &
       es1, &
! Imported Array Variables with intent (out):
       s3, &
       vg, &
       nadvmax, nindadv, xadv ! for Eulerian advection

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j3, &
       j6, &
       nf, &
       n, &
       nka, &
       nkt, &
       nkc, &
       nlev, &
       nrxn

  USE kpp_aer_Parameters, ONLY : & ! jjb rename NSPEC so that variables related to both mechanisms can be declared below
       nspec_a=>NSPEC
  USE kpp_tot_Parameters, ONLY : &
       nspec_t=>NSPEC

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  logical, intent(in) :: box
  integer, intent(in) :: n_bl

! Local scalars:
  integer :: ia,j,k,kc
  logical :: tr_ue
  real (kind=dp) :: x0
  ! sea salt initialisation
  real (kind=dp) :: xso42m,xhco3m,xno3m,xbrm,xclm,xim,xio3m,xiod

! Local arrays:
  real (kind=dp) :: &
       freep(nf), &
       is4(n,3), &
       xm(n), &
       x4(n), &
       x2(j1)

! Common blocks:
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3
  common /blck11/ rc(nkc,n)
  real (kind=dp) :: rc
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm
  common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
  real (kind=dp) :: conv2
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /blck78/ sa1(j2,nka)
  real (kind=dp) :: sa1
  common /budg/ bg(2,nrxn,nlev),il(nlev)
  real (kind=dp) :: bg
  integer :: il
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb63/ fcs(nka),xmol3(nka)
  real (kind=dp) :: fcs, xmol3
  common /kinv_i/ kinv
  integer :: kinv
  common /kpp_l1/ cloudt(nkc,n)
  logical :: cloudt
  common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss
  real (kind=dp) :: xcryssulf,xcrysss,xdelisulf,xdeliss
  common /kpp_laer/ henry_la(NSPEC_a,nf),xkmt_la(nf,nkc,NSPEC_a), &
       xkef_la(nf,nkc,NSPEC_a),xkeb_la(nf,nkc,NSPEC_a)
  real (kind=dp) :: henry_la, xkmt_la, xkef_la, xkeb_la
  common /kpp_ltot/ henry_lt(NSPEC_t,nf),xkmt_lt(nf,nkc,NSPEC_t), &
       xkef_lt(nf,nkc,NSPEC_t),xkeb_lt(nf,nkc,NSPEC_t)
  real (kind=dp) :: henry_lt, xkmt_lt, xkef_lt, xkeb_lt

! == End of declarations =======================================================


! print input concentrations and emission rates (user values):
! ------------------------------------------------------------
! initial mixing ratio of gas phase species in nmol mol-1 (=ppb)
  ! mixing ratio at ground
  write (jpfunprofc,6010)
6010 format (6x,'initial gas concentration at the surface [ppb]')
  write (jpfunprofc,6020) (s1_init_grd(j),j=1,j1)
  ! mixing ratio at top
  write (jpfunprofc,6012)
6012 format (6x,'initial gas concentration at the top [ppb]')
  write (jpfunprofc,6020) (s1_init_top(j),j=1,j1)
  ! emission rates of gas phase species in molecules/cm**2/s
  write (jpfunprofc,6025)
6025 format (6x,'emission rates [molecules/cm**2/s]')
  write (jpfunprofc,6020) (es1(j),j=1,j1)
6020 format (1x,10e12.5)


! Initialisation to 0 of some arrays, which are not necessarily updated for layers nf+1 to n
  cw(:,:) = 0._dp
  rc(:,:) = 0._dp
  cm(:,:) = 0._dp
  conv2(:,:) = 0._dp

! initiate all ions with marginal concentration to avoid computational problems
  sion1(:,:,:) = 0._dp
  sl1(:,:,:) = 0._dp

! jjb 25-10-2017
! initialise vg so that the first call to sedc/sedc_box see defined values
! It might be necessary to check the order the subroutines are called, and/or
! the calls hereafter: gasdrydep should probably be included as well
! but need to check
  vg(:) = 0._dp

! jjb 14/02/2017
!     initialise the reaction rate arrays
  henry_la(:,:) = 0._dp
  henry_lt(:,:) = 0._dp
  xkmt_la(:,:,:) = 0._dp
  xkmt_lt(:,:,:) = 0._dp
  xkef_la(:,:,:) = 0._dp
  xkef_lt(:,:,:) = 0._dp
  xkeb_la(:,:,:) = 0._dp
  xkeb_lt(:,:,:) = 0._dp

! conversion of gaseous species and air density
! air density: [rho]=kg/m^3
  do k=1,n
     cm3(k) = rho(k) * Avogadro / m_air * 1e-6   ! [air] in mlc/cm^3
     am3(k) = rho(k) / m_air                     ! [air] in mol/m^3
! conversion of gaseous species in ppb to mol/m**3(air)
     xm(k) = am3(k) * 1.e-9_dp                   ! ppb --> mol/m^3

! exp. decrease of concentrations from surface to free troposphere
     x4(k) = eta(k) / 1900._dp     ! jjb this may need to be adjusted
     x4(k) = min(1._dp, x4(k))
  end do

! Initialize arrays(j)
  do j=1,j1
     if(s1_init_grd(j).gt.0._dp) then
        ! avoid log(0) by adding a small value to s1_init_top
        x2(j) = -log(s1_init_grd(j)) + log(s1_init_top(j) + 1.e-10_dp)
     else
        x2(j) = 0._dp
! The interpolation method does not allow to calculate a gradient
        ! (see below: s1(:,:) = s1_init_grd(:)*... )
        ! with a top concentration > 0 and a ground concentration = 0
        ! Warn the user if this case arise
        !  jjb: this could be solved by simply adding a small value to groud conc=0, same as top conc above
        if(s1_init_top(j).gt.0._dp) then
           write (jpfunout,*)"Warning with gas species nb. ",ind_gas(j)
           write (jpfunout,*)"  Its top concentration is > 0 while its ground"
           write (jpfunout,*)"  concentration is = 0"
           write (jpfunout,*)"  The interpolation method does not allow this"
           write (jpfunout,*)"  See SR initc"
        end if
     end if
  end do

! ......................................................................
! Initialise gas concentration in the whole column
  s1(:,1) = 0._dp
  do k=2,n
     do j=1,j1

        ! halogen only in BL, but there no gradient
        if(gas_is_halo(j) .and. k<kinv .and. k>2 .and.trim(gas_name(j))/='HCl') then
           s1(j,k) = s1(j,k-1)
        else if(gas_is_halo(j) .and. k>=kinv .and. trim(gas_name(j))/='HCl' ) then
           s1(j,k) = 0._dp

        ! general case
        else
           s1(j,k) = s1_init_grd(j) * exp(x4(k) * x2(j)) * xm(k)
        end if
     end do
  end do
! ......................................................................

! initial radical concentrations in mol/m**3(air)
  do k=1,n
     do j=1,j5
        s3(j,k) = 0._dp
     end do
  end do
  xiod = 0._dp
  if (iod) xiod = 1._dp

! Eulerian configuration: read input file
  if (neula.eq.0) then
     open (jpfuneul,file=trim(cmechdir)//'euler_in.dat',status='old')
     ! First line = number of advected species
     read (jpfuneul,'(i3)') nadvmax
     allocate (nindadv(nadvmax), xadv(nadvmax))
     do j=1,nadvmax
        read (jpfuneul,5100) nindadv(j),xadv(j)
        ! Gas index = 0 is ignored
        ! If not zero, it must exist in the user list of gases, check:
        if (nindadv(j) /= 0 .and. ind_gas_rev(nindadv(j)) == 0) then
           write (jpfunerr,*) 'Error when reading euler_in.dat'
           write (jpfunerr,5101) '  gas index ',nindadv(j),' not found in gas list'
           write (jpfunerr,*) '  it will be ignored'
           nindadv(j) = 0 ! the gas will be ignored
        end if
     enddo
     close (jpfuneul)
  endif
5100 format (i3,d9.2)
5101 format (a,i3,a)


! define crystallization and deliquescene rel humidities (Seinfeld and Pandis,
! Fig 9.4, p. 519)
! usually several cloud cycles have been made therefore sulfate aerosol is
! usually shrinking; for sea salt crystallization point is used as well, because
! sea salt particle were produced as droplets and shrank, but did not get "dry"
! if rel hum is below crys rH in FT, and it's getting more humid then the deli rH
! has to be taken for reactivation of aerosol chemistry
  xcryssulf = 0.4_dp  ! crystallization humities
  xcrysss   = 0.42_dp
  xdelisulf = 0.7_dp  ! assumed based on mixing between different salts - should be
                      ! calculated explicitly if this starts to be critical (ie for
                      ! non-marine cases)
  xdeliss   = 0.75_dp
  write (jpfunout,*)'deliquescence rH  ',xdelisulf,xdeliss
  write (jpfunout,*)'crystallization rH',xcryssulf,xcrysss

! array cloudt: was there a cloud in this layer in previous timestep?
! initialize with "true" to assure that aerosol chemistry is on in FT for
! layers in which rH > xcrystallization
  do k=1,n
     do kc=1,nkc
        cloudt(kc,k) = .true.
     enddo
  enddo

! initial loading of aerosols with nh3,so4,fe(3),mn(2) (x0=mole/particle)
! watch out: sa1 is defined as sa1(j2,..) but addressed in j6 (=ion, sion1) terms
!            except for DOM which is in sl1 (therefore it is in j2)!!
  sa1(:,:) = 0._dp
  do ia=1,nka
     x0 = en(ia) * 1.e-3_dp * fcs(ia) / xmol3(ia)

!! ocean aerosol: particles with rn(ia)<.5 mum: 32% (NH4)2SO4, 64% NH4HSO4, 4% NH4NO3
! ocean aerosol: particles with rn(ia)<.5 mum: 34% (NH4)2SO4, 65.6% NH4HSO4, 0.4% NH4NO3
     if (iaertyp == 3) then
        if (rn(ia).lt.0.5_dp) then
           sa1(2,ia)  = x0 * 1.34_dp   !NH4+
           sa1(8,ia)  = x0 * 0.34_dp   !SO4=
!           sa1(13,ia) = x0 * 0.04_dp   !NO3-
           sa1(13,ia) = x0 * 0.004_dp  !NO3-
           sa1(19,ia) = x0 * 0.656_dp  !HSO4-
! larger particles: pure nacl
        else
! sea salt particle
! x0 = mol / particle
! all the xiii are scaled to the sum of all negative ions in seawater,
! Na+ is the sum of all positive ions; to get the correct molar ratios of
! Cl- or Br- to Na+, the lumped Na+ has to be multiplied by 0.806
           xso42m = 0.0485_dp
           xhco3m = 4.2e-3_dp
           xno3m  = 1.0e-7_dp
           xbrm   = 1.45e-3_dp
           xim    = 7.4e-8_dp / .545_dp * xiod
           xio3m  = 2.64e-7_dp / .545_dp * xiod
           xclm   = 1._dp - (xso42m+xhco3m+xno3m+xbrm+xim+xio3m)
           sa1( 8,ia) = xso42m * x0 ! SO4=
           sa1( 9,ia) = xhco3m * x0 ! HCO3-
           sa1(13,ia) = xno3m * x0  ! NO3-
           sa1(14,ia) = xclm * x0   ! Cl-
           sa1(20,ia) = x0          ! "Na+" -  eletronegativity
           sa1(24,ia) = xbrm * x0   ! Br-
           sa1(34,ia) = xim * x0    ! I-
           sa1(36,ia) = xio3m * x0  ! IO3-
           sa1(j2-j3+4,ia) = 0.27_dp * xbrm * x0 ! unspecified DOM
                                                 ! according to #2210: 0.27*[Br-]; enriched compared to ocean water ratio
        endif

     else if (iaertyp == 1) then
!! urban aerosol:
!! 2/3*xm = mole mass nh4no3; 1/3*xm = mole mass (nh4)2so4;
!! xm=en(ia)*fcs(ia)*1d-3: total soluble mole mass
!! --> 4/3*xm mole nh3, 2/3*xm mole no3, 1/3*xm mole so4
!        if (xmol3(ia).lt.130.) then
!           x0=en(ia)*1.d-03 * fcs(ia)/(3. * xmol3(ia))
!           sa1(3,ia)=x0 * 2.
!           sa1(4,ia)=x0 * 4.
!           sa1(6,ia)=x0
!        else

        ! Joyce et al 2014 case
        if (lpJoyce14bc) then
! soluble part of urban aerosol: pure H2SO4   PJ
! x0 = mol / particle
! molar ratios  "aeroPJ"  SR inic
           sa1(1,ia)  = x0 * 0.1868 * 2._dp !H+       PJ  (2*SO4=)
           sa1(2,ia)  = x0 * 0._dp          !NH4+     PJ
           sa1(8,ia)  = x0 * 0.1868_dp      !SO4=     PJ
           sa1(13,ia) = x0 * 0._dp          !NO3-     PJ
           if(rn(ia).le.0.5_dp)then         ! Cl-
!              sa1(14,ia)=x0 * 0.0356_dp     !Cl-      PJ (cl_b23c)
              sa1(14,ia) = x0 * 0.0227_dp   !Cl-      PJ (b25)
           else
              sa1(14,ia) = x0 * 0._dp       ! No Cl- in supermicron
           end if
           sa1(19,ia) = x0 * 0._dp          !HSO4-    PJ
!           sa1(j2-j3+4,ia) = x0 * 0.5757_dp !DOM      PJ  (obs DOM)
!           sa1(j2-j3+4,ia) = x0 * 0.7763_dp !DOM      PJ  (obs total PM2.5)
           sa1(j2-j3+4,ia) = x0 * 0.6642_dp !DOM      PJ  (b25)
        end if              ! Joyce
     end if                 ! iaertyp
  end do                    !ia


! print initial concentrations (continued)
  write (jpfunprofc,6030)
6030 format (6x,'sa1(4,nka)')
  write (jpfunprofc,6020) (sa1(4,ia),ia=1,nka)
  write (jpfunprofc,6040)
6040 format (6x,'sa1(6,nka)')
  write (jpfunprofc,6020) (sa1(6,ia),ia=1,nka)
  write (jpfunprofc,6050)
! 6050 format (6x,'sa1(j2-j3+4,nka)')
!  write (jpfunprofc,6020) (sa1(j2-j3+4,ia),ia=1,nka)
!  write (jpfunprofc,6060)
! 6060 format (6x,'sa1(j2-j3+5,nka)')
!  write (jpfunprofc,6020) (sa1(j2-j3+5,ia),ia=1,nka)
6050 format (6x,'sa1(14,nka)')
  write (jpfunprofc,6020) (sa1(14,ia),ia=1,nka)
  write (jpfunprofc,6060)
6060 format (6x,'sa1(24,nka)')
  write (jpfunprofc,6020) (sa1(24,ia),ia=1,nka)

! levels for  rate output
  il(1) =  5
  if (box) il(1)=n_bl
  il(2) = 15
  il(3) = 25
  il(4) = 35
  il(5) = 45
  il(6) = 55
  il(7) = 65
  il(8) = 75
  il(9) = 85
  il(10)= 95
  il(11)=105
  il(12)=115
  il(13)=125
  il(14)=135
  il(15)=145

! initial output for plotting
  do k=1,n
     is4(k,1) = am3(k) ! stay consistent with plot routine array size!
     is4(k,2) = 0._dp
     is4(k,3) = 0._dp
!     write (542,12) k,am3(k,1),am3(k,2)
!     write (543,12) k,cm3(k,1),cm3(k,2)
  enddo
! 12   format (i3,2d16.8)
  write (jpfunsg1) is4
  close (jpfunsg1)
  write (jpfunsr1) is4
  close (jpfunsr1)

! vertical grid analysis: find the layers to compute u10
  call aer_source_init

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!! Initialise vmean (constant factor calculation)
!      call v_mean_init
!! ... and make the first calculation
!      call v_mean (t(:nmax_chem_aer))

! Init calc of v_mean and henry
! initialize all variables also for box run
  call v_mean_a  (t,nf)
  call henry_a (t,nf)
  call st_coeff_a
  call v_mean_t  (t,nf)
  call henry_t (t,nf)
  call st_coeff_t

! initialize all levels also in box run
! free path length (lambda=freep):
  do k=1,nf
     freep(k) = 2.28e-5_dp * t(k) / p(k)
  enddo
  tr_ue = .false. !.true.
  call init_konc
  call fast_k_mt_a(freep,tr_ue,nf)
  call activ_init

end subroutine initc

!
!-----------------------------------------------------
!

subroutine liq_parm (xra,box,n_bl)

! Description :
! -----------
  ! parameter for liquid phase chemistry
  ! LWC, henry, k_mt needed for calculation of gas <--> liquid transfer
  ! (see cb kpp_l)
  ! nkc=4; 1: sulfate aerosol, 2: seasalt aerosol, 3: sulfate droplets
  !        4: seasalt droplets
  ! xra: aerodynamic resistence, needed for calculation of dry deposition velocities

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE constants, ONLY : &
! Imported Parameters:
       Avogadro, &
       m_air

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: xra ! aerodynamic resistence
  logical, intent(in) :: box
  integer, intent(in) :: n_bl

! Local scalars:
  integer :: iph3, iph4, itime, k, kc, nmin, nmin2, nmaxf, nmax
  logical :: update_now
! Local arrays:
  real (kind=dp) :: freep(nf)

! Common blocks:
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real(kind=dp) :: cw, cm
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /kpp_l1/ cloudt(nkc,n)
  logical :: cloudt

! == End of declarations =======================================================

  itime = int(time) ! jjb bugfix: previously lmin was used, but lmin doesn't change during
                    !   6 dd dub-timesteps... thus fast_k_mt_* routines were called 6 times
                    !   during a minute, then not during the next minute...
                    !   Now called once every 2 minutes, as expected.
  nmin  = 1
  nmin2 = 2
  nmaxf = nf
  nmax  = n
  if (box) then
     nmin  = n_bl
     nmaxf = n_bl
     nmax  = n_bl
  endif

! free path length (lambda=freep):
  do k=nmin,nmaxf
     freep(k) = 2.28e-5_dp * t(k) / p(k)
  enddo

! conversion of gaseous species and air density
! air density: [rho]=kg/m^3
  cm3(1:nmax) = rho(1:nmax) * Avogadro / m_air * 1e-6_dp ! [air] in mlc/cm^3
  am3(1:nmax) = rho(1:nmax) / m_air                      ! [air] in mol/m^3


! dry deposition velocities for gas phase
  call gasdrydep (xra,t,rho,freep)

  call cw_rc (nmaxf)

! Check if the tot mechanism will be called or not (iph3/4 = 1) and if a new bin is activated somewhere (update_now)
  iph3 = 0
  iph4 = 0
  update_now = .false.

  do k=nmin2,nmaxf
     do kc=1,nkc
        if (cm(kc,k).gt.0._dp) then
           if (.not.cloudt(kc,k)) update_now = .true.
           if (kc.eq.3) iph3 = 1
           if (kc.eq.4) iph4 = 1
        endif
     enddo
  enddo

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!! Compute the mean molecular speed (depends only on the temperature)
!    call v_mean (t(:nmax_chem_aer))

! Call all subroutines needed for aerosols (bins 1 & 2) (aer mechanism)
  call v_mean_a (t,nmaxf)      ! T varies slowly with time
  call henry_a (t,nmaxf)
  call st_coeff_a
  call equil_co_a (t,nmaxf)
  if (itime/120*120.eq.itime .or. update_now ) call fast_k_mt_a(freep,box,n_bl)

! If necessary, call all subroutines needed for droplets (bins 3 & 4) (tot mechanism)
  if (iph3.eq.1 .or. iph4.eq.1) then
     call v_mean_t (t,nmaxf)  ! T varies slowly with time
     call henry_t (t,nmaxf)
     call st_coeff_t
     call equil_co_t (t,nmaxf)
     if(itime/120*120.eq.itime .or. update_now) call fast_k_mt_t(freep,box,n_bl)
  endif

! calculate rate for surface reaction OH + Cl-
!   call gamma_surf (box,n_bl) ! jjb not used

! calculate rates for DRY heterogeneous reactions
  call dry_cw_rc (nmax)
  call dry_rates_g (t,p,nmax)
  call dry_rates_a (freep,nmaxf)
  call dry_rates_t (freep,nmaxf)

  call activ (box,n_bl)

end subroutine liq_parm


!
!--------------------------------------------------------------------------------
!

subroutine st_coeff_t

! Description :
! -----------
  ! sticking coefficients, needed for calculation of k_mt

  ! Reference: for instance Davidovits et al, Chem. Rev. 2006, 106, 1323-1354

! Modifications :
! -------------
  ! jjb: factorise some repetitive calculations

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       lpJoyce14bc

  USE constants, ONLY : &
! Imported Parameters:
       cal => cal15, &   ! calorie      [J]
       R => gas_const    ! gas constant [J/(mol*K)]

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'tot_Parameters.h' !additional common blocks and other definitions

  real (kind=dp), parameter :: CoR = cal/R

! Local scalars:
  integer :: j, k
  real (kind=dp) :: CoRT, RT, tcorr, zexp2
  real (kind=dp), external :: a_n2o5

! Common blocks:
  common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
  real (kind=dp) :: alpha, vmean
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho

! == End of declarations =======================================================


! Initialisation (default value = 0.1)
  alpha(:,:) = 0.1_dp

  do k=2,nf

     ! Compute a few factors used several times
     tcorr = 1._dp / t(k) - 1._dp / 298.15_dp
     RT    = R * t(k)
     CoRT  = cal/RT
     zexp2 = exp(2000._dp * tcorr)

! for some species a temperature dependence is calculated/estimated
! using:
! alpha = 1./(1.+1./(1./(1./alpha(T0)-1.)*EXP((-DeltaH/RGAS)*TCORR)))

! standard value
     alpha(ind_H2SO4,k) = 0.65_dp
     !alpha(ind_CH4,k)   = 0.1_dp ! default value
     !alpha(ind_HNO4,k)  = 0.1_dp
! experimentally determined values (from MOCCA)
     alpha(ind_O3P,k) = 1.0e-6_dp  ! JPL #15
     alpha(ind_O1D,k) = 1.0e-6_dp  ! JPL #15
     alpha(ind_O3,k)  = 2.0e-3_dp
     alpha(ind_O2,k)  = 1./(1._dp +1./(1./(1./1.0d-2 - 1._dp)*zexp2)) ! 06.04.00
!     alpha(ind_OH,k) = 4.0D-03
     alpha(ind_OH,k) = 1.0e-2_dp
     alpha(ind_HO2,k) = 2.0e-1_dp
!     alpha(ind_H2O2,k) =  9.1D-02
     alpha(ind_H2O2,k) = 1./(exp(-26.d3/RT+107.8456/R) + 1._dp)
     alpha(ind_NO,k) =  5.0D-05
     alpha(ind_NO2,k) =  1.5D-03

     if (lpJoyce14bc) then
        alpha(ind_NO3,k) =  2.5D-03 ! jjb: lower limit Thomas et al. (1989), Mihelcic et al. (1993) (cited by Finlayson Pitt & Pitt, Chemistry of the Upper and Lower Atmosphere, 2000
     else
!        alpha(ind_NO3,k) = 2.5D-03
        alpha(ind_NO3,k) = 4.0D-02
     end if
     if (lpJoyce14bc) then
        alpha(ind_N2O5,k) = a_n2o5(k,1)
     else
!        alpha(ind_N2O5,k) = 2.7D-02
!        alpha(ind_N2O5,k) = 1.0D-01 ! default value
     end if

     alpha(ind_HONO,k) = 4.0D-02
!     alpha(ind_HNO3,k) = 8.6D-02
     alpha(ind_HNO3,k) = 5.0D-01
     alpha(ind_NH3,k) = 6.0D-02
     alpha(ind_MO2,k) = 1./(1._dp +1./(1./(1./1.0d-2-1._dp)*zexp2))   ! 06.04.00
!     alpha(ind_ROOH,k) = 5.5D-03
!     alpha(ind_ROOH,k) = 0.01
     alpha(ind_ROOH,k)= 1./(exp(-6.5D3*CoRT+32.5*CoR)+1._dp)
     alpha(ind_HCHO,k) = 4.0D-02
!     alpha(ind_ACO2,k) = 1.8D-02  !HCOOH
     alpha(ind_ACO2,k)= 1./(exp(-7.9E3*CoRT+34.9*CoR)+1._dp)   !HCOOH
     alpha(ind_ACTA,k) = 6.7D-02  ! at 273 K, Jayne et al., 1991
     alpha(ind_CH3OH,k) = 5.6D-02  ! at 273 K, Jayne et al., 1991
     alpha(ind_C2H5OH,k) = 4.8D-02  ! at 273 K, Jayne et al., 1991
     alpha(ind_CO2,k) = 1./(1._dp +1./(1./(1./1.0d-2-1._dp)*zexp2)) !06.04.00
!     alpha(ind_HCl,k) = 7.2D-02
!     alpha(ind_HCl,k) = 0.1 ! default value
     alpha(ind_HCl,k) = 1./(exp(-3.072d3/t(k) + 1.283d1) + 1._dp) !T=290: 0.096, T=270: 0.190
!     alpha(ind_HOCl,k) = 7.2D-02
!     alpha(ind_HOCl,k) = 5.0D-01 see below
!     alpha(ind_ClNO3,k) = 1.0D-01 ! default value
!     alpha(ind_Cl2,k) = 5.5D-02
     alpha(ind_Cl2,k) = 1./(exp(-1.3d4*CoRT + 50.*CoR) + 1._dp)
!     alpha(ind_HBr,k) = 7.2D-02
!     alpha(ind_HBr,k) = 0.05
     alpha(ind_HBr,k) = 1./(exp(-3.94d3/t(k) + 1.664d1) + 1._dp) !T=290K: 0.017, T=270K: 0.130
     alpha(ind_HOBr,k) = 6.0D-01  ! #1077
     alpha(ind_HOCl,k) = alpha(ind_HOBr,k)
     alpha(ind_BrNO3,k) = 8.0D-01
!     alpha(ind_Br2,k) = 5.5D-02
     alpha(ind_Br2,k) = 1./(exp(-1.3D4*CoRT+50.*CoR)+1._dp)
!     alpha(ind_BrCl,k) = 5.5D-02
     alpha(ind_BrCl,k) = 0.33_dp !#840 alpha(ind_Cl2,k)
     alpha(ind_SO2,k) = 1.1D-01
!     alpha(ind_CH3SO3H,k) = 8.4D-02
     alpha(ind_CH3SO3H,k)= 1./(exp(-3.50D3*CoRT+16.7*CoR) + 1._dp) ! MSA, #955
     alpha(ind_DMS,k) = 1.0D-2 ! assumed
!     alpha(ind_DMSO,k) = 5.6D-02
     alpha(ind_DMSO,k) = 1./(exp(-5.12D3*CoRT+23.1*CoR) + 1._dp)  ! #955
     alpha(ind_DMSO2,k) = 1./(exp(-10.7D3*CoRT+43.0*CoR) + 1._dp) ! #955
     alpha(ind_CH3SO2H,k) = 2.0D-4 ! assumed #2123, MSIA
! no uptake of other DMS products like CH3SCH2OO, CH3S, CH3SO, CH3SO2, CH3SO3
     alpha(ind_INO3,k) = 1./(1._dp +1./(1./(1./1.0d-1 -1._dp)*zexp2))!06.04.00
!!     alpha(ind_HOI,k) = 7.2D-02
!!     alpha(ind_HOI,k) = 5.0D-01
     alpha(ind_HOI,k) = alpha(ind_HOBr,k)
!!     alpha(ind_HI,k) = 7.2D-02
     alpha(ind_HI,k) = 1./(exp(-4.13d3/t(k) + 1.715d1)+1._dp)
     alpha(ind_I2,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp) * zexp2))!06.04.00
     alpha(ind_IO,k) = 1./(1._dp +1./(1./(1./ 5.0d-1 -1._dp) * zexp2))!06.04.00
     alpha(ind_I2O2,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)* zexp2))!06.04.00
     alpha(ind_ICl,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp) * zexp2))!06.04.00
     alpha(ind_IBr,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp) * zexp2))!06.04.00
     alpha(ind_INO2,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)* zexp2))!06.04.00
!     alpha(ind_ClCHOk,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)*zexp2))!06.04.00
!     alpha(ind_BrCHOk,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)*zexp2))!06.04.00
!     alpha(ind_OIO,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp)*zexp2)) ! assumed, #980
     alpha(ind_OIO,k) = 1._dp
     alpha(ind_HIO3,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp)*zexp2))  ! assumed, #980
     alpha(ind_XOR,k) = 7.0d-2  ! same as bromoethanol, Jayne et al., 1991
!     alpha(ind_I2O,k) = ?
!     alpha(ind_I2O3,k) = ?

!     alpha(ind_Hg,k)    = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgO,k)   = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgCl,k)  = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgCl2,k) = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgBr,k)  = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgBr2,k) = 1.d-1       ! caution - wild guess, no information found!!

  enddo

! check that alpha <= 1 to avoid that the T-dependencies mess up the
!  numbers
  do k=2,nf
     do j=1,NSPEC
        alpha(j,k)=min(1.d0,alpha(j,k))
     end do
  end do

end subroutine st_coeff_t

!
!--------------------------------------------------------------------------------
!

subroutine st_coeff_a

! Description :
! -----------
  ! sticking coefficients, needed for calculation of k_mt

  ! Reference: for instance Davidovits et al, Chem. Rev. 2006, 106, 1323-1354

! Modifications :
! -------------
  ! jjb: factorise some repetitive calculations

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       lpJoyce14bc

  USE constants, ONLY : &
! Imported Parameters:
       cal => cal15, &   ! calorie      [J]
       R => gas_const    ! gas constant [J/(mol*K)]

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'aer_Parameters.h' !additional common blocks and other definitions

  real (kind=dp), parameter :: CoR = cal/R

! Local scalars:
  integer :: j, k
  real (kind=dp) :: CoRT, RT, tcorr, zexp2
  real (kind=dp), external :: a_n2o5

! Common blocks:
  common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
  real (kind=dp) :: alpha, vmean
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho

! == End of declarations =======================================================


! Initialisation (default value = 0.1)
  alpha(:,:) = 0.1_dp

  do k=2,nf

     ! Compute a few factors used several times
     tcorr = 1._dp / t(k) - 1._dp / 298.15_dp
     RT    = R * t(k)
     CoRT  = cal/RT
     zexp2 = exp(2000._dp * tcorr)

! for some species a temperature dependence is calculated/estimated
! using:
! alpha = 1./(1.+1./(1./(1./alpha(T0)-1.)*EXP((-DeltaH/RGAS)*TCORR)))

! standard value
     alpha(ind_H2SO4,k) = 0.65_dp
     !alpha(ind_CH4,k)   = 0.1_dp ! default value
     !alpha(ind_HNO4,k)  = 0.1_dp
! experimentally determined values (from MOCCA)
     alpha(ind_O3P,k) = 1.0e-6_dp  ! JPL #15
     alpha(ind_O1D,k) = 1.0e-6_dp  ! JPL #15
     alpha(ind_O3,k)  = 2.0e-3_dp
     alpha(ind_O2,k)  = 1./(1._dp +1./(1./(1./1.0d-2 - 1._dp)*zexp2)) ! 06.04.00
!     alpha(ind_OH,k) = 4.0D-03
     alpha(ind_OH,k) = 1.0e-2_dp
     alpha(ind_HO2,k) = 2.0e-1_dp
!     alpha(ind_H2O2,k) = 9.1D-02
     alpha(ind_H2O2,k) = 1./(exp(-26.d3/RT+107.8456/R) + 1._dp)
     alpha(ind_NO,k) = 5.0D-05
     alpha(ind_NO2,k) = 1.5D-03

     if (lpJoyce14bc) then
        alpha(ind_NO3,k) = 2.5D-03 ! jjb: lower limit Thomas et al. (1989), Mihelcic et al. (1993) (cited by Finlayson Pitt & Pitt, Chemistry of the Upper and Lower Atmosphere, 2000
     else
!        alpha(ind_NO3,k) = 2.5D-03
        alpha(ind_NO3,k) = 4.0D-02
     end if
     if (lpJoyce14bc) then
        alpha(ind_N2O5,k) = a_n2o5(k,1)
     else
!        alpha(ind_N2O5,k) = 2.7D-02
!        alpha(ind_N2O5,k) = 1.0D-01 ! default value
     end if

     alpha(ind_HONO,k) = 4.0D-02
!     alpha(ind_HNO3,k) = 8.6D-02
     alpha(ind_HNO3,k) = 5.0D-01
     alpha(ind_NH3,k) = 6.0D-02
     alpha(ind_MO2,k) = 1./(1._dp +1./(1./(1./1.0d-2-1._dp)*zexp2))   ! 06.04.00
!     alpha(ind_ROOH,k) = 5.5D-03
!     alpha(ind_ROOH,k) = 0.01
     alpha(ind_ROOH,k)= 1./(exp(-6.5D3*CoRT+32.5*CoR)+1._dp)
     alpha(ind_HCHO,k) = 4.0D-02
!     alpha(ind_ACO2,k) = 1.8D-02  !HCOOH
     alpha(ind_ACO2,k)= 1./(exp(-7.9E3*CoRT+34.9*CoR)+1._dp)   !HCOOH
     alpha(ind_ACTA,k) = 6.7D-02  ! at 273 K, Jayne et al., 1991
     alpha(ind_CH3OH,k) = 5.6D-02  ! at 273 K, Jayne et al., 1991
     alpha(ind_C2H5OH,k) = 4.8D-02  ! at 273 K, Jayne et al., 1991
     alpha(ind_CO2,k) = 1./(1._dp +1./(1./(1./1.0d-2-1._dp)*zexp2)) !06.04.00
!     alpha(ind_HCl,k) = 7.2D-02
!     alpha(ind_HCl,k) = 0.1 ! default value
     alpha(ind_HCl,k) = 1./(exp(-3.072d3/t(k) + 1.283d1) + 1._dp) !T=290: 0.096, T=270: 0.190
!     alpha(ind_HOCl,k) = 7.2D-02
!     alpha(ind_HOCl,k) = 5.0D-01 see below
!     alpha(ind_ClNO3,k) = 1.0D-01 ! default value
!     alpha(ind_Cl2,k) = 5.5D-02
     alpha(ind_Cl2,k) = 1./(exp(-1.3d4*CoRT + 50.*CoR) + 1._dp)
!     alpha(ind_HBr,k) = 7.2D-02
!     alpha(ind_HBr,k) = 0.05
     alpha(ind_HBr,k) = 1./(exp(-3.94d3/t(k) + 1.664d1) + 1._dp) !T=290K: 0.017, T=270K: 0.130
     alpha(ind_HOBr,k) = 6.0D-01  ! #1077
     alpha(ind_HOCl,k) = alpha(ind_HOBr,k)
     alpha(ind_BrNO3,k) = 8.0D-01
!     alpha(ind_Br2,k) = 5.5D-02
     alpha(ind_Br2,k) = 1./(exp(-1.3D4*CoRT+50.*CoR)+1._dp)
!     alpha(ind_BrCl,k) = 5.5D-02
     alpha(ind_BrCl,k) = 0.33_dp !#840 alpha(ind_Cl2,k)
     alpha(ind_SO2,k) = 1.1D-01
!     alpha(ind_CH3SO3H,k) = 8.4D-02
     alpha(ind_CH3SO3H,k)= 1./(exp(-3.50D3*CoRT+16.7*CoR) + 1._dp) ! MSA, #955
     alpha(ind_DMS,k) = 1.0D-2 ! assumed
!     alpha(ind_DMSO,k) = 5.6D-02
     alpha(ind_DMSO,k) = 1./(exp(-5.12D3*CoRT+23.1*CoR) + 1._dp)  ! #955
     alpha(ind_DMSO2,k) = 1./(exp(-10.7D3*CoRT+43.0*CoR) + 1._dp) ! #955
     alpha(ind_CH3SO2H,k) = 2.0D-4 ! assumed #2123, MSIA
! no uptake of other DMS products like CH3SCH2OO, CH3S, CH3SO, CH3SO2, CH3SO3
     alpha(ind_INO3,k) = 1./(1._dp +1./(1./(1./1.0d-1 -1._dp)*zexp2))!06.04.00
!!     alpha(ind_HOI,k) = 7.2D-02
!!     alpha(ind_HOI,k) = 5.0D-01
     alpha(ind_HOI,k) = alpha(ind_HOBr,k)
!!     alpha(ind_HI,k) = 7.2D-02
     alpha(ind_HI,k) = 1./(exp(-4.13d3/t(k) + 1.715d1)+1._dp)
     alpha(ind_I2,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp) * zexp2))!06.04.00
     alpha(ind_IO,k) = 1./(1._dp +1./(1./(1./ 5.0d-1 -1._dp) * zexp2))!06.04.00
     alpha(ind_I2O2,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)* zexp2))!06.04.00
     alpha(ind_ICl,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp) * zexp2))!06.04.00
     alpha(ind_IBr,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp) * zexp2))!06.04.00
     alpha(ind_INO2,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)* zexp2))!06.04.00
!     alpha(ind_ClCHOk,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)*zexp2))!06.04.00
!     alpha(ind_BrCHOk,k) = 1./(1._dp +1./(1./(1./ 1.0d-1 -1._dp)*zexp2))!06.04.00
!     alpha(ind_OIO,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp)*zexp2)) ! assumed, #980
     alpha(ind_OIO,k) = 1._dp
     alpha(ind_HIO3,k) = 1./(1._dp +1./(1./(1./ 1.0d-2 -1._dp)*zexp2))  ! assumed, #980
     alpha(ind_XOR,k) = 7.0d-2  ! same as bromoethanol, Jayne et al., 1991
!     alpha(ind_I2O,k) = ?
!     alpha(ind_I2O3,k) = ?

!     alpha(ind_Hg,k)    = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgO,k)   = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgCl,k)  = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgCl2,k) = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgBr,k)  = 1.d-1       ! caution - wild guess, no information found!!
!     alpha(ind_HgBr2,k) = 1.d-1       ! caution - wild guess, no information found!!

  enddo

! check that alpha <= 1 to avoid that the T-dependencies mess up the
!  numbers
  do k=2,nf
     do j=1,NSPEC
        alpha(j,k)=min(1.d0,alpha(j,k))
     end do
  end do

end subroutine st_coeff_a


!
!-----------------------------------------------------
!

      subroutine v_mean_init

! Description :
! -----------
!     Compute the mean molecular speed from Maxwell-Boltzmann distribution:
!     v_mean=sqrt(8*R_gas*T/(M*pi))      (M in kg/mol)

! Interface :
! ---------
!    SR v_mean_init is called during initialisation:
!      - by SR initc (no restart case)
!      - by SR str=main (restart case)

! Input :
! -----
!    - chemical species molar masses have been imported from the user defined files ('gas_species.csv')

! Output :
! ------
!    - v_mean_init computes the constant factor in vmean: vmean_init=sqrt(8*R_gas  /(M*pi))

! Externals :
! ---------
!    none

! Method :
! ------
!    Improve computing efficiency: split vmean calculation into a constant part, calculated here,
!    and a variable term (sqrt(T)) computed only once in SR v_mean

! Author :
! ------
!    Josue Bock


! Modifications :
! -------------
!
! 04-Jan-2017   Josue Bock   First version
!

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

      USE constants, ONLY : &
! Imported Parameters:
           gas_const, &
           pi

      USE gas_common, ONLY : &
! Imported Parameters:
           j1, &
           j5, &
           j4, &
! Imported Array Variables with intent (in):
           gas_mass, &
           rad_mass, &
           fix_mass, &
! Imported Array Variables with intent (out):
           vmean_init, &
           vmean

      USE global_params, ONLY : &
! Imported Parameters:
           nmax_chem_aer

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

! Local scalars:
      integer :: jtot
      integer :: jspec
      real (kind=dp) :: const_fact

! Local arrays:
      real (kind=dp) :: sqrt_mass (j1 + j5 + j4)

! == End of declarations =======================================================

      jtot = j1 + j5 + j4

      allocate ( vmean_init(jtot) )
      allocate ( vmean(jtot,nmax_chem_aer) )

      const_fact = sqrt(8.d0*gas_const/pi)

      do jspec = 1,j1
         sqrt_mass(jspec) = sqrt(gas_mass(jspec))
      end do
      do jspec = 1,j5
         sqrt_mass(j1+jspec) = sqrt(rad_mass(jspec))
      end do
      do jspec = 1,j4
         sqrt_mass(j1+j5+jspec) = sqrt(fix_mass(jspec))
      end do

      do jspec = 1,jtot
         vmean_init(jspec) = const_fact / sqrt_mass(jspec)
      end do

      end subroutine v_mean_init

!
!-----------------------------------------------------
!

     subroutine v_mean (temperature)

! Description :
! -----------
!     Compute the mean molecular speed from Maxwell-Boltzmann distribution:
!     v_mean=sqrt(8*R_gas*T/(M*pi))      (M in kg/mol)

! Interface :
! ---------
!    SR v_mean is called:
!      - during initialisation:
!        - by SR initc (no restart case)
!        - by SR str=main (restart case)
!      - during the run
!        - by SR liq_parm
!        - by SR box_update (which is call during time integration, but actually calls v_mean only during initialisation)

! Input :
! -----
!    - vmean_init has been computed by SR v_mean_init during initialisation

! Output :
! ------
!    - v_mean computes the mean molecular speed, using constant factor in vmean: vmean_init=sqrt(8*R_gas  /(M*pi))

! Externals :
! ---------
!    none

! Method :
! ------
!    Improve computing efficiency: split vmean calculation into a constant part, calculated here,
!    and a variable term (sqrt(T)) computed only once in SR v_mean

! Author :
! ------
!    Roland von Glasow


! Modifications :
! -------------
! 17-Jul-2015   Josue Bock   Commented unused /cb40/ (found in the code that vmean had been computed only between
!                               lcl and lct, but this had been changed, now from 1 to nmaxf
!
! 05-Mar-2016   Josue Bock   Forcheck errors 307E, 312E: commented lines related to CHBr2I, I2O, I2O3, I2O4, I2O5 and INO
!                               (undefined indexes)
!
! 17-Mar-2016   Josue Bock   Reindexed vmean array for computing efficiency: innermost is leftmost
!                               vmean(nf,NSPEC) -> vmean(NSPEC,nf)
!
! 22-Mar-2016   Josue Bock   Missing species added: SO3, H2SO2 and corrected values for HNO4 (63->79 g/mol), XOR (109->125 g/mol)
!
! 04-Jan-2017   Josue Bock   Major change to the code structure: split v_mean_init and v_mean (no longer v_mean_a and v_mean_t)
!                               This is ready in my sub-version jjb6v1.0*
!
! 11-Feb-2017   Josue Bock   Changed double-specific math function (dsqrt) into generic one (sqrt) in the whole file
!                               Also changed the exponent of hard-written mass values (e -> d) for a few species, for consistency

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

      USE gas_common, ONLY : &
! Imported Parameters:
           j1, &
           j5, &
           j4, &
! Imported Array Variables with intent (in):
           vmean_init, &
! Imported Array Variables with intent (out):
           vmean

      USE global_params, ONLY : &
! Imported Parameters:
           nmax_chem_aer

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      real (kind=dp), intent(in) :: temperature (nmax_chem_aer)
! Local scalars:
      integer :: jtot
      integer :: j,k
! Local arrays:
      real (kind=dp) :: sqrtt(nmax_chem_aer)

! == End of declarations =======================================================

      jtot = j1 + j5 + j4

      sqrtt = sqrt(temperature)

      do k=1,nmax_chem_aer
         do j=1,jtot
            vmean(j,k) = vmean_init(j) * sqrtt(k)
         end do
      end do

      end subroutine v_mean
!
!-----------------------------------------------------
!

subroutine v_mean_t (tt,nmaxf)

! Description :
! -----------
  ! mean molecular speed from Maxwell-Boltzmann distribution:
  ! v_mean=sqrt(8*R_gas*T/(M*pi))      (M in kg/mol)
  ! v_mean in m/s

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb vmean SO3 and HOSO2 were missing, HNO4 was wrong (HNO3 mass used instead of HNO4)

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'tot_Parameters.h' !additional common blocks and other definitions

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: tt(n)
  integer, intent(in) :: nmaxf

! Local scalars:
  integer :: k

! Common blocks:
  common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
  real (kind=dp) :: alpha, vmean

! Statement function
  real (kind=dp) :: func, a
  ! sqrt(8*R_gas/pi)=4.60138
  func(a,k)=sqrt(tt(k)/a)*4.60138

! == End of declarations =======================================================

! Initialisation of vmean
  vmean(:,:) = 0._dp

  do k=1,nmaxf
     vmean(ind_NO,k) = func(3.d-2,k)
     vmean(ind_NO2,k) = func(4.6d-2,k)
     vmean(ind_HNO3,k) = func(6.3d-2,k)
     vmean(ind_NH3,k) = func(1.7d-2,k)
     vmean(ind_SO2,k) = func(6.4d-2,k)
     vmean(ind_SO3,k) = func(8.0d-2,k)
     vmean(ind_HOSO2,k) = func(8.1d-2,k)
     vmean(ind_H2SO4,k) = func(9.8d-2,k)
     vmean(ind_O3,k) = func(4.8d-2,k)
     vmean(ind_CH4,k) = func(1.6d-2,k)
     vmean(ind_C2H6,k) = func(3.d-2,k)
!     vmean(ind_C3H8,k) = func(4.4d-2,k)
!     vmean(ind_ALKA,k) = func(,k)
     vmean(ind_ETHE,k) = func(2.8d-2,k)
!     vmean(ind_ALKE,k) = func(,k)
!     vmean(ind_AROM,k) = func(,k)
     vmean(ind_ACO2,k) = func(4.6d-2,k)
     vmean(ind_ACTA,k) = func(6.d-2,k)
     vmean(ind_HCHO,k) = func(3.d-2,k)
     vmean(ind_ALD2,k) = func(4.4d-2,k)  ! value for CH3CHO
     vmean(ind_H2O2,k) = func(3.4d-2,k)
     vmean(ind_ROOH,k) = func(4.8d-2,k) ! value for CH3OOH
     vmean(ind_HONO,k) = func(4.7d-2,k)
     vmean(ind_PAN,k) = func(1.21d-1,k)
!     vmean(ind_TPAN,k) = func(1.59d-1,k)
!     vmean(ind_KET,k) = func(,k)
!     vmean(ind_CRES,k) = func(1.08d-1,k)
!     vmean(ind_DIAL,k) = func(8.4d-2,k)
!     vmean(ind_GLYX,k) = func(5.8d-2,k)
!     vmean(ind_MGLY,k) = func(7.2d-2,k)
!     vmean(ind_NH4NO3,k) = func(8.1d-2,k)
     vmean(ind_HCl,k) = func(3.6d-2,k)
!     vmean(ind_R3N2,k) = func(,k)
!     vmean(ind_RAN1,k) = func(,k)
!     vmean(ind_RAN2,k) = func(,k)
     vmean(ind_N2O5,k) = func(1.08d-1,k)
     vmean(ind_HNO4,k) = func(7.9d-2,k)
     vmean(ind_NO3,k) = func(6.2d-2,k)
     vmean(ind_DMS,k) = func(6.2d-2,k)
     vmean(ind_HOCl,k) = func(5.2d-2,k)
     vmean(ind_ClNO2,k) = func(8.1d-2,k)
     vmean(ind_ClNO3,k) = func(9.7d-2,k)
     vmean(ind_Cl2,k) = func(7.1d-2,k)
     vmean(ind_Cl2O2,k) = func(1.029d-1,k)
     vmean(ind_HBr,k) = func(8.1d-2,k)
     vmean(ind_HOBr,k) = func(9.7d-2,k)
     vmean(ind_BrNO2,k) = func(1.26d-1,k)
     vmean(ind_BrNO3,k) = func(1.42d-1,k)
     vmean(ind_Br2,k) = func(1.6d-1,k)
     vmean(ind_BrCl,k) = func(1.15d-1,k)
     vmean(ind_HI,k) = func(1.28d-1,k)
     vmean(ind_HOI,k) = func(1.44d-1,k)
     vmean(ind_I2O2,k) = func(2.86d-1,k)
     vmean(ind_INO2,k) = func(1.73d-1,k)
     vmean(ind_INO3,k) = func(1.89d-1,k)
     vmean(ind_I2,k) = func(2.54d-1,k)
     vmean(ind_ICl,k) = func(1.62d-1,k)
     vmean(ind_IBr,k) = func(2.07d-1,k)
     vmean(ind_HIO3,k) = func(1.76d-1,k)
     vmean(ind_CH3I,k) = func(1.42d-1,k)
     vmean(ind_CH2I2,k) = func(2.68d-1,k)
     vmean(ind_CH2ClI,k) = func(1.76d-1,k)
     vmean(ind_C3H7I,k) = func(1.7d-1,k)
     vmean(ind_CH2BrI,k) = func(2.21d-1,k)
!     vmean(ind_CHBr2I,k) = func(3.d-1,k)
     vmean(ind_C2H5I,k) = func(1.56d-1,k)
     vmean(ind_DMS,k) = func(6.2d-2,k)
     vmean(ind_DMSO,k) = func(7.8d-2,k)
     vmean(ind_DMSO2,k) = func(9.4d-2,k)
     vmean(ind_DMOO,k) = func(9.3d-2,k) ! CH3SCH2OO
     vmean(ind_CH3S,k) = func(4.7d-2,k)
     vmean(ind_CH3SO,k) = func(6.3d-2,k)
     vmean(ind_CH3SO2,k) = func(7.9d-2,k)
     vmean(ind_CH3SO3,k) = func(9.5d-2,k)
     vmean(ind_CH3SO2H,k) = func(8.0d-2,k)   ! CH3S(O)OH, MSIA
     vmean(ind_CH3SO3H,k) = func(9.6d-2,k) ! CH3S(OO)OH, MSA
     vmean(ind_CO,k) = func(2.8d-2,k)
     vmean(ind_CO2,k) = func(4.4d-2,k)
!     vmean(ind_I2O,k) = func(2.70d-1,k)
!     vmean(ind_I2O3,k) = func(3.02d-1,k)
!     vmean(ind_I2O4,k) = func(3.18d-1,k)
!     vmean(ind_I2O5,k) = func(3.34d-1,k)
!     vmean(ind_INO,k) = func(1.57d-1,k)
     vmean(ind_Br2O,k) = func(1.76d-1,k)
     vmean(ind_ClONO,k) = func(8.15d-2,k)
     vmean(ind_ClO3,k) = func(8.35d-2,k)
     vmean(ind_Cl2O3,k) = func(1.19d-1,k)
     vmean(ind_CH3OH,k) = func(3.2d-2,k)
     vmean(ind_C2H5OH,k) = func(4.6d-2,k)
     vmean(ind_H2,k) = func(2.0d-3,k)
     vmean(ind_NHS,k) = func(5.8d-2,k)  ! C+N+S
     vmean(ind_RCl,k) = func(6.45d-2,k)  ! calculated using C2H5Cl
     vmean(ind_RBr,k) = func(1.27d-1,k)  ! calculated using CH3SBr
     vmean(ind_XOR,k) = func(1.09d-1,k)  ! calculated using bromoethanol
     vmean(ind_SOR,k) = func(9.4d-2,k)  ! calculated using CH3SCH2OOH
     vmean(ind_SPAN,k) = func(1.39d-1,k)  ! calculated using CH3SCH2OONO2
!     vmean(ind_Hg,k)    = func(2.00d-1,k)
!     vmean(ind_HgO,k)   = func(2.16d-1,k)
!     vmean(ind_HgCl,k)  = func(2.36d-1,k)
!     vmean(ind_HgCl2,k) = func(2.72d-1,k)
!     vmean(ind_HgBr,k)  = func(2.81d-1,k)
!     vmean(ind_HgBr2,k) = func(3.61d-1,k)

!#DEFRAD
     vmean(ind_OH,k) = func(1.7d-2,k)
     vmean(ind_HO2,k) = func(3.3d-2,k)
!     vmean(ind_AHO2,k) = func(6.3d-2,k)
     vmean(ind_MCO3,k) = func(7.5d-2,k)
     vmean(ind_MO2,k) = func(4.7d-2,k)
     vmean(ind_ETO2,k) = func(6.1d-2,k)
!     vmean(ind_KO2,k) = func(,k)
!     vmean(ind_R3O2,k) = func(,k)
!     vmean(ind_RAO2,k) = func(,k)
!     vmean(ind_TO2,k) = func(,k)
!     vmean(ind_TCO3,k) = func(1.15d-1,k)
!     vmean(ind_ZO2,k) = func(,k)
     vmean(ind_EO2,k) = func(7.7d-2,k)
!     vmean(ind_PO2,k) = func(,k)
     vmean(ind_CHO2,k) = func(4.6d-2,k)
!     vmean(ind_CRO2,k) = func(6.d-2,k)
!     vmean(ind_PRN1,k) = func(,k)
     vmean(ind_O1D,k) = func(1.6d-2,k)
     vmean(ind_Cl,k) = func(3.5d-2,k)
     vmean(ind_ClO,k) = func(5.1d-2,k)
     vmean(ind_OClO,k) = func(6.7d-2,k)
     vmean(ind_Br,k) = func(8.d-2,k)
     vmean(ind_BrO,k) = func(9.6d-2,k)
     vmean(ind_I,k) = func(1.27d-1,k)
     vmean(ind_IO,k) = func(1.43d-1,k)
     vmean(ind_OIO,k) = func(1.59d-1,k)
     vmean(ind_O3P,k) = func(1.6d-2,k)
     vmean(ind_ClRO2,k) = func(9.64d-2,k)  ! calculated using C2H5ClOO
     vmean(ind_BrRO2,k) = func(1.58d-1,k)  ! calculated using CH2SBrOO
     vmean(ind_IRO2,k) = func(1.73d-1,k)  ! calculated using CH2IOO

!#DEFFIX
     vmean(ind_O2,k) = func(3.2d-2,k)

  enddo

end subroutine v_mean_t


!
!-----------------------------------------------------
!

subroutine v_mean_a (tt,nmaxf)

! Description :
! -----------
  ! mean molecular speed from Maxwell-Boltzmann distribution:
  ! v_mean=sqrt(8*R_gas*T/(M*pi))      (M in kg/mol)
  ! v_mean in m/s

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb vmean SO3 and HOSO2 were missing, HNO4 was wrong (HNO3 mass used instead of HNO4)

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'aer_Parameters.h' !additional common blocks and other definitions

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: tt(n)
  integer, intent(in) :: nmaxf

! Local scalars:
  integer :: k

! Common blocks:
  common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
  real (kind=dp) :: alpha, vmean

! Statement function
  real (kind=dp) :: func, a
  ! sqrt(8*R_gas/pi)=4.60138
  func(a,k)=sqrt(tt(k)/a)*4.60138

! == End of declarations =======================================================

! Initialisation of vmean
  vmean(:,:) = 0._dp

  do k=1,nmaxf
     vmean(ind_NO,k) = func(3.d-2,k)
     vmean(ind_NO2,k) = func(4.6d-2,k)
     vmean(ind_HNO3,k) = func(6.3d-2,k)
     vmean(ind_NH3,k) = func(1.7d-2,k)
     vmean(ind_SO2,k) = func(6.4d-2,k)
     vmean(ind_SO3,k) = func(8.0d-2,k)
     vmean(ind_HOSO2,k) = func(8.1d-2,k)
     vmean(ind_H2SO4,k) = func(9.8d-2,k)
     vmean(ind_O3,k) = func(4.8d-2,k)
     vmean(ind_CH4,k) = func(1.6d-2,k)
     vmean(ind_C2H6,k) = func(3.d-2,k)
!     vmean(ind_C3H8,k) = func(4.4d-2,k)
!     vmean(ind_ALKA,k) = func(,k)
     vmean(ind_ETHE,k) = func(2.8d-2,k)
!     vmean(ind_ALKE,k) = func(,k)
!     vmean(ind_AROM,k) = func(,k)
     vmean(ind_ACO2,k) = func(4.6d-2,k)
     vmean(ind_ACTA,k) = func(6.d-2,k)
     vmean(ind_HCHO,k) = func(3.d-2,k)
     vmean(ind_ALD2,k) = func(4.4d-2,k)  ! value for CH3CHO
     vmean(ind_H2O2,k) = func(3.4d-2,k)
     vmean(ind_ROOH,k) = func(4.8d-2,k) ! value for CH3OOH
     vmean(ind_HONO,k) = func(4.7d-2,k)
     vmean(ind_PAN,k) = func(1.21d-1,k)
!     vmean(ind_TPAN,k) = func(1.59d-1,k)
!     vmean(ind_KET,k) = func(,k)
!     vmean(ind_CRES,k) = func(1.08d-1,k)
!     vmean(ind_DIAL,k) = func(8.4d-2,k)
!     vmean(ind_GLYX,k) = func(5.8d-2,k)
!     vmean(ind_MGLY,k) = func(7.2d-2,k)
!     vmean(ind_NH4NO3,k) = func(8.1d-2,k)
     vmean(ind_HCl,k) = func(3.6d-2,k)
!     vmean(ind_R3N2,k) = func(,k)
!     vmean(ind_RAN1,k) = func(,k)
!     vmean(ind_RAN2,k) = func(,k)
     vmean(ind_N2O5,k) = func(1.08d-1,k)
     vmean(ind_HNO4,k) = func(7.9d-2,k)
     vmean(ind_NO3,k) = func(6.2d-2,k)
     vmean(ind_DMS,k) = func(6.2d-2,k)
     vmean(ind_HOCl,k) = func(5.2d-2,k)
     vmean(ind_ClNO2,k) = func(8.1d-2,k)
     vmean(ind_ClNO3,k) = func(9.7d-2,k)
     vmean(ind_Cl2,k) = func(7.1d-2,k)
     vmean(ind_Cl2O2,k) = func(1.029d-1,k)
     vmean(ind_HBr,k) = func(8.1d-2,k)
     vmean(ind_HOBr,k) = func(9.7d-2,k)
     vmean(ind_BrNO2,k) = func(1.26d-1,k)
     vmean(ind_BrNO3,k) = func(1.42d-1,k)
     vmean(ind_Br2,k) = func(1.6d-1,k)
     vmean(ind_BrCl,k) = func(1.15d-1,k)
     vmean(ind_HI,k) = func(1.28d-1,k)
     vmean(ind_HOI,k) = func(1.44d-1,k)
     vmean(ind_I2O2,k) = func(2.86d-1,k)
     vmean(ind_INO2,k) = func(1.73d-1,k)
     vmean(ind_INO3,k) = func(1.89d-1,k)
     vmean(ind_I2,k) = func(2.54d-1,k)
     vmean(ind_ICl,k) = func(1.62d-1,k)
     vmean(ind_IBr,k) = func(2.07d-1,k)
     vmean(ind_HIO3,k) = func(1.76d-1,k)
     vmean(ind_CH3I,k) = func(1.42d-1,k)
     vmean(ind_CH2I2,k) = func(2.68d-1,k)
     vmean(ind_CH2ClI,k) = func(1.76d-1,k)
     vmean(ind_C3H7I,k) = func(1.7d-1,k)
     vmean(ind_CH2BrI,k) = func(2.21d-1,k)
!     vmean(ind_CHBr2I,k) = func(3.d-1,k)
     vmean(ind_C2H5I,k) = func(1.56d-1,k)
     vmean(ind_DMS,k) = func(6.2d-2,k)
     vmean(ind_DMSO,k) = func(7.8d-2,k)
     vmean(ind_DMSO2,k) = func(9.4d-2,k)
     vmean(ind_DMOO,k) = func(9.3d-2,k) ! CH3SCH2OO
     vmean(ind_CH3S,k) = func(4.7d-2,k)
     vmean(ind_CH3SO,k) = func(6.3d-2,k)
     vmean(ind_CH3SO2,k) = func(7.9d-2,k)
     vmean(ind_CH3SO3,k) = func(9.5d-2,k)
     vmean(ind_CH3SO2H,k) = func(8.0d-2,k)   ! CH3S(O)OH, MSIA
     vmean(ind_CH3SO3H,k) = func(9.6d-2,k) ! CH3S(OO)OH, MSA
     vmean(ind_CO,k) = func(2.8d-2,k)
     vmean(ind_CO2,k) = func(4.4d-2,k)
!     vmean(ind_I2O,k) = func(2.70d-1,k)
!     vmean(ind_I2O3,k) = func(3.02d-1,k)
!     vmean(ind_I2O4,k) = func(3.18d-1,k)
!     vmean(ind_I2O5,k) = func(3.34d-1,k)
!     vmean(ind_INO,k) = func(1.57d-1,k)
     vmean(ind_Br2O,k) = func(1.76d-1,k)
     vmean(ind_ClONO,k) = func(8.15d-2,k)
     vmean(ind_ClO3,k) = func(8.35d-2,k)
     vmean(ind_Cl2O3,k) = func(1.19d-1,k)
     vmean(ind_CH3OH,k) = func(3.2d-2,k)
     vmean(ind_C2H5OH,k) = func(4.6d-2,k)
     vmean(ind_H2,k) = func(2.0d-3,k)
     vmean(ind_NHS,k) = func(5.8d-2,k)  ! C+N+S
     vmean(ind_RCl,k) = func(6.45d-2,k)  ! calculated using C2H5Cl
     vmean(ind_RBr,k) = func(1.27d-1,k)  ! calculated using CH3SBr
     vmean(ind_XOR,k) = func(1.09d-1,k)  ! calculated using bromoethanol
     vmean(ind_SOR,k) = func(9.4d-2,k)  ! calculated using CH3SCH2OOH
     vmean(ind_SPAN,k) = func(1.39d-1,k)  ! calculated using CH3SCH2OONO2
!     vmean(ind_Hg,k)    = func(2.00d-1,k)
!     vmean(ind_HgO,k)   = func(2.16d-1,k)
!     vmean(ind_HgCl,k)  = func(2.36d-1,k)
!     vmean(ind_HgCl2,k) = func(2.72d-1,k)
!     vmean(ind_HgBr,k)  = func(2.81d-1,k)
!     vmean(ind_HgBr2,k) = func(3.61d-1,k)

!#DEFRAD
     vmean(ind_OH,k) = func(1.7d-2,k)
     vmean(ind_HO2,k) = func(3.3d-2,k)
!     vmean(ind_AHO2,k) = func(6.3d-2,k)
     vmean(ind_MCO3,k) = func(7.5d-2,k)
     vmean(ind_MO2,k) = func(4.7d-2,k)
     vmean(ind_ETO2,k) = func(6.1d-2,k)
!     vmean(ind_KO2,k) = func(,k)
!     vmean(ind_R3O2,k) = func(,k)
!     vmean(ind_RAO2,k) = func(,k)
!     vmean(ind_TO2,k) = func(,k)
!     vmean(ind_TCO3,k) = func(1.15d-1,k)
!     vmean(ind_ZO2,k) = func(,k)
     vmean(ind_EO2,k) = func(7.7d-2,k)
!     vmean(ind_PO2,k) = func(,k)
     vmean(ind_CHO2,k) = func(4.6d-2,k)
!     vmean(ind_CRO2,k) = func(6.d-2,k)
!     vmean(ind_PRN1,k) = func(,k)
     vmean(ind_O1D,k) = func(1.6d-2,k)
     vmean(ind_Cl,k) = func(3.5d-2,k)
     vmean(ind_ClO,k) = func(5.1d-2,k)
     vmean(ind_OClO,k) = func(6.7d-2,k)
     vmean(ind_Br,k) = func(8.d-2,k)
     vmean(ind_BrO,k) = func(9.6d-2,k)
     vmean(ind_I,k) = func(1.27d-1,k)
     vmean(ind_IO,k) = func(1.43d-1,k)
     vmean(ind_OIO,k) = func(1.59d-1,k)
     vmean(ind_O3P,k) = func(1.6d-2,k)
     vmean(ind_ClRO2,k) = func(9.64d-2,k)  ! calculated using C2H5ClOO
     vmean(ind_BrRO2,k) = func(1.58d-1,k)  ! calculated using CH2SBrOO
     vmean(ind_IRO2,k) = func(1.73d-1,k)  ! calculated using CH2IOO

!#DEFFIX
     vmean(ind_O2,k) = func(3.2d-2,k)

  enddo

end subroutine v_mean_a

!
!------------------------------------------------------
!

subroutine henry_t (tt,nmaxf)

! Description :
! -----------
  ! Temp dependent Henry constants
  ! inverse dimensionless Henry constant: k_(H,inv)^cc:=1/(k_H^cp*RT)
  ! in equilibrium: XXXaq = k_H^cp * LWC * XXXg

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
! jjb work done:
!     - removed pp (pressure) from the argument list: unused
!     - removed hard coded parameters, use modules instead
!     - improved computation efficiency: Tfact calculated only once per layer
!     - missing declarations and implicit none
!     - final conversion (inverse if henry /= 0.) optimised
!     - reindexed henry(NSPEC,nf) instead of (nf,NSPEC) for computing efficiency
!     - added initialisation of henry(:,:)

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'tot_Parameters.h'    !additional common blocks and other definitions

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: tt(n) ! temperature array
  integer       , intent(in) :: nmaxf ! max layer index where henry has to be computed

! Local scalars:
  integer :: j,k                  ! loop indexes
  real (kind=dp) :: func3, a0, b0 ! temperature dependency function, and its arguments
  real (kind=dp) :: FCT           ! conversion factor, see below
  real (kind=dp) :: Tfact         ! Tfact, local variable [K**-1]
!  real (kind=dp) :: xCO2, x3CO2, xNH3, xhp

! Common blocks:
  common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
       xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
  real (kind=dp) henry, xkmt, xkef, xkeb

! 0.082=8.3145*10^3/101325=R/p_0*10^3 : R includes conversion from M/atm --> mol/(m^3*Pa)
! so k_H^cp is taken as M/atm (the most common literature unit): SEE END OF SR

! func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))*0.082*tt(k0)
! func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))

! jjb 11/02/2017 slight improvement in calculation efficiency: compute Tfact only once per layer
!  Tfact = 1/T - 1/Tref, see below
  func3(a0,b0)=a0*exp(b0*Tfact)

! == End of declarations =======================================================

  henry(:,:) = 0._dp

  do k=1,nmaxf

!     henry(ind_H2SO4,k)=4.1d-12   !???
     henry(ind_H2SO4,k)=1.d+16  !NIST --> RS_Gmitro_Vermeulen
     henry(ind_CH4,k)=1.3d-3    !RS_Mackay
     henry(ind_C2H6,k)=2.0d-3   !RS_Mackay
!     henry(ind_C3H8,k)=1.4d-3   !RS_Mackay
!     henry(ind_ALKA,k)=9.5d-4   !ok
     henry(ind_ETHE,k)=4.9d-3   !ok
!     henry(ind_ALKE,k)=4.9d-3   !ok
!     henry(ind_AROM,k)=1.5d-1   !RS_Mackay: methylbenzene
!     henry(ind_KET,k)=10.d0     !RS_higher ketone
!     henry(ind_CRES,k)=8.2d2    !RS_Hine
!     henry(ind_DIAL,k)=10.0d0   !RS_Snider:propenal #?
!     henry(ind_GLYX,k)=3.6d5    !RS_Zhou
!     henry(ind_NH4NO3,k)=4.1d-12  !# was soll dieses Salz in Gasphase??
!     henry(ind_RAN2,k)=1.2d0    !RS_Kames: pentyl-nitrate
!     henry(ind_RAN1,k)=1.d0     !RS_hine
!     henry(ind_N2O5,k)=0.0      !RS_Sander  infinity!
!     henry(ind_ClNO2,k)=0.0     !  10^-2
!     henry(ind_ClNO3,k)=0.0     !  infinity
!     henry(ind_BrNO2,k)=0.0     !  10^-1
!     henry(ind_BrNO3,k)=0.0     !  infinity
     henry(ind_HI,k)=0.0d0        !   (dissociation)
     henry(ind_I2O2,k)=0.0d0      !
     henry(ind_INO2,k)=0.0d0      !
     henry(ind_INO3,k)=0.0d0      !
!     henry(ind_I2O,k)=0.0d0      !
!     henry(ind_I2O3,k)=0.0d0      !
!     henry(ind_OIO,k)  =  unknown but not needed as only surface reaction and no reversible uptake
!     henry(ind_HIO3,k) =                     -"-
     henry(ind_C3H7I,k)=1.1d-1  !RS_Hine
!     henry(ind_CH3SO2,k)= !RS_MOCCA
!     henry(ind_CH3SO3,k)= !RS_MOCCA
!     henry(ind_Hg,k)    = 1.3d-1                   ! #233 in #3127
!     henry(ind_HgO,k)   = 2.69d12                  ! #233 in #3127
!     henry(ind_HgCl,k)  = 2.75d6                   ! assumed HgCl2 (poss. higher??)
!     henry(ind_HgCl2,k) = 2.75d6                   ! #233 in #3127
!     henry(ind_HgBr,k)  = 2.75d6                   ! assumed HgBr2 (poss. higher??)
!     henry(ind_HgBr2,k) = 2.75d6                   ! #3127

! explicitly Temp dependent
     ! Tfact = 1/T - 1/Tref  with Tref = 298.15 K
     ! 1/298.15 = 3.3540d-3
     Tfact = 1.d0/tt(k) - 3.3540d-3

     henry(ind_NO,k)=func3(1.9d-03,1480.d0)   !RS_Lide
!     henry(ind_NO2,k)=func3(1.d-02,2500.d0)  !RS_Chameides
     henry(ind_NO2,k)=func3(6.4d-03,2500.d0)  !RS_MOCCA
!     henry(ind_HNO3,k)=func3(1.66d5,8694.d0) !RS_MOCCA
     henry(ind_HNO3,k)=func3(2.5d6/1.5d1,8694.d0) !RS_MOCCA_exakt
!     henry(ind_HNO4,k)=1.4d4     !Goetz, 1996, cited in Warneck, 1999, #695
     henry(ind_HNO4,k)=func3(1.2d4,6900.d0)   !06.04.00
     henry(ind_NH3,k)=func3(58.d0,4085.d0)    !RS_MOCCA
     henry(ind_SO2,k)=func3(1.2d0,3120.d0)    !RS_MOCCA
     henry(ind_O3,k)=func3(1.2d-02,2560.d0)   !RS_MOCCA
     henry(ind_ACO2,k)=func3(3.7d+03,5700.d0) !RS_MOCCA
     henry(ind_ACTA,k)=func3(4.1d+03,6300.d0) !RS_Johnson
     henry(ind_HCHO,k)=func3(7.0d+03,6425.d0) !RS_MOCCA
     henry(ind_ALD2,k)=func3(1.3d+01,5700.d0)  !RS_Benkelberg: acetaldehyde
     henry(ind_H2O2,k)=func3(1.d+05,6338.d0)   !RS_MOCCA
!     henry(ind_ROOH,k)=func3(7.45d+04,6620.d0) ! # aktualisieren
     henry(ind_ROOH,k)=func3(3.0d+02,5322.d0)  !RS_MOCCA
!     henry(ind_HONO,k)=func3(5.0d+01,4900.d0) !RS_Becker
     henry(ind_HONO,k)=func3(4.9d+01,4780.d0)  !RS_MOCCA
     henry(ind_PAN,k)=func3(2.8d0,6500.d0)     !RS_Kames
!     henry(ind_TPAN,k)=henry(22,k)
!     henry(ind_MGLY,k)=func3(3.7d+03,7553.d0) !RS_Betterton
!     henry(ind_HCl,k)=func3(1.17d0,9001.d0)   !RS_MOCCA
     henry(ind_HCl,k)=func3(2.d0/1.7d0,9001.d0) !RS_MOCCA_exakt
!     henry(ind_R3N2,k)=func3(1.d0,5450.d0)    !RS_kames: propyl nitrate
     henry(ind_NO3,k)=func3(2.d0,2000.d0)      !RS_MOCCA
     henry(ind_DMS,k)=func3(4.8d-1,3100.d0)    !RS_deBruyn
!     henry(ind_DMSO,k)=5.d4                   !RS_MOCCA
     henry(ind_DMSO,k)=func3(5.d4,6425.d0)     !RS_MOCCA 06.04.00
     henry(ind_DMSO2,k)= 1.d+16                !DMSO2=H2SO4, assumed
     henry(ind_CH3SO2H,k)=1.d+16               !MSIA=H2SO4, assumed
     henry(ind_CH3SO3H,k)=1.d+16               !MSA=H2SO4, assumed
     henry(ind_HOCl,k)=func3(6.7d2,5862.d0)    !RS_MOCCA
!     henry(ind_Cl2,k)=9.2d-2    !RS_MOCCA
     henry(ind_Cl2,k)=func3(9.1d-2,2500.d0)    !RS_MOCCA 06.04.00
     henry(ind_HBr,k)=func3(1.3d0,10239.d0)    !RS_MOCCA
     henry(ind_Br2,k)=func3(7.6d-1,4094.d0)    !RS_MOCCA
     henry(ind_BrCl,k)=func3(9.4d-1,5600.d0)   !RS_MOCCA
!     henry(ind_HOBr,k)=9.3d1    !RS_MOCCA
     henry(ind_HOBr,k)=func3(9.3d1,5862.d0)    !RS_MOCCA 06.04.00
!!     henry(ind_HI,k)=func3(2.5d9/K_a,9800.d0) !RS_Brimblecombe #K_a
     henry(ind_I2,k)=func3(3.d0,4431.d0) !RS_MOCCA
!!     henry(ind_HOI,k)=4.5d2     !RS_MOCCA
     henry(ind_HOI,k)=func3(4.5d2,5862.d0)     !RS_MOCCA 06.04.00
!!     henry(ind_ICl,k)=1.1d2     !RS_MOCCA
     henry(ind_ICl,k)=func3(1.1d2,5600.d0)     !RS_MOCCA 06.04.00
!!     henry(ind_IBr,k)=2.4d1     !RS_MOCCA
     henry(ind_IBr,k)=func3(2.4d1,5600.d0)     !RS_MOCCA 06.04.00
     henry(ind_CH3I,k)=func3(1.4d-1,4300.d0)   !RS_Moore
     henry(ind_CH2I2,k)=func3(2.3d0,5000.d0)   !RS_Moore
     henry(ind_CH2ClI,k)=func3(8.9d-1,4300.d0) !RS_Moore

! for radicals only OH, HO2, MO2 are defined ! #
!     henry(ind_OH,k)=25.d0      !RS_MOCCA
     henry(ind_OH,k)=func3(3.0d1,4300.d0)      !RS_MOCCA 06.04.00
!     henry(ind_HO2,k)=9.d3      !RS_MOCCA
     henry(ind_HO2,k)=func3(3.9d3,5900.d0)     !RS_MOCCA 06.04.00
     henry(ind_MO2,k)=func3(6.d0,5600.d0) !RS_Jacob
!     henry(ind_MO2,k)=6.d0      !RS_MOCCA
!!     henry(ind_IO,k)=4.5d2      !RS_MOCCA
     henry(ind_IO,k)=func3(4.5d2,5862.d0)      !RS_MOCCA 06.04.00
     henry(ind_CO2,k)=func3(3.1d-02,2423.d0)   !RS_MOCCA
     henry(ind_CO,k)=func3(9.9d-04,1300.d0)    !RS_MOCCA
     henry(ind_O2,k)=func3(1.3d-3,1500.d0)     !RS_Lide95 06.04.00
!     henry(ind_O2,k)=1.7d-3                   !RS_MOCCA
!     henry(ind_I2O4,k)=func3()
!     henry(ind_I2O5,k)=func3()
!     henry(ind_INO,k)=func3()
!     henry(ind_Br2O,k)=func3()
     henry(ind_ClONO,k)=4.6d-2  !RS_MOCCA, same as ClNO2
!     henry(ind_ClO3,k)=func3()
!     henry(ind_Cl2O3,k)=func3()
     henry(ind_CH3OH,k)=func3(1.6d2,5600.d0)   !RS_MOCCA
     henry(ind_C2H5OH,k)=func3(1.5d2,6400.d0)  !RS_MOCCA
     henry(ind_H2,k)=func3(7.8d-4,500.d0)      !RS_MOCCA
!     henry(ind_NHS,k)=func3()
!     henry(ind_RCl,k)=func3()
!     henry(ind_RBr,k)=func3()
     henry(ind_XOR,k)=func3(1.5d2,6400.d0) ! same as ethanol
     henry(ind_SOR,k)=func3(1.5d2,6400.d0) !  same as ethanol
!     henry(ind_SPAN,k)=func3()


! include HERE call to get effective henry coeffs (if needed),
! before inverse k_H are calculated
!     xCO2=func3(4.3d-7,-919.d0)
!     x3CO2=xCO2*func3(4.7d-11,-1787.d0)
!     xNH3=func3(1.71d-5,-4325.d0)/func3(1.d-14,-6716.d0)
!     if (sion1(1,1,k).gt.0.) then         !anhaengig von tropfenklasse! (--> pH)
!        xhp=(dmax(sion1(1,1,k),1.d-6))*1.d-3                   ! in mol/l !
!        henry(ind_CO2,k)=henry(ind_CO2,k)*(1+xCO2/xhp+x3CO2/(xhp**2))
!        henry(ind_NH3,k)=henry(ind_NH3,k)*(1+xNH3*xhp)
!     endif

  enddo

! unit: mol/(l*atm) --> mol(aq)/m3(aq) / mol(g)/m3(g) (i.e. dimensionless)
! FCT=1.d3*8.3145*T/p_0=0.082*T
! PLUS conversion to inverse Henry constant:  h_(H,inv)^cc = 1/k_H^cc
! i.e. mol(aq)/m3(aq) / mol(g)/m3(air) -->  mol(g)/m3(air) / mol(aq)/m3(aq)

  do k=1,nmaxf
     FCT = 0.0820577_dp * tt(k)
     do j=1,NSPEC
        if (henry(j,k).gt.0._dp) then
           henry(j,k) = 1._dp / (henry(j,k)*FCT)
! "else": henry=0 <=> k_H^cc=infinity
        end if
     end do
  end do

end subroutine henry_t


!
!------------------------------------------------------
!

subroutine henry_a (tt,nmaxf)

! Description :
! -----------
  ! Temp dependent Henry constants
  ! inverse dimensionless Henry constant: k_(H,inv)^cc:=1/(k_H^cp*RT)
  ! in equilibrium: XXXaq = k_H^cp * LWC * XXXg

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
! jjb work done:
!     - removed pp (pressure) from the argument list: unused
!     - removed hard coded parameters, use modules instead
!     - improved computation efficiency: Tfact calculated only once per layer
!     - missing declarations and implicit none
!     - final conversion (inverse if henry /= 0.) optimised
!     - reindexed henry(NSPEC,nf) instead of (nf,NSPEC) for computing efficiency
!     - added initialisation of henry(:,:)

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'aer_Parameters.h'    !additional common blocks and other definitions

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: tt(n) ! temperature array
  integer       , intent(in) :: nmaxf ! max layer index where henry has to be computed

! Local scalars:
  integer :: j,k                  ! loop indexes
  real (kind=dp) :: func3, a0, b0 ! temperature dependency function, and its arguments
  real (kind=dp) :: FCT           ! conversion factor, see below
  real (kind=dp) :: Tfact         ! Tfact, local variable [K**-1]
!  real (kind=dp) :: xCO2, x3CO2, xNH3, xhp

! Common blocks:
  common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
       xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
  real (kind=dp) henry, xkmt, xkef, xkeb

! 0.082=8.3145*10^3/101325=R/p_0*10^3 : R includes conversion from M/atm --> mol/(m^3*Pa)
! so k_H^cp is taken as M/atm (the most common literature unit): SEE END OF SR

! func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))*0.082*tt(k0)
! func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))

! jjb 11/02/2017 slight improvement in calculation efficiency: compute Tfact only once per layer
!  Tfact = 1/T - 1/Tref, see below
  func3(a0,b0)=a0*exp(b0*Tfact)

! == End of declarations =======================================================

  henry(:,:) = 0._dp

  do k=1,nmaxf

!     henry(ind_H2SO4,k)=4.1d-12   !???
     henry(ind_H2SO4,k)=1.d+16  !NIST --> RS_Gmitro_Vermeulen
     henry(ind_CH4,k)=1.3d-3    !RS_Mackay
     henry(ind_C2H6,k)=2.0d-3   !RS_Mackay
!     henry(ind_C3H8,k)=1.4d-3   !RS_Mackay
!     henry(ind_ALKA,k)=9.5d-4   !ok
     henry(ind_ETHE,k)=4.9d-3   !ok
!     henry(ind_ALKE,k)=4.9d-3   !ok
!     henry(ind_AROM,k)=1.5d-1   !RS_Mackay: methylbenzene
!     henry(ind_KET,k)=10.d0     !RS_higher ketone
!     henry(ind_CRES,k)=8.2d2    !RS_Hine
!     henry(ind_DIAL,k)=10.0d0   !RS_Snider:propenal #?
!     henry(ind_GLYX,k)=3.6d5    !RS_Zhou
!     henry(ind_NH4NO3,k)=4.1d-12  !# was soll dieses Salz in Gasphase??
!     henry(ind_RAN2,k)=1.2d0    !RS_Kames: pentyl-nitrate
!     henry(ind_RAN1,k)=1.d0     !RS_hine
!     henry(ind_N2O5,k)=0.0      !RS_Sander  infinity!
!     henry(ind_ClNO2,k)=0.0     !  10^-2
!     henry(ind_ClNO3,k)=0.0     !  infinity
!     henry(ind_BrNO2,k)=0.0     !  10^-1
!     henry(ind_BrNO3,k)=0.0     !  infinity
     henry(ind_HI,k)=0.0d0        !   (dissociation)
     henry(ind_I2O2,k)=0.0d0      !
     henry(ind_INO2,k)=0.0d0      !
     henry(ind_INO3,k)=0.0d0      !
!     henry(ind_I2O,k)=0.0d0      !
!     henry(ind_I2O3,k)=0.0d0      !
!     henry(ind_OIO,k)  =  unknown but not needed as only surface reaction and no reversible uptake
!     henry(ind_HIO3,k) =                     -"-
     henry(ind_C3H7I,k)=1.1d-1  !RS_Hine
!     henry(ind_CH3SO2,k)= !RS_MOCCA
!     henry(ind_CH3SO3,k)= !RS_MOCCA
!     henry(ind_Hg,k)    = 1.3d-1                   ! #233 in #3127
!     henry(ind_HgO,k)   = 2.69d12                  ! #233 in #3127
!     henry(ind_HgCl,k)  = 2.75d6                   ! assumed HgCl2 (poss. higher??)
!     henry(ind_HgCl2,k) = 2.75d6                   ! #233 in #3127
!     henry(ind_HgBr,k)  = 2.75d6                   ! assumed HgBr2 (poss. higher??)
!     henry(ind_HgBr2,k) = 2.75d6                   ! #3127

! explicitly Temp dependent
     ! Tfact = 1/T - 1/Tref  with Tref = 298.15 K
     ! 1/298.15 = 3.3540d-3
     Tfact = 1.d0/tt(k) - 3.3540d-3

     henry(ind_NO,k)=func3(1.9d-03,1480.d0)   !RS_Lide
!     henry(ind_NO2,k)=func3(1.d-02,2500.d0)  !RS_Chameides
     henry(ind_NO2,k)=func3(6.4d-03,2500.d0)  !RS_MOCCA
!     henry(ind_HNO3,k)=func3(1.66d5,8694.d0) !RS_MOCCA
     henry(ind_HNO3,k)=func3(2.5d6/1.5d1,8694.d0) !RS_MOCCA_exakt
!     henry(ind_HNO4,k)=1.4d4     !Goetz, 1996, cited in Warneck, 1999, #695
     henry(ind_HNO4,k)=func3(1.2d4,6900.d0)   !06.04.00
     henry(ind_NH3,k)=func3(58.d0,4085.d0)    !RS_MOCCA
     henry(ind_SO2,k)=func3(1.2d0,3120.d0)    !RS_MOCCA
     henry(ind_O3,k)=func3(1.2d-02,2560.d0)   !RS_MOCCA
     henry(ind_ACO2,k)=func3(3.7d+03,5700.d0) !RS_MOCCA
     henry(ind_ACTA,k)=func3(4.1d+03,6300.d0) !RS_Johnson
     henry(ind_HCHO,k)=func3(7.0d+03,6425.d0) !RS_MOCCA
     henry(ind_ALD2,k)=func3(1.3d+01,5700.d0)  !RS_Benkelberg: acetaldehyde
     henry(ind_H2O2,k)=func3(1.d+05,6338.d0)   !RS_MOCCA
!     henry(ind_ROOH,k)=func3(7.45d+04,6620.d0) ! # aktualisieren
     henry(ind_ROOH,k)=func3(3.0d+02,5322.d0)  !RS_MOCCA
!     henry(ind_HONO,k)=func3(5.0d+01,4900.d0) !RS_Becker
     henry(ind_HONO,k)=func3(4.9d+01,4780.d0)  !RS_MOCCA
     henry(ind_PAN,k)=func3(2.8d0,6500.d0)     !RS_Kames
!     henry(ind_TPAN,k)=henry(22,k)
!     henry(ind_MGLY,k)=func3(3.7d+03,7553.d0) !RS_Betterton
!     henry(ind_HCl,k)=func3(1.17d0,9001.d0)   !RS_MOCCA
     henry(ind_HCl,k)=func3(2.d0/1.7d0,9001.d0) !RS_MOCCA_exakt
!     henry(ind_R3N2,k)=func3(1.d0,5450.d0)    !RS_kames: propyl nitrate
     henry(ind_NO3,k)=func3(2.d0,2000.d0)      !RS_MOCCA
     henry(ind_DMS,k)=func3(4.8d-1,3100.d0)    !RS_deBruyn
!     henry(ind_DMSO,k)=5.d4                   !RS_MOCCA
     henry(ind_DMSO,k)=func3(5.d4,6425.d0)     !RS_MOCCA 06.04.00
     henry(ind_DMSO2,k)= 1.d+16                !DMSO2=H2SO4, assumed
     henry(ind_CH3SO2H,k)=1.d+16               !MSIA=H2SO4, assumed
     henry(ind_CH3SO3H,k)=1.d+16               !MSA=H2SO4, assumed
     henry(ind_HOCl,k)=func3(6.7d2,5862.d0)    !RS_MOCCA
!     henry(ind_Cl2,k)=9.2d-2    !RS_MOCCA
     henry(ind_Cl2,k)=func3(9.1d-2,2500.d0)    !RS_MOCCA 06.04.00
     henry(ind_HBr,k)=func3(1.3d0,10239.d0)    !RS_MOCCA
     henry(ind_Br2,k)=func3(7.6d-1,4094.d0)    !RS_MOCCA
     henry(ind_BrCl,k)=func3(9.4d-1,5600.d0)   !RS_MOCCA
!     henry(ind_HOBr,k)=9.3d1    !RS_MOCCA
     henry(ind_HOBr,k)=func3(9.3d1,5862.d0)    !RS_MOCCA 06.04.00
!!     henry(ind_HI,k)=func3(2.5d9/K_a,9800.d0) !RS_Brimblecombe #K_a
     henry(ind_I2,k)=func3(3.d0,4431.d0) !RS_MOCCA
!!     henry(ind_HOI,k)=4.5d2     !RS_MOCCA
     henry(ind_HOI,k)=func3(4.5d2,5862.d0)     !RS_MOCCA 06.04.00
!!     henry(ind_ICl,k)=1.1d2     !RS_MOCCA
     henry(ind_ICl,k)=func3(1.1d2,5600.d0)     !RS_MOCCA 06.04.00
!!     henry(ind_IBr,k)=2.4d1     !RS_MOCCA
     henry(ind_IBr,k)=func3(2.4d1,5600.d0)     !RS_MOCCA 06.04.00
     henry(ind_CH3I,k)=func3(1.4d-1,4300.d0)   !RS_Moore
     henry(ind_CH2I2,k)=func3(2.3d0,5000.d0)   !RS_Moore
     henry(ind_CH2ClI,k)=func3(8.9d-1,4300.d0) !RS_Moore

! for radicals only OH, HO2, MO2 are defined ! #
!     henry(ind_OH,k)=25.d0      !RS_MOCCA
     henry(ind_OH,k)=func3(3.0d1,4300.d0)      !RS_MOCCA 06.04.00
!     henry(ind_HO2,k)=9.d3      !RS_MOCCA
     henry(ind_HO2,k)=func3(3.9d3,5900.d0)     !RS_MOCCA 06.04.00
     henry(ind_MO2,k)=func3(6.d0,5600.d0) !RS_Jacob
!     henry(ind_MO2,k)=6.d0      !RS_MOCCA
!!     henry(ind_IO,k)=4.5d2      !RS_MOCCA
     henry(ind_IO,k)=func3(4.5d2,5862.d0)      !RS_MOCCA 06.04.00
     henry(ind_CO2,k)=func3(3.1d-02,2423.d0)   !RS_MOCCA
     henry(ind_CO,k)=func3(9.9d-04,1300.d0)    !RS_MOCCA
     henry(ind_O2,k)=func3(1.3d-3,1500.d0)     !RS_Lide95 06.04.00
!     henry(ind_O2,k)=1.7d-3                   !RS_MOCCA
!     henry(ind_I2O4,k)=func3()
!     henry(ind_I2O5,k)=func3()
!     henry(ind_INO,k)=func3()
!     henry(ind_Br2O,k)=func3()
     henry(ind_ClONO,k)=4.6d-2  !RS_MOCCA, same as ClNO2
!     henry(ind_ClO3,k)=func3()
!     henry(ind_Cl2O3,k)=func3()
     henry(ind_CH3OH,k)=func3(1.6d2,5600.d0)   !RS_MOCCA
     henry(ind_C2H5OH,k)=func3(1.5d2,6400.d0)  !RS_MOCCA
     henry(ind_H2,k)=func3(7.8d-4,500.d0)      !RS_MOCCA
!     henry(ind_NHS,k)=func3()
!     henry(ind_RCl,k)=func3()
!     henry(ind_RBr,k)=func3()
     henry(ind_XOR,k)=func3(1.5d2,6400.d0) ! same as ethanol
     henry(ind_SOR,k)=func3(1.5d2,6400.d0) !  same as ethanol
!     henry(ind_SPAN,k)=func3()


! include HERE call to get effective henry coeffs (if needed),
! before inverse k_H are calculated
!     xCO2=func3(4.3d-7,-919.d0)
!     x3CO2=xCO2*func3(4.7d-11,-1787.d0)
!     xNH3=func3(1.71d-5,-4325.d0)/func3(1.d-14,-6716.d0)
!     if (sion1(1,1,k).gt.0.) then         !anhaengig von tropfenklasse! (--> pH)
!        xhp=(dmax(sion1(1,1,k),1.d-6))*1.d-3                   ! in mol/l !
!        henry(ind_CO2,k)=henry(ind_CO2,k)*(1+xCO2/xhp+x3CO2/(xhp**2))
!        henry(ind_NH3,k)=henry(ind_NH3,k)*(1+xNH3*xhp)
!     endif

  enddo

! unit: mol/(l*atm) --> mol(aq)/m3(aq) / mol(g)/m3(g) (i.e. dimensionless)
! FCT=1.d3*8.3145*T/p_0=0.082*T
! PLUS conversion to inverse Henry constant:  h_(H,inv)^cc = 1/k_H^cc
! i.e. mol(aq)/m3(aq) / mol(g)/m3(air) -->  mol(g)/m3(air) / mol(aq)/m3(aq)

  do k=1,nmaxf
     FCT = 0.0820577_dp * tt(k)
     do j=1,NSPEC
        if (henry(j,k).gt.0._dp) then
           henry(j,k) = 1._dp / (henry(j,k)*FCT)
! "else": henry=0 <=> k_H^cc=infinity
        end if
     end do
  end do

end subroutine henry_a


!
!------------------------------------------------------
!

subroutine cw_rc (nmaxf)

! Description :
! -----------
  ! mean radius and LWC for "chemical" particles size bins

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
! jjb work done
!     cleaning
!     modules, instead of hard coded parameters
!     rc has to be always defined
!     useless test removed to set cm (case rH>rH_deliq : whatever cloud value is ok to define cm)
!     inconsistency between .lt. and .gt. : case "==" was thus missing. corrected.
!     implicit none, and all missing declarations
!     ial definition only once, out of the main do loop for computing efficiency

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       ifeed

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer, intent(in) :: nmaxf ! max index for calculation  ! jjb might add a test on it, error if > n, warning if > nf

! Local parameters:
  real (kind=dp), parameter :: xpi = 4._dp / 3._dp * pi

!  real (kind=dp), parameter :: cwm=1.d-05  ! jjb old value, already commented in v741
  real (kind=dp), parameter :: cwm = 1.d-1 ! Threshold for switching chemistry in bins 1 & 2 ("aerosols")
  real (kind=dp), parameter :: cwmd = 1.d2 ! Threshold for switching chemistry in bins 3 & 4 ("droplets")
                                           !   here "d" stands for droplets, not dry !

! Local scalars:
  real (kind=dp) :: cm1, cm2, cm3, cm4
  real (kind=dp) :: cw1, cw2, cw3, cw4
  real (kind=dp) :: rc1, rc2, rc3, rc4
  real (kind=dp) :: x0, x1
  integer :: ia, ial, jt ! loop indexes for 2D particle grid
  integer :: k  ! index for grid number in do loops

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka

  common /blck11/ rc(nkc,n)
  real (kind=dp) :: rc

  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm

  common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
  real (kind=dp) :: conv2

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), & ! only rq and e are used
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw, ew, rn, rw, en, e, dew, rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n) ! only ff is used
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n) ! only feu is used
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

  common /kpp_l1/ cloud(nkc,n)
  logical :: cloud

  common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss
  real (kind=dp) :: xcryssulf,xcrysss,xdelisulf,xdeliss

  common /kinv_i/ kinv
  integer :: kinv

! == End of declarations =======================================================

! rc(nkc,n): mean radius of droplets in m
! cw(nkc,n): LWC in given radius range in m^3(aq)g/m^3(air)

! Define lower bound depending on nucleation settings
  if (ifeed.eq.2) then
     ial = 2
  else
     ial = 1
  endif

! Main loop
  do k=2,nmaxf

! OLD options to skip some layers
!         if (k.lt.lcl.or.k.gt.lct) go to 1000
!         if (feu(k).lt.xcryssulf.and.feu(k).lt.xcrysss) go to 1000
! see SR kon: if rH>0.7 microphysics is calculated, this
! is also start for liquid aerosol chemistry

! Initialisation
     rc1=0._dp ; rc2=0._dp ; rc3=0._dp ; rc4=0._dp
     cw1=0._dp ; cw2=0._dp ; cw3=0._dp ; cw4=0._dp
     cm1=0._dp ; cm2=0._dp ; cm3=0._dp ; cm4=0._dp

! here TOTAL particle volume is used for calculating LWC
! this is correct only for completely soluble aerosol
!
! the cw variables and therefore also cvv?/conv2 use
! the volume of solution (cvv? converts to mol/l_solution)
! for the calculation of the activity coefficients the
! molality is needed (mol/kg_solvent) ==> cm is only water
! volume (m^3 water/m^3 air)


! small aerosol
     do ia=ial,ka
        do jt=1,kw(ia)
           x0  = ff(jt,ia,k) * xpi * rq(jt,ia)**3
           cw1 = cw1 + x0
           rc1 = rc1 + x0 * rq(jt,ia)
           x1  = ff(jt,ia,k) * e(jt)
           cm1 = cm1 + x1
        enddo
! small droplets
        do jt=kw(ia)+1,nkt
           x0  = ff(jt,ia,k) * xpi * rq(jt,ia)**3
           cw3 = cw3 + x0
           rc3 = rc3 + x0 * rq(jt,ia)
           x1  = ff(jt,ia,k) * e(jt)
           cm3 = cm3 + x1
        enddo
     enddo
! large aerosol
     do ia=ka+1,nka
        do jt=1,kw(ia)
           x0  = ff(jt,ia,k) * xpi * rq(jt,ia)**3
           cw2 = cw2 + x0
           rc2 = rc2 + x0 * rq(jt,ia)
           x1  = ff(jt,ia,k) * e(jt)
           cm2 = cm2 + x1
        enddo
! large droplets
        do jt=kw(ia)+1,nkt
           x0  = ff(jt,ia,k) * xpi * rq(jt,ia)**3
           cw4 = cw4 + x0
           rc4 = rc4 + x0 * rq(jt,ia)
           x1  = ff(jt,ia,k) * e(jt)
           cm4 = cm4 + x1
        enddo
     enddo

! conversion: um^3/cm^3 --> m^3(aq)/m^3(air):10^-12
!           : um        --> m               :10^-6
! aerosol: cw must be greater than 1.d-13  => cwm
! droplet: cw must be greater than 1.d-10  => cwmd
! conversion for cm: mg(water)/cm^3(air)
!                       --> m^3(wat)/m^3(air):10^-3

! define rc for all rH because rc is used in calculation of particle
! chemistry sedimentation
     if (cw1.gt.0._dp) then
        rc(1,k) = rc1 / cw1 * 1.d-6
     else
        rc(1,k) = 0._dp
     end if
     if (cw2.gt.0._dp) then
        rc(2,k) = rc2 / cw2 * 1.d-6
     else
        rc(2,k) = 0._dp
     end if
     if (cw3.gt.0._dp) then
        rc(3,k)=rc3 / cw3 * 1.d-6
     else
        rc(3,k) = 0._dp
     end if
     if (cw4.gt.0._dp) then
        rc(4,k) = rc4 / cw4 * 1.d-6
     else
        rc(4,k) = 0._dp
     end if

     cw(1,k) = cw1 * 1.d-12
     cw(2,k) = cw2 * 1.d-12
     cw(3,k) = cw3 * 1.d-12
     cw(4,k) = cw4 * 1.d-12


! cm (old:cw) is used as switch for aerosol chemistry therefore:
!     off: below crys point
!     on : only if above threshold value and crys rH (if cloud was true) or
!          deli rH (if cloud was false)

     if (feu(k).lt.min(xcryssulf,xcrysss)) then
        if (k.le.kinv) print*,k,feu(k),' below both crystal. points'
        cm(:,k)    = 0._dp
        conv2(:,k) = 0._dp

     else

        if ( (cw1.ge.cwm) .and. &                        ! cw1 has to be > cwm, and...
             ((cloud(1,k).and.feu(k).ge.xcryssulf).or. & !  ... if cloud is true, rH >= rH_crys is enough
             (feu(k).ge.xdelisulf)) ) then               !  ... else (if cloud is false), rH >= rH_deliq is necessary
                                                         !         (cloud can be true as well in the later case)
           cm(1,k)    = cm1 * 1.d-3
           conv2(1,k) = 1.d9 / cw1 ! 1.d-3/cw(1,k)
        else
           cm(1,k)    = 0._dp
           conv2(1,k) = 0._dp
        endif

        if ( (cw2.ge.cwm) .and. &
             ((cloud(2,k).and.feu(k).ge.xcrysss).or. &   ! same comments
             (feu(k).ge.xdeliss)) ) then
           cm(2,k)    = cm2 * 1.d-3
           conv2(2,k) = 1.d9 / cw2
        else
           cm(2,k)    = 0._dp
           conv2(2,k) = 0._dp
        endif

        if (cw3.ge.cwmd) then
           cm(3,k)    = cm3 * 1.d-3
           conv2(3,k) = 1.d9 / cw3
        else
           cm(3,k)    = 0._dp
           conv2(3,k) = 0._dp
        endif

        if (cw4.ge.cwmd) then
           cm(4,k)    = cm4 * 1.d-3
           conv2(4,k) = 1.d9 / cw4
        else
           cm(4,k)    = 0._dp
           conv2(4,k) = 0._dp
        endif

     end if

  end do

end subroutine cw_rc


!
!--------------------------------------------------------
!

subroutine fast_k_mt_t (freep,box,n_bl)  !_1D

! Description :
! -----------
  ! transfer coefficient after Schwarz, 1986 (see Sander & Crutzen '96, JGR, 9127)
  ! but no mean values used (like in SR k_mt_a/t) but integrated values

  ! transfer coefficients (variable xkmt) computed only if LWC above threshold for each layer and each chemical bin (with cm used as switch)
  ! sedimentation velocity (variable vt) computed whatever LWC

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb: added all missing declarations and implicit none
  !      shorten the code by introducing lmax and llchem, to avoid duplicated part of code
  !       for vt calculation ("dry" case, meaning cm=0)

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       ifeed

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nka, &
       nkt, &
       nkc

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!      USE gas_common, ONLY : &
!! Imported Parameters:
!           j1, &
!           j5, &
!! Imported Array Variables with intent (in):
!           gas_k2m_t, &
!           rad_k2m_t, &
!           vm=>vmean
!
!      USE kpp_tot_Global, ONLY : &
!           SPC_NAMES

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'tot_Parameters.h' !additional common blocks and other definitions

! Subroutine arguments:
  real(kind=dp), intent(in) :: freep(nf)
  logical, intent(in) :: box
  integer, intent(in) :: n_bl

! Local parameters:
  integer, parameter :: nx=50                      ! number of exchanged chemical species
  real (kind=dp), parameter :: z4pi3 = 4._dp * pi / 3._dp
! Local scalars:
  integer :: ia, iia_0, iia_e, jt, jjt_0, jjt_e    ! running indexes and limits for ia, jt
  integer :: k, kc, l, lmax, nmin, nmax            ! running indexes, and limits for k
  logical :: llchem                                ! switch to compute xkmt (cm>0) or not (cm=0)
  real (kind=dp) :: rqq, x1, xk1, xx1, x2, xvs
  real (kind=dp), external :: vterm                ! terminal velocity (m/s) computed as a function of r, t, p
  ! Local arrays:
  integer :: lex(nx)
  real (kind=dp) :: rqm(nkt,nka)

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real(kind=dp) :: cw, cm

!  common /cb40/ time,lday,lst,lmin,it,lcl,lct ! jjb only time is used, for potential error message
!  real (kind=dp) :: time
!  integer :: lday, lst, lmin, it, lcl, lct

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw, ew, rn, rw, en, e, dew, rq
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
  real(kind=dp) :: alpha, vmean
  common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
       xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
  real(kind=dp) :: henry, xkmt, xkef, xkeb
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real(kind=dp) :: vt, vd, vdm

! == End of declarations =======================================================

! Define the list of species indexes that are exchanged
  data lex /ind_NO2,ind_HNO3,ind_NH3,ind_SO2,ind_H2SO4,ind_O3, &
          ind_ACO2,ind_HCHO,ind_H2O2,ind_HONO,ind_HCl,ind_N2O5,ind_HNO4, &
          ind_NO3,ind_OH,ind_HO2,ind_MO2,ind_CO2,ind_O2,ind_ROOH, &
          ind_HOCl,ind_Cl2,ind_HBr,ind_HOBr,ind_Br2,ind_BrCl,ind_DMSO, &
          ind_ClNO3,ind_BrNO3,ind_CH3SO3H,ind_DMS,ind_CH3SO2H,ind_DMSO2, &
          ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr,ind_OIO,ind_INO2, &
          ind_INO3,ind_HI,ind_I2O2,ind_HIO3,ind_NO,ind_ACTA,ind_CH3OH, &
          ind_C2H5OH,ind_XOR,ind_SOR/
!    &    1,2,3,4,5/           !ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr/
! ind_Hg,ind_HgO,ind_HgCl,ind_HgCl2,ind_HgBr,ind_HgBr2; also change/check nx

! change rq in um to rqm in m
  rqm(:,:) = rq(:,:) * 1.d-6

! Define min/max layer index where transfer coefficients are computed
  if (box) then
     nmin=n_bl
     nmax=n_bl
  else
     nmin=2
     nmax=nf     
  endif

! loop over vertical grid
  do k=nmin,nmax
! loop over the nkc different chemical bins
     do kc=1,nkc

        ! Compute xkmt and vt (if cm > 0) or vt only (if cm = 0)
        if (cm(kc,k).gt.0._dp) then  ! switch changed from cw
           lmax = nx
           llchem = .true.
        else
           lmax = 1
           llchem = .false.
        end if

! loop over the species to be exchanged between gas and aqueous phase---
        do l=1,lmax

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!! jjb search for desired species index, from lex table, in gas_m2k_t table
!               jspec = 1
!               ! search in non radical gas list
!               do while(gas_k2m_t(jspec) /= lex(l) .and. jspec+1 <= j1)
!                  jspec = jspec+1
!               end do
!
!               ! search in radical gas list
!               if(gas_k2m_t(jspec) /= lex(l)) then
!               jspec = 1
!               do while(rad_k2m_t(jspec) /= lex(l) .and. jspec+1 <= j5)
!                  jspec = jspec+1
!               end do
!
!               if(rad_k2m_t(jspec) /= lex(l)) then
!                  if(trim(SPC_NAMES(lex(l))) == 'O2') then
!                     jspec = j1+j5+1
!                  else
!                     if(time<121.d0 .and. k==nmin.and.kc==1) then
!                        print*,"in fast_k_mt_t, error"
!                        print*,spc_names(lex(l))," kpp index ",lex(l)
!                     end if
!                     cycle
!                     !stop "stopped by SR fast_k_mt_t"
!                  end if
!               else
!                  jspec = jspec + j1 ! add offset in vmean array (radical case)
!               end if
!
!               end if
!! end jjb

! define summation limits (1) ---
           if (kc.eq.1.or.kc.eq.3) then
              if (ifeed.eq.2) then
                 iia_0 = 2
              else
                 iia_0 = 1
              endif
              iia_e = ka
           endif
           if (kc.eq.2.or.kc.eq.4) then
              iia_0 = ka+1
              iia_e = nka
           endif

! fast version without logarithmic integration
           ! Initialisation of local variables
           if (llchem) then
              x1  = 0._dp
              xk1 = 0._dp
           end if
           if (l.eq.1) xx1 = 0._dp
           if (alpha(lex(l),k).gt.0._dp) x1=4./(3.*alpha(lex(l),k))

           do ia=iia_0,iia_e
! define summation limits (2)
              if (kc.eq.1.or.kc.eq.2) then
                 jjt_0 = 1
                 jjt_e = kw(ia)
              endif
              if (kc.eq.3.or.kc.eq.4) then
                 jjt_0 = kw(ia)+1
                 jjt_e = nkt
              endif

              do jt=jjt_0,jjt_e
! conversion: um      --> m               :10^-6
                 rqq = rqm(jt,ia)
! kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
!     with D_g=lambda*vmean/3.
! here a volume weighted value is calculated, therefore weighting with r^3:
! kmt=4/3*pi/L*sum(a)*sum(r){r^3*N*kt}

! < jjb 21-12-2016
                 if (llchem) then
                    x2 = vmean(lex(l),k) / (rqq/freep(k)+x1) ![1/s]
!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!                     x2=vm(jspec,k)/(rqq/freep(k)+x1) ![1/s]
! jjb 21-12-2016 >

! conversion: 1/cm^3 --> 1/m^3(air):10^6
                    xk1 = xk1 + x2 * rqq*rqq * ff(jt,ia,k) * 1.d6
                 end if
! LWC weighted sedimentation velocity
                 if (l.eq.1) then
                    xvs = vterm(rqq,t(k),p(k))
                    xx1 = xx1 + rqq*rqq*rqq * xvs * ff(jt,ia,k) * 1.e6_dp
                 endif
              enddo !jt
           enddo !ia

! k_mt=4*pi/(3*LWC)*sum
           if (cw(kc,k).gt.0._dp) then
              if (llchem) xkmt(k,kc,lex(l)) = z4pi3 / cw(kc,k) * xk1 ![1/s]
! sedimentation velocity:
              if (l.eq.1) vt(kc,k) = z4pi3 / cw(kc,k) * xx1
           end if

        enddo ! l

     enddo    ! kc
  end do      ! k

end subroutine fast_k_mt_t


!
!--------------------------------------------------------
!

subroutine fast_k_mt_a (freep,box,n_bl)  !_1D

! Description :
! -----------
  ! transfer coefficient after Schwarz, 1986 (see Sander & Crutzen '96, JGR, 9127)
  ! but no mean values used (like in SR k_mt_a/t) but integrated values

  ! transfer coefficients (variable xkmt) computed only if LWC above threshold for each layer and each chemical bin (with cm used as switch)
  ! sedimentation velocity (variable vt) computed whatever LWC

! Author :
! ------
  ! Roland von Glasow

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb: added all missing declarations and implicit none
  !      shorten the code by introducing lmax and llchem, to avoid duplicated part of code
  !       for vt calculation ("dry" case, meaning cm=0)
  ! jjb: bugfix with kc, must compute vt up to nkc_l (for use in SR sedl), not only 1,2

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       ifeed,        &
       nkc_l

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nka, &
       nkt, &
       nkc

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!      USE gas_common, ONLY : &
!! Imported Parameters:
!           j1, &
!           j5, &
!! Imported Array Variables with intent (in):
!           gas_k2m_a, &
!           rad_k2m_a, &
!           vm=>vmean
!
!      USE kpp_aer_Global, ONLY : &
!           SPC_NAMES

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'aer_Parameters.h' !additional common blocks and other definitions

! Subroutine arguments:
  real(kind=dp), intent(in) :: freep(nf)
  logical, intent(in) :: box
  integer, intent(in) :: n_bl

! Local parameters:
  integer, parameter :: nx=50                      ! number of exchanged chemical species
  real (kind=dp), parameter :: z4pi3 = 4._dp * pi / 3._dp
! Local scalars:
  integer :: ia, iia_0, iia_e, jt, jjt_0, jjt_e    ! running indexes and limits for ia, jt
  integer :: k, kc, l, lmax, nmin, nmax            ! running indexes, and limits for k
  logical :: llchem                                ! switch to compute xkmt (cm>0) or not (cm=0)
  real (kind=dp) :: rqq, x1, xk1, xx1, x2, xvs
  real (kind=dp), external :: vterm                ! terminal velocity (m/s) computed as a function of r, t, p
  ! Local arrays:
  integer :: lex(nx)
  real (kind=dp) :: rqm(nkt,nka)

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real(kind=dp) :: cw, cm

!  common /cb40/ time,lday,lst,lmin,it,lcl,lct ! jjb only time is used, for potential error message
!  real (kind=dp) :: time
!  integer :: lday, lst, lmin, it, lcl, lct

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw, ew, rn, rw, en, e, dew, rq
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
  real(kind=dp) :: alpha, vmean
  common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
       xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
  real(kind=dp) :: henry, xkmt, xkef, xkeb
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real(kind=dp) :: vt, vd, vdm

! == End of declarations =======================================================

! Define the list of species indexes that are exchanged
  data lex /ind_NO2,ind_HNO3,ind_NH3,ind_SO2,ind_H2SO4,ind_O3, &
          ind_ACO2,ind_HCHO,ind_H2O2,ind_HONO,ind_HCl,ind_N2O5,ind_HNO4, &
          ind_NO3,ind_OH,ind_HO2,ind_MO2,ind_CO2,ind_O2,ind_ROOH, &
          ind_HOCl,ind_Cl2,ind_HBr,ind_HOBr,ind_Br2,ind_BrCl,ind_DMSO, &
          ind_ClNO3,ind_BrNO3,ind_CH3SO3H,ind_DMS,ind_CH3SO2H,ind_DMSO2, &
          ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr,ind_OIO,ind_INO2, &
          ind_INO3,ind_HI,ind_I2O2,ind_HIO3,ind_NO,ind_ACTA,ind_CH3OH, &
          ind_C2H5OH,ind_XOR,ind_SOR/
!    &    1,2,3,4,5/           !ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr/
! ind_Hg,ind_HgO,ind_HgCl,ind_HgCl2,ind_HgBr,ind_HgBr2; also change/check nx

! change rq in um to rqm in m
  rqm(:,:) = rq(:,:) * 1.d-6

! Define min/max layer index where transfer coefficients are computed
  if (box) then
     nmin=n_bl
     nmax=n_bl
  else
     nmin=2
     nmax=nf     
  endif

! loop over vertical grid
  do k=nmin,nmax
! loop over the nkc different chemical bins
     do kc=1,nkc_l

        ! Compute xkmt and vt (if cm > 0) or vt only (if cm = 0)
        if (cm(kc,k).gt.0._dp) then  ! switch changed from cw
           lmax = nx
           llchem = .true.
        else
           lmax = 1
           llchem = .false.
        end if

! loop over the species to be exchanged between gas and aqueous phase---
        do l=1,lmax

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!! jjb search for desired species index, from lex table, in gas_m2k_a table
!               !print*,spc_names(lex(l)),j1,j5
!               jspec = 1
!               ! search in non radical gas list
!               do while(gas_k2m_a(jspec) /= lex(l) .and. jspec+1 <= j1)
!                  jspec = jspec+1
!               end do
!
!               ! search in radical gas list
!               if(gas_k2m_a(jspec) /= lex(l)) then
!               jspec = 1
!               do while(rad_k2m_a(jspec) /= lex(l) .and. jspec+1 <= j5)
!                  jspec = jspec+1
!               end do
!
!               if(rad_k2m_a(jspec) /= lex(l)) then
!
!                  if(trim(SPC_NAMES(lex(l))) == 'O2') then
!                     jspec = j1+j5+1
!                  else
!                     if(time<121.d0 .and. k==nmin.and.kc==1) then
!                        print*,"in fast_k_mt_t, error"
!                        print*,spc_names(lex(l))," kpp index ",lex(l)
!                     end if
!                     cycle
!                     !stop "stopped by SR fast_k_mt_a"
!                  end if
!               else
!                  jspec = jspec + j1 ! add offset in vmean array (radical case)
!               end if
!
!               end if
!               !print*,"final jspec =",jspec
!! end jjb

! define summation limits (1) ---
           if (kc.eq.1.or.kc.eq.3) then
              if (ifeed.eq.2) then
                 iia_0 = 2
              else
                 iia_0 = 1
              endif
              iia_e = ka
           endif
           if (kc.eq.2.or.kc.eq.4) then
              iia_0 = ka+1
              iia_e = nka
           endif

! fast version without logarithmic integration
           ! Initialisation of local variables
           if (llchem) then
              x1  = 0._dp
              xk1 = 0._dp
           end if
           if (l.eq.1) xx1 = 0._dp
           if (alpha(lex(l),k).gt.0._dp) x1=4./(3.*alpha(lex(l),k))

           do ia=iia_0,iia_e
! define summation limits (2)
              if (kc.eq.1.or.kc.eq.2) then
                 jjt_0 = 1
                 jjt_e = kw(ia)
              endif
              if (kc.eq.3.or.kc.eq.4) then
                 jjt_0 = kw(ia)+1
                 jjt_e = nkt
              endif

              do jt=jjt_0,jjt_e
! conversion: um      --> m               :10^-6
                 rqq = rqm(jt,ia)
! kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
!     with D_g=lambda*vmean/3.
! here a volume weighted value is calculated, therefore weighting with r^3:
! kmt=4/3*pi/L*sum(a)*sum(r){r^3*N*kt}

! < jjb 21-12-2016
                 if (llchem) then
                    x2 = vmean(lex(l),k) / (rqq/freep(k)+x1) ![1/s]
!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!                     x2=vm(jspec,k)/(rqq/freep(k)+x1) ![1/s]
! jjb 21-12-2016 >

! conversion: 1/cm^3 --> 1/m^3(air):10^6
                    xk1 = xk1 + x2 * rqq*rqq * ff(jt,ia,k) * 1.d6
                 end if
! LWC weighted sedimentation velocity
                 if (l.eq.1) then
                    xvs = vterm(rqq,t(k),p(k))
                    xx1 = xx1 + rqq*rqq*rqq * xvs * ff(jt,ia,k) * 1.e6_dp
                 endif
              enddo !jt
           enddo !ia

! k_mt=4*pi/(3*LWC)*sum
           if (cw(kc,k).gt.0._dp) then
              if (llchem) xkmt(k,kc,lex(l)) = z4pi3 / cw(kc,k) * xk1 ![1/s]
! sedimentation velocity:
              if (l.eq.1) vt(kc,k) = z4pi3 / cw(kc,k) * xx1
           end if

        enddo ! l

     enddo    ! kc
  end do      ! k

end subroutine fast_k_mt_a


!
!------------------------------------------------------
!

subroutine equil_co_t (tt,nmaxf)

! Description :
! -----------
  ! equilibrium constant (see MOCCA)
  ! xkef: forward
  ! xkeb: backward
  ! activity coefficients included via xgamma
  ! acidity constants Ka = XXXaf/XXXab [mol/l];
  ! and other equilibrium constants (f=forward, b=backward reaction);
  ! converted to [mol/m3(air)];
  ! absolute values are chosen arbitrarily to ensure fast equilibration;

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb: added all missing declarations and implicit none

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       j6, &
       nf, &
       n, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'tot_Parameters.h' !additional common blocks and other definitions

! Subroutine arguments:
  real(kind=dp), intent(in) :: tt(n)
  integer, intent(in) :: nmaxf

! Local scalars:
  integer :: k, kc
  real (kind=dp) :: cv2, xsw

! Statement function
  real (kind=dp) :: a0, b0, funa
  funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))

! Common blocks:
  common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
  real (kind=dp) :: conv2
  common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
       xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
  real(kind=dp) :: henry, xkmt, xkef, xkeb
  common /kpp_mol/ xgamma(nf,j6,nkc)
  real (kind=dp) :: xgamma

! == End of declarations =======================================================


  do k=2,nmaxf
     do kc=1,nkc
        cv2 = conv2(kc,k)
        xsw = 1._dp
        if (cv2.eq.0.) xsw = 0._dp ! jjb no need to compute all the xkef & xkeb which will be 0. !
!        if (cv2.gt.0.) then       ! but still better than initialising all NSPEC list !

! absolute values of back- and forward reactions are choosen arbitrarily to
! ensure a quick (but numerically stable) equilibrium
! if illogical error messages arise, try fiddling with the absolute values
! of equilibrium reactions

        xkef(k,kc,ind_H2O) = xsw * funa(1.0d-5,-6716.d0,k) !schneller wird schlechter!
        xkeb(k,kc,ind_H2O) = 1.0D9 * cv2 * xgamma(k,1,kc) * xgamma(k,3,kc)

!        xkef(k,kc,ind_HO2) = xsw * 2.d5  !Chameides '84
        xkef(k,kc,ind_HO2) = xsw * 1.6d5  !Weinstein-loyd and Schwartz, '91
        xkeb(k,kc,ind_HO2) = 1.d10 * cv2 * xgamma(k,1,kc) * xgamma(k,11,kc)

        xkef(k,kc,ind_ACO2) = xsw * 1.8D0
        xkeb(k,kc,ind_ACO2) = 1.0D4 * cv2 * xgamma(k,1,kc) * xgamma(k,16,kc)

!        xkef(k,kc,ind_CO2) = xsw * funa(4.3d3,-913.d0,k)
!        xkeb(k,kc,ind_CO2) = 1.0D10 * cv2
        xkef(k,kc,ind_CO2) = xsw * funa(4.3d-2,-913.d0,k)
        xkeb(k,kc,ind_CO2) = 1.0D5 * cv2 * xgamma(k,1,kc) * xgamma(k,9,kc)

        xkef(k,kc,ind_HONO) = xsw * funa(5.1d+3,-1260.d0,k)
        xkeb(k,kc,ind_HONO) = 1.0D7 * cv2 * xgamma(k,1,kc) * xgamma(k,12,kc)

        xkef(k,kc,ind_HNO3) = xsw * funa(1.54d+10,8700.d0,k)
        xkeb(k,kc,ind_HNO3) = 1.0D9 * cv2 * xgamma(k,1,kc) * xgamma(k,13,kc)

        xkef(k,kc,ind_HNO4) = xsw * 2.0D3
        xkeb(k,kc,ind_HNO4) = 2.0D8 * cv2

        xkef(k,kc,ind_NH3) = xsw * funa(1.7d5,-4325.d0,k)
        xkeb(k,kc,ind_NH3) = 1.0D10 *  cv2 * xgamma(k,3,kc) * xgamma(k,2,kc)

        xkef(k,kc,ind_HSO3ml1) = xsw * funa(6.0d2,1120.d0,k) * xgamma(k,5,kc)
        xkeb(k,kc,ind_HSO3ml1) = 1.0D10 *  cv2 * xgamma(k,1,kc) * xgamma(k,6,kc)

        xkef(k,kc,ind_H2SO4) = xsw * 1.0d12 !Seinfeld, Pandis (1998), p.391
        xkeb(k,kc,ind_H2SO4) = 1.0d9 * cv2 * xgamma(k,1,kc) * xgamma(k,19,kc)

        xkef(k,kc,ind_HSO4ml1) = xsw * funa(1.02d+6,2720.d0,k) * xgamma(k,19,kc)
        xkeb(k,kc,ind_HSO4ml1) = 1.0D8 * cv2 * xgamma(k,1,kc) * xgamma(k,8,kc)

        xkef(k,kc,ind_SO2) = xsw * funa(1.7d8,2090.d0,k)
        xkeb(k,kc,ind_SO2) = 1.0D10 * cv2 * xgamma(k,1,kc) * xgamma(k,5,kc)

        xkef(k,kc,ind_HCHO) = 1.d10 * cv2  !Chameides '84   !mechanism changed
        xkeb(k,kc,ind_HCHO) = xsw * 1.d5                    !mechanism changed

        xkef(k,kc,ind_HCl) = xsw * funa(1.7d10,6896.d0,k)
        xkeb(k,kc,ind_HCl) = 1.0D4 * cv2 * xgamma(k,1,kc) * xgamma(k,14,kc)

        xkef(k,kc,ind_Cl2ml1) = xsw * 5.2d4 * xgamma(k,15,kc)    !Chameides '84
        xkeb(k,kc,ind_Cl2ml1) = 1.d10 * cv2 * xgamma(k,14,kc)

!        xkef(k,kc,ind_Cl2) = 1.1D-3
!        xkeb(k,kc,ind_Cl2) = 2.1D2 * cv2

        xkef(k,kc,ind_HOCl) = xsw * 3.2D2
        xkeb(k,kc,ind_HOCl) = 1.0D10 * cv2 * xgamma(k,1,kc) * xgamma(k,22,kc)

        xkef(k,kc,ind_HBr) = xsw * 1.0D13
        xkeb(k,kc,ind_HBr) = 1.0D4 *  cv2 * xgamma(k,1,kc) * xgamma(k,24,kc)

        xkef(k,kc,ind_Br2) = xsw * funa(2.95d4,-4068.d0,k) * xgamma(k,25,kc)  !Liu et al, 2002, #2109
        xkeb(k,kc,ind_Br2) = funa(1.17d10,-1812.d0,k) *  cv2 * xgamma(k,24,kc)

        xkef(k,kc,ind_HOBr) = xsw * funa(2.3d1,-3091.d0,k)
        xkeb(k,kc,ind_HOBr) = 1.0D10 * cv2 * xgamma(k,1,kc) * xgamma(k,26,kc)

        xkef(k,kc,ind_BrCl2ml1) = funa(5.d9,1143.d0,k) * cv2 * xgamma(k,14,kc) !#894
        xkeb(k,kc,ind_BrCl2ml1) = xsw * 1.3D9 * xgamma(k,28,kc)
! no activities used, they "cancel out"
        xkef(k,kc,ind_Br2Clml1) = 5.d9 * cv2 !*xgamma(k,24,kc) !if original equilibrium
        xkeb(k,kc,ind_Br2Clml1) = xsw * 2.8D5!*xgamma(k,29,kc) !       -"-
        xkef(k,kc,ind_Br2l1) = 5.d9 * cv2  !*xgamma(k,14,kc)    !if original equilibrium
        xkeb(k,kc,ind_Br2l1) = xsw * 3.85D9!*xgamma(k,29,kc)    !       -"-
        xkef(k,kc,ind_ICl) = 1.0D11  * cv2 * xgamma(k,14,kc)
        xkeb(k,kc,ind_ICl) = xsw * 1.3D9 * xgamma(k,37,kc)
        xkef(k,kc,ind_IBr) = 1.0D11 * cv2 * xgamma(k,24,kc)
        xkeb(k,kc,ind_IBr) = xsw * 3.5D8 * xgamma(k,38,kc)
! new (speculative) ICl <--> IBr equilibria, assumed to yield the same
! ICl/IBr ratio as in the BrCl <--> Br2 eqilibria (BrCl/Br2)
! no activities used, they are not known!! 11.04.01
        xkef(k,kc,ind_IClBrml1) = 5.d9 * cv2 !*xgamma(k,24,kc)
        xkeb(k,kc,ind_IClBrml1) = xsw * 2.8D5!*xgamma(k,29,kc)
        xkef(k,kc,ind_I2) = 5.d9 * cv2  !*xgamma(k,14,kc)
        xkeb(k,kc,ind_I2) = xsw * 3.85D9!*xgamma(k,29,kc)
!        xkef(k,kc,ind_HIO2) = xsw * 2.0D3
!        xkeb(k,kc,ind_HIO2) = 2.0D9 * cv2
        xkef(k,kc,ind_HIO3) = xsw * 1.57D4
        xkeb(k,kc,ind_HIO3) = 1.0D5 * cv2

! mercury; no activity coefficients; forward not faster than 1.d14
!
!        xkef(k,kc,ind_HgOHpl1)    = 4.27d14 * cv2           ! #493,  K_eq = 4.27d10
!        xkeb(k,kc,ind_HgOHpl1)    = xsw * 1.d4
!        xkef(k,kc,ind_HgOH2l1)    = 2.6d14 * cv2            ! #4174, K_eq = 2.6d11
!        xkeb(k,kc,ind_HgOH2l1)    = xsw * 1.d3
!        xkef(k,kc,ind_HgSO3l1)    = 5.01d14 * cv2           ! #493,  K_eq = 5.01d12
!        xkeb(k,kc,ind_HgSO3l1)    = xsw * 1.d2
!        xkef(k,kc,ind_HgSO322ml1) = 2.5d14 * cv2            ! #4158, K_eq = 2.5d11
!        xkeb(k,kc,ind_HgSO322ml1) = xsw * 1.d3
!        xkef(k,kc,ind_HgOHCll1)   = 2.69d14 * cv2            ! #4158, K_eq = 2.69d7
!        xkeb(k,kc,ind_HgOHCll1)   = xsw * 1.d7
!        xkef(k,kc,ind_HgClpl1)    = 2.0d14 * cv2             ! #4174, K_eq = 2.0d7
!        xkeb(k,kc,ind_HgClpl1)    = xsw * 1.d7
!        xkef(k,kc,ind_HgCl2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d6
!        xkeb(k,kc,ind_HgCl2l1)    = xsw * 1.d8
!        xkef(k,kc,ind_HgCl3ml1)   = 6.7d8 * cv2             ! #4174, K_eq = 6.7d0
!        xkeb(k,kc,ind_HgCl3ml1)   = xsw * 1.d8
!        xkef(k,kc,ind_HgCl42ml1)  = 1.3d9 * cv2             ! #4174, K_eq = 1.3d1
!        xkeb(k,kc,ind_HgCl42ml1)  = xsw * 1.d8
!        xkef(k,kc,ind_HgOHBrl1)   = 2.69d14 * cv2            ! #4158, assumed
!        xkeb(k,kc,ind_HgOHBrl1)   = xsw * 1.d7
!        xkef(k,kc,ind_HgBrpl1)    = 1.1d14 * cv2             ! #4174, K_eq = 1.1d9
!        xkeb(k,kc,ind_HgBrpl1)    = xsw * 1.d5
!        xkef(k,kc,ind_HgBr2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d8
!        xkeb(k,kc,ind_HgBr2l1)    = xsw * 1.d6
!        xkef(k,kc,ind_HgBr3ml1)   = 1.5d10 * cv2             ! #4174, K_eq = 1.5d2
!        xkeb(k,kc,ind_HgBr3ml1)   = xsw * 1.d8
!        xkef(k,kc,ind_HgBr42ml1)  = 2.3d10 * cv2             ! #4174, K_eq = 2.3d1
!        xkeb(k,kc,ind_HgBr42ml1)  = xsw * 1.d9

!     else ! conv2(kc,k) = 0
!        xkef(k,kc,:) = 0. ! jjb xsw was set to 0. in this case, leading to xkef = 0.
!        xkeb(k,kc,:) = 0. ! jjb cv2 was tested equal to 0 in this case, leading to xkeb = 0.
!     end if

     enddo
  enddo

end subroutine equil_co_t


!
!------------------------------------------------------
!

subroutine equil_co_a (tt,nmaxf)

! Description :
! -----------
  ! equilibrium constant (see MOCCA)
  ! xkef: forward
  ! xkeb: backward
  ! activity coefficients included via xgamma
  ! acidity constants Ka = XXXaf/XXXab [mol/l];
  ! and other equilibrium constants (f=forward, b=backward reaction);
  ! converted to [mol/m3(air)];
  ! absolute values are chosen arbitrarily to ensure fast equilibration;

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb: added all missing declarations and implicit none

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       j6, &
       nf, &
       n, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  include 'aer_Parameters.h' !additional common blocks and other definitions

! Subroutine arguments:
  real(kind=dp), intent(in) :: tt(n)
  integer, intent(in) :: nmaxf

! Local scalars:
  integer :: k, kc
  real (kind=dp) :: cv2, xsw

! Statement function
  real (kind=dp) :: a0, b0, funa
  funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))

! Common blocks:
  common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
  real (kind=dp) :: conv2
  common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
       xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
  real(kind=dp) :: henry, xkmt, xkef, xkeb
  common /kpp_mol/ xgamma(nf,j6,nkc)
  real (kind=dp) :: xgamma

! == End of declarations =======================================================


  do k=2,nmaxf
     do kc=1,2 !nkc
        cv2 = conv2(kc,k)
        xsw = 1._dp
        if (cv2.eq.0.) xsw = 0._dp ! jjb no need to compute all the xkef & xkeb which will be 0. !
!        if (cv2.gt.0.) then     ! but still better than initialising all NSPEC list !

! absolute values of back- and forward reactions are choosen arbitrarily to
! ensure a quick (but numerically stable) equilibrium
! if illogical error messages arise, try fiddling with the absolute values
! of equilibrium reactions

        xkef(k,kc,ind_H2O) = xsw * funa(1.0d-5,-6716.d0,k) !schneller wird schlechter!
        xkeb(k,kc,ind_H2O) = 1.0D9 * cv2 * xgamma(k,1,kc) * xgamma(k,3,kc)

!        xkef(k,kc,ind_HO2) = xsw * 2.d5  !Chameides '84
        xkef(k,kc,ind_HO2) = xsw * 1.6d5  !Weinstein-loyd and Schwartz, '91
        xkeb(k,kc,ind_HO2) = 1.d10 * cv2 * xgamma(k,1,kc) * xgamma(k,11,kc)

        xkef(k,kc,ind_ACO2) = xsw * 1.8D0
        xkeb(k,kc,ind_ACO2) = 1.0D4 * cv2 * xgamma(k,1,kc) * xgamma(k,16,kc)

!        xkef(k,kc,ind_CO2) = xsw * funa(4.3d3,-913.d0,k)
!        xkeb(k,kc,ind_CO2) = 1.0D10 * cv2
        xkef(k,kc,ind_CO2) = xsw * funa(4.3d-2,-913.d0,k)
        xkeb(k,kc,ind_CO2) = 1.0D5 * cv2 * xgamma(k,1,kc) * xgamma(k,9,kc)

        xkef(k,kc,ind_HONO) = xsw * funa(5.1d+3,-1260.d0,k)
        xkeb(k,kc,ind_HONO) = 1.0D7 * cv2 * xgamma(k,1,kc) * xgamma(k,12,kc)

        xkef(k,kc,ind_HNO3) = xsw * funa(1.54d+10,8700.d0,k)
        xkeb(k,kc,ind_HNO3) = 1.0D9 * cv2 * xgamma(k,1,kc) * xgamma(k,13,kc)

        xkef(k,kc,ind_HNO4) = xsw * 2.0D3
        xkeb(k,kc,ind_HNO4) = 2.0D8 * cv2

        xkef(k,kc,ind_NH3) = xsw * funa(1.7d5,-4325.d0,k)
        xkeb(k,kc,ind_NH3) = 1.0D10 *  cv2 * xgamma(k,3,kc) * xgamma(k,2,kc)

        xkef(k,kc,ind_HSO3ml1) = xsw * funa(6.0d2,1120.d0,k) * xgamma(k,5,kc)
        xkeb(k,kc,ind_HSO3ml1) = 1.0D10 * cv2 * xgamma(k,1,kc) * xgamma(k,6,kc)

        xkef(k,kc,ind_H2SO4) = xsw * 1.0d12 !Seinfeld, Pandis (1998), p.391
        xkeb(k,kc,ind_H2SO4) = 1.0d9 * cv2 * xgamma(k,1,kc) * xgamma(k,19,kc)

        xkef(k,kc,ind_HSO4ml1) = xsw * funa(1.02d+6,2720.d0,k) * xgamma(k,19,kc)
        xkeb(k,kc,ind_HSO4ml1) = 1.0D8 * cv2 * xgamma(k,1,kc) * xgamma(k,8,kc)

        xkef(k,kc,ind_SO2) = xsw * funa(1.7d8,2090.d0,k)
        xkeb(k,kc,ind_SO2) = 1.0D10 * cv2 * xgamma(k,1,kc) * xgamma(k,5,kc)

        xkef(k,kc,ind_HCHO) = 1.d10 * cv2  !Chameides '84   !mechanism changed
        xkeb(k,kc,ind_HCHO) = xsw * 1.d5                    !mechanism changed

        xkef(k,kc,ind_HCl) = xsw * funa(1.7d10,6896.d0,k)
        xkeb(k,kc,ind_HCl) = 1.0D4 * cv2 * xgamma(k,1,kc) * xgamma(k,14,kc)

        xkef(k,kc,ind_Cl2ml1) = xsw * 5.2d4 * xgamma(k,15,kc)    !Chameides '84
        xkeb(k,kc,ind_Cl2ml1) = 1.d10 * cv2 * xgamma(k,14,kc)

!        xkef(k,kc,ind_Cl2) = 1.1D-3
!        xkeb(k,kc,ind_Cl2) = 2.1D2 * cv2

        xkef(k,kc,ind_HOCl) = xsw * 3.2D2
        xkeb(k,kc,ind_HOCl) = 1.0D10 * cv2 * xgamma(k,1,kc) * xgamma(k,22,kc)

        xkef(k,kc,ind_HBr) = xsw * 1.0D13
        xkeb(k,kc,ind_HBr) = 1.0D4 * cv2 * xgamma(k,1,kc) * xgamma(k,24,kc)

        xkef(k,kc,ind_Br2) = xsw * funa(2.95d4,-4068.d0,k) * xgamma(k,25,kc)  !Liu et al, 2002, #2109
        xkeb(k,kc,ind_Br2) = funa(1.17d10,-1812.d0,k) * cv2 * xgamma(k,24,kc)

        xkef(k,kc,ind_HOBr) = xsw * funa(2.3d1,-3091.d0,k)
        xkeb(k,kc,ind_HOBr) = 1.0D10 * cv2 * xgamma(k,1,kc) * xgamma(k,26,kc)

        xkef(k,kc,ind_BrCl2ml1) = funa(5.d9,1143.d0,k) * cv2 * xgamma(k,14,kc) !#894
        xkeb(k,kc,ind_BrCl2ml1) = xsw * 1.3D9 * xgamma(k,28,kc)
! no activities used, they "cancel out"
        xkef(k,kc,ind_Br2Clml1) = 5.d9 * cv2 !*xgamma(k,24,kc)  !if original equilibrium
        xkeb(k,kc,ind_Br2Clml1) = xsw * 2.8D5!*xgamma(k,29,kc)  !       -"-
        xkef(k,kc,ind_Br2l1) = 5.d9 * cv2  !*xgamma(k,14,kc)    !if original equilibrium
        xkeb(k,kc,ind_Br2l1) = xsw * 3.85D9!*xgamma(k,29,kc)    !       -"-
        xkef(k,kc,ind_ICl) = 1.0D11 * cv2 * xgamma(k,14,kc)
        xkeb(k,kc,ind_ICl) = xsw * 1.3D9 * xgamma(k,37,kc)
        xkef(k,kc,ind_IBr) = 1.0D11 * cv2 *xgamma(k,24,kc)
        xkeb(k,kc,ind_IBr) = xsw * 3.5D8 * xgamma(k,38,kc)
! new (speculative) ICl <--> IBr equilibria, assumed to yield the same
! ICl/IBr ratio as in the BrCl <--> Br2 eqilibria (BrCl/Br2)
! no activities used, they are not known!! 11.04.01
        xkef(k,kc,ind_IClBrml1) = 5.d9 * cv2 !*xgamma(k,24,kc)
        xkeb(k,kc,ind_IClBrml1) = xsw * 2.8D5!*xgamma(k,29,kc)
        xkef(k,kc,ind_I2) = 5.d9 * cv2  !*xgamma(k,14,kc)
        xkeb(k,kc,ind_I2) = xsw * 3.85D9!*xgamma(k,29,kc)
!        xkef(k,kc,ind_HIO2) = xsw * 2.0D3
!        xkeb(k,kc,ind_HIO2) = 2.0D9 * cv2
        xkef(k,kc,ind_HIO3) = xsw * 1.57D4
        xkeb(k,kc,ind_HIO3) = 1.0D5 * cv2

! mercury; no activity coefficients; forward not faster than 1.d14
!
!        xkef(k,kc,ind_HgOHpl1)    = 4.27d14 * cv2           ! #493,  K_eq = 4.27d10
!        xkeb(k,kc,ind_HgOHpl1)    = xsw * 1.d4
!        xkef(k,kc,ind_HgOH2l1)    = 2.6d14 * cv2            ! #4174, K_eq = 2.6d11
!        xkeb(k,kc,ind_HgOH2l1)    = xsw * 1.d3
!        xkef(k,kc,ind_HgSO3l1)    = 5.01d14 * cv2           ! #493,  K_eq = 5.01d12
!        xkeb(k,kc,ind_HgSO3l1)    = xsw * 1.d2
!        xkef(k,kc,ind_HgSO322ml1) = 2.5d14 * cv2            ! #4158, K_eq = 2.5d11
!        xkeb(k,kc,ind_HgSO322ml1) = xsw * 1.d3
!        xkef(k,kc,ind_HgOHCll1)   = 2.69d14 * cv2            ! #4158, K_eq = 2.69d7
!        xkeb(k,kc,ind_HgOHCll1)   = xsw * 1.d7
!        xkef(k,kc,ind_HgClpl1)    = 2.0d14 * cv2             ! #4174, K_eq = 2.0d7
!        xkeb(k,kc,ind_HgClpl1)    = xsw * 1.d7
!        xkef(k,kc,ind_HgCl2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d6
!        xkeb(k,kc,ind_HgCl2l1)    = xsw * 1.d8
!        xkef(k,kc,ind_HgCl3ml1)   = 6.7d8 * cv2             ! #4174, K_eq = 6.7d0
!        xkeb(k,kc,ind_HgCl3ml1)   = xsw * 1.d8
!        xkef(k,kc,ind_HgCl42ml1)  = 1.3d9 * cv2             ! #4174, K_eq = 1.3d1
!        xkeb(k,kc,ind_HgCl42ml1)  = xsw * 1.d8
!        xkef(k,kc,ind_HgOHBrl1)   = 2.69d14 * cv2            ! #4158, assumed
!        xkeb(k,kc,ind_HgOHBrl1)   = xsw * 1.d7
!        xkef(k,kc,ind_HgBrpl1)    = 1.1d14 * cv2             ! #4174, K_eq = 1.1d9
!        xkeb(k,kc,ind_HgBrpl1)    = xsw * 1.d5
!        xkef(k,kc,ind_HgBr2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d8
!        xkeb(k,kc,ind_HgBr2l1)    = xsw * 1.d6
!        xkef(k,kc,ind_HgBr3ml1)   = 1.5d10 * cv2             ! #4174, K_eq = 1.5d2
!        xkeb(k,kc,ind_HgBr3ml1)   = xsw * 1.d8
!        xkef(k,kc,ind_HgBr42ml1)  = 2.3d10 * cv2             ! #4174, K_eq = 2.3d1
!        xkeb(k,kc,ind_HgBr42ml1)  = xsw * 1.d9

!     else ! conv2(kc,k) = 0
!        xkef(k,kc,:) = 0. ! jjb xsw was set to 0. in this case, leading to xkef = 0.
!        xkeb(k,kc,:) = 0. ! jjb cv2 was tested equal to 0 in this case, leading to xkeb = 0.
!     end if

     enddo
  enddo

end subroutine equil_co_a


!
!------------------------------------------------------------------------
!

subroutine konc

! Description :
! -----------
  ! new concentrations of chemical species in liquid phases
  ! due to changes in particle size distribution
  !
  ! all changes are calculated via changes of the volume of
  ! the nkc size bins

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb : vol1 now split into vol1_a for aerosol and vol1_d for droplets, so that arrays are filled
  !       (previously vol1(2/4,1:ka) and vol1(1/3,ka+1:nka) were empty)
  ! jjb: added all missing declarations and implicit none, minor rewritting

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:
  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nf, &
       n, &
       nka, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: ia, ii, jj, k, l
  real (kind=dp) :: del, delta, dp_1, dp_2, dp_3, dp_4, xs

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck07/ part_o_a(nka,nf+1),part_o_d(nka,nf+1), &
                  part_n_a(nka,nf+1),part_n_d(nka,nf+1),pntot(nkc,nf+1)
  real (kind=dp) :: part_o_a, part_o_d, part_n_a, part_n_d, pntot
  common /blck08/ vol1_a(nka,nf+1),vol1_d(nka,nf+1),vol2(nkc,nf+1)
  real (kind=dp) :: vol1_a, vol1_d, vol2
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1

! == End of declarations =======================================================

! vol2 and vol1_a/d old liquid volume in bin (1-nkc) and row/class (1-nka)
! (um^3/cm^3), used in SR konc to shift moles from aerosol to drop
! and vice versa.
! part_o and part_n old and new part. conc. in row/class (1-nkc,1-nka) (cm^-3)

  do k=2,nf

! small (dry) particles
     do ia=1,ka
        dp_1 = part_o_a(ia,k) - part_n_a(ia,k)
        dp_3 = part_o_d(ia,k) - part_n_d(ia,k)
! Change of particle number should cancel out. If not, arise a warning message
        if (abs(dp_1+dp_3).gt.1.d-10) then
           write (jpfunout,6000) 'Warning SR konc dp_1 > dp_3',k,ia,dp_1,dp_3
        end if
! search bin that loses moles:
        if (dp_1.ge.1.d-10) then
           ii = 1
        else ! (dp_3.ge.1.d-10)
           ii = 3
        end if
! no change if number of growing particles is too small:
        if (abs(dp_1).lt.1.d-10) then
           xs = 0._dp
        else
           xs = 1._dp
        end if
! mol/m^3(air) -1-> mol/l(aq) -2-> mol/m^3(air,row) -3-> mol/part(row)
! -4-> mol/m^3(air,change,row)
! mol/m^3(air) -1-> mol/l(aq):
!     mol/m^3(air)/(vol2*1.d-12) = mol/m^3(aq)
!     mol/m^3(aq)*1.d-3 --> mol/l(aq)
! mol/l(aq) -2-> mol/m^3(air,row):
!     mol/l(aq) *vol1_@(ia,k)*1.d-9 = mol/m^3(air,row)
! mol/m^3(air,row) -3-> mol/part(row):
!     mol/m^3(air,row)*1.d-6 cm^3(air,row)/part = mol/part(row)
! mol/part(row) -4-> mol/m^3(air,change,row):
!     mol/part(row)*1.d6*part/cm^3(air,row,change) = mol/m^3(air,row,change)
        if (ii.eq.1) then
           jj=3
           if (vol2(ii,k).gt.0..and.part_o_a(ia,k).gt.0.) then
              delta = vol1_a(ia,k)/vol2(ii,k)*dp_1/part_o_a(ia,k)*xs
           else
              delta = 0._dp
           end if
        else ! ii == 3
           jj=1
           if (vol2(ii,k).gt.0..and.part_o_d(ia,k).gt.0.) then
              delta = vol1_d(ia,k)/vol2(ii,k)*dp_3/part_o_d(ia,k)*xs
           else
              delta = 0._dp
           end if
        endif
! Check delta is within 0 and 1, and change concentrations only if >0
        if (delta.lt.0._dp) then
           write (jpfunout,6001) k,ia,delta,'Warning SR konc s <'
        else if (delta.gt.1._dp) then
           write (jpfunout,6001) k,ia,delta,'Warning SR konc s >'
        else if (delta.gt.0._dp) then
           do l=1,j2
              del = sl1(l,ii,k) * delta
              sl1(l,ii,k) = max(0._dp, sl1(l,ii,k)-del)
              sl1(l,jj,k) = max(0._dp, sl1(l,jj,k)+del)
           enddo
           do l=1,j6
              del = sion1(l,ii,k) * delta
              sion1(l,ii,k) = max(0._dp, sion1(l,ii,k)-del)
              sion1(l,jj,k) = max(0._dp, sion1(l,jj,k)+del)
           enddo
        !else      ! the only remaining case here is delta == 0
        !   cycle  ! leave concentrations unchanged in this case
        end if
     enddo

! large (dry) particles (see above for comments)
     do ia=ka+1,nka
        dp_2 = part_o_a(ia,k) - part_n_a(ia,k)
        dp_4 = part_o_d(ia,k) - part_n_d(ia,k)
        if (abs(dp_2+dp_4).gt.1.d-10) then
           write (jpfunout,6000) 'Warning SR konc dp_2 > dp_4',k,ia,dp_2,dp_4
        end if
        if (dp_2.ge.1.d-10) then
           ii=2
        else ! (dp_4.ge.1.d-10)
           ii=4
        end if
        if (abs(dp_2).lt.1.d-10) then
           xs = 0._dp
        else
           xs = 1._dp
        end if
        if (ii.eq.2) then
           jj=4
           if (vol2(ii,k).gt.0..and.part_o_a(ia,k).gt.0.) then
              delta = vol1_a(ia,k)/vol2(ii,k)*dp_2/part_o_a(ia,k)*xs
           else
              delta = 0._dp
           end if
        else ! ii == 4
           jj=2
           if (vol2(ii,k).gt.0..and.part_o_d(ia,k).gt.0.) then
              delta = vol1_d(ia,k)/vol2(ii,k)*dp_4/part_o_d(ia,k)*xs
           else
              delta = 0._dp
           end if
        endif

        if (delta.lt.0._dp) then
           write (jpfunout,6001) k,ia,delta,'SR konc l <'
        else if (delta.gt.1._dp) then
           write (jpfunout,6001) k,ia,delta,'SR konc l >'
        else if (delta.gt.0._dp) then
           do l=1,j2
              del = sl1(l,ii,k)*delta
              sl1(l,ii,k) = max(0._dp,sl1(l,ii,k)-del)
              sl1(l,jj,k) = max(0._dp,sl1(l,jj,k)+del)
           enddo
           do l=1,j6
              del = sion1(l,ii,k)*delta
              sion1(l,ii,k) = max(0._dp,sion1(l,ii,k)-del)
              sion1(l,jj,k) = max(0._dp,sion1(l,jj,k)+del)
           enddo
        end if
     enddo

! Transfer all remaining species from droplet bins (#3 or #4) to the corresponding aerosol bin
! (#1 and #2, respectively) if too little particle remains in the droplet bins     
     do l=1,j2
        if (pntot(3,k).lt.1.d-7) then
           sl1(l,1,k) = sl1(l,1,k) + max(0._dp, sl1(l,3,k))
           sl1(l,3,k) = 0._dp
        endif
        if (pntot(4,k).lt.1.d-7) then
           sl1(l,2,k) = sl1(l,2,k) + max(0._dp, sl1(l,4,k))
           sl1(l,4,k) = 0._dp
        endif
     enddo

     do l=1,j6
        if (pntot(3,k).lt.1.d-7) then
           sion1(l,1,k) = sion1(l,1,k) + max(0._dp, sion1(l,3,k))
           sion1(l,3,k) = 0._dp
        endif
        if (pntot(4,k).lt.1.d-7) then
           sion1(l,2,k) = sion1(l,2,k) + max(0._dp, sion1(l,4,k))
           sion1(l,4,k) = 0._dp
        endif
     enddo

  enddo ! k

6000 format (a,2i3,2es10.3)
6001 format (2i3,es10.3,a)

end subroutine konc

!
!-----------------------------------------------------------------------------------
!

subroutine init_konc

! Description :
! -----------
  ! Init loading of aerosols with ions.
  ! If rH>0.7 aerosols get activated and liquid chemistry is started.
  ! The ions are transported vertically by turbulence and sedimentation
  ! independent of activation.

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! jjb work done: turned ap into one dimension only (previously ap(n,nka))
  !                reindexed sa1 for computing efficiency
  !                introduced lpJoyce14bc special case
  !                cleaned, missing declarations and implicit none

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:
  USE config, ONLY : &
       ifeed, &
       lpJoyce14bc

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j3, &
       j6, &
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: ia, ial, jt, k
  real (kind=dp) :: ap(nka)

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /blck78/ sa1(j2,nka)
  real(kind=dp) :: sa1
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

! == End of declarations =======================================================

! sa1 in mole/particle; sl1, sion1 in mole m**-3
! ap total number of aerosols in cm**-3
! of sulfate and sea salt aerosol regime

  do k=2,n-1

     if (ifeed.eq.2) then
        ial = 2
     else
        ial = 1
     endif
     do ia=ial,nka
        ap(ia) = 0._dp
        do jt=1,nkt
           ap(ia) = ap(ia)+ff(jt,ia,k)
        enddo
        if (ap(ia).eq.0._dp) write(jpfunout,6002) ia,k,'init(ap)=0'
     enddo
! no3-, nh4+ and so4= due to nucleation scavenging and evaporation
     do ia=1,ka
        if (lpJoyce14bc) &
             sion1( 1,1,k)=sion1( 1,1,k)+ap(ia)*sa1( 1,ia)*1.d6 !H+ PJ
        sion1( 2,1,k)=sion1( 2,1,k)+ap(ia)*sa1( 2,ia)*1.d6
        sion1( 8,1,k)=sion1( 8,1,k)+ap(ia)*sa1( 8,ia)*1.d6
        sion1(13,1,k)=sion1(13,1,k)+ap(ia)*sa1(13,ia)*1.d6
        if (lpJoyce14bc) &
             sion1(14,1,k)=sion1(14,1,k)+ap(ia)*sa1(14,ia)*1.d6 !Cl- (cl_b23c) PJ
        sion1(17,1,k)=sion1(19,1,k)                          !mixing tracer
        sion1(19,1,k)=sion1(19,1,k)+ap(ia)*sa1(19,ia)*1.d6
        if (lpJoyce14bc) &
             sl1(j2-j3+4,1,k)=sl1(j2-j3+4,1,k)+ap(ia)*sa1(j2-j3+4,ia)*1.d6  !DOM    PJ
! NO3-, NH4-, SO4=, HSO4- : j6 used in sion1, sa1
! Br-, HCO3-, I-, IO3-, Cl-: j6 used in sion1, sa1
     enddo
     if (sion1(8,1,k).eq.0.) write (jpfunout,6003) k,'init(so4=)=0'

     do ia=ka+1,nka
        if (lpJoyce14bc) then
           sion1( 1,2,k)=sion1( 1,2,k)+ap(ia)*sa1( 1,ia)*1.d6 !H+ PJ
           sion1( 2,2,k)=sion1( 2,2,k)+ap(ia)*sa1( 2,ia)*1.d6 !NH4+ PJ
        end if
        sion1( 8,2,k)=sion1( 8,2,k)+ap(ia)*sa1( 8,ia)*1.d6 !SO4=
        sion1( 9,2,k)=sion1( 9,2,k)+ap(ia)*sa1( 9,ia)*1.d6 !HCO3-
        sion1(13,2,k)=sion1(13,2,k)+ap(ia)*sa1(13,ia)*1.d6 !NO3-
        sion1(14,2,k)=sion1(14,2,k)+ap(ia)*sa1(14,ia)*1.d6 !Cl-
        sion1(17,2,k)=sion1(14,2,k)                          !mixing tracer
        if (lpJoyce14bc) then
           sion1(19,2,k)=sion1(19,2,k)+ap(ia)*sa1(19,ia)*1.d6 !HSO4- PJ
        end if
        sion1(20,2,k)=sion1(20,2,k)+ap(ia)*sa1(20,ia)*1.d6 !Na+; inert
        sion1(24,2,k)=sion1(24,2,k)+ap(ia)*sa1(24,ia)*1.d6 !Br-
        sion1(34,2,k)=sion1(34,2,k)+ap(ia)*sa1(34,ia)*1.d6 !I-
        sion1(36,2,k)=sion1(36,2,k)+ap(ia)*sa1(36,ia)*1.d6 !IO3-
        sl1(j2-j3+4,2,k)=sl1(j2-j3+4,2,k)+ap(ia)*sa1(j2-j3+4,ia)*1.d6 !DOM
     enddo
     if (sion1(14,2,k).eq.0.) write (jpfunout,6003) k,'init(cl-)=0'
  enddo

6002 format (2i3,a)
6003 format (i3,a)

end subroutine init_konc


!
!-----------------------------------------------------------------------------------
!

subroutine aer_source_init

! Description :
! -----------
  ! find indexes and weights of layers around 10m heigh for accurate u10 calculation

! Method :
! -----------
  ! find k10m, k10p such that eta(k10m) <= 10 < eta(k10p)
  ! linear interpolation to compute the weights

! Author :
! ------
  ! Josué Bock

! Modifications :
! -------------
  ! 05-Jun-2021  Josué Bock  First version

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:
  USE config, ONLY : &
! Imported Parameters:
       detamin, &
       lpmona, lpsmith

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  integer :: k
  real (kind=dp) :: zza, zzb, zzr
  common /blck09/ k10m, k10p, w10m, w10p
  integer :: k10m, k10p
  real (kind=dp) :: w10m, w10p
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

! == End of declarations =======================================================

  k=1
  do while (eta(k+1).lt.10._dp)
     k = k+1
  end do
  k10m = k
  k10p = k+1

  ! linear interpolation:
  ! a*k     + b = eta(k)
  ! a*(k+1) + b = eta(k+1) = eta(k) + a
  !  thus a=detamin and b=eta(k)-k*detamin
  ! a*r     + b = 10 with k <= r < k+1
  zza = detamin
  zzb = eta(k)-k*detamin
  zzr = (10._dp - zzb) / zza

  w10p = zzr - real(k,dp)
  w10m = 1._dp - w10p

  write (jpfunout,6000)
  write (jpfunout,6010) k10m,w10m
  write (jpfunout,6011) k10p,w10p
  write (jpfunout,6012) lpmona, lpsmith
6000 format ('---------------------------------------------------------------------')
6010 format ('u10 will be computed from wind at index',i3,' with weight',f5.2)
6011 format ('and index',i3,' with weight',f5.2)
6012 format ('emission parameterisation: Monahan',l2,' Smith',l2)

end subroutine aer_source_init


!
!-----------------------------------------------------------------------------------
!

subroutine aer_source (box,dd,z_mbl,n_bl)

! Description :
! -----------
  ! calculation of sea salt aerosol source

! Author :
! ------
  ! Roland von Glasow

! Modifications :
! -------------
  ! work done:
  ! 28-Apr-2021  Josué Bock  Propagate rho3 from constants module, instead of hard coded values
  ! 04-Jun-2021  Josué Bock  

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:
  USE config, ONLY : &
! Imported Parameters:
       lpJoyce14bc,  &          ! special emission scheme: Joyce et al, ACP 2014 (base case)
       lpmona,       &          ! Monahan et al., 1986
       lpsmith                  ! Smith et al, 1993

  USE constants, ONLY : &
! Imported Parameters:
       pi, &
       rho3, &                  ! aerosol density
       rhow                     ! water density

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j3, &
       j6, &
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  logical, intent(in) ::  box
  integer, intent(in) :: n_bl
  real (kind=dp), intent(in) :: dd, z_mbl

! Local parameters:
  ! optimisation: define parameters that will be computed only once
  real (kind=dp), parameter :: zrho_frac = rho3 / rhow
  real (kind=dp), parameter :: z4pi3 = 4.e-09_dp * pi / 3._dp

! Local scalars:
  integer :: ia, iamin, iamax, ikc, jt, jt_low, jtt, k_in
  real (kind=dp) :: a1, a2, f1, f2, r01, r02
  real (kind=dp) :: bb, df, df1, df2
  real (kind=dp) :: a0, b0, rg, eg, rr
  real (kind=dp) :: d_z, u10

  real (kind=dp), external :: rgl

! Common blocks
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck09/ k10m, k10p, w10m, w10p
  integer :: k10m, k10p
  real (kind=dp) :: w10m, w10p
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real(kind=dp) :: sl1, sion1
  common /blck78/ sa1(j2,nka)
  real(kind=dp) :: sa1
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb44/ a0m,b0m(nka)
  real (kind=dp) :: a0m,b0m
  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw, ew, rn, rw, en, e, dew, rq
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
  common /ff_0/ ff_0(nka)
  real(kind=dp) :: ff_0
  common /sss/ brsss,clsss,xnasss
  real (kind=dp) :: brsss, clsss, xnasss

! == End of declarations =======================================================

! Set the index and the height of the layer where new particles are injected
  if (box) then
     k_in = n_bl
     d_z  = z_mbl
  else
     k_in = 2
     d_z  = deta(2)
  endif
! this is also defined for box:
  u10 = w10m*sqrt(u(k10m)**2+v(k10m)**2) + w10p*sqrt(u(k10p)**2+v(k10p)**2)

  if (lpsmith) then
     a1  = 10**(0.0676_dp * u10 + 2.43_dp)
     a2  = 10**(0.959_dp * u10**0.5 - 1.476_dp)
     f1  = 3.1_dp
     f2  = 3.3_dp
     r01 = 2.1_dp
     r02 = 9.2_dp
  endif

! Define radius bounds for aerosol source
  if (lpJoyce14bc) then
     ! sfc source, "quick and dirty"   PJ
     iamin = 1
     iamax = nka
  else
     ! general case: only large aerosols are emitted
     iamin = ka+1
     iamax = nka
  end if

  do ia=iamin, iamax

! compare SR equil in str.f
! aerosols in equilibrium with ambient rH:
     feu(k_in) = min(feu(k_in), 0.99999_dp)
! equilibrium radius
     a0 = a0m / t(k_in)
     b0 = b0m(ia) * zrho_frac
! b0=b0m*rho3/rhow; rho3=2000; rhow=1000
     rg = rgl(rn(ia),a0,b0,feu(k_in))
     eg = z4pi3 * (rg**3 - rn(ia)**3)

     jtclass: do jt=1,nkt
! find jt that corresponds to particle size at ambient rH
!        if (rq(jt,ia).ge.rg) then
        if (eg.le.ew(jt)) then
!           rr=rq(jt,ia)
! source functions are for rH=80%
           rr = rgl(rn(ia),a0,b0,0.8_dp)
! find jt-index that corresponds to rr for this dry aerosol radius
           jt_low = -1
           do jtt=1,nkt
              if (rq(jtt,ia).le.rr) jt_low = jtt
           enddo
!           print *,ia,jt_low,rr,rq(jt_low,ia)

! =============              =============
!               GENERAL CASE
! =============              =============
           if (.not.lpJoyce14bc) then
              if (lpmona) then
! aerosol source after Monahan et al. 86 cited in Gong et al, 97, JGR, 102, p3805-3818
!  "log" is log10
!                 bb=(0.380_dp - log(rr)) / 0.65
                 bb=(0.380_dp - log10(rr)) / 0.65
!                 print *,'bb',bb
                 df1=1.373*u10**3.41*rr**(-3)*(1._dp + 0.057 *rr**1.05) * 10**(1.19*exp(-bb**2))
!                 if (rr.gt.10..and.rr.lt.75.) df22=8.6d-6*exp(2.08*u10)*rr**(-2)
!                 if (rr.gt.75..and.rr.lt.100.) df23=4.83d-2*exp(2.08*u10)*rr**(-4)
!                 if (rr.gt.100.) df24=8.6d6*exp(2.08*u10)*rr**(-8)
                 df=df1      !+df22+df23+df24
              endif
              if (lpsmith) then
! aerosol source after Smith et al., 93, QJRMS,119,809-824
                 df1=a1*exp(-f1*(log(rr/r01)**2))
                 df2=a2*exp(-f2*(log(rr/r02)**2))
                 df=df1+df2
              endif
!              print *,df,rr,rq(jt,ia)-rq(jt,ia-1)
! df in m^-2 um^-1 s^-1, convert to: cm^-3 s^-1
! m^-2 --> m^-3: 1/dz, m^-3 --> cm^-3: 10-6, integrate over r
!              if (jt.eq.1) then
              if (jt_low.eq.1) then
!                 df=df*(rq(jt+1,ia)-rq(jt,ia))/d_z*1.d-6
                 df=df*(rq(jt_low+1,ia)-rq(jt_low,ia))/d_z*1.d-6
              else
!                 df=df*(rw(jt,ia)-rw(jt-1,ia))/d_z*1.d-6
                 df=df*(rw(jt_low,ia)-rw(jt_low-1,ia))/d_z*1.d-6
              endif
!              print *,ia,df,df*86400
! =============              =============
!               SPECIAL CASE
! =============              =============
           else if (lpJoyce14bc) then
              ! "aeroPJ"  SR aer_source, 2 hour source, 2 hr delay
              if(lday.eq.0.AND.lst.ge.2.AND.lst.le.3)then
                 df=ff_0(ia)/86400. ! magnitude: replenished in one day
                 df=df*15.03        ! modifer to match DEC_FBX constraint (twk3)
!                 df=df*35.26        ! modifer to match STD_URB constraint
              else
                 df=0._dp
              endif
           end if

! Select the correct kc bin
           if (.not.lpJoyce14bc) then
              ! general case: iamin = ka+1, only large aerosol (kc=2)
              ikc = 2
           else if (lpJoyce14bc) then
              if (ia.le.ka) then
                 ikc = 1
              else
                 ikc = 2
              end if
           end if


! new distribution
! change number conc and aerosol conc (Br-, I-, IO3-, Cl-, Na+, DOM)
           ff(jt,ia,k_in)=ff(jt,ia,k_in)+df*dd
           if (lpJoyce14bc) then
              sion1( 1,ikc,k_in)=sion1( 1,ikc,k_in)+df*dd*sa1( 1,ia)*1.d6 !H+
           end if
           sion1( 8,ikc,k_in)=sion1( 8,ikc,k_in)+df*dd*sa1( 8,ia)*1.d6 !SO4=
           if (.not.lpJoyce14bc) then
              sion1( 9,ikc,k_in)=sion1( 9,ikc,k_in)+df*dd*sa1( 9,ia)*1.d6 !HCO3-
              sion1(13,ikc,k_in)=sion1(13,ikc,k_in)+df*dd*sa1(13,ia)*1.d6 !NO3-
           end if
           sion1(14,ikc,k_in)=sion1(14,ikc,k_in)+df*dd*sa1(14,ia)*1.d6 !Cl-
           if (.not.lpJoyce14bc) then
              sion1(20,ikc,k_in)=sion1(20,ikc,k_in)+df*dd*sa1(20,ia)*1.d6 !Na+, ion balance
              sion1(24,ikc,k_in)=sion1(24,ikc,k_in)+df*dd*sa1(24,ia)*1.d6 !Br-
              sion1(34,ikc,k_in)=sion1(34,ikc,k_in)+df*dd*sa1(34,ia)*1.d6 !I-
              sion1(36,ikc,k_in)=sion1(36,ikc,k_in)+df*dd*sa1(36,ia)*1.d6 !IO3-
           end if
           sl1(j2-j3+4,ikc,k_in)=sl1(j2-j3+4,ikc,k_in)+df*dd*sa1(j2-j3+4,ia)*1.d6 !DOM
           brsss=brsss+df*dd*sa1(24,ia)*1.d6
           clsss=clsss+df*dd*sa1(14,ia)*1.d6
           xnasss=xnasss+df*dd*sa1(20,ia)*1.d6
!           print *,jt,f(2,ia,jt),df
           ! once the appropriate jt has been found and refilled, exit jt loop
           exit jtclass
        endif                 ! eg.le.ew(jt)
     enddo jtclass            ! jt

  enddo ! ia

end subroutine aer_source


!
!-----------------------------------------------------------------------
!

      subroutine plo_ppH
! print NH3,NH4+,HCl,Cl-,HNO3,NO3-,SO2,H2SO4,HSO4-,SO4=
! to calculate the potential pH
! use for parameterisation development for global models

      USE gas_common, ONLY : &
! Imported Array Variables with intent (in):
           s1

      USE global_params, ONLY : &
! Imported Parameters:
           j2, &
           j6, &
           nf, &
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

! Local scalars:
      character (len=8),parameter :: fname = 'ppHa.out'
      integer :: k

! Local arrays:
      real (kind=dp) :: i0(n,25)

! Common blocks:
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      real (kind=dp) :: sl1, sion1

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct
!- End of header ---------------------------------------------------------------

      do k=2,nf
            i0(k,1)=s1(4,k)
            i0(k,2)=s1(30,k)
            i0(k,3)=s1(3,k)
            i0(k,4)=s1(5,k)
            i0(k,5)=s1(6,k)

            i0(k,6)=sl1(4,1,k)
            i0(k,7)=sl1(30,1,k)
            i0(k,8)=sl1(3,1,k)
            i0(k,9)=sl1(5,1,k)
            i0(k,10)=sl1(6,1,k)
            i0(k,11)=sl1(4,2,k)
            i0(k,12)=sl1(30,2,k)
            i0(k,13)=sl1(3,2,k)
            i0(k,14)=sl1(5,2,k)
            i0(k,15)=sl1(6,2,k)

            i0(k,16)=sion1(2,1,k)
            i0(k,17)=sion1(14,1,k)
            i0(k,18)=sion1(13,1,k)
            i0(k,19)=sion1(19,1,k)
            i0(k,20)=sion1(8,1,k)
            i0(k,21)=sion1(2,2,k)
            i0(k,22)=sion1(14,2,k)
            i0(k,23)=sion1(13,2,k)
            i0(k,24)=sion1(19,2,k)
            i0(k,25)=sion1(8,2,k)
      enddo

 3000 continue
      open (73, file=fname,status='unknown',form='unformatted', &
       position='append',err=3000)
      write (73) lday,lst,lmin,i0
      close (73)

      end subroutine plo_ppH


!
!-----------------------------------------------------------
!

      subroutine kpp_driver (box,dd_ch,n_bl)
! interface between MISTRA and the KPP gas phase chemistry

      USE config, ONLY : &
           halo, &
           iod, &
           neula

      USE constants, ONLY : &
! Imported Parameters:
       pi

      USE gas_common, ONLY : &
! Imported Array Variables with intent (inout):
           s1, &
           s3, &
           ind_gas_rev, nadvmax, nindadv, xadv ! for Eulerian advection


      USE global_params, ONLY : &
! Imported Parameters:
           j2, &
           j6, &
           nf, &
           n, &
           nkc, &
           nphrxn

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit double precision (a-h,o-z)
!     logical halo,iod,box,new_a,new_c,cloud,short,ros3 ! jjb ros3 removed, thus new_a, new_c & short as well
      logical box,cloud

      common /blck01/ am3(n),cm3(n)

      common /blck12/ cw(nkc,n),cm(nkc,n)
      real(kind=dp) :: cw, cm
      common /blck13/ conv2(nkc,n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      real (kind=dp) :: sl1, sion1
      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
      real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
      common /band_rat/ photol_j(nphrxn,n)
      common /kinv_i/ kinv
      integer :: kinv
      common /cb_1/ air_cc,te,h2oppm,pk
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_l1/ cloud(nkc,n)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated
!     common /ph_r/ ph_rat(47) ! jjb ! jjb test include in _a_g_t common blocks
      dimension ph_rat(nphrxn)

      airmolec=6.022d+20/18.0
! calculate u0 for prelim photolysis
      rlat=alat*1.745329d-02
      rdec=declin*1.745329d-02
      zeit=lst*3600.+dfloat(lmin-1)*60.
!      zeit=zeit+dtg  !quite exact, but photolysis rates are calculated in
! greater intervals and variable dtg is known only in the old chemical module
!      horang=7.272205d-05*zeit-3.1415927
      horang=7.272205d-05*zeit-pi
      u00=dcos(rdec)*dcos(rlat)*dcos(horang)+dsin(rdec)*dsin(rlat)
      ru0=6371.*u00
      u0=8./(dsqrt(ru0**2+102000.)-ru0)

      call mass_ch
! loop over all vertical layers
!      n_max=n
      n_max=n-1
      n_min=2
      if (box) then
         n_max=n_bl
         n_min=n_bl
      endif

! eliminate negative values
      where (s1 < 0.d0) s1=0.d0 ! jjb new using where construct for array s1
      where (s3 < 0.d0) s3=0.d0 ! jjb new using where construct for array s3
      !where (sl1  < 0.d0)   sl1=0.d0
      !where (sion1 < 0.d0) sion1=0.d0 ! at the moment, done in each mechanism

      do k=n_min,n_max !c aer#
!         if (k.eq.2) short=.false. !c aer# ! jjb ros2 only
! cloud#      do k=n_max,2,-1 !c cloud#
! cloud#         if (k.eq.n_max) short=.false. !c cloud#

! define temp, H2O, air, ..
! air, h2o in mlc/cm^3 are used for 3rd order reactions that are
! formulated as quasi-2nd order
         te=t(k)
         air_cc=cm3(k)                      ![air] in mlc/cm^3
         air=am3(k)                         ![air] in mol/m^3
!        co=am3(k,2)                        ![co]  in mol/m^3 ! jjb removed, unused
         h2o=xm1(k)*rho(k)/1.8d-2           ![h2o] in mol/m^3
         h2o_cc=xm1(k)*airmolec*rho(k)      ![h2o] in mlc/cm^3
         h2oppm=h2o_cc*1.d6/air_cc        ![h2o] in ppm=umol/mol
         pk=p(k)
         dt_ch=dd_ch
         !tkpp=time
         tkpp=1.d0
         !tkpp=0.d0 ! jjb change here. kpp has not to know the real time, and having large TIN would reduce minimum timestep allowed in KPP
!         if (k.lt.nf) then ! conv2 is dimension n now, and initialised to 0
            cvv1=conv2(1,k)
            cvv2=conv2(2,k)
            cvv3=conv2(3,k)
            cvv4=conv2(4,k)
!         else ! jjb check that there is no pb with index nf
!            cvv1=0.
!            cvv2=0.
!            cvv3=0.
!            cvv4=0.
!         endif

! photolysis rates
         if (u0.ge.3.48d-2) then
            do i=1,nphrxn
               ph_rat(i)=( photol_j(i,k)+photol_j(i,k-1) ) / 2 ! jjb photol_j defined for levels, not for layers
            enddo
! end of prelim j rates
         else
            do i=1,nphrxn
               ph_rat(i)=0.
            enddo
         end if

! set halogen rates to zero (if wanted)
         xhal=1.
         xiod=1.
         if (.not.halo) then
            xhal=0.
            xiod=0.
         endif
         if (.not.iod) xiod=0.
! set liquid rates to zero
         xliq1=1.
         xliq2=1.
         xliq3=1.
         xliq4=1.
! cm is switch now (not cw)
! avoid index out of bounds in cm:
         if (k.ge.nf) then
            xliq1=0.
            xliq2=0.
            xliq3=0.
            xliq4=0.
         else
            if (cm(1,k).eq.0.) xliq1=0.
            if (cm(2,k).eq.0.) xliq2=0.
            if (cm(3,k).eq.0.) xliq3=0.
            if (cm(4,k).eq.0.) xliq4=0.
         endif

         if (xliq1.eq.1..and..not.cloud(1,k)) then
            cloud(1,k)=.true.
            print *,'new aerosol layer,l1 : ',k
         endif
         if (xliq2.eq.1..and..not.cloud(2,k)) then
            cloud(2,k)=.true.
            print *,'new aerosol layer,l2 : ',k
         endif
         if (xliq3.eq.1..and..not.cloud(3,k)) then
            cloud(3,k)=.true.
            print *,'new cloud layer,l3 : ',k
         endif
         if (xliq4.eq.1..and..not.cloud(4,k)) then
            cloud(4,k)=.true.
            print *,'new cloud layer,l4 : ',k
         endif
         if (xliq1.eq.0..and.cloud(1,k))  cloud(1,k)=.false.
         if (xliq2.eq.0..and.cloud(2,k))  cloud(2,k)=.false.
         if (xliq3.eq.0..and.cloud(3,k))  cloud(3,k)=.false.
         if (xliq4.eq.0..and.cloud(4,k))  cloud(4,k)=.false.

! set heterogeneous rates to zero
         xhet1=1.
         xhet2=1.
         if (xliq1.eq.1.) xhet1=0.
         if (xliq2.eq.1.) xhet2=0.

! advection if eulerian view (neula=0), xadv in mol/(mol*day)

         if (neula.eq.0) then
            if (k.le.kinv) then
               do j=1,nadvmax
                  if (nindadv(j).ne.0) &
         s1(ind_gas_rev(nindadv(j)),k)=s1(ind_gas_rev(nindadv(j)),k)+xadv(j)*dt_ch*air/86400.
               enddo
            endif
         endif

!         print('(i4,6f4.1)'),k,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2
! call kpp chemical integrator
            if (xliq1.eq.1..or.xliq2.eq.1.) then
               if (xliq3.eq.1..or.xliq4.eq.1.) then
! gas + aerosol + cloud
!                  print *,k," tot ",tkpp,dt_ch
!                 call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4 &
!                      ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,air,co,h2o) ! jjb
!                 call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4 &
!          ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2,air,co,h2o) ! jjb working version
!                  call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4 &
!           ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2,air,co,h2o &
!           ,ph_rat)                                                   ! jjb test version to handle ph_rat another way
                  call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4 &
           ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2,air,h2o &
           ,ph_rat)                                                   ! jjb co removed
!         print *,'layer nb',k,'aqueous chem? TOT mech'
!         print *,'xliq1-4=',xliq1,xliq2,xliq3,xliq4
!         print *,'xhet1-2=',xhet1,xhet2
!         print *,'cloud(1-4,k)=',cloud(1,k),cloud(2,k),cloud(3,k) &
!        ,cloud(4,k)

               else
! gas + aerosol
!                  print *,k," aer ",tkpp,dt_ch
!                 call aer_drive (tkpp,dt_ch,k,cvv1,cvv2,xhal, &
!                      xiod,xliq1,xliq2,xhet1,xhet2,air,co,h2o) ! jjb working version
!                  call aer_drive (tkpp,dt_ch,k,cvv1,cvv2,xhal, &
!                       xiod,xliq1,xliq2,xhet1,xhet2,air,co,h2o &
!           ,ph_rat)                                             ! jjb test version to handle ph_rat another way
                  call aer_drive (tkpp,dt_ch,k,cvv1,cvv2,xhal, &
                       xiod,xliq1,xliq2,xhet1,xhet2,air,h2o &
           ,ph_rat)                                             ! jjb co removed
!         print *,'layer nb',k,'aqueous chem? AER mech'
!         print *,'xliq1-4=',xliq1,xliq2,xliq3,xliq4
!         print *,'xhet1-2=',xhet1,xhet2
!         print *,'cloud(1-4,k)=',cloud(1,k),cloud(2,k),cloud(3,k) &
!        ,cloud(4,k)
               endif
            else
! gas only
! het reactions on dry aerosol always on!
               xhet1 = 1.
               xhet2 = 1.
!               print *,k," gas ",tkpp,dt_ch
!              call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2,conv1, & ! jjb working version
!                   air,co,h2o)
!               call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2,conv1, & ! jjb test version to handle ph_rat another way
!                    air,co,h2o,ph_rat)
!              call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2,conv1, & ! jjb co removed
!                   air,h2o,ph_rat)
!               write(1717,*)'Gas layer',time,k
               call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2, &      ! jjb conv1 removed
                    air,h2o,ph_rat)
!         print *,'layer nb',k,'aqueous chem? GAS mech'
!         print *,'xliq1-4=',xliq1,xliq2,xliq3,xliq4
!         print *,'xhet1-2=',xhet1,xhet2
!         print *,'cloud(1-4,k)=',cloud(1,k),cloud(2,k),cloud(3,k) &
!       ,cloud(4,k)
            endif

      enddo ! k

! eliminate negative values
      where (s1 < 0.d0) s1=0.d0 ! jjb new using where construct for array s1
      where (s3 < 0.d0) s3=0.d0 ! jjb new using where construct for array s3

      sl1(:,:,:) = max(0.d0,sl1(:,:,:))
      sion1(:,:,:) = max(0.d0,sion1(:,:,:))

      if (lmin.eq.10) call ionbalance (box,n_bl)

      end subroutine kpp_driver

!
!-----------------------------------------------------
!

      subroutine ionbalance (box,n_bl)
! check ion balance

      USE global_params, ONLY : &
! Imported Parameters:
           j2, &
           j6, &
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      logical, intent(in) :: box
      integer, intent(in) :: n_bl

      integer :: k, kc
      integer :: n_min, n_max

      real (kind=dp) :: xpos(nkc),xneg(nkc)

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      real (kind=dp) :: sl1, sion1


      write (103,*) lday,lst,lmin,' aerosol'
      write (104,*) lday,lst,lmin,' droplet'

      n_min=1
      n_max=n
      if (box) then
         n_min=n_bl
         n_max=n_bl
      endif

      do k=n_min,n_max
         do kc=1,nkc
            xpos(kc)=sion1(1,kc,k)+sion1(2,kc,k)+sion1(20,kc,k)
            xneg(kc)=sion1(3,kc,k)+sion1(4,kc,k)+sion1(5,kc,k)+ &
                 2*sion1(6,kc,k)+sion1(7,kc,k)+2*sion1(8,kc,k) &
                 +sion1(9,kc,k)+sion1(10,kc,k)+sion1(11,kc,k) &
                 +sion1(12,kc,k)+sion1(13,kc,k)+sion1(14,kc,k) &
                 +sion1(15,kc,k)+sion1(16,kc,k)+sion1(19,kc,k) &
                 +sion1(21,kc,k) &
                 +sion1(22,kc,k)+sion1(23,kc,k)+sion1(24,kc,k) &
                 +sion1(25,kc,k)+sion1(26,kc,k)+sion1(27,kc,k) &
                 +sion1(28,kc,k)+sion1(29,kc,k)+sion1(30,kc,k) &
                 +sion1(31,kc,k)+sion1(32,kc,k)+sion1(33,kc,k) &
                 +sion1(34,kc,k)+sion1(35,kc,k)+sion1(36,kc,k) &
                 +sion1(37,kc,k)+sion1(38,kc,k)+sion1(39,kc,k)
         enddo
         write (103,101) k,xpos(1),xneg(1),xpos(1)-xneg(1), &
                         xpos(2),xneg(2),xpos(2)-xneg(2)
         write (104,101) k,xpos(3),xneg(3),xpos(3)-xneg(3), &
                         xpos(4),xneg(4),xpos(4)-xneg(4)
      enddo
 101  format (i3,6d16.8)

      end subroutine ionbalance
!
!------------------------------------------------------------
!

      subroutine dry_cw_rc (nmax)
! calculates LWC and mean radius for "dry" aerosol

      USE config, ONLY : &
           ifeed

      USE constants, ONLY : &
! Imported Parameters:
       pi

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           n, &
           nka, &
           nkt, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: nmax

! Local parameters:
      real (kind=dp), parameter :: xpi = 4./3.*pi

! Local scalars:
      real (kind=dp) :: cwd1, cwd2
      real (kind=dp) :: rcd1, rcd2
      real (kind=dp) :: x0
      integer :: ia, ial, jt, k

! Common blocks:
      common /blck06/ kw(nka),ka
      integer :: kw, ka

      common /blck11/ rcd(nkc,n)
      real (kind=dp) :: rcd

      common /blck12/ cwd(nkc,n),cm(nkc,n)
      real (kind=dp) :: cwd, cm

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), & ! only rq is used
                    e(nkt),dew(nkt),rq(nkt,nka)
      real (kind=dp) :: enw, ew, rn, rw, en, e, dew, rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n) ! only ff is used
      real (kind=dp) :: ff, fsum
      integer :: nar

!- End of header ---------------------------------------------------------------

      if (ifeed.eq.2) then
         ial = 2
      else
         ial = 1
      endif

! rcd(n,2): mean radius of aerosols in m -------
      do k=nf+1,nmax
! calculate rc for all rel. humidities drier than crystallization
         rcd1=0.
         rcd2=0.
         cwd1=0.
         cwd2=0.

! here TOTAL particle volume is used for calculating LWC
! this is correct only for completely soluble aerosol
! small aerosol

         do ia=ial,ka
            do jt=1,kw(ia)
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cwd1=cwd1+x0
               rcd1=rcd1+x0*rq(jt,ia)
            enddo
         enddo
! large aerosol
         do ia=ka+1,nka
            do jt=1,kw(ia)
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cwd2=cwd2+x0
               rcd2=rcd2+x0*rq(jt,ia)
            enddo
         enddo
! conversion: um^3/cm^3 --> m^3(aq)/m^3(air):10^-12
!           : um        --> m               :10^-6
         if (cwd1.gt.0.d0) then
            rcd(1,k)=rcd1/cwd1*1.d-6
         else
            rcd(1,k) = 0.d0
         end if

         if (cwd2.gt.0.d0) then
            rcd(2,k)=rcd2/cwd2*1.d-6
         else
            rcd(2,k) = 0.d0
         end if

         cwd(1,k)=cwd1*1.d-12
         cwd(2,k)=cwd2*1.d-12

      end do

      end subroutine dry_cw_rc


!
!------------------------------------------------------------
!

      subroutine dry_rates_g (tt,pp,nmax)
! calculates kmt for heterogeneous reactions on dry aerosol
! for ALL levels without real aerosol chemistry
! 1: sulfate aerosol, 2: seasalt aerosol

! to be included in gas-phase mechanism

      USE global_params, ONLY : &
! Imported Parameters:
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      include 'gas_Parameters.h'     !additional common blocks and other definitions

      real (kind=dp), intent(in) :: tt(n) ! Temperature
      real (kind=dp), intent(in) :: pp(n) ! Pressure
      integer, intent(in) :: nmax           ! Maximum layer for calculations

      integer,parameter :: ndr=4 ! number of species
      integer :: idr(ndr)
      integer :: k,kc,l
      real (kind=dp) :: x1, FCT
      real (kind=dp) :: xgamma(NSPEC,2),freep(n),vmean(NSPEC,n)
      common /blck11/ rcd(nkc,n) ! average particle radius in chemical bin
      real (kind=dp) :: rcd      !  note that in the context of this subroutine, this
                                 !  radius is considered as a "dry" radius, thus labeled rcd

      common /kpp_dryg/ xkmtd(n,2,NSPEC),henry(n,NSPEC),xeq(n,NSPEC)
      real (kind=dp) :: xkmtd, henry, xeq

      real (kind=dp) :: func, funa, func3 ! local function
      real (kind=dp) :: a,a0,b0           !   and their variables
      integer :: k0                       !   ...

!      data idr/ind_HNO3,ind_N2O5,ind_BrNO3,ind_ClNO3,ind_HOBr/
!     data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4,ind_HCl/     ! jjb HCl is not handled here (but it should be, probably !)
      data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4/           ! jjb removed at the moment

      func(a,k)=dsqrt(tt(k)/a)*4.60138
      funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))
      func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))

! calculate freep for all heights:
      do k=2,nmax
         freep(k)=2.28d-5 * tt(k) / pp(k)
      enddo

! define gamma's for all species ----
      xgamma(ind_HNO3,1)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_HNO3,2)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_N2O5,1)  = 0.02       !estimated from JPL 2003
      xgamma(ind_N2O5,2)  = 0.02       !estimated from JPL 2003
      xgamma(ind_NH3,1)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_NH3,2)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_H2SO4,1) = 0.1        !estimated from alpha=0.65
      xgamma(ind_H2SO4,2) = 0.1        !estimated from alpha=0.65
!      xgamma(ind_HCl,1)   =
!      xgamma(ind_HCl,2)   =
!      xgamma(ind_BrNO3,1) = 0.3
!      xgamma(ind_BrNO3,2) = 0.3
!      xgamma(ind_ClNO3,1) = 0.3
!      xgamma(ind_ClNO3,2) = 0.3
!      xgamma(ind_HOBr,1)  = 0.2
!      xgamma(ind_HOBr,2)  = 0.2
!      xgamma(ind_HOCl,1)  = 0.
!      xgamma(ind_HOCl,2)  = 0.
!      xgamma(ind_HOI,1)   = 0.
!      xgamma(ind_HOI,2)   = 0.

! the following are needed to calculate dry rates w/ assumption of Henry's
! law equilibrium

! define equilibrium constant
! obviously without Pitzer coefficients
      do k=2,nmax
         xeq(k,ind_HNO3) = funa(1.54d+1,8700.d0,k) ! jjb note this is different in SR equil_co_*
      enddo

! define inverse Henry's constant
      do k=2,nmax
         henry(k,ind_HNO3)=func3(2.5d6/xeq(k,ind_HNO3),8694.d0,k) ! jjb note this is different in SR equil_co_*
         FCT=0.0820577*tt(k)
         do l=1,ndr
            if (henry(k,idr(l)).ne.0.d0) henry(k,idr(l))=1./ & ! jjb .d0
                 (henry(k,idr(l))*FCT)
! "else": henry=0 <=> k_H^cc=infinity
         enddo
      enddo

! mean molecular speed (see SR v_mean_* for details)
! v_mean in m/s; sqrt(8*R_gas/pi)=4.60138
      do k=2,nmax
         vmean(ind_HNO3,k)  = func(6.3d-2,k)
         vmean(ind_N2O5,k)  = func(1.08d-1,k)
         vmean(ind_NH3,k)   = func(1.7d-2,k)
         vmean(ind_H2SO4,k) = func(9.8d-2,k)
!         vmean(ind_HCl,k)   = func(3.6d-2,k)
!         vmean(ind_BrNO3,k) = func(1.42d-1,k)
!         vmean(ind_ClNO3,k) = func(9.7d-2,k)
!         vmean(ind_HOBr,k)  = func(9.7d-2,k)
!         vmean(ind_HOCl,k)  = func(5.2d-2,k)
!         vmean(ind_HOI,k)   = func(1.44d-1,k)
      enddo

! calculate kmt ----

! k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha))
! if : D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
! k_mt=v_mean/(r**2/lambda+4*r/(3*alpha))

! lambda=freep : mean free path
!     x1=0. ! jjb useless
      do k=2,nmax
         do l=1,ndr
            do kc=1,2
               if (xgamma(idr(l),kc).gt.0..and.rcd(kc,k).gt.0.) then
                  x1=1./(rcd(kc,k)*(rcd(kc,k)/freep(k)+4./(3.*xgamma &
                       (idr(l),kc))))
               else
                  x1=0.
               endif
               xkmtd(k,kc,idr(l))=vmean(idr(l),k)*x1
!                  print *,k,idr(l),kc
!                  print *,xkmtd(k,kc,idr(l)),x1,vmean(idr(l),k)
            enddo
         enddo
      enddo

!      do k=2,n
!            write (432, 1001) k,xkmtd(k,1,ind_HNO3),1./xkmtd(k,1,ind_HNO3), &
!              xkmtd(k,2,ind_HNO3),1./xkmtd(k,2,ind_HNO3), &
!              xkmtd(k,1,ind_HNO3)*cwd(k,1),xkmtd(k,2,ind_HNO3)*cwd(k,2)
!      enddo
! 1001 format(i4, 6d16.8)

      end subroutine dry_rates_g


!
!------------------------------------------------------------
!

      subroutine dry_rates_a (freep,nmaxf)
! calculates kmt for heterogeneous reactions on dry aerosol
! 1: sulfate aerosol, 2: seasalt aerosol

! this SR calculates the gammas for the case where one aerosol class is
! allready above its crystallization point, but the other one is not yet
! so for the first class the complete aerosol chemistry is active, for
! the second only the reactions on dry aerosol (regulated via xliq1/2 and
! xhet1/2 in SR kpp_driver)

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h'     !additional common blocks and other definitions

      parameter (ndr=4)
      common /blck11/ rcd(nkc,n)
!     common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &     ! jjb  none of the objects of the common block is used
!                   e(nkt),dew(nkt),rq(nkt,nka)         !   (after commenting rqm = rq line below)
!      double precision enw, ew, rn, rw, en, e, dew, rq
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_drya/ xkmtd(nf,2,NSPEC),xeq(nf,NSPEC)
!      common /kpp_dryp/ rcd(n,2),cwd(n,2)
!     dimension xgamma(NSPEC,2),freep(nf),idr(ndr),rqm(nkt,nka) ! jjb rqm is unreferenced
      dimension xgamma(NSPEC,2),freep(nf),idr(ndr)              ! jjb thus removed


!      data idr/ind_HNO3,ind_N2O5,ind_BrNO3,ind_ClNO3,ind_HOBr/
!     data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4,ind_HCl/     ! jjb HCl is not handled here (but it should be, probably !)
      data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4/           ! jjb removed at the moment

!     funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3)) ! jjb another copy-paste mistake from SR dry_rates_g
      funa(a0,b0,k)=a0*exp(b0*(1/t(k)-3.354d-3))  ! jjb here the temperature is passed through a CB, not as a parameter

! change rq in um to rqm in m            ! jjb rqm is not used in this SR
!     do ia=1,nka                        !   ( see commented block at the end)
!        do jt=1,nkt
!           rqm(jt,ia)=rq(jt,ia)*1.d-6
!        enddo
!     enddo

! define gamma's for all species ----
      xgamma(ind_HNO3,1)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_HNO3,2)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_N2O5,1)  = 0.02       !estimated from JPL 2003
      xgamma(ind_N2O5,2)  = 0.02       !estimated from JPL 2003
      xgamma(ind_NH3,1)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_NH3,2)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_H2SO4,1) = 0.1        !estimated from alpha=0.65
      xgamma(ind_H2SO4,2) = 0.1        !estimated from alpha=0.65
!      xgamma(ind_BrNO3,1)=0.3
!      xgamma(ind_BrNO3,2)=0.3
!      xgamma(ind_ClNO3,1)=0.3
!      xgamma(ind_ClNO3,2)=0.3
!      xgamma(ind_HOBr,1)=0.2
!      xgamma(ind_HOBr,2)=0.2
!      xgamma(ind_HOCl,1)=0.
!      xgamma(ind_HOCl,2)=0.
!      xgamma(ind_HOI,1)=0.
!      xgamma(ind_HOI,2)=0.


! the following are needed to calculate dry rates w/ assumption of Henry's
! law equilibrium

! define equilibrium constant
! obviously without Pitzer coefficients
!     do k=2,nmax  ! jjb nmax is not defined (copy-paste mistake from SR dry_rates_g)
      do k=2,nmaxf ! jjb nmaxf is the correct index here
         xeq(k,ind_HNO3) = funa(1.54d+1,8700.d0,k) ! jjb note this is different in SR equil_co_*
      enddo

! calculate kmt ----

! k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha))
! if : D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
! k_mt=v_mean/(r**2/lambda+4*r/(3*alpha))

! lambda=freep : mean free path
!     x1=0. ! jjb useless
      do k=2,nmaxf
         do l=1,ndr
            do kc=1,2
               if (xgamma(idr(l),kc).gt.0..and.rcd(kc,k).gt.0.) then
                  x1=1./(rcd(kc,k)*(rcd(kc,k)/freep(k)+4./(3.*xgamma &
                       (idr(l),kc))))
               else
                  x1=0.
               endif
               xkmtd(k,kc,idr(l))=vmean(idr(l),k)*x1
!                  print *,k,idr(l),kc
!                  print *,xkmtd(k,kc,idr(l)),x1,vmean(idr(l),k)
            enddo
         enddo
      enddo

!      write (434,1100)
!      do k=2,n
!            write (434, 1101) k,xkmtd(k,1,ind_HNO3),1./xkmtd(k,1,ind_HNO3), &
!              xkmtd(k,2,ind_HNO3),1./xkmtd(k,2,ind_HNO3), &
!              xkmtd(k,1,ind_HNO3)*cwd(k,1),xkmtd(k,2,ind_HNO3)*cwd(k,2)
!      enddo
! 1100 format ('bulk values')
! 1101 format(i4, 6d16.8)
! 1102 format ('integrated values')


! now the same but with integrated values - - - - - - - - - - - - - - - -
! if working: include smart way of saving CPU time
! this leads to kmt that are about 10x greater than the bulk approach,
! nevertheless use bulk as this is what gamma's are determined for, plus
! it's faster (CPU-wise)


!! loop over vertical grid
!      do  k=2,nmaxf
!! loop over the nkc different chemical bins
!         do kc=1,2!nkc
!! loop over the species to be exchanged between gas and aqueous phase---
!            do l=1,ndr
!               xkmtd(k,kc,idr(l))=0.
!! define summation limits (1) ---
!               if (kc.eq.1) then
!                  iia_0=1
!                  iia_e=ka
!               endif
!               if (kc.eq.2) then
!                  iia_0=ka+1
!                  iia_e=nka
!               endif
!               ! kc eq. 3/4
!! fast version without logarithmic integration
!               x1=0.
!               xk1=0.
!               if (l.eq.1) xx1=0.
!               if (alpha(k,idr(l)).gt.0.) x1=4./(3.*alpha(k,idr(l)))
!               do ia=iia_0,iia_e
!! define summation limits (2)
!                  if (kc.eq.1) then
!                     jjt_0=1
!                     jjt_e=kw(ia)
!                  endif
!                  if (kc.eq.2) then
!                     jjt_0=1
!                     jjt_e=kw(ia)
!                  endif
!                  ! kc eq. 3/4
!                  do jt=jjt_0,jjt_e
!! conversion: um      --> m               :10^-6
!                     rqq=rqm(jt,ia)
!! kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
!!     with D_g=lambda*vmean/3.
!! here a volume weighted value is calculated, therefore weighting with r^3:
!! kmt=4/3*pi/L*sum(a)*sum(r){r^3*N*kt}
!                     x2=vmean(idr(l),k)/(rqq/freep(k)+x1) ![1/s]
!! conversion: 1/cm^3 --> 1/m^3(air):10^6
!                     xk1=xk1+x2*rqq*rqq*ff(jt,ia,k)*1.d6
!                  enddo !jt
!               enddo !ia
!! k_mt=4*pi/(3*LWC)*sum
!!               xkmtd(k,kc,idr(l))=4.*3.1415927/(3.*cwd(k,kc))*xk1 ![1/s]
!               xkmtd(k,kc,idr(l))=4.*pi/(3.*cwd(k,kc))*xk1 ![1/s]
!            enddo !l
!         enddo                  ! kc
!      enddo                     !k

!      write (434,1102)
!      do k=2,n
!            write (434, 1101) k,xkmtd(k,1,ind_HNO3),1./xkmtd(k,1,ind_HNO3), &
!              xkmtd(k,2,ind_HNO3),1./xkmtd(k,2,ind_HNO3), &
!              xkmtd(k,1,ind_HNO3)*cwd(k,1),xkmtd(k,2,ind_HNO3)*cwd(k,2)
!      enddo

      end subroutine dry_rates_a

!
!------------------------------------------------------------
!

      subroutine dry_rates_t (freep,nmaxf)
! calculates kmt for heterogeneous reactions on dry aerosol
! 1: sulfate aerosol, 2: seasalt aerosol

! this SR calculates the gammas for the case where one aerosol class is
! allready above its crystallization point, but the other one is not yet
! so for the first class the complete aerosol chemistry is active, for
! the second only the reactions on dry aerosol (regulated via xliq1/2 and
! xhet1/2 in SR kpp_driver)

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h'     !additional common blocks and other definitions

      parameter (ndr=4)
      common /blck11/ rcd(nkc,n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho
      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_dryt/ xkmtd(nf,2,NSPEC),xeq(nf,NSPEC)
!     common /kpp_dryp/ rcd(n,2),cwd(n,2)
!     dimension xgamma(NSPEC,2),freep(nf),idr(ndr),rqm(nkt,nka) ! jjb rqm is unreferenced
      dimension xgamma(NSPEC,2),freep(nf),idr(ndr)              ! jjb thus removed


!      data idr/ind_HNO3,ind_N2O5,ind_BrNO3,ind_ClNO3,ind_HOBr/
!     data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4,ind_HCl/     ! jjb HCl is not handled here (but it should be, probably !)
      data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4/           ! jjb removed at the moment

!     funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3)) ! jjb another copy-paste mistake from SR dry_rates_g
      funa(a0,b0,k)=a0*exp(b0*(1/t(k)-3.354d-3))  ! jjb here the temperature is passed through a CB, not as a parameter

! change rq in um to rqm in m            ! jjb rqm is not used in this SR
!     do ia=1,nka                        !   ( see commented block at the end)
!        do jt=1,nkt
!           rqm(jt,ia)=rq(jt,ia)*1.d-6
!        enddo
!     enddo

! define gamma's for all species ----
      xgamma(ind_HNO3,1)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_HNO3,2)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_N2O5,1)  = 0.02       !estimated from JPL 2003
      xgamma(ind_N2O5,2)  = 0.02       !estimated from JPL 2003
      xgamma(ind_NH3,1)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_NH3,2)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_H2SO4,1) = 0.1        !estimated from alpha=0.65
      xgamma(ind_H2SO4,2) = 0.1        !estimated from alpha=0.65
!      xgamma(ind_BrNO3,1)=0.3
!      xgamma(ind_BrNO3,2)=0.3
!      xgamma(ind_ClNO3,1)=0.3
!      xgamma(ind_ClNO3,2)=0.3
!      xgamma(ind_HOBr,1)=0.2
!      xgamma(ind_HOBr,2)=0.2
!      xgamma(ind_HOCl,1)=0.
!      xgamma(ind_HOCl,2)=0.
!      xgamma(ind_HOI,1)=0.
!      xgamma(ind_HOI,2)=0.


! the following are needed to calculate dry rates w/ assumption of Henry's
! law equilibrium

! define equilibrium constant
! obviously without Pitzer coefficients
!     do k=2,nmax  ! jjb nmax is not defined (copy-paste mistake from SR dry_rates_g)
      do k=2,nmaxf ! jjb nmaxf is the correct index here
         xeq(k,ind_HNO3) = funa(1.54d+1,8700.d0,k) ! jjb note this is different in SR equil_co_*
      enddo

! calculate kmt ----

! k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha))
! if : D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
! k_mt=v_mean/(r**2/lambda+4*r/(3*alpha))

! lambda=freep : mean free path
!     x1=0. ! jjb useless
      do k=2,nmaxf
         do l=1,ndr
            do kc=1,2
               if (xgamma(idr(l),kc).gt.0..and.rcd(kc,k).gt.0.) then
                  x1=1./(rcd(kc,k)*(rcd(kc,k)/freep(k)+4./(3.*xgamma &
                       (idr(l),kc))))
               else
                  x1=0.
               endif
               xkmtd(k,kc,idr(l))=vmean(idr(l),k)*x1
            enddo
         enddo
      enddo

      end subroutine dry_rates_t

!
!------------------------------------------------------------
!

      subroutine activ (box,n_bl)
! front end for Beiping Luo's model to calculate activity coefficients

! jjb work done:
!     cleaned, declared all variables, implicit none, removed old, commented stuff
!     cleaned final write (+ avoid useless write when cm < cm_min)
!     use parameter for min cm value
!     (probable) bugfix to be done: loop kc = 1,nkc: use nkc_l instead of nkc, if nkc_l is also
!               used (future development) to restrict tot mechanism
!     use xip (not xit) to check is pitzer will be called or not

! 03-Mar-2017  <J. Bock> When Pitzer is not called, set gamma=1 (default value) to avoid previous
!                        values to be still used

! 04-Mar-2017  <J. Bock> Introduce "xip2" to check the lower ionic strengh value allowed in Piter.
!                        See SR gammann: 1./xo4 where xo4=(omega1*I2)**4 and I2=sqrt(I) (I=xip here)
!                        min(omega(i,j)) = 0.43 thus xip2 = 0.01*xip**2

      USE global_params, ONLY : &
! Imported Parameters:
           j2, &
           j6, &
           nf, &
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      logical, intent(in) :: box
      integer, intent(in) :: n_bl

! Local parameters:
      integer, parameter :: nk = 10  ! Number of vertical levels outputed
      integer, parameter :: kk(nk) = (/2,12,22,32,42,52,62,72,82,92/) ! Indexes of levels in output

! Local scalars:
      integer :: jg, k, kc   ! loop indexes for liq. species, vertical layers, liq. bins
      integer :: nmin, nmax  ! lower and upper bounds for calculations

      real (kind=dp) :: wact
      real (kind=dp) :: xip2

! Local arrays:
      real (kind=dp) :: wa(nkc,nf)
      real (kind=dp) :: xip(nkc,nf), xit(nkc,nf)

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
      real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

      common /blck12/ cw(nkc,n),cm(nkc,n)
      real (kind=dp) :: cw, cm

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      real (kind=dp) :: sl1, sion1

      common /kpp_mol/ xgamma(nf,j6,nkc)
      real (kind=dp) :: xgamma


! Initialise lower and upper bounds (in the vertical grid) to do the calculations
      nmin=2
      nmax=nf
      if (box) then
         nmin=n_bl
         nmax=n_bl
      endif

! Initialise gammas
      xgamma(:,:,:) = 1.d0

      do k=nmin,nmax
         do kc=1,nkc

            ! Skip if cm too small
            if (cm(kc,k) <= tiny(0.d0)) cycle

! calculate ionic strength of species accounted for in pitzer module
!       ionic strength: I=0.5*sum(molality*charge^2)
! concentrations --> molality:
!  mol/m^3_air * m^3_air/m^3_solvent * 10^-3 * 1 dm^3/kg (density of water) = mol/kg_solvent
!   sion1      * cm^-1               * 1.d-3
            xip(kc,k)=0.5*(sion1(1,kc,k)+sion1(2,kc,k)+sion1(20,kc,k)+ &
                 sion1(19,kc,k)+4.*sion1(8,kc,k)+sion1(13,kc,k)+ &
                 sion1(14,kc,k))
            xip(kc,k)=xip(kc,k)*1.d-3/cm(kc,k)
! calculate ionic strength of all species
            xit(kc,k)=0.5*(sion1(1,kc,k)+sion1(2,kc,k)+sion1(3,kc,k)+ &
                 sion1(4,kc,k)+sion1(5,kc,k)+4.*sion1(6,kc,k)+ &
                 sion1(7,kc,k)+4.*sion1(8,kc,k)+sion1(9,kc,k)+ &
                 sion1(10,kc,k)+sion1(11,kc,k)+sion1(12,kc,k)+ &
                 sion1(13,kc,k)+sion1(14,kc,k)+sion1(15,kc,k)+ &
                 sion1(16,kc,k)+sion1(19,kc,k)+sion1(20,kc,k)+ &
                 sion1(21,kc,k)+sion1(22,kc,k)+sion1(23,kc,k)+ &
                 sion1(24,kc,k)+sion1(25,kc,k)+sion1(26,kc,k)+ &
                 sion1(27,kc,k)+sion1(28,kc,k)+sion1(29,kc,k)+ &
                 sion1(30,kc,k)+sion1(31,kc,k)+sion1(32,kc,k)+ &
                 sion1(33,kc,k)+sion1(34,kc,k)+sion1(35,kc,k)+ &
                 sion1(36,kc,k)+sion1(37,kc,k)+sion1(38,kc,k)+ &
                 sion1(39,kc,k)) &
                 *1.d-3/cm(kc,k)


! check that ionic strength remains in a "reasonable" range; if too high there
! will also be a "floating overflow" in SR pitzer
            xip2 = 0.01*xip(kc,k)**2 ! <jjb> see explanations in header. 0.01 < 0.43**4
            !if (xit(kc,k).le.80..and.xit(kc,k).gt.0.) then
            if (xip(kc,k).le.80..and.xip2.gt.tiny(0.)) then
! calculate pitzer coefficients for layer k and liquid size bin kc
               call pitzer (k,kc,wact) ! for sulfate data can be used for ionic strengths
                                       ! up to 40 M, for sea salt only up to 6 M
               wa(kc,k)=wact
            else
               print *,'ionic strength > 80. or <=0.',k,kc
               print *,'I=',xip(kc,k),xit(kc,k)
               cycle
            endif

! gamma's calculated by SR pitzer are for molalities, used here are molarities ==> conversion factor needed:
! mol/kg(solvent) --> mol/l(solution): cm/cw (assuming unity density for the solvent water)


! don't apply to all xgamma's only to those that are <> 1:
            !if (cw(kc,k).gt.0.d0) then ! jjb: test not needed, already excluded by cm above
               do jg=1,j6
                  if (xgamma(k,jg,kc).ne.1.) xgamma(k,jg,kc) = &
                       xgamma(k,jg,kc) * cm(kc,k)/cw(kc,k)
               enddo
            !end if

! define gamma's for species that are not included in pitzer module
! L+J: Liang and Jacobson, 1999, JGR, 104, 13749, #554
! C+S: Chameides and Stelson, 1992, JGR, 97,20565,#470
            xgamma(k,3,kc)=xgamma(k,13,kc) !OH-   = NO3- (assumed, Luo, pers comm 2000)
            xgamma(k,5,kc)=xgamma(k,19,kc) !HSO3- = HSO4- (L+J)
            xgamma(k,6,kc)=xgamma(k,8,kc)  !SO3=  = SO4=  (L+J)
            xgamma(k,7,kc)=xgamma(k,19,kc) !SO4-  = HSO4- (L+J)
            xgamma(k,9,kc)=xgamma(k,5,kc)  !HCO3- = HSO3- (C+S)
            xgamma(k,11,kc)=xgamma(k,5,kc) !O2-   = Cl2- = HSO3- (C+S)
            xgamma(k,12,kc)=xgamma(k,13,kc)!NO2-  = NO3-  (L+J)
            xgamma(k,15,kc)=xgamma(k,5,kc) !Cl2-  = HSO3- (C+S)
            xgamma(k,16,kc)=xgamma(k,5,kc)  !HCOO- = HSO3- (assumed)????
!            xgamma(k,21,kc)=xgamma(k,,kc) !NO4-  = ??   (assumed)
            xgamma(k,22,kc)=xgamma(k,14,kc)!ClO-  = Cl-   (assumed)
            xgamma(k,24,kc)=xgamma(k,14,kc)!Br-   = Cl-   (assumed)
            xgamma(k,25,kc)=xgamma(k,5,kc) !Br2-  = HSO3- (C+S)
            xgamma(k,26,kc)=xgamma(k,24,kc)!BrO-  = Br-   (assumed)
!            xgamma(k,28,kc)=xgamma(k,,kc)  !BrCl2-= ? not used in SR equil_co*
!            xgamma(k,29,kc)=xgamma(k,,kc)  !Br2Cl-= ?       -"-
            xgamma(k,37,kc)=xgamma(k,5,kc) !ICl2- = HSO3- (assumed)
            xgamma(k,38,kc)=xgamma(k,5,kc) !IBr2- = HSO3- (assumed)

         enddo  !kc
      enddo     !k

! output
      if (lmin/30*30.eq.lmin) then
 3000    continue
         open (109,file='gam.out',status='unknown',position='append', &
              err=3000)
         do k=1,nk
            do kc=1,nkc
               write (109,20) kk(k),kc,lday,lst,lmin,feu(kk(k)), &
                              cm(kc,kk(k)),cw(kc,kk(k)),xm2(kk(k))*1.d-3

               ! Skip next 'write' instructions if cm too small
               if (cm(kc,kk(k)) <= tiny(0.d0)) then
                  write(109,*)' -> cm too small, xgamma set equal to 1.'
                  cycle
               end if

               write (109,13) xip(kc,kk(k)),xit(kc,kk(k)),wa(kc,kk(k))
               write (109,10) (xgamma(kk(k),jg,kc),jg=1,j6)
               write (109,10) (sion1(jg,kc,kk(k))*1.d-3/cm(kc,kk(k)), &
                               jg=1,j6)
            enddo
         enddo
         close (109)
      endif
 10   format (8d14.6)
 13   format (3d14.6)
 20   format ('layer ',i3,' bin ',i1,' day ',2i3,':',i2,' rH ',d12.4, &
              ' LWC ',3d12.4)

      end subroutine activ

!
!------------------------------------------------------------
!

      subroutine activ_init
! initialise activity coefficients (gammas)

      USE global_params, ONLY : &
! Imported Parameters:
           j6, &
           nf, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      common /kpp_mol/ xgamma(nf,j6,nkc)
      real (kind=dp) xgamma

      xgamma(:,:,:) = 1.d0

      open (109,file='gam.out',status='unknown',form='formatted')
      write(109,*)'File opened by SR activ_init, written by SR activ'
      write(109,*)'  For each outputed layer:'
      write(109,*)'  partial ionic strength, total ionic str, water act'
      write(109,*)'  gamma (all ionic species)'
      write(109,*)'  molalities (all ionic species)'
      write(109,*)'  --------------------------------------------------'
      close (109)

      end subroutine activ_init

!
!-------------------------------------------------------------
!

      subroutine gasdrydep (xra,tt,rho,freep)
! calculate dry deposition velocities for gas phase after Sehmel cited in
! Seinfeld and Pandis, 1999


! jjb IMPORTANT NOTE
!
! Some points need to be addressed in this subroutine.

! 1) use mk tables to convert mistra <--> KPP indexes
! 2) change the if test in the middle (after hs(i) and "FCT" definitions) : if hs .gt. 0 (instead of .ne.)
!    and replace the test at the end : if hs == -1 instead of if hs == -1/FCT
! 3) correction of one formula, double check

! jjb 27-05-2021
!   In this routine, arrays vm, hs and f0 (for gas) have dimension ind_gas(j1) (= j1_fake)
!   This allows the user tu use the same indexing as used in gas_species, instead of the
!   "compressed" indexing in j1
!   The output of this routine, vg, is converted into compressed dimension through ind_gas
!   indexing of arrays vm, hs and f0, which pick up only relevant values.

      USE config, ONLY : &
           lpJoyce14bc

      USE data_surface, ONLY : &
           ustern               ! frictional velocity

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           n, &
           nkc

      USE gas_common, ONLY : &
! Imported Parameters:
           j1, &
! Imported Array Variables with intent(in):
           ind_gas, &
!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!           vmean, &
           vg

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      include 'aer_Parameters.h'     !additional common blocks and other definitions

! Subroutine arguments
      real (kind=dp), intent(in) :: xra,tt(n),rho(n),freep(nf)
! Local variables
      integer :: i, k
      real (kind=dp) :: FCT, rb_fact
      real (kind=dp) :: sac, xeta, xnu, xr_ac, xr_cl, xr_clO, xr_clS, &
           xr_dc, xr_gs, xr_gsO, xr_gsS, xrb, xrc
      real (kind=dp) :: vm(ind_gas(j1)),hs(ind_gas(j1)),f0(ind_gas(j1))

! Common blocks
      common /cb48/ sk,sl,dtrad(n),dtcon(n) ! sk is used in Joyce config
      real (kind=dp) :: sk, sl, dtrad, dtcon
      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
           xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      real (kind=dp) :: henry, xkmt, xkef, xkeb
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      real (kind=dp) :: alpha, vmean


! function used in calculation of Hstar
      real (kind=dp) :: funa, a0, b0
      funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))

! gas phase dry deposition velocity:v_d=1/(ra + rb + rc)
!    ra=1/(kappa ustar) (ln (z/z0) +Phi) ;where z=height of surface (constant flux) layer
!                                         Phi takes stratification into account
!      see SR partdep, ra=xra
!    rb=5.*Sc^(2/3)/ustar ;Sc=nu/D, nu=kinematic viscosity of air, D diffusivity of gas
!      ustar=ustern
!      gas phase diffusivity is calculated for each species using the mean molecular speed vmean:
!        D=lambda*vmean/3.  (lambda=freep)
! rc after Sehmel, 1980
!    rc=2.54d+4/(Hstar * T * ustar) ; Hstar=effective Henry constant
! rc after Wesely, 1989
!     rc=1./(H*/(10^5*r_gsS)+f_0/r_gsO)  ; r_gsS=1., r_gsO=2000., f_0: standard 0.1
!                                                                      low reactivity:0., high:1.

!      Hstar calculated from henry and sea water pH=8.1 (Riley and Skirrow, 1965)

      k=2 ! only in lowest model layer
      xeta=1.8325d-5*(416.16/(tt(k)+120.))*((tt(k)/296.16)**1.5) !dynamic viscosity of air, Jacobsen p. 92
      xnu=xeta/rho(k)  !kinematic viscosity of air

! get vmean for j1-list (no deposition or radicals or fixed)
! get vmean for j1-order from KPP names
      vm(1)=vmean(ind_NO,k)
      vm(2)=vmean(ind_NO2,k)
      vm(3)=vmean(ind_HNO3,k)
      vm(4)=vmean(ind_NH3,k)
      vm(5)=vmean(ind_SO2,k)
      vm(6)=vmean(ind_H2SO4,k)
      vm(7)=vmean(ind_O3,k)
      vm(8)=vmean(ind_CH4,k)
      vm(9)=vmean(ind_C2H6,k)
!      vm(10)=vmean(ind_C3H8,k)
!      vm(11)=vmean(ind_ALKA,k)
      vm(12)=vmean(ind_ETHE,k)
!      vm(13)=vmean(ind_ALKE,k)
!      vm(14)=vmean(ind_AROM,k)
      vm(15)=vmean(ind_ACO2,k)
      vm(16)=vmean(ind_ACTA,k)
      vm(17)=vmean(ind_HCHO,k)
      vm(18)=vmean(ind_ALD2,k)
      vm(19)=vmean(ind_H2O2,k)
      vm(20)=vmean(ind_ROOH,k)
      vm(21)=vmean(ind_HONO,k)
      vm(22)=vmean(ind_PAN,k)
!      vm(23)=vmean(ind_TPAN,k)
!      vm(24)=vmean(ind_KET,k)
!      vm(25)=vmean(ind_CRES,k)
!      vm(26)=vmean(ind_DIAL,k)
!      vm(27)=vmean(ind_GLYX,k)
!      vm(28)=vmean(ind_MGLY,k)
!      vm(29)=vmean(ind_NH4NO3,k)
      vm(30)=vmean(ind_HCl,k)
!      vm(31)=vmean(ind_R3N2,k)
!      vm(32)=vmean(ind_RAN2,k)
!      vm(33)=vmean(ind_RAN1,k)
      vm(34)=vmean(ind_N2O5,k)
      vm(35)=vmean(ind_HNO4,k)
      vm(36)=vmean(ind_NO3,k)
      vm(37)=vmean(ind_DMS,k)
      vm(38)=vmean(ind_HOCl,k)
      vm(39)=vmean(ind_ClNO2,k)
      vm(40)=vmean(ind_ClNO3,k)
      vm(41)=vmean(ind_Cl2,k)
      vm(42)=vmean(ind_HBr,k)
      vm(43)=vmean(ind_HOBr,k)
      vm(44)=vmean(ind_BrNO2,k)
      vm(45)=vmean(ind_BrNO3,k)
      vm(46)=vmean(ind_Br2,k)
      vm(47)=vmean(ind_BrCl,k)
      vm(48)=vmean(ind_HI,k)
      vm(49)=vmean(ind_HOI,k)
      vm(50)=vmean(ind_I2O2,k)
      vm(51)=vmean(ind_INO2,k)
      vm(52)=vmean(ind_INO3,k)
      vm(53)=vmean(ind_I2,k)
      vm(54)=vmean(ind_ICl,k)
      vm(55)=vmean(ind_IBr,k)
      vm(56)=vmean(ind_CH3I,k)
      vm(57)=vmean(ind_CH2I2,k)
      vm(58)=vmean(ind_CH2ClI,k)
      vm(59)=vmean(ind_C3H7I,k)
      vm(60)=vmean(ind_DMSO,k)
      vm(61)=vmean(ind_CH3SO2,k)
      vm(62)=vmean(ind_CH3SO3,k)
      vm(63)=vmean(ind_CH3SO3H,k) !MSA=CH3S(OO)OH
      vm(64)=vmean(ind_CO,k)
      vm(65)=vmean(ind_Cl2O2,k)
      vm(66)=vmean(ind_DMOO,k) ! CH3SCH2OO
      vm(67)=vmean(ind_CH3S,k)
      vm(68)=vmean(ind_CH3SO,k)
      vm(69)=vmean(ind_CH3SO2H,k) ! MSIA=CH3S(O)OH
      vm(70)=vmean(ind_DMSO2,k)
      vm(71)=vmean(ind_CH2BrI,k)
!      vm(72)=vmean(ind_CHBr2I,k)
      vm(73)=vmean(ind_C2H5I,k)
      vm(74)=vmean(ind_HIO3,k)
!      vm(75)=vmean(ind_NUCV,k)
      vm(76)=vmean(ind_SO3,k)
      vm(77)=vmean(ind_HOSO2,k)
      vm(78)=vmean(ind_CO2,k)
!      vm(79)=vmean(ind_I2O,k)
!      vm(80)=vmean(ind_I2O3,k)
!      vm(81)=vmean(ind_I2O4,k)
!      vm(82)=vmean(ind_I2O5,k)
!      vm(83)=vmean(ind_INO,k)
      vm(84)=vmean(ind_Br2O,k)
      vm(85)=vmean(ind_ClONO,k)
      vm(86)=vmean(ind_ClO3,k)
      vm(87)=vmean(ind_Cl2O3,k)
      vm(88)=vmean(ind_CH3OH,k)
      vm(89)=vmean(ind_C2H5OH,k)
      vm(90)=vmean(ind_H2,k)
      vm(91)=vmean(ind_NHS,k)
      vm(92)=vmean(ind_RCl,k)
      vm(93)=vmean(ind_RBr,k)
      vm(94)=vmean(ind_XOR,k)
      vm(95)=vmean(ind_SOR,k)
      vm(96)=vmean(ind_SPAN,k)
!      vm(97)=vmean(ind_Hg,k)
!      vm(98)=vmean(ind_HgO,k)
!      vm(99)=vmean(ind_HgCl,k)
!      vm(100)=vmean(ind_HgCl2,k)
!      vm(101)=vmean(ind_HgBr,k)
!      vm(102)=vmean(ind_HgBr2,k)

!     vm(radical 29)=vmean(ind_OIO,k)

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
      !!do jspec = 1,j1
      !   vm(:) = vmean(:j1,k)
      !!end do

! jjb add hs explicit initialisation
      hs(:) = 0._dp

! get henry constant for j1-order from KPP names and calculate Hstar
      sac=10.**(-8.1d0) ![H+]=10^(- pH), pH=8.1
      hs(1)=henry(ind_NO,k)
      hs(2)=henry(ind_NO2,k)
      hs(3)=henry(ind_HNO3,k)
      hs(4)=henry(ind_NH3,k)
      hs(5)=henry(ind_SO2,k)
      hs(6)=henry(ind_H2SO4,k)
      hs(7)=henry(ind_O3,k)
      hs(8)=henry(ind_CH4,k)
      hs(9)=henry(ind_C2H6,k)
!      hs(10)=henry(ind_C3H8,k)
!      hs(11)=henry(ind_ALKA,k)
      hs(12)=henry(ind_ETHE,k)
!      hs(13)=henry(ind_ALKE,k)
!      hs(14)=henry(ind_AROM,k)
      hs(15)=henry(ind_ACO2,k)
      hs(16)=henry(ind_ACTA,k)
      hs(17)=henry(ind_HCHO,k)
      hs(18)=henry(ind_ALD2,k)
      hs(19)=henry(ind_H2O2,k)
      hs(20)=henry(ind_ROOH,k)
      hs(21)=henry(ind_HONO,k)
      hs(22)=henry(ind_PAN,k)
!      hs(23)=henry(ind_TPAN,k)
!      hs(24)=henry(ind_KET,k)
!      hs(25)=henry(ind_CRES,k)
!      hs(26)=henry(ind_DIAL,k)
!      hs(27)=henry(ind_GLYX,k)
!      hs(28)=henry(ind_MGLY,k)
!      hs(29)=henry(ind_NH4NO3,k)
      hs(30)=henry(ind_HCl,k)
!      hs(31)=henry(ind_R3N2,k)
!      hs(32)=henry(ind_RAN2,k)
!      hs(33)=henry(ind_RAN1,k)
      if (.not.lpJoyce14bc) hs(34)=-1.!henry(ind_N2O5,k)
      hs(35)=henry(ind_HNO4,k)
      hs(36)=henry(ind_NO3,k)
      hs(37)=henry(ind_DMS,k)
      hs(38)=henry(ind_HOCl,k)
      hs(39)=henry(ind_ClNO2,k)
      if (.not.lpJoyce14bc) hs(40)=-1.!henry(ind_ClNO3,k)
      hs(41)=henry(ind_Cl2,k)
      hs(42)=henry(ind_HBr,k)
      hs(43)=henry(ind_HOBr,k)
      hs(44)=henry(ind_BrNO2,k)
      if (.not.lpJoyce14bc) hs(45)=-1.!henry(ind_BrNO3,k)
      hs(46)=henry(ind_Br2,k)
      hs(47)=henry(ind_BrCl,k)
      if (.not.lpJoyce14bc) hs(48)=-1.!henry(ind_HI,k)!*(1. + /sac)
      hs(49)=henry(ind_HOI,k)!*(1. + /sac)
      hs(50)=henry(ind_I2O2,k)
      hs(51)=henry(ind_INO2,k)
      if (.not.lpJoyce14bc) hs(52)=-1.!henry(ind_INO3,k)
      hs(53)=henry(ind_I2,k)
      hs(54)=henry(ind_ICl,k)
      hs(55)=henry(ind_IBr,k)
      hs(56)=henry(ind_CH3I,k)
      hs(57)=henry(ind_CH2I2,k)
      hs(58)=henry(ind_CH2ClI,k)
      hs(59)=henry(ind_C3H7I,k)
      hs(60)=henry(ind_DMSO,k)
      hs(61)=henry(ind_CH3SO2,k)
      hs(62)=henry(ind_CH3SO3,k)
      hs(63)=henry(ind_CH3SO3H,k)! MSA=CH3S(OO)OH
      hs(64)=henry(ind_CO,k)
      hs(65)=henry(ind_Cl2O2,k)
      hs(66)=henry(ind_DMOO,k) ! CH3SCH2OO
      hs(67)=henry(ind_CH3S,k)
      hs(68)=henry(ind_CH3SO,k)
      hs(69)=henry(ind_CH3SO2H,k) ! MSIA=CH3S(O)OH
      hs(70)=henry(ind_DMSO2,k)
      hs(71)=henry(ind_CH2BrI,k)
!      hs(72)=henry(ind_CHBr2I,k)
      hs(73)=henry(ind_C2H5I,k)
      hs(74)=henry(ind_HIO3,k)
!      hs(75)=henry(ind_NUCV,k)
      hs(76)=henry(ind_SO3,k)
      hs(77)=henry(ind_HOSO2,k)
      hs(78)=henry(ind_CO2,k)
!      hs(79)=henry(ind_I2O,k)
!      hs(80)=henry(ind_I2O3,k)
!      hs(81)=henry(ind_I2O4,k)
!      hs(82)=henry(ind_I2O5,k)
!      hs(83)=henry(ind_INO,k)
      hs(84)=henry(ind_Br2O,k)
      hs(85)=henry(ind_ClONO,k)
      hs(86)=henry(ind_ClO3,k)
      hs(87)=henry(ind_Cl2O3,k)
      hs(88)=henry(ind_CH3OH,k)
      hs(89)=henry(ind_C2H5OH,k)
      hs(90)=henry(ind_H2,k)
      hs(91)=henry(ind_NHS,k)
      hs(92)=henry(ind_RCl,k)
      hs(93)=henry(ind_RBr,k)
      hs(94)=henry(ind_XOR,k)
      hs(95)=henry(ind_SOR,k)
      hs(96)=henry(ind_SPAN,k)
!      hs(97)=henry(ind_Hg,k)
!      hs(98)=henry(ind_HgO,k)
!      hs(99)=henry(ind_HgCl,k)
!      hs(100)=henry(ind_HgCl2,k)
!      hs(101)=henry(ind_HgBr,k)
!      hs(102)=henry(ind_HgBr2,k)

! henry constant were already transformed into dimensionless values
      FCT=0.0820577*tt(k)
      do i=1,j1
         if (hs(ind_gas(i)).ne.0.)hs(ind_gas(i))=1./(hs(ind_gas(i))*FCT)
         f0(ind_gas(i))=0.1
      enddo
! some species dissociate:
      hs(3)=hs(3)*(1. +funa(1.54d+1,8700.d0,k)/sac)
      hs(4)=hs(4)*(1. + funa(1.7d-5,-4325.d0,k)* &
           sac/funa(1.d-14,-6710.d0,k))
      hs(5)=hs(5)*(1.+ funa(1.7d-2,2090.d0,k)/sac + &
           funa(1.7d-2,2090.d0,k)*funa(6.0d-8,1120.d0,k) /sac**2)
      hs(6)=hs(6)*(1. + 1.0d+3/sac + &
           1.0d+3*funa(1.02d-2,2720.d0,k)/sac**2)
      hs(30)=hs(30)*(1. + funa(1.7d6,6896.d0,k)/sac)
      hs(38)=hs(38)*(1. + 3.2d-8/sac)
      hs(42)=hs(42)*(1. + 1.d9/sac)
      hs(43)=hs(43)*(1. + funa(2.3d-9,-3091.d0,k)/sac)
! define some f_0 that deviate from standard (Pandis and Seinfeld, Table 19.3):
      f0(1) =0.  !NO
      f0(3) =0.  !HNO3
      f0(4) =0.  !NH3
      f0(5) =0.  !SO2
      f0(7) =1.  !O3
      f0(8) =0.  !CH4
      f0(9) =0.  !C2H6
      f0(10)=0.  !C3H8
      f0(11)=0.  !ALKA
      f0(14)=0.  !AROM
      f0(15)=0.  !HCOOH, formic acid
      f0(16)=0.  !CH3COOH, acetic acid
      f0(17)=0.  !HCHO
      f0(19)=1.  !H2O2
      f0(20)=1.  !ROOH
      f0(30)=0.  !HCl
!      f0(34)=1.  !N2O5
      f0(35)=0.  !HNO4
      f0(36)=1.  !NO3
      f0(42)=0.  !HBr
!      f0(43)=1. !HOBr - didn't find a reference for that, so I commented it again
!      f0(49)=1. !HOI
!      f0() =0.  !CH3CHO, acetaldehyde
!      f0() =0.  !CH3CH2CHO, propionaldehyde
!      f0() =0.3 !CH3OOH, methylhydroperoxide

! calculate vg from rb and rc for each species
!     rb_fact=5./ustern*(xnu*freep(k)/3.)**(2./3.) ! jjb mistake ?
      rb_fact=5./ustern*(3.*xnu/freep(k))**(2./3.)
!      rc_fact=2.54d+4/(tt(k)*ustern)  !Sehmel, 1980
!      print *,rb_fact,rc_fact
      do i=1,j1
       if (vm(ind_gas(i)).eq.0.) then
         vg(i)=0.
       else
        if (.not.lpJoyce14bc) then
! ======= GENERAL CASE =========
         if (hs(ind_gas(i)).ne.0.) then
!           vg(i)=1./(xra+(rb_fact/(vm(ind_gas(i))**(2./3.)))+ &
!                 (rc_fact/hs(ind_gas(i)))) !Sehmel, 1980
            vg(i)=1./(xra+(rb_fact/(vm(ind_gas(i))**(2./3.)))+ &
                  1./(hs(ind_gas(i))*1.d-5+f0(ind_gas(i))/2000.))  !Wesely, 1989

!! jjb 27-05-2021: under development, not validated yet. Use the old v_mean_a/t routines for now
!            vg(i)=1./(xra+(rb_fact/(vm(i)**(2./3.)))+1./(hs(ind_gas(i))* &
!                 1.d-5+f0(ind_gas(i))/2000.))  !Wesely, 1989

         else
!            vg(i)= 1./(xra+(rb_fact/(vm(ind_gas(i))**(2./3.)))+ &
!                   (rc_fact/1.d-7)) !hs=0 => 1/hs -> infinity,
                                                              !here just a value for hstar is chosen
                                                              !to get a small v_dd
            if (f0(ind_gas(i)) > 0.) then
!               print*,xra,rb_fact,vm(i),f0(ind_gas(i))
               vg(i)=1./(xra+(rb_fact/(vm(ind_gas(i))**(2./3.))) &
                     +1./(f0(ind_gas(i))/2000.))
            else
               vg(i)=0.
            end if
         endif
         if (hs(ind_gas(i)).eq.(-1./FCT)) &
              vg(i)=1./(xra+(rb_fact/(vm(ind_gas(i))**(2./3.)))+.1)     ! "infinite solubility"
! ======= SPECIAL CASE =========
        else if (lpJoyce14bc) then
! no need to check if hs=0 as long as not both hs=0 and f0=0
          xrb = rb_fact/(vm(ind_gas(i))**(2./3.))
! rc = 1/ ((1/r_dc+r_cl) + 1/(r_ac + r_gs)); (19.50)
! land use type: 6, mixed forest, wet-lands
! season: 4 - winter, snow, subfreezing
          xr_clS =  400.d0
          xr_clO =  600.d0
          xr_ac  = 1500.d0
          xr_gsS =  100.d0
          xr_gsO = 3500.d0
! the following equations are 19.54 - 19.56
          xr_dc = 100.d0*(1.d0+1000.d0/(sk+10.d0)) ! assume flat terrain, so last term in 19.54 = 1.
          xr_cl = 1.d0/(hs(ind_gas(i))*1.d-5/xr_clS &
               + f0(ind_gas(i))/xr_clO)
          xr_gs = 1.d0/(hs(ind_gas(i))*1.d-5/xr_gsS &
               + f0(ind_gas(i))/xr_gsO)
          xrc = 1.d0/(1.d0/(xr_dc + xr_cl) + 1.d0/(xr_ac + xr_gs))
          vg(i)=1.d0/(xra+xrb+xrc)

        end if ! general/special case (lpJoyce14bc)

       endif ! vm == 0
      enddo

      end subroutine gasdrydep


!
!-----------------------------------------------------------------------------
!

      subroutine mass_ch
! calculates total number of Br and Cl atoms in all phases [mol/m2]
! including deposited atoms
! emitted atoms (see SR aer_source) are subtracted

      USE config, ONLY: &
           nkc_l

      USE gas_common, ONLY: &
           s1, &
           j1_br, ind_gas_br, &
           j1_cl, ind_gas_cl, &
           s3, &
           j5_br, ind_rad_br, &
           j5_cl, ind_rad_cl

      USE global_params, ONLY : &
! Imported Parameters:
           j2, &
           j3, &
           j6, &
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit double precision (a-h,o-z)
      logical cl_out

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      real (kind=dp) :: detw, deta, eta, etw

      common /sss/ brsss,clsss,xnasss

      if ((lday.eq.0.and.lst.eq.0.and.lmin.le.1).or.lmin/30*30.eq.lmin) &
           then
      cl_out=.false.
      if (nkc_l.gt.2) cl_out=.true.
! gas phase----------------------------
! atmosphere
      brg=0.
      clg=0.
      brgd=0.
      clgd=0.
      do k=2,n
         ! Brominated species
         brgp=0.
         do j=1,j1_br
            brgp=brgp+ind_gas_br(1,j)*s1(ind_gas_br(2,j),k)
         end do
         do j=1,j5_br
            brgp=brgp+ind_rad_br(1,j)*s3(ind_rad_br(2,j),k)
         end do
         brg=brg+brgp*detw(k) ! to get mol/m^2

         ! ditto for chlorinated species
         clgp=0.
         do j=1,j1_cl
            clgp=clgp+ind_gas_cl(1,j)*s1(ind_gas_cl(2,j),k)
         end do
         do j=1,j5_cl
            clgp=clgp+ind_rad_cl(1,j)*s3(ind_rad_cl(2,j),k)
         end do
         clg=clg+clgp*detw(k)
      enddo

!  + deposited ! already in mol/m^2
      do j=1,j1_br
         brgd=brgd+ind_gas_br(1,j)*s1(ind_gas_br(2,j),1)
      end do
      do j=1,j5_br
         brgd=brgd+ind_rad_br(1,j)*s3(ind_rad_br(2,j),1)
      end do

      do j=1,j1_cl
         clgd=clgd+ind_gas_cl(1,j)*s1(ind_gas_cl(2,j),1)
      end do
      do j=1,j5_cl
         clgd=clgd+ind_rad_cl(1,j)*s3(ind_rad_cl(2,j),1)
      end do

! aqueous phase-------------------------
! atmosphere
      bra1=0.
      bra2=0.
      bra3=0.
      bra4=0.
      cla1=0.
      cla2=0.
      cla3=0.
      cla4=0.
      xnaa2=0.
      bra1d=0.
      bra2d=0.
      bra3d=0.
      bra4d=0.
      cla1d=0.
      cla2d=0.
      cla3d=0.
      cla4d=0.
      xnaa2d=0.
      do k=2,n
         bra1=bra1+(sl1(42,1,k)+sl1(43,1,k)+sl1(44,1,k)+sl1(45,1,k)+ &
              2.*sl1(46,1,k)+sl1(47,1,k)+sl1(55,1,k)+sl1(j2-j3+9,1,k) &
              +sion1(24,1,k)+2.*sion1(25,1,k)+sion1(26,1,k)+ &
          sion1(27,1,k)+sion1(28,1,k)+2.*sion1(29,1,k)+2.*sion1(38,1,k)) &
              *detw(k)
         bra2=bra2+(sl1(42,2,k)+sl1(43,2,k)+sl1(44,2,k)+sl1(45,2,k)+ &
              2.*sl1(46,2,k)+sl1(47,2,k)+sl1(55,2,k)+sl1(j2-j3+9,2,k) &
              +sion1(24,2,k)+2.*sion1(25,2,k)+sion1(26,2,k)+ &
          sion1(27,2,k)+sion1(28,2,k)+2.*sion1(29,2,k)+2.*sion1(38,2,k)) &
              *detw(k)
         if (cl_out) then
         bra3=bra3+(sl1(42,3,k)+sl1(43,3,k)+sl1(44,3,k)+sl1(45,3,k)+ &
              2.*sl1(46,3,k)+sl1(47,3,k)+sl1(55,3,k)+sl1(j2-j3+9,3,k) &
              +sion1(24,3,k)+2.*sion1(25,3,k)+sion1(26,3,k)+ &
          sion1(27,3,k)+sion1(28,3,k)+2.*sion1(29,3,k)+2.*sion1(38,3,k)) &
              *detw(k)
         bra4=bra4+(sl1(42,4,k)+sl1(43,4,k)+sl1(44,4,k)+sl1(45,4,k)+ &
              2.*sl1(46,4,k)+sl1(47,4,k)+sl1(55,4,k)+sl1(j2-j3+9,4,k) &
              +sion1(24,4,k)+2.*sion1(25,4,k)+sion1(26,4,k)+ &
          sion1(27,4,k)+sion1(28,4,k)+2.*sion1(29,4,k)+2.*sion1(38,4,k)) &
              *detw(k)
         endif

         cla1=cla1+(sl1(30,1,k)+sl1(38,1,k)+sl1(39,1,k)+sl1(40,1,k)+ &
             2.*sl1(41,1,k)+sl1(47,1,k)+sl1(54,1,k)+sl1(58,1,k) &
             +sl1(j2-j3+8,1,k) &
             +sion1(14,1,k)+2.*sion1(15,1,k)+sion1(22,1,k)+sion1(23,1,k) &
             +2.*sion1(28,1,k)+sion1(29,1,k)+2.*sion1(37,1,k))*detw(k)
         cla2=cla2+(sl1(30,2,k)+sl1(38,2,k)+sl1(39,2,k)+sl1(40,2,k)+ &
             2.*sl1(41,2,k)+sl1(47,2,k)+sl1(54,2,k)+sl1(58,2,k) &
             +sl1(j2-j3+8,2,k) &
             +sion1(14,2,k)+2.*sion1(15,2,k)+sion1(22,2,k)+sion1(23,2,k) &
             +2.*sion1(28,2,k)+sion1(29,2,k)+2.*sion1(37,2,k))*detw(k)
        if (cl_out) then
         cla3=cla3+(sl1(30,3,k)+sl1(38,3,k)+sl1(39,3,k)+sl1(40,3,k)+ &
             2.*sl1(41,3,k)+sl1(47,3,k)+sl1(54,3,k)+sl1(58,3,k) &
             +sl1(j2-j3+8,3,k) &
             +sion1(14,3,k)+2.*sion1(15,3,k)+sion1(22,3,k)+sion1(23,3,k) &
             +2.*sion1(28,3,k)+sion1(29,3,k)+2.*sion1(37,3,k))*detw(k)
         cla4=cla4+(sl1(30,4,k)+sl1(38,4,k)+sl1(39,4,k)+sl1(40,4,k)+ &
             2.*sl1(41,4,k)+sl1(47,4,k)+sl1(54,4,k)+sl1(58,4,k) &
             +sl1(j2-j3+8,4,k) &
             +sion1(14,4,k)+2.*sion1(15,4,k)+sion1(22,4,k)+sion1(23,4,k) &
             +2.*sion1(28,4,k)+sion1(29,4,k)+2.*sion1(37,4,k))*detw(k)
        endif

         xnaa2=xnaa2+sion1(20,2,k)*detw(k)
      enddo

! + deposited
      bra1d=sl1(42,1,1)+sl1(43,1,1)+sl1(44,1,1)+sl1(45,1,1)+ &
           2.*sl1(46,1,1)+sl1(47,1,1)+sl1(55,1,1)+sl1(j2-j3+9,1,1) &
           +sion1(24,1,1)+2.*sion1(25,1,1)+sion1(26,1,1)+ &
           sion1(27,1,1)+sion1(28,1,1)+2.*sion1(29,1,1)+ &
           2.*sion1(38,1,1)
      bra2d=sl1(42,2,1)+sl1(43,2,1)+sl1(44,2,1)+sl1(45,2,1)+ &
           2.*sl1(46,2,1)+sl1(47,2,1)+sl1(55,2,1)+sl1(j2-j3+9,2,1) &
           +sion1(24,2,1)+2.*sion1(25,2,1)+sion1(26,2,1)+ &
           sion1(27,2,1)+sion1(28,2,1)+2.*sion1(29,2,1)+ &
           2.*sion1(38,2,1)
       if (cl_out) then
      bra3d=sl1(42,3,1)+sl1(43,3,1)+sl1(44,3,1)+sl1(45,3,1)+ &
           2.*sl1(46,3,1)+sl1(47,3,1)+sl1(55,3,1)+sl1(j2-j3+9,3,1) &
           +sion1(24,3,1)+2.*sion1(25,3,1)+sion1(26,3,1)+ &
           sion1(27,3,1)+sion1(28,3,1)+2.*sion1(29,3,1)+ &
           2.*sion1(38,3,1)
      bra4d=sl1(42,4,1)+sl1(43,4,1)+sl1(44,4,1)+sl1(45,4,1)+ &
           2.*sl1(46,4,1)+sl1(47,4,1)+sl1(55,4,1)+sl1(j2-j3+9,4,1) &
           +sion1(24,4,1)+2.*sion1(25,4,1)+sion1(26,4,1)+ &
           sion1(27,4,1)+sion1(28,4,1)+2.*sion1(29,4,1)+ &
           2.*sion1(38,4,1)
        endif

      cla1d=sl1(30,1,1)+sl1(38,1,1)+sl1(39,1,1)+sl1(40,1,1)+ &
           2.*sl1(41,1,1)+sl1(47,1,1)+sl1(54,1,1)+sl1(58,1,1) &
           +sl1(j2-j3+8,1,1) &
           +sion1(14,1,1)+2.*sion1(15,1,1)+sion1(22,1,1)+sion1(23,1,1) &
           +2.*sion1(28,1,1)+sion1(29,1,1)+2.*sion1(37,1,1)
      cla2d=sl1(30,2,1)+sl1(38,2,1)+sl1(39,2,1)+sl1(40,2,1)+ &
           2.*sl1(41,2,1)+sl1(47,2,1)+sl1(54,2,1)+sl1(58,2,1) &
           +sl1(j2-j3+8,2,1) &
           +sion1(14,2,1)+2.*sion1(15,2,1)+sion1(22,2,1)+sion1(23,2,1) &
           +2.*sion1(28,2,1)+sion1(29,2,1)+2.*sion1(37,2,1)
       if (cl_out) then
      cla3d=sl1(30,3,1)+sl1(38,3,1)+sl1(39,3,1)+sl1(40,3,1)+ &
           2.*sl1(41,3,1)+sl1(47,3,1)+sl1(54,3,1)+sl1(58,3,1) &
           +sl1(j2-j3+8,3,1) &
           +sion1(14,3,1)+2.*sion1(15,3,1)+sion1(22,3,1)+sion1(23,3,1) &
           +2.*sion1(28,3,1)+sion1(29,3,1)+2.*sion1(37,3,1)
      cla4d=sl1(30,4,1)+sl1(38,4,1)+sl1(39,4,1)+sl1(40,4,1)+ &
           2.*sl1(41,4,1)+sl1(47,4,1)+sl1(54,4,1)+sl1(58,4,1) &
           +sl1(j2-j3+8,4,1) &
           +sion1(14,4,1)+2.*sion1(15,4,1)+sion1(22,4,1)+sion1(23,4,1) &
           +2.*sion1(28,4,1)+sion1(29,4,1)+2.*sion1(37,4,1)
      endif

      xnaa2d=sion1(20,1,1)+sion1(20,2,1)+sion1(20,3,1)+sion1(20,4,1)

! output


 100  continue
      open (74,file='mass.out',status='unknown', position='append' &
           ,err=100)
      write (74,10) lday,lst,lmin
       if (cl_out) then
        write (74,20) brg, brgd, bra1, bra2, bra1d,bra2d, brsss*detw(1), &
           brg+bra1+bra2+brgd+bra1d+bra2d-brsss*detw(1)
        write (74,20) brg, brgd, bra3, bra4,bra3d, bra4d, brsss*detw(1), &
           brg+bra3+bra4+brgd+bra3d+bra4d-brsss*detw(1),brg+bra1+bra2+ &
           brgd+bra1d+bra2d+bra3+bra4+brgd+bra3d+bra4d-brsss*detw(1)
        write (74,25) clg, clgd, cla1, cla2,cla1d, cla2d, clsss*detw(1), &
           clg+cla1+cla2+clgd+cla1d+cla2d-clsss*detw(1)
        write (74,25) clg, clgd, cla3, cla4,cla3d, cla4d, clsss*detw(1), &
           clg+cla3+cla4+clgd+cla3d+cla4d-clsss*detw(1),clg+cla1+cla2+ &
           clgd+cla1d+cla2d+cla3+cla4+clgd+cla3d+cla4d-clsss*detw(1)
        write (74,30) xnaa2, xnaa2d,xnasss*detw(1),xnaa2+xnaa2d-xnasss* &
           detw(1)
       else
        write (74,20) brg, brgd, bra1, bra2, bra1d,bra2d, brsss*detw(1), &
           brg+bra1+bra2+brgd+bra1d+bra2d-brsss*detw(1)
        write (74,25) clg, clgd, cla1, cla2, cla1d,cla2d, clsss*detw(1), &
           clg+cla1+cla2+clgd+cla1d+cla2d-clsss*detw(1)
        write (74,30) xnaa2, xnaa2d,  xnasss*detw(1), xnaa2+xnaa2d &
           -xnasss*detw(1)
       endif
      close (74)

 10   format(3i3)
 20   format('br: ',8d15.7)
 25   format('cl: ',8d15.7)
 30   format('na: ',45x,d15.7,15x,3d15.7)

      endif

      end subroutine mass_ch



!
!----------------------------------------------------------------
!


      subroutine ave_parms (n_bl,nz_box)
! average cw and rc over BL in box model if BL_box=.true.
!      implicit double precision (a-h,o-z)

      USE global_params, ONLY : &
! Imported Parameters:
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit double precision (a-h,o-z)

      include 'gas_Parameters.h' !additional common blocks and other definitions

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      real(kind=dp) :: cw, cm
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_dryg/ xkmtd(n,2,NSPEC),henry(n,NSPEC),xeq(n,NSPEC)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb unused now

! note: activity coeff. are being correctly calculated with the averaged
! concentrations in layer n_bl

! start averaging one level above "working level"
! - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3

      do k=nstart,nz_box
         cwsum  = 0.
         rcsum = 0.
         cmsum  = 0.
!         cwmsum = 0.
!         acmsum = 0. ! jjb removed
!         convsum= 0. ! jjb unused below, now removed
         do kc=1,nkc
            cwsum  = cwsum  + cw(kc,k)
            rcsum = rcsum + rc(kc,k)
            cmsum  = cmsum  + cm(kc,k)
!            acmsum = acmsum + acm(k,kc) ! jjb removed
!            convsum= convsum+ conv2(k,kc) ! jjb unused
         enddo
         cw(kc,n_bl)  = cwsum / (nz_box-nstart+1)
         rc(kc,n_bl) = rcsum / (nz_box-nstart+1)
         cm(kc,n_bl)  = cmsum / (nz_box-nstart+1)
!         cwm(kc,n_bl) = cwmsum / (nz_box-nstart+1)
!         acm(n_bl,kc) = acmsum / (nz_box-nstart+1) ! jjb removed
      enddo


      am3s = 0.
!      am32s = 0.
      cm3s = 0.
!      cm32s = 0.
      ps    = 0.
      rhos  = 0.
      ts    = 0.

      do k=nstart,nz_box
         am3s = am3s + am3(k)
!         am32s = am32s + am3(k,2)
         cm3s = cm3s + cm3(k)
!         cm32s = cm32s + cm3(k,2)
         ps    = ps    + p(k)
         rhos  = rhos  + rho(k)
         ts    = ts    + t(k)
      enddo
      am3(n_bl) = am3s / (nz_box-nstart+1)
!      am3(n_bl,2) = am32s / (nz_box-nstart+1)
      cm3(n_bl) = cm3s / (nz_box-nstart+1)
!      cm3(n_bl,2) = cm32s / (nz_box-nstart+1)

      p(n_bl)     = ps   / (nz_box-nstart+1)
      rho(n_bl)   = rhos / (nz_box-nstart+1)
      t(n_bl)     = ts   / (nz_box-nstart+1)


      do j=1,NSPEC
         do kc=1,2
            xkmtds = 0.
            do k=nstart,nz_box
               xkmtds = xkmtds + xkmtd(k,kc,j)
            enddo
            xkmtd(n_bl,kc,j) = xkmtds / (nz_box-nstart+1)
         enddo
      enddo

      end subroutine ave_parms


!
!-------------------------------------------------------
!

      subroutine ave_j (nz_box,n_bl)
! calculate average photolysis rates for BL and put these values on level n_bl

      USE global_params, ONLY : &
! Imported Parameters:
           n, &
           nphrxn

      implicit double precision (a-h,o-z)

      common /band_rat/ photol_j(nphrxn,n)
      double precision photol_j

! start averaging one level above "working level"
! - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3

      do j=1,nphrxn ! jjb
         xjsum = 0.
         do k=nstart,nz_box
            xjsum = xjsum + photol_j(j,k)
         enddo
         photol_j(j,n_bl) = xjsum/(nz_box-nstart+1)
      enddo

      end subroutine ave_j

!
!----------------------------------------------------------------
!

      subroutine ave_aer (n_bl,nz_box)
! average k_mt over BL in box model if BL_box=.true.

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           nkc

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions

      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
           xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
!     dimension fs(n,nka),ffsum(nka) ! jjb both not used here

! start averaging one level above "working level"
! - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3


      do j=1,NSPEC
         do kc=1,2
            xkmsum = 0.
            xkfsum = 0.
            xkbsum = 0.
            do k=nstart,nz_box
               xkmsum = xkmsum + xkmt(k,kc,j)
               xkfsum = xkfsum + xkef(k,kc,j)
               xkbsum = xkbsum + xkeb(k,kc,j)
            enddo
            xkmt(n_bl,kc,j) = xkmsum / (nz_box-nstart+1)
            xkef(n_bl,kc,j) = xkfsum / (nz_box-nstart+1)
            xkeb(n_bl,kc,j) = xkbsum / (nz_box-nstart+1)
         enddo
      enddo

      do j=1,NSPEC
         xhensum = 0.
         xalsum  = 0.
         xvmsum  = 0.
         do k=nstart,nz_box
            xhensum = xhensum + henry(j,k)
            xalsum  = xalsum  + alpha(j,k)
            xvmsum  = xvmsum  + vmean(j,k)
         enddo
         henry(j,n_bl)  = xhensum / (nz_box-nstart+1)
         alpha(j,n_bl)  = xalsum  / (nz_box-nstart+1)
         vmean(j,n_bl)  = xvmsum  / (nz_box-nstart+1)
      enddo

! particle size distribution
! some logical mistake - leeds to very low LWC after next call to SR cw_rc
!      do ia=1,nka
!         ffsum(ia) = 0.
!      enddo
!! "dry" the aerosol and sum up
!      do k=nstart,nz_box
!         do ia=1,nka
!            do jt=1,nkt
!               fs(k,ia)=fs(k,ia)+ff(jt,ia,k)
!            enddo
!            ffsum(ia)=ffsum(ia)+fs(k,ia)
!         enddo
!      enddo
!! average
!      do ia=1,nka
!         ff(1,ia,n_bl)=ffsum(ia)  / (nz_box-nstart+1)
!         do jt=2,nkt
!            ff(jt,ia,n_bl)=0.
!         enddo
!      enddo

      end subroutine ave_aer

!
!----------------------------------------------------------------
!

      subroutine ave_tot (n_bl,nz_box)
! average k_mt over BL in box model if BL_box=.true.

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h' !additional common blocks and other definitions

      common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
           xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      real(kind=dp) :: henry, xkmt, xkef, xkeb
      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)

! start averaging one level above "working level"
! - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3


      do j=1,NSPEC
         do kc=1,nkc
            xkmsum = 0.
            xkfsum = 0.
            xkbsum = 0.
            do k=nstart,nz_box
               xkmsum = xkmsum + xkmt(k,kc,j)
               xkfsum = xkfsum + xkef(k,kc,j)
               xkbsum = xkbsum + xkeb(k,kc,j)
            enddo
            xkmt(n_bl,kc,j) = xkmsum / (nz_box-nstart+1)
            xkef(n_bl,kc,j) = xkfsum / (nz_box-nstart+1)
            xkeb(n_bl,kc,j) = xkbsum / (nz_box-nstart+1)
         enddo
      enddo

      do j=1,NSPEC
         xhensum = 0.
         xalsum  = 0.
         xvmsum  = 0.
         do k=nstart,nz_box
            xhensum = xhensum + henry(j,k)
            xalsum  = xalsum  + alpha(j,k)
            xvmsum  = xvmsum  + vmean(j,k)
         enddo
         henry(j,n_bl)  = xhensum / (nz_box-nstart+1)
         alpha(j,n_bl)  = xalsum  / (nz_box-nstart+1)
         vmean(j,n_bl)  = xvmsum  / (nz_box-nstart+1)
      enddo

      end subroutine ave_tot

!
!----------------------------------------------------------------
!

      subroutine print_k_mt_a


      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           nkc

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions

      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
           xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)

      print *,'print aerosol kmt'
      do k=1,nf
         print *,k,xkmt(k,1,ind_O3),xkmt(k,2,ind_O3)
      enddo

      end subroutine print_k_mt_a



!
!-------------------------------------------------------
!

      subroutine set_box_gas (nlevbox,n_bl)
!     pick the values from the designated level: nlevbox

      USE global_params, ONLY : &
! Imported Parameters:
           n

      implicit double precision (a-h,o-z)

      include 'gas_Parameters.h' !additional common blocks and other definitions

      common /kpp_dryg/ xkmtd(n,2,NSPEC),henry(n,NSPEC),xeq(n,NSPEC)

      do kc=1,2!nkc
         do j=1,NSPEC
            xkmtd(n_bl,kc,j) = xkmtd(nlevbox,kc,j)
         enddo
      enddo

      end subroutine set_box_gas

!
!-------------------------------------------------------
!

      subroutine set_box_lev_a (nlevbox,n_bl)
!     pick the values from the designated level: nlevbox

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           n, &
           nka, &
           nkt, &
           nkc, &
           mb

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions

      common /cb11/ totrad (mb,n)
      double precision totrad

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      real (kind=dp) :: ff, fsum
      integer :: nar

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      real(kind=dp) :: cw, cm
      common /blck13/ conv2(nkc,n)
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
           xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_drya/ xkmtd(nf,2,NSPEC),xeq(nf,NSPEC)

      do kc=1,2!nkc
!         acm(n_bl,kc)  = acm(nlevbox,kc) ! jjb removed
!        conv2(n_bl,kc)= conv2(nlevbox,kc)! jjb updated
         conv2(kc,n_bl) = conv2(kc,nlevbox)
         cw(kc,n_bl)   = cw(kc,nlevbox)
         rc(kc,n_bl)  = rc(kc,nlevbox)
         do j=1,NSPEC
            xkmt(n_bl,kc,j) = xkmt(nlevbox,kc,j)
            xkef(n_bl,kc,j) = xkef(nlevbox,kc,j)
            xkeb(n_bl,kc,j) = xkeb(nlevbox,kc,j)
            xkmtd(n_bl,kc,j) = xkmtd(nlevbox,kc,j)
         enddo
      enddo

      am3(n_bl) = am3(nlevbox)
!      am3(n_bl,2) = am3(nlevbox,2)
      cm3(n_bl) = cm3(nlevbox)
!      cm3(n_bl,2) = cm3(nlevbox,2)

      p(n_bl)     = p(nlevbox)
      rho(n_bl)   = rho(nlevbox)
      t(n_bl)     = t(nlevbox)
      totrad(1,n_bl) = totrad(1,nlevbox) ! jjb check this. Why is only the first wavelength band picked up?

      do j=1,NSPEC
         alpha(j,n_bl) = alpha(j,nlevbox)
         henry(j,n_bl) = henry(j,nlevbox)
         vmean(j,n_bl) = vmean(j,nlevbox)
      enddo

! SR cw_rc is still being called, therefore also init f
      do ia=1,nka
         do jt=1,nkt
            ff(jt,ia,n_bl) = ff(jt,ia,nlevbox)
         enddo
      enddo

      end subroutine set_box_lev_a

!
!-------------------------------------------------------
!


      subroutine set_box_lev_t (nlevbox,n_bl)
!     pick the values from the designated level: nlevbox

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           n, &
           nkc, &
           mb

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h' !additional common blocks and other definitions

      common /cb11/ totrad (mb,n)
      double precision totrad

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      real(kind=dp) :: cw, cm
      common /blck13/ conv2(nkc,n)
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
           xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      real(kind=dp) :: henry, xkmt, xkef, xkeb
      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated

      xph3=0.
      xph4=0.
      if (cm(n_bl,3).gt.0.) xph3 = 1.
      if (cm(n_bl,4).gt.0.) xph4 = 1.

      if (xph3.eq.1..and.xph4.eq.1.) return

      do kc=1,nkc
!         acm(n_bl,kc)  = acm(nlevbox,kc) ! jjb removed
!        conv2(n_bl,kc)= conv2(nlevbox,kc) ! jjb updated
         conv2(kc,n_bl) = conv2(kc,nlevbox)
         cw(kc,n_bl) = cw(kc,nlevbox)
         rc(kc,n_bl) = rc(kc,nlevbox)
         do j=1,NSPEC
            xkmt(n_bl,kc,j) = xkmt(nlevbox,kc,j)
            xkef(n_bl,kc,j) = xkef(nlevbox,kc,j)
            xkeb(n_bl,kc,j) = xkeb(nlevbox,kc,j)
         enddo
      enddo

      am3(n_bl) = am3(nlevbox)
!      am3(n_bl,2) = am3(nlevbox,2)
      cm3(n_bl) = cm3(nlevbox)
!      cm3(n_bl,2) = cm3(nlevbox,2)

      p(n_bl)     = p(nlevbox)
      rho(n_bl)   = rho(nlevbox)
      t(n_bl)     = t(nlevbox)
      totrad(1,n_bl) = totrad(1,nlevbox) ! jjb check this. Why is only the first wavelength band picked up?

      do j=1,NSPEC
         alpha(j,n_bl) = alpha(j,nlevbox)
         henry(j,n_bl) = henry(j,nlevbox)
         vmean(j,n_bl) = vmean(j,nlevbox)
      enddo

      end subroutine set_box_lev_t


!
!-------------------------------------------------------
!

      subroutine print_vals (nlevbox,n_bl)
!     test output

      USE global_params, ONLY : &
! Imported Parameters:
           nf, &
           n, &
           nkc, &
           mb

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions

      common /cb11/ totrad (mb,n)
      double precision totrad

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      real(kind=dp) :: cw, cm
      common /blck13/ conv2(nkc,n)
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC), &
           xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)


      print *,'output info', lst,lmin

      k=n_bl

      print *,'k = ',k
      do kc=1,2!nkc
         print *,'kc = ',kc
!        print *,acm(k,kc),conv2(k,kc),cw(k,kc) ! jjb acm removed
         print *,conv2(kc,k),cw(kc,k)
         print *,rc(kc,k),cm(kc,k)
         do j=1,NSPEC
            print *,xkmt(k,kc,j),xkef(k,kc,j),xkeb(k,kc,j)
         enddo
      enddo

      print *,'parameters'
!     print *,am3(k,1),am3(k,2),p(k) ! jjb am3 CO removed
      print *,am3(k),p(k)
!     print *,cm3(k,1),cm3(k,2),rho(k) ! jjb cm3 CO removed
      print *,cm3(k),rho(k)
      print *,t(k), totrad(1,k)

      print *,'alpha,henry,vmean'
      do j=1,NSPEC
         print *,alpha(j,k), henry(j,k), vmean(j,k)
      enddo

      k=nlevbox

      print *,'crap'

      print *,'k = ',k
      do kc=1,2!nkc
         print *,'kc = ',kc
!        print *,acm(k,kc),conv2(k,kc),cw(k,kc) ! jjb acm removed
         print *,conv2(kc,k),cw(kc,k)
         print *,rc(kc,k),cm(kc,k)
         do j=1,NSPEC
            print *,xkmt(k,kc,j),xkef(k,kc,j),xkeb(k,kc,j)
         enddo
      enddo

      print *,'parameters'
!     print *,am3(k,1),am3(k,2),p(k) ! jjb am3 CO removed
      print *,am3(k),p(k)
!     print *,cm3(k,1),cm3(k,2),rho(k) ! jjb cm3 CO removed
      print *,cm3(k),rho(k)
      print *,t(k), totrad(1,k)

      print *,'alpha,henry,vmean'
      do j=1,NSPEC
         print *,alpha(j,k), henry(j,k), vmean(j,k)
      enddo

      print *,' used for k=n_bl'
      k=n_bl
      do kc=1,2                 !nkc
         print *,k,kc,cw(kc,k)
         print *,xkmt(k,kc,ind_O3),xkmt(k,kc,ind_HCl), &
                 xkmt(k,kc,ind_HNO3)
      enddo

      print *,'now column vals'
      do k=1,nf
         do kc=1,2!nkc
            print *,k,kc,cw(kc,k)
            print *,xkmt(k,kc,ind_O3),xkmt(k,kc,ind_HCl), &
                 xkmt(k,kc,ind_HNO3)
         enddo
      enddo

      end subroutine print_vals

!
!-------------------------------------------------------
!
! jjb commented as unused 24/03/2016

!$$$      subroutine gamma_surf (box,n_bl)
!$$$c     calculation of reaction rate coefficients for surface reactions
!$$$c     note that the value that is used here for the accommodation
!$$$c     coefficient is in fact a reaction probability (=gamma)
!$$$c     after Knipping and Dabdub 2002 (#1692) and von Glasow (2006)
!$$$
!$$$! transfer coefficient after Schwarz, 1986 (see Sander & Crutzen '96, JGR, 9127)
!$$$! but no mean values used (like in SR k_mt_a/t) but integrated values
!$$$
!$$$      USE global_params, ONLY : &
!$$$! Imported Parameters:
!$$$           j2,
!$$$           j6,
!$$$           nf,
!$$$           n,
!$$$           nka,
!$$$           nkt,
!$$$           nkc
!$$$
!$$$      USE precision, ONLY : &
!$$$! Imported Parameters:
!$$$           dp
!$$$
!$$$      implicit double precision (a-h,o-z)
!$$$      logical box
!$$$
!$$$      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
!$$$                    e(nkt),dew(nkt),rq(nkt,nka)
!$$$      double precision enw,ew,rn,rw,en,e,dew,rq
!$$$      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
!$$$      real (kind=dp) :: ff, fsum
!$$$      integer :: nar
!$$$
!$$$      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
!$$$      real(kind=dp) :: theta, thetl, t, talt, p, rho
!$$$      common /blck06/ kw(nka),ka
!$$$      integer :: kw, ka
!$$$      common /blck12/ cw(nkc,n),cm(nkc,n)
!$$$      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
!$$$      real (kind=dp) :: sl1, sion1
!$$$!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
!$$$!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o &
!$$$!           (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
!$$$!     common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc) ! jjb none of the objects of the common block is used
!$$$!     real (kind=dp) :: vt, vd, vdm
!$$$!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated
!$$$      common /k_surf/ xkmt_OHClm(nf,nkc)
!$$$
!$$$      dimension freep(nf),rqm(nkt,nka),xkmt_surf(nf,nkc)
!$$$
!$$$      func(a,k)=dsqrt(t(k)/a)*4.60138
!$$$
!$$$! change rq in um to rqm in m
!$$$      do ia=1,nka
!$$$         do jt=1,nkt
!$$$            rqm(jt,ia)=rq(jt,ia)*1.d-6
!$$$         enddo
!$$$      enddo
!$$$
!$$$      nmin=2
!$$$      nmax=nf
!$$$      if (box) then
!$$$        nmin=n_bl
!$$$        nmax=n_bl
!$$$      endif
!$$$
!$$$! free path length (lambda=freep):
!$$$      do k=nmin,nmax
!$$$         freep(k)=2.28d-5 * t(k) / p(k)
!$$$      enddo
!$$$! loop over vertical grid
!$$$      do 1000 k=nmin,nmax
!$$$         vmean_OH= func(1.7d-2,k)
!$$$! loop over the nkc different chemical bins
!$$$         do kc=1,4
!$$$            if (cm(kc,k).eq.0.) goto 1001 ! switch changed from cw
!$$$! define summation limits (1) ---
!$$$            if (kc.eq.1) then
!$$$               iia_0=1
!$$$               iia_e=ka
!$$$            endif
!$$$            if (kc.eq.2) then
!$$$               iia_0=ka+1
!$$$               iia_e=nka
!$$$            endif
!$$$            if (kc.eq.3) then
!$$$               iia_0=1
!$$$               iia_e=ka
!$$$            endif
!$$$            if (kc.eq.4) then
!$$$               iia_0=ka+1
!$$$               iia_e=nka
!$$$            endif
!$$$! fast version without logarithmic integration
!$$$            x1=0.
!$$$            xk1=0.
!$$$! --- OH + Cl- --> 0.5 Cl2 + OH- ---
!$$$! "best guess": use alpha from Knipping and Dabdub and gas phase limitation:
!$$$! alpha = 0.02 * gamma_s * [Cl-], where [Cl-] is in mol/l, gamma_s=2
!$$$! C(ind_Clmlx) / LWC *1.d-3:  in mol/m3_air * m3_air/m3_aq * m3_aq/l_aq
!$$$            if (cw(kc,k).gt.0.) then
!$$$               gamma = min(1.,4.d-5*sion1(14,kc,k)/cw(kc,k))
!$$$            else
!$$$               gamma = 0.d0
!$$$            end if
!$$$            if (gamma.gt.0.) x1= 4./(3. * gamma)
!$$$! -------------------------------------------
!$$$            do ia=iia_0,iia_e
!$$$! define summation limits (2)
!$$$               if (kc.eq.1) then
!$$$                  jjt_0=1
!$$$                  jjt_e=kw(ia)
!$$$               endif
!$$$               if (kc.eq.2) then
!$$$                  jjt_0=1
!$$$                  jjt_e=kw(ia)
!$$$               endif
!$$$               if (kc.eq.3) then
!$$$                  jjt_0=kw(ia)+1
!$$$                  jjt_e=nkt
!$$$               endif
!$$$               if (kc.eq.4) then
!$$$                  jjt_0=kw(ia)+1
!$$$                  jjt_e=nkt
!$$$               endif
!$$$               do jt=jjt_0,jjt_e
!$$$! conversion: um      --> m               :10^-6
!$$$                  rqq=rqm(jt,ia)
!$$$! kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
!$$$!     with D_g=lambda*vmean/3.
!$$$! here a volume weighted value is calculated, therefore weighting with r^3:
!$$$! kmt=4/3*piL*/sum(a)*sum(r){r^3*N*kt}
!$$$! --- OH + Cl- --> 0.5 Cl2 + OH- ----
!$$$                  x2=vmean_OH/(rqq/freep(k)+x1) ![1/s]
!$$$! conversion: 1/cm^3 --> 1/m^3(air):10^6
!$$$                  xk1=xk1+x2*rqq*rqq*ff(jt,ia,k)*1.d6
!$$$               enddo            !jt
!$$$            enddo               !ia
!$$$! k_mt=4*pi/(3*LWC)*sum
!$$$            if (cw(kc,k).gt.0.d0) then
!$$$!               xkmt_surf(k,kc)=4.*3.1415927/(3.*cw(kc,k))*xk1 ![1/s]
!$$$               xkmt_surf(k,kc)=4.*pi/(3.*cw(kc,k))*xk1 ![1/s]
!$$$            end if
!$$$
!$$$! kmt is supposed to be a first-order rate coefficient but in KPP kmt will
!$$$! be multiplied by [Cl-]/[Br-], so divide by [X-] here:
!$$$
!$$$            if (sion1(14,kc,k).gt.0.) &
!$$$                 xkmt_OHClm(k,kc)=xkmt_surf(k,kc)/sion1(14,kc,k)
!$$$
!$$$!            print *,'SR gamma_surf'
!$$$!            print *,lday,lst,lmin
!$$$!            print *,k,kc,cw(k,kc)!,sion1(14,kc,k)
!$$$!            print *,xkmt_OHClm(k,kc)
!$$$! multiply with conc of X- to get "real" kmt:
!$$$!            print *,xkmt_OHClm(k,kc)*sion1(14,kc,k),xkmt_OHClm(k,kc)* &
!$$$!                 sion1(24,kc,k)
!$$$
!$$$ 1001       continue
!$$$         enddo                  ! kc
!$$$
!$$$ 1000 continue                  !k
!$$$
!$$$      end subroutine gamma_surf



!
!-------------------------------------------------------------
!

! functions for KPP

! third body reaction rate constants are formulated as pseudo
! second order: cm^3/s, so [aircc]=cm^3/s has to be used

      double precision function farr (a,b)
! Arrhenius relation

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a
      integer :: b

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      farr=a*exp(b/te)
      end function farr

!----------------------------------------------------------------

      double precision function farr_sp (a,b,c,d)
! Arrhenius relation (with complex temperature dependence)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a,c
      integer :: b,d

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      farr_sp=a*((te/b)**c)*exp(d/te)
      end function farr_sp

!----------------------------------------------------------------

      double precision function ATK_3 (a1,a2,b1,b2,fc)
! calculate third body reactions according to Atkinson '92

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,a2,b1,b2,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2
      x2=fc

      atk_3=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      end function ATK_3

!----------------------------------------------------------------

      double precision function ATK_3a (a1,a2,b1,b2,tfc)
! calculate third body reactions according to Atkinson '92

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,a2,b1,b2,tfc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2
      x2=exp(-te/tfc)
      atk_3a=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      end function ATK_3a

!----------------------------------------------------------------

      double precision function ATK_3c (a1,b1,fc)
! calculate third body reactions according to Atkinson '92

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,b1,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*exp(-10000./te)*aircc
      b0=b1*exp(-10900./te)
      x2=fc
      if (fc.eq.0.) x2=exp(-te/250.)+exp(-1050./te)
      atk_3c=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      end function ATK_3c

!----------------------------------------------------------------

      double precision function ATK_3d (a1,b1,fc)
! calculate third body reactions according to IUPAC 2004

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,b1,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*exp(-8000./te)*aircc
      b0=b1*exp(-8820./te)
      x2=fc
      atk_3d=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      end function ATK_3d

!----------------------------------------------------------------

      double precision function ATK_3e (a1,a2,b1,b2,fc)
! calculate third body reactions according to Atkinson '92 (used for OH + OIO)
!  NOT USED: Plane et al., 2006 reported a pressure dependent rate coefficient
!  from a theoretical study, however they suggested an Arrhenius relationship
!  for use in atmospheric modelling

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,a2,b1,b2,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2*exp(46./te)
      x2=fc

      atk_3e=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      end function ATK_3e

!----------------------------------------------------------------

      double precision function ATK_3f (a1,a2,b1,b2,fc)
! calculate third body reactions according to Atkinson '92 (used for OCLO + O3P)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,a2,b1,b2,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*aircc*(te/298.)**a2
      b0=b1*(te/298.)**b2
      x2=fc

      atk_3f=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      end function ATK_3f

!----------------------------------------------------------------

      double precision function sHNO3 (a1,b1,a2,b2,a3,b3)
! calculate special rate function for OH + HNO3 (JPL 2003)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,a1,a2,a3,func,tte
      integer :: b0,b1,b2,b3

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      func(a0,b0)=a0*exp(b0*tte)
      tte=1./te
      sHNO3=func(a1,b1)+(func(a3,b3)*aircc/ &
           (1+func(a3,b3)*aircc/func(a2,b2)))
      end function sHNO3

!----------------------------------------------------------------

      double precision function fbck (a1,a2,b1,b2,fc,ak,bk)
! calculate thermal decomposition rate from forward and
! equilibrium rate

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,a2,b1,b2,fc,x1,x2,ak,bk

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.d0)**a2
      b0=b1*(te/300.d0)**b2
      x2=fc

      x1=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      fbck=x1/(ak*dexp(bk/te))
      end function fbck

!----------------------------------------------------------------

      double precision function fbckJ (a1,a2,b1,b2,ak,bk)
! calculate thermal decomposition rate from forward and
! equilibrium rate

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,a2,b1,b2,x1,x2,ak,bk

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.d0)**a2
      b0=b1*(te/300.d0)**b2
      x2= 0.6d0

      x1=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      fbckJ=x1/(ak*dexp(bk/te))
      end function fbckJ

!----------------------------------------------------------------

      double precision function fbck2 (a1,a2,b1,b2,fc,ck)
! calculate thermal decomposition rate from forward and
! equilibrium rate (used for BrNO3 decomposition)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,b0,a1,a2,b1,b2,fc,x1,x2,ak,bk,ck

      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

! parameters to calculate K_eq in atm-1 (Orlando and Tyndall, 1996)
      ak=5.44d-9
      bk=14192.d0
      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2
      x2=fc

      x1=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)* &
           dlog10(a0/b0))))
      fbck2 = 0.d0
      if (ck.ne.0.d0) fbck2=x1/(ak*dexp(bk/te)*8.314/101325.*te/ck)
      end function fbck2


!----------------------------------------------------------------

      double precision function sp_17_old (a1)
! calculate special rate function for rxn 17 (OH + CO)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a1
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      sp_17_old=a1*(1+0.6*pk/101325.) !(pressure in atm)
      end function sp_17_old

!----------------------------------------------------------------

      double precision function sp_17 (a,b)
! calculate special rate function for rxn 17 (OH + CO)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a,b
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      sp_17=a*(1.d0+aircc/b) ! simpler IUPAC parametrization
!$$$      sp_17=a1*(1+0.6*pk/101325.) !(pressure in atm)
      end function sp_17

!----------------------------------------------------------------

      double precision function sp_23 (a1,b1,a2,b2,a3,b3)
! calculate special rate function for rxn 23 (HO2+HO2 - including H2O correction)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0,a1,a2,a3,func,tte
      integer :: b0,b1,b2,b3
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      func(a0,b0)=a0*exp(b0*tte)
      tte=1./te
      sp_23=(func(a1,b1)+func(a2*aircc,b2))* &
           (1+(func(a3*aircc*h2oppm*1.0d-6,b3)))
      end function sp_23

!----------------------------------------------------------------

      double precision function sp_29 (a1,b1,a2,b2,c)
! calculate special rate function for rxn 29

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a1,b1,a2,b2,c,fun1,fun2,fun3,num,den,z, &
       a0,b0,c0,d0
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      fun1(a0,b0,c0)=a0*(b0**c0)
      fun2(a0,b0)=1./(1.+(dlog10(a0/b0))*(dlog10(a0/b0)))
      fun3(a0,b0,c0,d0)=a0/(1.+a0/b0)*(c0**d0)

      num=aircc*fun1(a1,te,b1)
      den=    fun1(a2,te,b2)
      z=fun2(num,den)
      sp_29=fun3(num,den,c,z)
      end function sp_29

!----------------------------------------------------------------

      double precision function fcn (x1)
! rate constant for thermal decomposition of ClNO3 (#2730)

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: x1,x2,xmg
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      x2=8.314*te
      xmg=pk/x2
      fcn=10**(-6.16)*dexp(-90.7d3/x2)*xmg*x1
      end function fcn

!----------------------------------------------------------------

      double precision function farr2 (a0,b0)
! Arrhenius function but with b0 as the value for T=298K

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      real (kind=dp) :: a0
      integer :: b0
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk
! 1/298.=3.3557d-3
      farr2=a0*exp(dble(b0)*(1.d0/te-3.3557d-3))
      end function farr2

!----------------------------------------------------------------

      double precision function fhet_t (a0,b0,c0)
! heterogeneous rate function
! ClFCT     = 5.0D2                ; factor for H02/H01, i.e Cl-/H2O
! BrFCT     = 3.0D5                ; factor for H03/H01, i.e Br-/H2O
! a0=1..4  liquid size class
! b0=1..3  branch of het reaction: H2O, Cl-, Br-
! c0=1..3  gas phase reactant:     N2O5, ClNO3, BrNO3


      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      INCLUDE 'tot_Parameters.h'
      INCLUDE 'tot_Global.h'
      integer, intent(in) :: a0,b0,c0

      real (kind=dp) :: h2oa, hetT
      real (kind=dp) :: xbr, xtr

      if (a0.eq.1) then
         h2oa=FIX(indf_H2Ol1)
         hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)

      else if (a0.eq.2) then
         h2oa=FIX(indf_H2Ol2)
         hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)

      else if (a0.eq.3) then
         h2oa=FIX(indf_H2Ol3)
         hetT=h2oa + 5.0D2*C(ind_Clml3) + 3.0D5*C(ind_Brml3)

      else if (a0.eq.4) then
         h2oa=FIX(indf_H2Ol4)
         hetT=h2oa + 5.0D2*C(ind_Clml4) + 3.0D5*C(ind_Brml4)

      else ! undefined case, shouldn't happend
         stop 'Wrong a0 index in function fhet_t'
      endif


      if (b0.eq.1) then
         xbr=h2oa

      else if (b0.eq.2) then
         xbr=5.0D2

      else if (b0.eq.3) then
         xbr=3.0D5

      else ! undefined case, shouldn't happend
         stop 'Wrong b0 index in function fhet_t'
      endif


      if (c0.eq.1) then
         xtr=yxkmt(a0,ind_N2O5)
      else if (c0.eq.2) then
         xtr=yxkmt(a0,ind_ClNO3)
      else if (c0.eq.3) then
         xtr=yxkmt(a0,ind_BrNO3)
      else ! undefined case, shouldn't happend
         stop 'Wrong c0 index in function fhet_t'
      endif


      if (hetT.gt.0.d0) then
         fhet_t=xtr * ycw(a0) * xbr /hetT
      else
         fhet_t=0.
      endif

      end function fhet_t

!----------------------------------------------------------------

      double precision function fliq_60 (a1,b1,c,d)
! calculate special rate function for rxn 60

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a1,c,d
      integer :: b1
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

      if (d.gt.0.d0) then
!         fliq_60=farr2(a1,b1)*c/(c+0.1/d)
         fliq_60=a1*dexp(dble(b1)*(1.d0/te-3.3557d-3))*c/(c+0.1d0/d)
      else
         fliq_60=0.
      endif
      end

!----------------------------------------------------------------

      double precision function dmin2 (a)
! confine rate constant to upper limit (diffusion control)
! a=k; dclim=upper limit due to diffusion-control

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a,dclim

      dclim = 1.d10

      dmin2 = dmin1(a,dclim )

      end

!----------------------------------------------------------------

      double precision function dmin3 (a)
! confine rate constant to upper limit (diffusion control)
! a=k; dclim=upper limit due to diffusion-control
! factor 2.d0 is to account for larger upper limit for
!    2nd order reactions between differtly-charged ions

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a,dclim

      dclim = 1.d10

      dmin3 = dmin1(a,dclim*2.d0 )

      end

!----------------------------------------------------------------

      double precision function flsc (a,b,c,d)
! calculate special rate function
! after #s_364, Schmitz (1999), eq.(4) / #s_650, Schmitz (2000)
!: dio3/dt = k1*[IO3-][H+]^2[I-]^2 + k2*[IO3-][H+]^2[I-]
! a=k1, b=H+, c=I-, d=cvvz

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a,b,c,d

      if (d.gt.0.d0) then
         flsc=( a*b**2*d**4 + 1.2d3*b**2/c*d**3 )
      else
         flsc=0.
      endif
      end

!----------------------------------------------------------------

      double precision function flsc4 (a,b,c)
! calculate special rate function
! after #s_650, Schmitz (2000)
! a=k4, b=H+, c=cvvz

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a,b,c

      if (c.gt.0.d0) then
         flsc4=( a*b*c**3 )
      else
         flsc4=0.
      endif
      end

!----------------------------------------------------------------

      double precision function flsc5 (a,b,c)
! calculate special rate function
! after #s_650, Schmitz (2000)
! a=k5, b=H+, c=cvvz

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a,b,c

      if (c.gt.0.d0) then
        flsc5=( a*b**2*c**4 )
      else
        flsc5=0.
      endif
      end

!----------------------------------------------------------------

      function flsc6 (a,b)

! Description :
! -----------
!    calculate special rate function
!    after G. Schmitz (pers.comm.)
!    a=k6, b=H+

! Author :
! ------
!    Roland von Glasow

! Modifications :
! -------------
!  17-Oct-2016  Josue Bock  implicit none (u8.5)
!
!  05-Mar-2017  Josue Bock  introduced a lower limit for [H+] to prevent flsc6
!                           to reach very high values if [H+] is very small
!                           The chosen value might need to be further adjusted
!
!  19-Oct-2017  Josue Bock  - added a test to warn only if negative values (avoid warning
!                             messages to be displayed just after initialisation)
!                           - cosmetic improvements: header, precision from module, ...
! End modifications
!-----------------------------------------------------------------------------------------------------------------------


! Declarations:
! ------------
! Modules used:
      USE precision, ONLY: &
! Imported Type Definitions:
           dp                   ! kind double precision real

      implicit none

! Function result
      real(kind=dp) :: flsc6

! Function arguments
! Scalar arguments with intent(in):
      real(kind=dp), intent(in) :: a
      real(kind=dp), intent(in) :: b

! Local parameters:
      real(kind=dp), parameter :: min_Hp = 1.e-15_dp ! <jjb> introduce a lower limit for [H+] to avoid overflow flsc6

!- End of header ------------------------------------------------------------

!      if (b.gt.0.d0) then
      if (b.gt.min_Hp) then
         flsc6= a/b
      else
         flsc6=0._dp
         if (b.lt.0._dp) print*,"Warning: flsc6 encountered [H+]<0",b
      endif

      end function flsc6

!----------------------------------------------------------------

      double precision function uplim (a,b,c,d)

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a,b,c,d,dclim
! dclim = upper limit for diffusion-controlled reactions
! a=k-; b=k+; alpha=b/dclim; c=H+; d=cvvz;

      dclim = 1.d10
      if (d.gt.0.d0) then
         if (c<0.) print*,"Warning uplim encountered [H+]<0" ! Track unexpected case
!        uplim = ( a/(1.+b/dclim*c*d) )
        uplim = ( a/(1.+b/dclim*max(c,0.d0)*d) ) ! <jjb> avoid unexpected values if [H+]<0
      else
        uplim = 0.
      endif
      end

!----------------------------------------------------------------

      double precision function uparm (a0,b0,c,d,e)
! Arrhenius function but with b0 as the value for T=298K
! dclim = upper limit for diffusion-controlled reactions
! c=k+; alpha=c/dclim; d=[H+]; e=cvvz

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a0,c,d,e,dclim
      integer :: b0
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk
! 1/298.=3.3557d-3

      dclim = 1.d10
      if (d.gt.0.d0) then
        uparm=a0*exp(dble(b0)*(1.d0/te-3.3557d-3))/(1.+c/dclim*d*e)
      else
        uparm=0.
      endif
      end

!----------------------------------------------------------------

      double precision function uplip (a,b,c)

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a,b,c,dclim
! dclim = upper limit for diffusion-controlled 3rd order reactions
! with H+ as reactant
! a=k+; b=[H+]; c=cvvz; alpha=a/dclim

      dclim = 1.d10
      if (c.gt.0.d0) then
         if (b<0.) print*,"Warning uplip encountered [H+]<0" ! Track unexpected case
!         uplip = ( a/(1.+a/dclim*b*c)*c**2 )
         uplip = ( a/(1.+a/dclim*max(b,0.d0)*c)*c**2 ) ! <jjb> avoid unexpected values if [H+]<0
      else
        uplip = 0.
      endif
      end

!----------------------------------------------------------------

      double precision function uparp (a0,b0,c,d)
! Arrhenius function but with b0 as the value for T=298K
! dclim = upper limit for diffusion-controlled reactions
! c=H+; d=cvvz; alpha=f(a0,b0)/dclim

      USE precision, ONLY : &
! Imported Parameters:
           dp


      implicit none

      real (kind=dp) :: a0,c,d,dclim
      integer :: b0
      common /cb_1/ aircc,te,h2oppm,pk
      real (kind=dp) :: aircc,te,h2oppm,pk

! 1/298.=3.3557d-3

      dclim = 1.d10
      if (d.gt.0.d0) then
        uparp=a0*exp(dble(b0)*(1.d0/te-3.3557d-3))/ &
             (1.+(a0*exp(dble(b0)*(1.d0/te-3.3557d-3)))/dclim*c*d)*d**2
      else
        uparp=0.
      endif
      end

!----------------------------------------------------------------

! jjb: old function for heterogeneous rates, kept for legacy

      double precision function fdhet_a (a0,b0,c0)
! heterogeneous rate function

! a0=1..2  liquid size class
! b0=1     gas phase reactant:     HNO3
! c0=2..3  branch of het reaction: Cl-, Br-

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none

      INCLUDE 'aer_Parameters.h'
      INCLUDE 'aer_Global.h'

      integer, intent(in) :: a0,b0,c0
      real (kind=dp) :: hetT, xbr, xkt

! calculate het_total
      if (a0.eq.1) then
         hetT=C(ind_Clml1) + C(ind_Brml1)
!        branching ratio = 1, because f(X) is assumed to be 1
!        in KPP the rate is multiplied with [Cl-/Br-] so the branching
!        expression k1=k*cl-/hetT is correctly implemented
         if (c0.eq.2) xbr=1.
         if (c0.eq.3) xbr=1.
      else if (a0.eq.2) then
         hetT=C(ind_Clml2) + C(ind_Brml2)
!        branching ratio
         if (c0.eq.2) xbr=1.
         if (c0.eq.3) xbr=1.
      endif
! mass transfer coefficient
      if (b0.eq.1) xkt=yxkmtd(a0,ind_HNO3)

! kmt in (m^3_air/(m^3_aq*s)) therefore multiplication with LWC (m^3_aq/m^3_air)
! to get k in 1/s
      if (hetT.gt.0.d0) then
         fdhet_a=xkt * ycwd(a0) * xbr /hetT
      else
         fdhet_a=0.
      endif

      end function fdhet_a

!----------------------------------------------------------------

      double precision function fhet_da (xliq,xhet,a0,b0,c0)
! heterogeneous rate function
! ClFCT     = 5.0D2                ; factor for H02/H01, i.e Cl-/H2O
! BrFCT     = 3.0D5                ; factor for H03/H01, i.e Br-/H2O
! a0=1..2  liquid size class
! b0=1..3  branch of het reaction: H2O, Cl-, Br-
! c0=1..3  gas phase reactant:     N2O5, ClNO3, BrNO3

      implicit double precision (a-h,o-z)

      INCLUDE 'aer_Parameters.h'
      INCLUDE 'aer_Global.h'
      integer a0,b0,c0


      if (xhet.eq.0.) then
         if (c0.eq.1) xtr=yxkmt(a0,ind_N2O5)
         if (c0.eq.2) xtr=yxkmt(a0,ind_ClNO3)
         if (c0.eq.3) xtr=yxkmt(a0,ind_BrNO3)
         if (a0.eq.1) then
            h2oa=FIX(indf_H2Ol1)
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycw(a0)
         else if (a0.eq.2) then
            h2oa=FIX(indf_H2Ol2)
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycw(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=FIX(indf_H2Ol1)
            if (a0.eq.2) hetT=FIX(indf_H2Ol2)
         endif
      else
!        if (c0.eq.2) xtr=yxkmtd(a0,ind_N2O5)  ! jjb wrong index
!        if (c0.eq.3) xtr=yxkmtd(a0,ind_BrNO3) ! jjb wrong index
!        if (c0.eq.4) xtr=yxkmtd(a0,ind_ClNO3) ! jjb wrong index
         if (c0.eq.1) xtr=yxkmtd(a0,ind_N2O5)  ! jjb corrected
         if (c0.eq.2) xtr=yxkmtd(a0,ind_BrNO3) ! jjb corrected
         if (c0.eq.3) xtr=yxkmtd(a0,ind_ClNO3) ! jjb corrected
!         print*,xliq,a0,c0,xtr
         if (a0.eq.1) then
            h2oa=55.55*ycwd(1)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycwd(a0)
         else if (a0.eq.2) then
            h2oa=55.55*ycwd(2)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycwd(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=55.55*ycwd(1)*1.d+3
            if (a0.eq.2) hetT=55.55*ycwd(2)*1.d+3
         endif
      endif

      if (b0.eq.1) xbr=h2oa
      if (b0.eq.2) xbr=5.0D2
      if (b0.eq.3) xbr=3.0D5

      if (hetT.gt.0.d0) then
         fhet_da=xtr * yw * xbr /hetT
      else
         fhet_da=0.
      endif
      if ((c0.eq.2.or.c0.eq.3.or.b0.eq.2.or.b0.eq.3).and.xhal.eq.0.) &
           fhet_da=0.
      if (xliq.eq.0.) fhet_da=0.
!      print*,xliq,a0,c0,fhet_da
      end function fhet_da

!----------------------------------------------------------

      double precision function fhet_dt (xliq,xhet,a0,b0,c0)
! heterogeneous rate function
! ClFCT     = 5.0D2                ; factor for H02/H01, i.e Cl-/H2O
! BrFCT     = 3.0D5                ; factor for H03/H01, i.e Br-/H2O
! a0=1..2  liquid size class
! b0=1..3  branch of het reaction: H2O, Cl-, Br-
! c0=1..3  gas phase reactant:     N2O5, ClNO3, BrNO3

      implicit double precision (a-h,o-z)

      INCLUDE 'tot_Parameters.h'
      INCLUDE 'tot_Global.h'
      integer a0,b0,c0


      if (xhet.eq.0.) then
         if (c0.eq.1) xtr=yxkmt(a0,ind_N2O5)
         if (c0.eq.2) xtr=yxkmt(a0,ind_ClNO3)
         if (c0.eq.3) xtr=yxkmt(a0,ind_BrNO3)
         if (a0.eq.1) then
            h2oa=FIX(indf_H2Ol1)
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycw(a0)
         endif
         if (a0.eq.2) then
            h2oa=FIX(indf_H2Ol2)
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycw(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=FIX(indf_H2Ol1)
            if (a0.eq.2) hetT=FIX(indf_H2Ol2)
         endif
      else
!        if (c0.eq.2) xtr=yxkmtd(a0,ind_N2O5)  ! jjb wrong index
!        if (c0.eq.3) xtr=yxkmtd(a0,ind_BrNO3) ! jjb wrong index
!        if (c0.eq.4) xtr=yxkmtd(a0,ind_ClNO3) ! jjb wrong index
         if (c0.eq.1) xtr=yxkmtd(a0,ind_N2O5)  ! jjb corrected
         if (c0.eq.2) xtr=yxkmtd(a0,ind_BrNO3) ! jjb corrected
         if (c0.eq.3) xtr=yxkmtd(a0,ind_ClNO3) ! jjb corrected
         if (a0.eq.1) then
            h2oa=55.55*ycwd(1)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycwd(a0)
         endif
         if (a0.eq.2) then
            h2oa=55.55*ycwd(2)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycwd(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=55.55*ycwd(1)*1.d+3
            if (a0.eq.2) hetT=55.55*ycwd(2)*1.d+3
         endif
      endif

      if (b0.eq.1) xbr=h2oa
      if (b0.eq.2) xbr=5.0D2
      if (b0.eq.3) xbr=3.0D5

      if (hetT.gt.0.d0) then
         fhet_dt=xtr * yw * xbr /hetT
      else
         fhet_dt=0.
      endif
      if ((c0.eq.2.or.c0.eq.3.or.b0.eq.2.or.b0.eq.3).and.xhal.eq.0.) &
           fhet_dt=0.
      if (xliq.eq.0.) fhet_dt=0.

      end function fhet_dt

!----------------------------------------------------------

      double precision function fdhetg (na,nb)
! heterogeneous rate function
! a0=1..2  liquid size class

! na is bin number
! nb is reaction:
!   1 HNO3
!   2 N2O5
!   3 NH3
!   4 H2SO4

      implicit double precision (a-h,o-z)

      INCLUDE 'gas_Parameters.h'
      INCLUDE 'gas_Global.h'
      integer na,nb
! net mass transfer coefficient including Henry's law equilibrium for HNO3

      if (nb.eq.1) then
! not limited by Henry's law:
!         xkt=yxkmtd(na,ind_HNO3) * ycwd(na)

! limited by Henry's law:
! see Diss RvG (3.10): dcg/dt = ... + kmt(LWC*Cg - Ca/H), this term is implemented here
! note that Cg is multiplied to rate in KPP, therefore the heterogeneous reaction rate is:
!     kmt(LWC - Ca/(Cg*H)) = x1 + x2
!     in x2 the "aqueous" concentration of HNO3 on the dry aerosol is calculated via
!     Henry's law and a HARDCODED particle pH = 2

         x1 = yxkmtd(na,ind_HNO3) * ycwd(na)
! index out of bounds in C(ind_NO3mlz) as this is not known in gas_Parameters.h
!         if (na.eq.1) caq=((C(ind_HNO3l1)+C(ind_NO3ml1))*1.d-2)/ &
!              (yxeq(ind_HNO3) + 1.d-2)
!         if (na.eq.2) caq=((C(ind_HNO3l2)+C(ind_NO3ml2))*1.d-2)/ &
!              (yxeq(ind_HNO3) + 1.d-2)
! if pH=2, the fractionation between HNO3 and NO3- can be calculated with equil. const:
! [NO3-]=Kq [HNO3] 1/[H+] = 1500 [HNO3] at pH=2
! this is all VERY rough and should be replaced!!
         if (na.eq.1) caq=((C(ind_HNO3l1)*1.5d3)*1.d-2)/ &
              (yxeq(ind_HNO3) + 1.d-2)
         if (na.eq.2) caq=((C(ind_HNO3l2)*1.5d3)*1.d-2)/ &
              (yxeq(ind_HNO3) + 1.d-2)
         x2 = 0.d0
          if (C(ind_HNO3).ne.0.d0.and.yhenry(ind_HNO3).ne.0.d0) &
              x2=-yxkmtd(na,ind_HNO3)/(C(ind_HNO3)*yhenry(ind_HNO3))*caq
         xkt = max(0.d0,(x1 + x2))
      endif
! kmt in (m^3_air/(m^3_aq*s)) therefore multiplication with LWC (m^3_aq/m^3_air)
! to get k in 1/s:
      if (nb.eq.2) xkt=yxkmtd(na,ind_N2O5) * ycwd(na)
      if (nb.eq.3) xkt=yxkmtd(na,ind_NH3) * ycwd(na)
      if (nb.eq.4) xkt=yxkmtd(na,ind_H2SO4) * ycwd(na)


      fdhetg=xkt

!      print *,k,na,nb
!      print *,fdhetg,xkt,ycwd(na)
!      if (nb.eq.1) write (440, 1001) k,na,nb,fdhetg,yxkmtd(na,ind_HNO3) &
!           * ycwd(na)
! 1001 format(3i4, 6d16.8)

      end function fdhetg

!----------------------------------------------------------------

      double precision function fdheta (na,nb)
! heterogeneous rate function
! a0=1..2  liquid size class

      implicit double precision (a-h,o-z)

      INCLUDE 'aer_Parameters.h'
      INCLUDE 'aer_Global.h'
      integer na,nb
! see explanation in FCN fdhetg

      if (nb.eq.1) then
         x1 = yxkmtd(na,ind_HNO3) * ycwd(na)
         caq = 0.d0
         if ((yxeq(ind_HNO3)+1.d-2).ne.0.d0) then
            if (na.eq.1) caq=((C(ind_HNO3l1)+C(ind_NO3ml1))*1.d-2)/ &
                 (yxeq(ind_HNO3) + 1.d-2)
            if (na.eq.2) caq=((C(ind_HNO3l2)+C(ind_NO3ml2))*1.d-2)/ &
                 (yxeq(ind_HNO3) + 1.d-2)
         endif
         x2 = 0.d0
          if (C(ind_HNO3).ne.0.d0.and.yhenry(ind_HNO3).ne.0.d0) &
              x2=-yxkmtd(na,ind_HNO3)/(C(ind_HNO3)*yhenry(ind_HNO3))*caq
         xkt = max(0.d0,(x1 + x2))
      endif
      if (nb.eq.2) xkt=yxkmtd(na,ind_N2O5) * ycwd(na)
      if (nb.eq.3) xkt=yxkmtd(na,ind_NH3) * ycwd(na)
      if (nb.eq.4) xkt=yxkmtd(na,ind_H2SO4) * ycwd(na)

      fdheta=xkt

      end function fdheta

!-----------------------------------------------------------------------------


      double precision function fdhett (na,nb)
! heterogeneous rate function
! a0=1..2  liquid size class

      implicit double precision (a-h,o-z)

      INCLUDE 'tot_Parameters.h'
      INCLUDE 'tot_Global.h'
      integer na,nb
! see explanation in FCN fdhetg

      if (nb.eq.1) then
         x1 = yxkmtd(na,ind_HNO3) * ycwd(na)
         caq = 0.d0
         if ((yxeq(ind_HNO3)+1.d-2).ne.0.d0) then
            if (na.eq.1) caq=((C(ind_HNO3l1)+C(ind_NO3ml1))*1.d-2)/ &
                 (yxeq(ind_HNO3) + 1.d-2)
            if (na.eq.2) caq=((C(ind_HNO3l2)+C(ind_NO3ml2))*1.d-2)/ &
                 (yxeq(ind_HNO3) + 1.d-2)
         endif
         x2 = 0.d0
          if (C(ind_HNO3).ne.0.d0.and.yhenry(ind_HNO3).ne.0.d0) &
              x2=-yxkmtd(na,ind_HNO3)/(C(ind_HNO3)*yhenry(ind_HNO3))*caq
         xkt = max(0.d0,(x1 + x2))
      endif
      if (nb.eq.2) xkt=yxkmtd(na,ind_N2O5) * ycwd(na)
      if (nb.eq.3) xkt=yxkmtd(na,ind_NH3) * ycwd(na)
      if (nb.eq.4) xkt=yxkmtd(na,ind_H2SO4) * ycwd(na)

      fdhett=xkt

      end function fdhett

!-----------------------------------------------------------------------------
!     double precision function DMS_add (c) ! jjb argument not used
      double precision function DMS_add ()  ! jjb removed (also in mech/master_gas.eqn)
! calculate special rate function for DMS + OH addition; IUPAC 10/06
! k(298K)=2.2d-12 cm3/(mlc s)
      implicit none

      double precision o2,tte

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      o2=0.21*aircc
      tte=1./te
      DMS_add=9.5d-39*exp(5270.*tte)*o2/(1.+7.5d-29*exp(5610.*tte)*o2)
      end function DMS_add

!----------------------------------------------------------------------------

      function a_n2o5(k,kc)
!     reactive uptake coeff of N2O5 (Bertram and Thornton, 2009)  PJ

! jjb: adapt for Mistra v9, implicit none
      ! cleaning: remove /cb52a/ f_sum

      USE global_params, ONLY : &
! Imported Parameters:
           j2, &
           j6, &
           n, &
           nkc

      USE precision, ONLY : &
! Imported Parameters:
           dp

      implicit none
      integer, intent (in) :: k, kc
      real(kind=dp) :: a_n2o5
      real(kind=dp) :: denom, xclm, xh2o, xk2f, xno3m
      real(kind=dp), parameter :: ppsmall = 1e-25_dp

      common /blck12/ cw(nkc,n),cm(nkc,n)
      real(kind=dp) :: cw, cm
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      real(kind=dp) :: sl1, sion1
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho

      if (cw(kc,k).gt.0.d0) then
!     convert from (mol m-3) to (mol L-1), species 13: NO3-
!     convert from (mol m-3) to (mol L-1), species 14: Cl-
         xno3m = sion1(13,kc,k)/cw(kc,k) * 1e-3_dp
         xclm  = sion1(14,kc,k)/cw(kc,k) * 1e-3_dp
      else
         xno3m = 0._dp
         xclm = 0._dp
      end if

!     initialize water to (mol L-1)
      xh2o = 0._dp
      if (cm(1,k).gt.0._dp) then
      if (cw(1,k).gt.0._dp) then
         xh2o = 55.55_dp * (cm(1,k) / cw(1,k))
      end if
      end if

!     big honking reactive uptake coefficient parameterization
      xk2f = 1.15e6_dp - 1.15e6_dp * exp(-0.13_dp * xh2o)

      denom = 1._dp
      if (xno3m.gt.0.d0) &
!      if (xno3m.gt.ppsmall) &
           denom = 1._dp + 6.e-2_dp*xh2o/xno3m + 29._dp*xclm/xno3m
      a_n2o5 = 3.2e-8_dp * xk2f * (1._dp - (1._dp/denom))

!      if (k.eq.2) print *,'gamma(N2O5),k=2',k,a_n2o5,xh2o,xno3m &
!           ,xclm,denom,t(k),p(k)
!      if (k.eq.10) print *,'gamma(N2O5),k=10',k,a_n2o5,xh2o,xno3m &
!           ,xclm,denom,t(k),p(k)

      end function a_n2o5

!-----------------------------------------------------------------------------

!      double precision function xkHgBr (x1)
!! rate coefficient for recombination Hg+Br --> HgBr (Donohoue et al., 2006, #4161
!      double precision aircc,te,h2oppm,pk,x1,x2
!      common /cb_1/ aircc,te,h2oppm,pk
!
!! reaction is 3rd order, multiply with conversion factors here (instead of in
!! master_gas.eqn) in order to have an argument for the function
!! note that the fit that they give does not exactly reproduce the measured values in their Tables 1 and 2)
!      xkHgBr = 1.46d-32 * (te/298.)**(-1.86) * x1 * x1
!      end function xkHgBr

!-----------------------------------------------------------------------------

!      double precision function xkHgBrBr (x1)
!! rate coefficient for recombination HgBr+Br --> HgBr2 (Goodsite et al., 2004, #3244
!      double precision aircc,te,h2oppm,pk,x1
!      common /cb_1/ aircc,te,h2oppm,pk
!
!! reaction is 2nd order, multiply with conversion factor here (instead of in
!! master_gas.eqn) in order to have an argument for the function
!      xkHgBrBr = 2.5d-10*(te/298.d0)**(-.57d0) * x1
!      end function xkHgBrBr

!-----------------------------------------------------------------------------

!      double precision function xkGood (x1)
!! rate coefficient for HgBr dissociation, use Goodsite et al., 2004
!! (#3244) but scaled with ratio of Donohoue and Goodsite as
!! suggested in Seigneur and Lohmann, 2008 (#4136)
!      double precision aircc,te,h2oppm,pk,x1,x2,xkHgBr
!      common /cb_1/ aircc,te,h2oppm,pk
!
!      if (te.eq.0.d0) then
!         xkGood = 0.d0
!      else
!! argument for kHgBr is unity to avoid scaling with conversion factors
!         x2=xkHgBr(1.d0) / (1.1d-12 * (te/298.d0) **(-2.37))
!! multiply x2 with air density as Donohoue is 3rd order and Goodsite 2nd order
!! [M in mol/m3] = p / RT; convert to molec/cm3 : * conv1
!         x2 = x2 * pk / (8.314 * te) * x1
!         xkGood = 1.2d10 * exp(-8357.d0 / te) * x2
!      endif
!
!      end function xkGood

!-----------------------------------------------------------------------------



