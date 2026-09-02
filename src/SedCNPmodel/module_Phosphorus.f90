module Phosphorus 

   use CNPfunctions
   use module_SedCNPvariables
!=====||__WHQ5403q__||=====!
!
   use SedCNP_config,         only: SedCNPmodel
!
!=====||__WHQ5403q__||=====!
   use CNPparams

   implicit none

   contains

   subroutine SoilPhosphorusCycle(i, j, lcover, itime, TSOIL, THETA, &
      THETA_S, THETA_F, So_WStorage, Quptk, Qsurf, Qperc, Qintf)

      implicit none

      ! --- Input arguments
      integer,               intent(in) :: i, j, lcover, itime
      real(8), dimension(4), intent(in) :: TSOIL, THETA, So_WStorage, Quptk
      real(8),               intent(in) :: THETA_S, THETA_F, Qsurf  !BK20260318
      real(8), dimension(4), intent(in) :: Qperc
      real(8),               intent(in) :: Qintf

      ! --- Local variables
      integer               :: isl, status, imonth, st, ips
      character(len=2)      :: smonth
      real(8)               :: dt, rate0, KTBMBM, KTHETA, KTC, beta, RK4SOL, dPI
      real(8), allocatable, dimension(:) :: Qintf_layer, Weight  !BK20260612
      real(8)               :: TotalWeight  !BK20251225
      real(8)               :: dx    ! grid cell size [m]  !BK20260301
      real(8)               :: q_out, w_strg  !(m3) !BK20260301      
      integer               :: k  !BK20251225

      ! Decomposition and transport fluxes
      real(8) :: dLPOPLITD, dLPOPRESD, dLPOPEXCD, dLPOPMAND
      real(8) :: dRPOPD, dLDOPD, dRDOPD, dMBMPD
      real(8) :: dLPOPLITsurf, dLPOPRESsurf, dLPOPEXCsurf, dLPOPMANsurf
      real(8), allocatable :: dPIPAdmix(:), dPIPSdmix(:)
      real(8), allocatable :: dPIPAsurf_flux(:), dPIPSsurf_flux(:)
      real(8), allocatable :: dSSA(:), dDmixSoilWT(:,:)

      ! Initial State and Flux Scaling Variables
      real(8) :: So_LPOPLIT_init, So_LPOPRES_init, So_LPOPEXC_init, &
                 So_LPOPMAN_init
      real(8) :: So_RPOP_init, So_LDOP_init, So_RDOP_init, So_MBMP_init
      real(8) :: So_PO4_init, So_PIPA_init, So_PIPS_init
      real(8) :: TotalFlux_LPOPLIT, TotalFlux_LPOPRES, TotalFlux_LPOPEXC, &
                 TotalFlux_LPOPMAN
      real(8) :: TotalFlux_RPOP, TotalFlux_LDOP, TotalFlux_RDOP, &
                 TotalFlux_MBMP
      real(8) :: TotalFlux_PO4, TotalFlux_PIPA, TotalFlux_PIPS
      real(8) :: Scale_LPOPLIT, Scale_LPOPRES, Scale_LPOPEXC, Scale_LPOPMAN
      real(8) :: Scale_RPOP, Scale_LDOP, Scale_RDOP, Scale_MBMP
      real(8) :: Scale_PO4, Scale_PIPA, Scale_PIPS

      allocate(dPIPAdmix(nps), dPIPSdmix(nps), dPIPAsurf_flux(nps), &
               dPIPSsurf_flux(nps), dSSA(nst), dDmixSoilWT(nst,nps))

      ! ====================================================================
      ! I. INITIALIZE & UPDATE STORAGES WITH EXTERNAL/INTERNAL INPUTS
      ! ====================================================================
      
      dt = real(SedCNPmodel%SedCNP_timestep)
      dx = SedCNP_hydro%dx(1)    !BK20260301
      smonth = trim(dateSedCNP%olddate(6:7))
      read(smonth,*) imonth

      allocate(Qintf_layer(domain%nsl), Weight(domain%nsl)) !BK20260612
      
      ! Apportion total interflow Qintf to layers based on water content above 
      ! field capacity  !BK20251225
      TotalWeight = 0.0
      DO k = 1, domain%nsl
         Weight(k) = max(0.0, THETA(k) - THETA_F) * DSOIL(k)
         TotalWeight = TotalWeight + Weight(k)
      END DO
      if (TotalWeight > 0.0) then
         DO k = 1, domain%nsl
            Qintf_layer(k) = Qintf * Weight(k) / TotalWeight
         END DO
      else
         DO k = 1, domain%nsl
            Qintf_layer(k) = 0.0
         END DO
      endif
      
      ! =================================================================
      ! II. CALCULATE-SCALE-UPDATE LOGIC
      ! =================================================================
      
      DO isl = 1, domain%nsl
         
         ! -----------------------------------------------------------------
         ! A. Store initial values for simultaneous flux calculation
         ! (previous timestep state + current timestep inputs)         
         ! -----------------------------------------------------------------

         if (isl == 1) then
            So_LPOPLITin(i,j) = So_LPOPLITinso(lcover,imonth) * dt
            So_LPOPEXCin(i,j) = So_LPOPEXCinso(lcover,imonth) * dt
            So_LPOPMANin(i,j) = So_LPOPMANinso(lcover,imonth) * dt
            So_LPOPRESin(i,isl,j) = So_LPOPRESinso(lcover,imonth) * dt * FROOT(lcover, isl)
            So_PO4fxin(i,j)   = So_PO4fert(lcover,imonth) * dt

            So_LPOPLIT_init = So_LPOPLIT0(i,j) + So_LPOPLITin(i,j)
            So_LPOPEXC_init = So_LPOPEXC0(i,j) + So_LPOPEXCin(i,j)
            So_LPOPMAN_init = So_LPOPMAN0(i,j) + So_LPOPMANin(i,j)
            So_LPOPRES_init = So_LPOPRES0(i,isl,j) + So_LPOPRESin(i,isl,j)
            So_RPOP_init    = So_RPOP0(i,isl,j)
            So_LDOP_init    = So_LDOP0(i,isl,j) + So_LDOPgwso0(i,j)  !BK20260316
            So_RDOP_init    = So_RDOP0(i,isl,j) + So_RDOPgwso0(i,j)  !BK20260316
            So_MBMP_init    = So_MBMP0(i,isl,j)
            So_PO4_init     = So_PO40(i,isl,j)  + So_PO4gwso0(i,j) + So_PO4fxin(i,j)  !BK20260316
            So_PIPA_init    = So_PIPA0(i,isl,j)
            So_PIPS_init    = So_PIPS0(i,isl,j)
         else
            So_LPOPRESin(i,isl,j) = So_LPOPRESinso(lcover,imonth) * dt * FROOT(lcover, isl)

            So_LPOPLIT_init = 0.0  ! Does not exist in lower layers
            So_LPOPEXC_init = 0.0  ! Does not exist in lower layers
            So_LPOPMAN_init = 0.0  ! Does not exist in lower layers
            So_LPOPRES_init = So_LPOPRES0(i,isl,j) + So_LPOPRESin(i,isl,j)
            So_RPOP_init    = So_RPOP0(i,isl,j)
            So_LDOP_init    = So_LDOP0(i,isl,j) + So_LDOPperc0(i,isl-1,j)
            So_RDOP_init    = So_RDOP0(i,isl,j) + So_RDOPperc0(i,isl-1,j)
            So_MBMP_init    = So_MBMP0(i,isl,j)
            So_PO4_init     = So_PO40(i,isl,j)  + So_PO4perc0(i,isl-1,j)
            So_PIPA_init    = So_PIPA0(i,isl,j)
            So_PIPS_init    = So_PIPS0(i,isl,j)
         endif

         ! -----------------------------------------------------------------
         ! B. Calculate all transformation and transport fluxes
         ! -----------------------------------------------------------------

         KTBMBM = calc_KTB (TSOIL(isl), TMBM, MMBM)
         KTHETA = calc_KTHETA (THETA(isl), THETA_S, THETA_F)
         KTC    = EXP (0.115 * TSOIL(isl) - 2.88)

         ! --- Organic P Decomposition ---
         if (isl == 1) then
            rate0 = KTBMBM * KTHETA * KLIT(lcover)
            RK4SOL = calc_RK4SOL(-rate0, So_LPOPLIT_init)
            dLPOPLITD = So_LPOPLIT_init - RK4SOL
            
            rate0 = KTBMBM * KTHETA * KEXC(lcover)
            RK4SOL = calc_RK4SOL(-rate0, So_LPOPEXC_init)
            dLPOPEXCD = So_LPOPEXC_init - RK4SOL
            
            rate0 = KTBMBM * KTHETA * KMAN(lcover)
            RK4SOL = calc_RK4SOL(-rate0, So_LPOPMAN_init)
            dLPOPMAND = So_LPOPMAN_init - RK4SOL
         endif
         
         rate0 = KTBMBM * KTHETA * KRES(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_LPOPRES_init)
         dLPOPRESD = So_LPOPRES_init - RK4SOL
         
         rate0 = KTBMBM * KTHETA * KRPOM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_RPOP_init)
         dRPOPD = So_RPOP_init - RK4SOL
         
         rate0 = KTBMBM * KTHETA * KLDOM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_LDOP_init)
         dLDOPD = So_LDOP_init - RK4SOL
         
         rate0 = KTBMBM * KTHETA * KRDOM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_RDOP_init)
         dRDOPD = So_RDOP_init - RK4SOL
         
         rate0 = KTBMBM * KTHETA * KMBM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_MBMP_init)
         dMBMPD = So_MBMP_init - RK4SOL
         
         ! --- Inorganic P Transformations & Plant Uptake ---
         if ((Quptk(isl) > 0.0) .and. (So_WStorage(isl) > 0.0)) then
            So_PO4uptk(i,isl,j) = (1.0 - exp(-Quptk(isl) * domain%areaxy / &
               1000.0 * domain%DT / So_WStorage(isl))) * So_PO4_init
         else
            So_PO4uptk(i,isl,j) = 0.0
         endif
         
         if (So_PO4_init .ge. SIGMAPO4 * So_PIPA_init) then
            if (So_PIPA_init > 0.0) then
               beta = So_PO4_init / (SIGMAPO4 * So_PIPA_init)
            else
               beta = 100.0
            endif
            So_PO4ADS(i,isl,j) = KTC * KTHETA * KADSPO4 * &
               (1.0 - exp(-beta)) * &
               (So_PO4_init - SIGMAPO4 * So_PIPA_init) * dt
         else
            if (So_PO4_init > 0.0) then
               beta = (SIGMAPO4 * So_PIPA_init) / So_PO4_init
            else
               beta = 100.0
            endif
            So_PO4ADS(i,isl,j) = 0.1 * KTC * KTHETA * KADSPO4 * &
               (1.0 - exp(-beta)) * &
               (So_PO4_init - SIGMAPO4 * So_PIPA_init) * dt
         endif
         
         if (So_PIPA_init .ge. SIGMAPIPA * So_PIPS_init) then
            if (So_PIPS_init > 0.0) then
               beta = So_PIPA_init / (SIGMAPIPA * So_PIPS_init)
            else
               beta = 100.0
            endif
            So_PIPAADS(i,isl,j) = KTC * KTHETA * KADSPIPA * &
               (1.0 - exp(-beta)) * &
               (So_PIPA_init - SIGMAPIPA * So_PIPS_init) * dt
         else
            if (So_PIPA_init > 0.0) then
               beta = (SIGMAPIPA * So_PIPS_init) / So_PIPA_init
            else
               beta = 100.0
            endif
            So_PIPAADS(i,isl,j) = 0.1 * KTC * KTHETA * KADSPIPA * &
               (1.0 - exp(-beta)) * &
               (So_PIPA_init - SIGMAPIPA * So_PIPS_init) * dt
         endif

         ! --- Transport by surface runoff
         if (isl == 1) then
            q_out = Qsurf / 1000. * dx * dx * dt ! (m3)
            w_strg = THETA(1) * Dmix * dx * dx   ! (m3)  !BK20260302
            !POM
            dPI = calc_dPI(etaPOMsurf, q_out, w_strg)  !BK20260302
            dLPOPLITsurf = So_LPOPLIT_init * (Dmix / DSOIL(isl)) * dPI
            dLPOPEXCsurf = So_LPOPEXC_init * (Dmix / DSOIL(isl)) * dPI
            dLPOPMANsurf = So_LPOPMAN_init * (Dmix / DSOIL(isl)) * dPI
            dLPOPRESsurf = So_LPOPRES_init * (Dmix / DSOIL(isl)) * dPI
            So_LPOPsurf(i,j) = dLPOPLITsurf + dLPOPRESsurf + &
                               dLPOPEXCsurf + dLPOPMANsurf
            So_RPOPsurf(i,j) = So_RPOP_init * (Dmix / DSOIL(isl)) * dPI
            So_MBMPsurf(i,j) = So_MBMP_init * (Dmix / DSOIL(isl)) * dPI
            !DOM
            dPI = calc_dPI(etaDOMsurf, q_out, w_strg)  !BK20260302
            So_LDOPsurf(i,j) = So_LDOP_init * (Dmix / DSOIL(isl)) * dPI
            So_RDOPsurf(i,j) = So_RDOP_init * (Dmix / DSOIL(isl)) * dPI
            !PO4
            dPI = calc_dPI(etaPO4surf, q_out, w_strg)  !BK20260302
            So_PO4surf(i,j) = So_PO4_init * (Dmix / DSOIL(isl)) * dPI
            
            if (SedCNP_hydro%sltyp(i,j) /= SedCNP_hydro%sltyp(i,j)) then
               st = 14  ! HARDCODED for 'water'
            else
               st = int(SedCNP_hydro%sltyp(i,j))
            endif
            dSSA(st) = 0.0
            do ips = 1, nps
               dSSA(st) = dSSA(st) + SSA(ips) * overSed%SOILPSF(st,ips)
               dDmixSoilWT(st,ips) = overSed%SOILPSF(st,ips) * &
                  domain%areaxy * Dmix * (1.0 - THETA_S) * 2650.0
            enddo
            
            So_PIPAsurf(i,j) = 0.0
            So_PIPSsurf(i,j) = 0.0
            do ips = 1, nps
               !BK20260213 bug fixed
               dPIPAdmix(ips) = So_PIPA_init * (Dmix / DSOIL(isl)) * & 
                  (SSA(ips) * overSed%SOILPSF(st,ips) / dSSA(st))
               dPIPSdmix(ips) = So_PIPS_init * (Dmix / DSOIL(isl)) * & 
                  (SSA(ips) * overSed%SOILPSF(st,ips) / dSSA(st))

               if (overSed%Ssurf(ips,i,j) < 1.e-20 .or. &
                   overSed%Ssurf(ips,i,j) /= overSed%Ssurf(ips,i,j)) then
                  dPIPAsurf_flux(ips)=0.0
                  dPIPSsurf_flux(ips)=0.0
               else
                  dPIPAsurf_flux(ips) = dPIPAdmix(ips) * &
                     overSed%Ssurf(ips,i,j) / dDmixSoilWT(st,ips)
                  dPIPSsurf_flux(ips) = dPIPSdmix(ips) * &
                     overSed%Ssurf(ips,i,j) / dDmixSoilWT(st,ips)
               endif
               So_PIPAsurf(i,j) = So_PIPAsurf(i,j) + dPIPAsurf_flux(ips) * dt
               So_PIPSsurf(i,j) = So_PIPSsurf(i,j) + dPIPSsurf_flux(ips) * dt
            enddo
         endif
         
         ! --- Transport by Interflow
         q_out = Qintf_layer(isl) / 1000. * dx * dx * dt  ! (m3)
         w_strg = THETA(isl) * dx * dx    ! (m3)  !BK20260302   
         !DOM         
         dPI = calc_dPI(etaDOMintf, q_out, w_strg)  !BK20260302
         So_LDOPintf(i,isl,j) = So_LDOP_init * dPI
         So_RDOPintf(i,isl,j) = So_RDOP_init * dPI
         !PO4      
         dPI = calc_dPI(etaPO4intf, q_out, w_strg)  !BK20260302
         So_PO4intf(i,isl,j) = So_PO4_init * dPI

         ! --- Transport by Percolation         
         q_out = Qperc(isl) / 1000. * dx * dx * dt  ! (m3)
         w_strg = THETA(isl) * dx * dx    ! (m3)  !BK20260302   
         !DOM
         dPI = calc_dPI(etaDOMperc, q_out, w_strg)  !BK20260302
         So_LDOPperc(i,isl,j) = So_LDOP_init * dPI
         So_RDOPperc(i,isl,j) = So_RDOP_init * dPI
         !PO4
         dPI = calc_dPI(etaPO4perc, q_out, w_strg)  !BK20260302
         So_PO4perc(i,isl,j) = So_PO4_init * dPI

         ! -----------------------------------------------------------------
         ! C. Scale fluxes to ensure mass conservation
         ! -----------------------------------------------------------------

         if (isl == 1) then
            TotalFlux_LPOPLIT = dLPOPLITD + dLPOPLITsurf
            if (TotalFlux_LPOPLIT > So_LPOPLIT_init .and. &
                  TotalFlux_LPOPLIT > 0.0) then
               Scale_LPOPLIT = So_LPOPLIT_init / TotalFlux_LPOPLIT
               dLPOPLITD = dLPOPLITD * Scale_LPOPLIT
               dLPOPLITsurf = dLPOPLITsurf * Scale_LPOPLIT
            endif
            
            TotalFlux_LPOPEXC = dLPOPEXCD + dLPOPEXCsurf
            if (TotalFlux_LPOPEXC > So_LPOPEXC_init .and. &
                  TotalFlux_LPOPEXC > 0.0) then
               Scale_LPOPEXC = So_LPOPEXC_init / TotalFlux_LPOPEXC
               dLPOPEXCD = dLPOPEXCD * Scale_LPOPEXC
               dLPOPEXCsurf = dLPOPEXCsurf * Scale_LPOPEXC
            endif
            
            TotalFlux_LPOPMAN = dLPOPMAND + dLPOPMANsurf
            if (TotalFlux_LPOPMAN > So_LPOPMAN_init .and. &
                  TotalFlux_LPOPMAN > 0.0) then
               Scale_LPOPMAN = So_LPOPMAN_init / TotalFlux_LPOPMAN
               dLPOPMAND = dLPOPMAND * Scale_LPOPMAN
               dLPOPMANsurf = dLPOPMANsurf * Scale_LPOPMAN
            endif
         endif
         
         if (isl == 1) then
            TotalFlux_LPOPRES = dLPOPRESD + dLPOPRESsurf
         else
            TotalFlux_LPOPRES = dLPOPRESD
         endif
         if (TotalFlux_LPOPRES > So_LPOPRES_init .and. &
               TotalFlux_LPOPRES > 0.0) then
            Scale_LPOPRES = So_LPOPRES_init / TotalFlux_LPOPRES
            dLPOPRESD = dLPOPRESD * Scale_LPOPRES
            if(isl == 1) dLPOPRESsurf = dLPOPRESsurf * Scale_LPOPRES
         endif
         
         if (isl == 1) then
             TotalFlux_RPOP = dRPOPD + So_RPOPsurf(i,j)
         else
             TotalFlux_RPOP = dRPOPD
         endif
         if (TotalFlux_RPOP > So_RPOP_init .and. TotalFlux_RPOP > 0.0) then
            Scale_RPOP = So_RPOP_init / TotalFlux_RPOP
            dRPOPD = dRPOPD * Scale_RPOP
            if(isl == 1) So_RPOPsurf(i,j) = So_RPOPsurf(i,j) * Scale_RPOP
         endif

         if (isl == 1) then
            TotalFlux_LDOP = dLDOPD + So_LDOPsurf(i,j) + &
               So_LDOPintf(i,isl,j) + So_LDOPperc(i,isl,j)
         else
            TotalFlux_LDOP = dLDOPD + So_LDOPintf(i,isl,j) + &
               So_LDOPperc(i,isl,j)
         endif
         if(TotalFlux_LDOP > So_LDOP_init .and. TotalFlux_LDOP > 0.0) then
            Scale_LDOP = So_LDOP_init / TotalFlux_LDOP
            dLDOPD = dLDOPD * Scale_LDOP
            if(isl == 1) So_LDOPsurf(i,j) = So_LDOPsurf(i,j) * Scale_LDOP
            So_LDOPintf(i,isl,j) = So_LDOPintf(i,isl,j) * Scale_LDOP
            So_LDOPperc(i,isl,j) = So_LDOPperc(i,isl,j) * Scale_LDOP
         endif
         
         if (isl == 1) then
            TotalFlux_RDOP = dRDOPD + So_RDOPsurf(i,j) + &
               So_RDOPintf(i,isl,j) + So_RDOPperc(i,isl,j)
         else
            TotalFlux_RDOP = dRDOPD + So_RDOPintf(i,isl,j) + &
               So_RDOPperc(i,isl,j)
         endif
         if(TotalFlux_RDOP > So_RDOP_init .and. TotalFlux_RDOP > 0.0) then
            Scale_RDOP = So_RDOP_init / TotalFlux_RDOP
            dRDOPD = dRDOPD * Scale_RDOP
            if(isl == 1) So_RDOPsurf(i,j) = So_RDOPsurf(i,j) * Scale_RDOP
            So_RDOPintf(i,isl,j) = So_RDOPintf(i,isl,j) * Scale_RDOP
            So_RDOPperc(i,isl,j) = So_RDOPperc(i,isl,j) * Scale_RDOP
         endif

         if (isl == 1) then
            TotalFlux_MBMP = dMBMPD + So_MBMPsurf(i,j)
         else
            TotalFlux_MBMP = dMBMPD
         endif
         if (TotalFlux_MBMP > So_MBMP_init .and. TotalFlux_MBMP > 0.0) then
            Scale_MBMP = So_MBMP_init / TotalFlux_MBMP
            dMBMPD = dMBMPD * Scale_MBMP
            if(isl == 1) So_MBMPsurf(i,j) = So_MBMPsurf(i,j) * Scale_MBMP
         endif
         
         if (So_PO4ADS(i,isl,j) < 0.0) then
            So_PO4ADS(i,isl,j) = &
               min(abs(So_PO4ADS(i,isl,j)), So_PIPA_init) * (-1.0)
         endif
         if (So_PIPAADS(i,isl,j) < 0.0) then
            So_PIPAADS(i,isl,j) = &
               min(abs(So_PIPAADS(i,isl,j)), So_PIPS_init) * (-1.0)
         endif
         
         if (isl == 1) then
            TotalFlux_PO4 = So_PO4uptk(i,isl,j) + So_PO4surf(i,j) + &
               So_PO4intf(i,isl,j) + So_PO4perc(i,isl,j)
         else
            TotalFlux_PO4 = So_PO4uptk(i,isl,j) + So_PO4intf(i,isl,j) + &
               So_PO4perc(i,isl,j)
         endif
         if(So_PO4ADS(i,isl,j) > 0.0) &
            TotalFlux_PO4 = TotalFlux_PO4 + So_PO4ADS(i,isl,j)
         
         if (TotalFlux_PO4 > So_PO4_init .and. TotalFlux_PO4 > 0.0) then
            Scale_PO4 = So_PO4_init / TotalFlux_PO4
            So_PO4uptk(i,isl,j) = So_PO4uptk(i,isl,j) * Scale_PO4
            if(So_PO4ADS(i,isl,j) > 0.0) &
               So_PO4ADS(i,isl,j) = So_PO4ADS(i,isl,j) * Scale_PO4
            if(isl == 1) So_PO4surf(i,j) = So_PO4surf(i,j) * Scale_PO4
            So_PO4intf(i,isl,j) = So_PO4intf(i,isl,j) * Scale_PO4
            So_PO4perc(i,isl,j) = So_PO4perc(i,isl,j) * Scale_PO4
         endif
         
         if (isl == 1) then
            TotalFlux_PIPA = So_PIPAsurf(i,j)
         else
            TotalFlux_PIPA = 0.0
         endif
         if(So_PIPAADS(i,isl,j) > 0.0) &
            TotalFlux_PIPA = TotalFlux_PIPA + So_PIPAADS(i,isl,j)

         if (TotalFlux_PIPA > So_PIPA_init .and. TotalFlux_PIPA > 0.0) then
            Scale_PIPA = So_PIPA_init / TotalFlux_PIPA
            if(isl==1) So_PIPAsurf(i,j) = So_PIPAsurf(i,j) * Scale_PIPA
            if(So_PIPAADS(i,isl,j) > 0.0) &
               So_PIPAADS(i,isl,j) = So_PIPAADS(i,isl,j) * Scale_PIPA
         endif
         
         if (isl == 1) then
            TotalFlux_PIPS = So_PIPSsurf(i,j)
         else
            TotalFlux_PIPS = 0.0
         endif
         if (TotalFlux_PIPS > So_PIPS_init .and. TotalFlux_PIPS > 0.0) then
            Scale_PIPS = So_PIPS_init / TotalFlux_PIPS
            if(isl==1) So_PIPSsurf(i,j) = So_PIPSsurf(i,j) * Scale_PIPS
         endif

         ! Recalculate So_LPOPsurf with potentially scaled component fluxes 
         if (isl == 1) then
            So_LPOPsurf(i,j) = dLPOPLITsurf + dLPOPRESsurf + &
                                 dLPOPEXCsurf + dLPOPMANsurf
         endif

         ! -----------------------------------------------------------------
         ! D. Update final storages using scaled fluxes
         ! -----------------------------------------------------------------

         if (isl == 1) then
            So_LPOPLIT(i,j) = So_LPOPLIT_init - dLPOPLITD - dLPOPLITsurf
            So_LPOPEXC(i,j) = So_LPOPEXC_init - dLPOPEXCD - dLPOPEXCsurf
            So_LPOPMAN(i,j) = So_LPOPMAN_init - dLPOPMAND - dLPOPMANsurf
            So_LPOPRES(i,isl,j) = So_LPOPRES_init - dLPOPRESD - dLPOPRESsurf
            
            So_RPOP(i,isl,j) = So_RPOP_init + FDEC * FREF * (dLPOPLITD + &
               dLPOPRESD + dLPOPEXCD + dLPOPMAND) - dRPOPD - So_RPOPsurf(i,j)
            So_LDOP(i,isl,j) = So_LDOP_init + FDEC * (1. - FREF) * &
               (dLPOPLITD + dLPOPRESD + dLPOPEXCD + dLPOPMAND) + &
               FDEC * dMBMPD - dLDOPD - So_LDOPsurf(i,j) - &
               So_LDOPintf(i,isl,j) - So_LDOPperc(i,isl,j)
            So_RDOP(i,isl,j) = So_RDOP_init + FDEC * (dLDOPD + dRPOPD + &
               dRDOPD)- dRDOPD - So_RDOPsurf(i,j) - So_RDOPintf(i,isl,j) - &
               So_RDOPperc(i,isl,j)
            So_MBMP(i,isl,j) = So_MBMP_init + FBIO * (dLPOPLITD + dLPOPRESD + &
               dLPOPEXCD + dLPOPMAND + dRPOPD + dLDOPD + dRDOPD + dMBMPD) - &
               dMBMPD - So_MBMPsurf(i,j)
            So_PO4(i,isl,j) = So_PO4_init + FMET * (dLPOPLITD + dLPOPRESD + &
               dLPOPEXCD + dLPOPMAND + dRPOPD + dLDOPD + dRDOPD + dMBMPD) - &
               So_PO4ADS(i,isl,j) - So_PO4uptk(i,isl,j) - So_PO4surf(i,j) - &
               So_PO4intf(i,isl,j) - So_PO4perc(i,isl,j)
            So_PIPA(i,isl,j) = So_PIPA_init + So_PO4ADS(i,isl,j) - &
               So_PIPAADS(i,isl,j) - So_PIPAsurf(i,j)
            So_PIPS(i,isl,j) = So_PIPS_init + So_PIPAADS(i,isl,j) - &
               So_PIPSsurf(i,j)
         else
            So_LPOPRES(i,isl,j) = So_LPOPRES_init - dLPOPRESD
            So_RPOP(i,isl,j) = So_RPOP_init + FDEC * FREF * dLPOPRESD - dRPOPD
            So_LDOP(i,isl,j) = So_LDOP_init + FDEC * (1. - FREF) * &
               dLPOPRESD + FDEC*dMBMPD - dLDOPD - So_LDOPintf(i,isl,j) - &
               So_LDOPperc(i,isl,j)
            So_RDOP(i,isl,j) = So_RDOP_init + FDEC * (dLDOPD + dRPOPD + &
               dRDOPD) - dRDOPD - So_RDOPintf(i,isl,j) - So_RDOPperc(i,isl,j)
            So_MBMP(i,isl,j) = So_MBMP_init + FBIO * (dLPOPRESD + dRPOPD + &
               dLDOPD + dRDOPD + dMBMPD) - dMBMPD
            So_PO4(i,isl,j)  = So_PO4_init + FMET * (dLPOPRESD + dRPOPD + &
               dLDOPD + dRDOPD + dMBMPD) - So_PO4ADS(i,isl,j) - &
               So_PO4uptk(i,isl,j) - So_PO4intf(i,isl,j) - So_PO4perc(i,isl,j)
            So_PIPA(i,isl,j) = So_PIPA_init + So_PO4ADS(i,isl,j) - &
               So_PIPAADS(i,isl,j)
            So_PIPS(i,isl,j) = So_PIPS_init + So_PIPAADS(i,isl,j)
         endif
      END DO

      ! Enforce non-negative soil storages
      if (So_LPOPLIT(i,j) < 1.0e-20) So_LPOPLIT(i,j) = 0.0
      if (So_LPOPEXC(i,j) < 1.0e-20) So_LPOPEXC(i,j) = 0.0
      if (So_LPOPMAN(i,j) < 1.0e-20) So_LPOPMAN(i,j) = 0.0
      DO isl = 1, domain%nsl
        if (So_LPOPRES(i,isl,j) < 1.0e-20) So_LPOPRES(i,isl,j) = 0.0
        if (So_RPOP(i,isl,j) < 1.0e-20) So_RPOP(i,isl,j) = 0.0
        if (So_LDOP(i,isl,j) < 1.0e-20) So_LDOP(i,isl,j) = 0.0
        if (So_RDOP(i,isl,j) < 1.0e-20) So_RDOP(i,isl,j) = 0.0
        if (So_MBMP(i,isl,j) < 1.0e-20) So_MBMP(i,isl,j) = 0.0
        if (So_PO4(i,isl,j) < 1.0e-20) So_PO4(i,isl,j) = 0.0
        if (So_PIPA(i,isl,j) < 1.0e-20) So_PIPA(i,isl,j) = 0.0
        if (So_PIPS(i,isl,j) < 1.0e-20) So_PIPS(i,isl,j) = 0.0
      END DO
      
      ! ------------------------------------------------------------------------
      ! III. UPDATE OUTPUTS TO OTHER COMPONENTS
      ! ------------------------------------------------------------------------

      So_LDOPsogw(i,j) = So_LDOPperc(i,domain%nsl,j)  !BK20260318
      So_RDOPsogw(i,j) = So_RDOPperc(i,domain%nsl,j)  !BK20260318
      So_PO4sogw(i,j)  = So_PO4perc(i,domain%nsl,j)  !BK20260318
      
   end subroutine SoilPhosphorusCycle
   

   ! subroutine AquiferPhosphorusCycle(i,j,Qaqso,Aq_WStorage,Qaqgw)
   !    implicit none
      
   !    ! --- Input arguments
   !    integer, intent(in) :: i, j
   !    real(8), intent(in) :: Qaqso, Aq_WStorage, Qaqgw

   !    ! --- Local Variables
   !    real(8) :: KTC, KTHETA, beta
   !    real(8) :: q_out, dx    !BK20260302
   !    real(8) :: dPI, dt      
      
   !    ! --- Variables for calculate-scale-update logic
   !    real(8) :: Aq_LDOP_init, Aq_RDOP_init, Aq_PO4_init, Aq_PIPA_init
   !    real(8) :: Aq_PIPS_init
   !    real(8) :: TotalFlux_LDOP, TotalFlux_RDOP, TotalFlux_PO4
   !    real(8) :: TotalFlux_PIPA, TotalFlux_PIPS
   !    real(8) :: Scale_LDOP, Scale_RDOP, Scale_PO4, Scale_PIPA, Scale_PIPS

   !    dt = real(SedCNPmodel%SedCNP_timestep)
   !    dx = SedCNP_hydro%dx(1)    !BK20260301

   !    ! ====================================================================
   !    ! I. INITIALIZE STATE AND ADD INFLOWS
   !    ! ====================================================================

   !    Aq_LDOP_init = Aq_LDOP0(i,j) + So_LDOPsoaq0(i,j)
   !    Aq_RDOP_init = Aq_RDOP0(i,j) + So_RDOPsoaq0(i,j)
   !    Aq_PO4_init  = Aq_PO40(i,j)  + So_PO4soaq0(i,j)
   !    Aq_PIPA_init = Aq_PIPA0(i,j)
   !    Aq_PIPS_init = Aq_PIPS0(i,j)

   !    ! ====================================================================
   !    ! II. CALCULATE TRANSFORMATIONS & TRANSPORT FLUXES
   !    ! ====================================================================
      
   !    ! --- Transformations (Sorption only) ---
   !    KTC = EXP(0.115 * TAQUIFER - 2.88)
   !    KTHETA = 0.6 ! Assumed constant for saturated conditions
      
   !    ! adorption/desorption between PO4 and Active PIP (PIPA)
   !    if (Aq_PO4_init .ge. SIGMAPO4 * Aq_PIPA_init) then
   !       if (Aq_PIPA_init > 0.0) then
   !          beta = Aq_PO4_init / (SIGMAPO4 * Aq_PIPA_init)
   !       else
   !          beta = 100.0
   !       endif
   !       Aq_PO4ADS(i,j) = KTC * KTHETA * KADSPO4 * (1.0 - exp(-beta)) * &
   !          (Aq_PO4_init - SIGMAPO4 * Aq_PIPA_init) * dt
   !    else
   !       if (Aq_PO4_init > 0.0) then
   !          beta = (SIGMAPO4*Aq_PIPA_init) / Aq_PO4_init
   !       else
   !          beta = 100.0
   !       endif
   !       Aq_PO4ADS(i,j) = 0.1 * KTC * KTHETA * KADSPO4 * &
   !          (1.0 - exp(-beta)) * (Aq_PO4_init - SIGMAPO4 * Aq_PIPA_init) * dt
   !    endif
      
   !    ! adsorption/desorption between Active PIP (PIPA) and Stable PIP (PIPS)
   !    if (Aq_PIPA_init .ge. SIGMAPIPA * Aq_PIPS_init) then
   !       if (Aq_PIPS_init > 0.0) then
   !          beta = Aq_PIPA_init / (SIGMAPIPA * Aq_PIPS_init)
   !       else
   !          beta = 100.0
   !       endif
   !       Aq_PIPAADS(i,j) = KTC * KTHETA * KADSPIPA * &
   !          (1.0 - exp(-beta)) * (Aq_PIPA_init - SIGMAPIPA * Aq_PIPS_init) * dt
   !    else
   !       if(Aq_PIPA_init > 0.0) then
   !          beta = (SIGMAPIPA * Aq_PIPS_init) / Aq_PIPA_init
   !       else
   !          beta = 100.0
   !       endif
   !       Aq_PIPAADS(i,j) = 0.1 * KTC * KTHETA * KADSPIPA * &
   !          (1.0 - exp(-beta)) * (Aq_PIPA_init - SIGMAPIPA * Aq_PIPS_init) * dt
   !    endif
      
   !    ! --- Transport ---
   !    ! Transport for irrigation abstraction
   !    if (Aq_WStorage > (Qaqso * domain%areaxy / 1000.0 * dt)) then
   !       Aq_LDOPaqso(i,j) = (Qaqso * domain%areaxy / 1000.0 * dt) * &
   !          Aq_LDOP_init / Aq_WStorage
   !       Aq_RDOPaqso(i,j) = (Qaqso * domain%areaxy / 1000.0 * dt) * &
   !          Aq_RDOP_init / Aq_WStorage
   !       Aq_PO4aqso(i,j)  = (Qaqso * domain%areaxy / 1000.0 * dt) * &
   !          Aq_PO4_init / Aq_WStorage
   !    else
   !       Aq_LDOPaqso(i,j) = Aq_LDOP_init
   !       Aq_RDOPaqso(i,j) = Aq_RDOP_init
   !       Aq_PO4aqso(i,j)  = Aq_PO4_init
   !    endif
      
   !    ! Transport to groundwater bucket
   !    q_out = Qaqgw / 1000. * dx * dx * dt  ! (m3)
   !    !DOM
   !    dPI = calc_dPI(etaDOMaqgw, q_out, Aq_Wstorage)  !BK20260302
   !    Aq_LDOPaqgw(i,j) = Aq_LDOP_init*dPI
   !    Aq_RDOPaqgw(i,j) = Aq_RDOP_init*dPI
   !    !PO4
   !    dPI = calc_dPI(etaPO4aqgw, q_out, Aq_Wstorage)  !BK20260302
   !    Aq_PO4aqgw(i,j) = Aq_PO4_init*dPI

   !    ! ====================================================================
   !    ! III. SCALE FLUXES FOR MASS CONSERVATION
   !    ! ====================================================================

   !    TotalFlux_LDOP = Aq_LDOPaqso(i,j) + Aq_LDOPaqgw(i,j)
   !    if(TotalFlux_LDOP > Aq_LDOP_init .and. TotalFlux_LDOP > 0.0) then
   !       Scale_LDOP=Aq_LDOP_init/TotalFlux_LDOP
   !       Aq_LDOPaqso(i,j) = Aq_LDOPaqso(i,j) * Scale_LDOP
   !       Aq_LDOPaqgw(i,j) = Aq_LDOPaqgw(i,j) * Scale_LDOP
   !    endif
      
   !    TotalFlux_RDOP = Aq_RDOPaqso(i,j) + Aq_RDOPaqgw(i,j)
   !    if(TotalFlux_RDOP > Aq_RDOP_init .and. TotalFlux_RDOP > 0.0) then
   !       Scale_RDOP=Aq_RDOP_init/TotalFlux_RDOP
   !       Aq_RDOPaqso(i,j) = Aq_RDOPaqso(i,j) * Scale_RDOP
   !       Aq_RDOPaqgw(i,j) = Aq_RDOPaqgw(i,j) * Scale_RDOP
   !    endif
      
   !    ! Sum and scale PO4 sinks
   !    TotalFlux_PO4 = Aq_PO4aqso(i,j) + Aq_PO4aqgw(i,j)
   !    if(Aq_PO4ADS(i,j) > 0.0) &
   !       TotalFlux_PO4 = TotalFlux_PO4 + Aq_PO4ADS(i,j)
   !    if(TotalFlux_PO4 > Aq_PO4_init .and. TotalFlux_PO4 > 0.0) then
   !       Scale_PO4 = Aq_PO4_init / TotalFlux_PO4
   !       Aq_PO4aqso(i,j) = Aq_PO4aqso(i,j) * Scale_PO4
   !       Aq_PO4aqgw(i,j) = Aq_PO4aqgw(i,j) * Scale_PO4
   !       if(Aq_PO4ADS(i,j) > 0.0) &
   !          Aq_PO4ADS(i,j) = Aq_PO4ADS(i,j) * Scale_PO4
   !    endif
      
   !    ! Sum and scale PIPA sinks
   !    TotalFlux_PIPA = 0.0
   !    if(Aq_PO4ADS(i,j) < 0.0) &
   !       TotalFlux_PIPA = TotalFlux_PIPA - Aq_PO4ADS(i,j)
   !    if(Aq_PIPAADS(i,j) > 0.0) &
   !       TotalFlux_PIPA = TotalFlux_PIPA + Aq_PIPAADS(i,j)
   !    if(TotalFlux_PIPA > Aq_PIPA_init .and. TotalFlux_PIPA > 0.0) then
   !       Scale_PIPA = Aq_PIPA_init / TotalFlux_PIPA
   !       if(Aq_PO4ADS(i,j) < 0.0) &
   !          Aq_PO4ADS(i,j) = Aq_PO4ADS(i,j) * Scale_PIPA
   !       if(Aq_PIPAADS(i,j) > 0.0) &
   !          Aq_PIPAADS(i,j) = Aq_PIPAADS(i,j) * Scale_PIPA
   !    endif

   !    ! Sum and scale PIPS sinks
   !    TotalFlux_PIPS = 0.0
   !    if(Aq_PIPAADS(i,j) < 0.0) &
   !       TotalFlux_PIPS = TotalFlux_PIPS - Aq_PIPAADS(i,j)
   !    if(TotalFlux_PIPS > Aq_PIPS_init .and. TotalFlux_PIPS > 0.0) then
   !       Scale_PIPS = Aq_PIPS_init / TotalFlux_PIPS
   !       if(Aq_PIPAADS(i,j) < 0.0) &
   !          Aq_PIPAADS(i,j) = Aq_PIPAADS(i,j) * Scale_PIPS
   !    endif

   !    ! ====================================================================
   !    ! IV. UPDATE FINAL STORAGES
   !    ! ====================================================================

   !    Aq_LDOP(i,j) = Aq_LDOP_init - Aq_LDOPaqso(i,j) - Aq_LDOPaqgw(i,j)
   !    Aq_RDOP(i,j) = Aq_RDOP_init - Aq_RDOPaqso(i,j) - Aq_RDOPaqgw(i,j)
   !    Aq_PO4(i,j)  = Aq_PO4_init  - Aq_PO4aqso(i,j) - Aq_PO4aqgw(i,j) - &
   !                   Aq_PO4ADS(i,j)
   !    Aq_PIPA(i,j) = Aq_PIPA_init + Aq_PO4ADS(i,j) - Aq_PIPAADS(i,j)
   !    Aq_PIPS(i,j) = Aq_PIPS_init + Aq_PIPAADS(i,j)

   !    ! Enforce non-negative aquifer storages
   !    if (Aq_LDOP(i,j) < 1.0e-20) Aq_LDOP(i,j) = 0.0
   !    if (Aq_RDOP(i,j) < 1.0e-20) Aq_RDOP(i,j) = 0.0
   !    if (Aq_PO4(i,j) < 1.0e-20) Aq_PO4(i,j) = 0.0
   !    if (Aq_PIPA(i,j) < 1.0e-20) Aq_PIPA(i,j) = 0.0
   !    if (Aq_PIPS(i,j) < 1.0e-20) Aq_PIPS(i,j) = 0.0

   ! end subroutine AquiferPhosphorusCycle


   subroutine GroundwaterPhosphorusCycle(gwid, Gw_WStorage, Qgwch)

      implicit none
      
      ! --- Input arguments
      integer, intent(in)  :: gwid
      real(8), intent(in)  :: Gw_WStorage    ![m3]       !BK20260320
      real(8), intent(in)  :: Qgwch          ![mm s-1]   !BK20260320

      ! --- Local variables
      integer              :: i, j, istat
      real(8)              :: q_out          ![m3]       !BK20260302
      real(8)              :: dPI, dx, dt                !BK20260302
      real(8), allocatable :: Qgwso(:,:)     ![mm s-1]   !BK20260320
      real(8), allocatable :: Qgwso_bas(:)   ![m3]       !BK20260320
      
      ! --- Variables for calculate-scale-update logic
      real(8) :: Gw_LDOP_init, Gw_RDOP_init, Gw_PO4_init
      real(8) :: TotalFlux_LDOP, TotalFlux_RDOP, TotalFlux_PO4
      real(8) :: Scale_LDOP, Scale_RDOP, Scale_PO4
      real(8) :: Scale_gwso

      dt = real(SedCNPmodel%SedCNP_timestep)
      dx = SedCNP_hydro%dx(1)    !BK20260301

      allocate(Qgwso(domain%ix,domain%jx),stat=istat)    !BK20260320
      allocate(Qgwso_bas(domain%nbasin),stat=istat)      !BK20260320

      ! ====================================================================
      ! I. AGGREGATE INFLOWS FROM SOIL
      ! ====================================================================

      Gw_LDOPsogw0(gwid) = 0.0  !BK20260316
      Gw_RDOPsogw0(gwid) = 0.0  !BK20260316
      Gw_PO4sogw0(gwid)  = 0.0  !BK20260316
      
      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               Gw_LDOPsogw0(gwid) = Gw_LDOPsogw0(gwid) + So_LDOPsogw0(i,j)  !BK20260316
               Gw_RDOPsogw0(gwid) = Gw_RDOPsogw0(gwid) + So_RDOPsogw0(i,j)  !BK20260316
               Gw_PO4sogw0(gwid)  = Gw_PO4sogw0(gwid) + So_PO4sogw0(i,j)    !BK20260316
            END IF
         END DO
      END DO

      ! ====================================================================
      ! II. INITIALIZE STATE AND ADD INFLOWS
      ! ====================================================================

      Gw_LDOP_init = Gw_LDOP0(gwid) + Gw_LDOPsogw0(gwid)  !BK20260316
      Gw_RDOP_init = Gw_RDOP0(gwid) + Gw_RDOPsogw0(gwid)  !BK20260316
      Gw_PO4_init  = Gw_PO40(gwid)  + Gw_PO4sogw0(gwid)  !BK20260316

      ! ====================================================================
      ! III. CALCULATE TRANSPORT FLUXES 
      ! ====================================================================

      ! --- groundwater abstraction for irrigation !BK20260320
      Qgwso_bas(gwid) = 0.0
      Gw_LDOPgwso(gwid) = 0.0
      Gw_RDOPgwso(gwid) = 0.0
      Gw_PO4gwso(gwid) = 0.0

      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               ! Clamp irridepth to >= 0 to exclude fill values (e.g. -9999)
               ! Qgwso(i,j) = SedCNP_hydro%irridepth(i,j) / dt  ! [mm s-1]
               Qgwso(i,j) = max(0.0d0, SedCNP_hydro%irridepth(i,j)) / dt  ! [mm s-1]  !BK20260508
               Qgwso_bas(gwid) = Qgwso_bas(gwid) + Qgwso(i,j) / 1000.0 &
                     * domain%areaxy * dt  ! [m3]
            END IF
         END DO
      END DO

      ! if (Gw_WStorage > Qgwso_bas(gwid)) then
      if (Gw_WStorage > 0.0d0 .and. Gw_WStorage > Qgwso_bas(gwid)) then  !BK20260508
         Gw_LDOPgwso(gwid) = Gw_LDOP_init * Qgwso_bas(gwid) / Gw_WStorage
         Gw_RDOPgwso(gwid) = Gw_RDOP_init * Qgwso_bas(gwid) / Gw_WStorage
         Gw_PO4gwso(gwid) = Gw_PO4_init * Qgwso_bas(gwid) / Gw_WStorage
      ! else
      else if (Gw_WStorage > 0.0d0) then  !BK20260508
         Gw_LDOPgwso(gwid) = Gw_LDOP_init
         Gw_RDOPgwso(gwid) = Gw_RDOP_init
         Gw_PO4gwso(gwid) = Gw_PO4_init
      else                         !---BK20260508
         Gw_LDOPgwso(gwid) = 0.0d0
         Gw_RDOPgwso(gwid) = 0.0d0
         Gw_PO4gwso(gwid) = 0.0d0  !---BK20260508
      endif

      ! --- groundwater discharge to channel      
      q_out = Qgwch / 1000. * dx * dx * dt  ! (m3)
      !DOM
      dPI = calc_dPI(etaDOMgwch, q_out, Gw_WStorage)  !BK20260302
      Gw_LDOPgwch(gwid) = Gw_LDOP_init * dPI
      Gw_RDOPgwch(gwid) = Gw_RDOP_init * dPI
      !PO4
      dPI = calc_dPI(etaPO4gwch, q_out, Gw_WStorage)  !BK20260302
      Gw_PO4gwch(gwid) = Gw_PO4_init * dPI

      ! ====================================================================
      ! IV. SCALE FLUXES TO ENSURE MASS CONSERVATION
      ! ====================================================================

      TotalFlux_LDOP = Gw_LDOPgwso(gwid) + Gw_LDOPgwch(gwid)
      if (TotalFlux_LDOP > Gw_LDOP_init .and. TotalFlux_LDOP > 0.0) then
         Scale_LDOP = Gw_LDOP_init / TotalFlux_LDOP
         Gw_LDOPgwso(gwid) = Gw_LDOPgwso(gwid) * Scale_LDOP
         Gw_LDOPgwch(gwid) = Gw_LDOPgwch(gwid) * Scale_LDOP
      endif
      
      TotalFlux_RDOP = Gw_RDOPgwso(gwid) + Gw_RDOPgwch(gwid)
      if (TotalFlux_RDOP > Gw_RDOP_init .and. TotalFlux_RDOP > 0.0) then
         Scale_RDOP = Gw_RDOP_init / TotalFlux_RDOP
         Gw_RDOPgwso(gwid) = Gw_RDOPgwso(gwid) * Scale_RDOP
         Gw_RDOPgwch(gwid) = Gw_RDOPgwch(gwid) * Scale_RDOP
      endif
     
      TotalFlux_PO4 = Gw_PO4gwso(gwid) + Gw_PO4gwch(gwid)
      if (TotalFlux_PO4 > Gw_PO4_init .and. TotalFlux_PO4 > 0.0) then
         Scale_PO4 = Gw_PO4_init / TotalFlux_PO4
         Gw_PO4gwso(gwid) = Gw_PO4gwso(gwid) * Scale_PO4
         Gw_PO4gwch(gwid) = Gw_PO4gwch(gwid) * Scale_PO4
      endif

      ! --- scale flux *gwso(i,j) at the grid-level  !BK20260316
      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               if (Qgwso_bas(gwid) > 1.0e-9) then
                  Scale_gwso = (Qgwso(i,j) / 1000.0 * domain%areaxy * dt) / Qgwso_bas(gwid)
               else
                  Scale_gwso = 0.0
               endif
               So_LDOPgwso(i,j) = Gw_LDOPgwso(gwid) * Scale_gwso
               So_RDOPgwso(i,j) = Gw_RDOPgwso(gwid) * Scale_gwso
               So_PO4gwso(i,j) = Gw_PO4gwso(gwid) * Scale_gwso
            END IF
         END DO
      END DO

      ! ====================================================================
      ! V. UPDATE FINAL STORAGES
      ! ====================================================================
      Gw_LDOP(gwid) = Gw_LDOP_init - Gw_LDOPgwso(gwid) - Gw_LDOPgwch(gwid)
      Gw_RDOP(gwid) = Gw_RDOP_init - Gw_RDOPgwso(gwid) - Gw_RDOPgwch(gwid)
      Gw_PO4(gwid)  = Gw_PO4_init  - Gw_PO4gwso(gwid) - Gw_PO4gwch(gwid)

      ! Enforce non-negative groundwater storages
      if (Gw_LDOP(gwid) < 1.0e-20) Gw_LDOP(gwid) = 0.0
      if (Gw_RDOP(gwid) < 1.0e-20) Gw_RDOP(gwid) = 0.0
      if (Gw_PO4(gwid) < 1.0e-20) Gw_PO4(gwid) = 0.0

   end subroutine GroundwaterPhosphorusCycle
   

   subroutine ChannelPhosphorusCycle(ich,itime,Qdsch,WtopWdth,WStorage)
    
      implicit none
      
      ! --- Input arguments
      integer, intent(in) :: ich, itime
      real(8), intent(in)  :: Qdsch, WtopWdth, WStorage
      
      ! --- Local variables
      integer :: i, j, isl, ips, ia, im, month, stat, chid, gwid
      real(8) :: rate0, dW, dZ, dXi1, dXi2, dLambda_L, dLambda_N, dLambda_P
      real(8) :: dPHI_N, dPHI_P, dLPOPD, dLPOPZ, dLPOPS, dLDOPD, dRPOPD
      real(8) :: dRPOPS, dRDOPD, dMBMPD, dAZP, dALGCconc, dTSSconc
      real(8) :: dTA(ntm), dMA(ntm), dt, abs_ratio
      real(8) :: RK4SOL, KTBMBM, KTBZOO, KTBALG, KTHETA, ALPHA, dPI, KTC, beta
      real(8) :: KZG(domain%nch), KZR(domain%nch), KZM(domain%nch)
      real(8) :: KAG(domain%nch,nalg), KAR(domain%nch,nalg)
      real(8) :: TotalSuspendedSSA  !BK20260213
      real(8) :: KAE(domain%nch,nalg), KAM(domain%nch,nalg)
      real(8), PARAMETER :: e = 2.718281828459
      character(len=2) :: monthstr

      ! --- Variables for calculate-scale-update logic
      real(8) :: Ch_ALGP_init(nalg), Ch_ZOOP_init, Ch_LPOP_init
      real(8) :: Ch_LDOP_init, Ch_RPOP_init, Ch_RDOP_init, Ch_MBMP_init
      real(8) :: Ch_PO4_init, Ch_PIPA_init, Ch_PIPS_init
      real(8) :: TotalFlux_ALGP(nalg), Scale_ALGP(nalg)
      real(8) :: TotalFlux_ZOOP, Scale_ZOOP
      real(8) :: TotalFlux_LPOP, Scale_LPOP
      real(8) :: TotalFlux_LDOP, Scale_LDOP
      real(8) :: TotalFlux_RPOP, Scale_RPOP
      real(8) :: TotalFlux_RDOP, Scale_RDOP
      real(8) :: TotalFlux_MBMP, Scale_MBMP
      real(8) :: TotalFlux_PO4, Scale_PO4
      real(8) :: TotalFlux_PIPA, Scale_PIPA
      real(8) :: TotalFlux_PIPS, Scale_PIPS

      dt = real(SedCNPmodel%SedCNP_timestep)
      chid = SedCNP_hydro%linkID(ich)
      gwid = domain%gwid_ch(chid)

      ! ========================================================================
      ! I. INFLOWS
      ! ========================================================================
      
      ! --- Inflow from upstream channels
      Ch_ALGP(ich,1:nalg) = Ch_ALGP0(ich,1:nalg) + Ch_ALGPusch0(ich,1:nalg)
      Ch_ZOOP(ich) = Ch_ZOOP0(ich) + Ch_ZOOPusch0(ich)
      Ch_MBMP(ich) = Ch_MBMP0(ich) + Ch_MBMPusch0(ich)
      Ch_LPOP(ich) = Ch_LPOP0(ich) + Ch_LPOPusch0(ich)
      Ch_RPOP(ich) = Ch_RPOP0(ich) + Ch_RPOPusch0(ich)
      Ch_LDOP(ich) = Ch_LDOP0(ich) + Ch_LDOPusch0(ich)
      Ch_RDOP(ich) = Ch_RDOP0(ich) + Ch_RDOPusch0(ich)
      Ch_PO4(ich)  = Ch_PO40(ich)  + Ch_PO4usch0(ich)
      Ch_PIPA(ich) = Ch_PIPA0(ich) + Ch_PIPAusch0(ich)
      Ch_PIPS(ich) = Ch_PIPS0(ich) + Ch_PIPSusch0(ich)

      ! --- Lateral inflows (surface, interflow, gw) ---
      IF (gwid >= 0) THEN
         ! Initialisation
         Ch_LPOPsurf0(ich) = 0.0
         Ch_RPOPsurf0(ich) = 0.0
         Ch_LDOPsurf0(ich) = 0.0
         Ch_RDOPsurf0(ich) = 0.0
         Ch_MBMPsurf0(ich) = 0.0
         Ch_PO4surf0(ich)  = 0.0
         Ch_PIPAsurf0(ich) = 0.0
         Ch_PIPSsurf0(ich) = 0.0
         Ch_LDOPintf0(ich) = 0.0
         Ch_RDOPintf0(ich) = 0.0
         Ch_PO4intf0(ich)  = 0.0
         
         DO i = 1, domain%ix
            DO j = 1, domain%jx
               IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN

                  ! surface runoff
                  Ch_LPOPsurf0(ich) = Ch_LPOPsurf0(ich) + So_LPOPsurf0(i,j)
                  Ch_RPOPsurf0(ich) = Ch_RPOPsurf0(ich) + So_RPOPsurf0(i,j)
                  Ch_LDOPsurf0(ich) = Ch_LDOPsurf0(ich) + So_LDOPsurf0(i,j)
                  Ch_RDOPsurf0(ich) = Ch_RDOPsurf0(ich) + So_RDOPsurf0(i,j)
                  Ch_MBMPsurf0(ich) = Ch_MBMPsurf0(ich) + So_MBMPsurf0(i,j)
                  Ch_PO4surf0(ich)  = Ch_PO4surf0(ich)  + So_PO4surf0(i,j)
                  Ch_PIPAsurf0(ich) = Ch_PIPAsurf0(ich) + So_PIPAsurf0(i,j)
                  Ch_PIPSsurf0(ich) = Ch_PIPSsurf0(ich) + So_PIPSsurf0(i,j)

                  ! interflow
                  DO isl = 1, domain%nsl
                     Ch_LDOPintf0(ich) = Ch_LDOPintf0(ich) + &
                        So_LDOPintf0(i,isl,j)
                     Ch_RDOPintf0(ich) = Ch_RDOPintf0(ich) + &
                        So_RDOPintf0(i,isl,j)
                     Ch_PO4intf0(ich)  = Ch_PO4intf0(ich) + &
                        So_PO4intf0(i,isl,j)
                  END DO

               END IF
            END DO
         END DO
         
         ! groundwater inflow
         Ch_LDOPgwch0(ich) = Gw_LDOPgwch0(gwid)
         Ch_RDOPgwch0(ich) = Gw_RDOPgwch0(gwid)
         Ch_PO4gwch0(ich)  = Gw_PO4gwch0(gwid)

         ! update storage
         Ch_LPOP(ich) = Ch_LPOP(ich) + Ch_LPOPsurf0(ich)
         Ch_RPOP(ich) = Ch_RPOP(ich) + Ch_RPOPsurf0(ich)
         Ch_LDOP(ich) = Ch_LDOP(ich) + Ch_LDOPsurf0(ich) + Ch_LDOPintf0(ich) + &
            Ch_LDOPgwch0(ich)
         Ch_RDOP(ich) = Ch_RDOP(ich) + Ch_RDOPsurf0(ich) + Ch_RDOPintf0(ich) + &
            Ch_RDOPgwch0(ich)
         Ch_MBMP(ich) = Ch_MBMP(ich) + Ch_MBMPsurf0(ich)
         Ch_PO4(ich)  = Ch_PO4(ich)  + Ch_PO4surf0(ich)  + Ch_PO4intf0(ich) + &
            Ch_PO4gwch0(ich)
         Ch_PIPA(ich) = Ch_PIPA(ich) + Ch_PIPAsurf0(ich)
         Ch_PIPS(ich) = Ch_PIPS(ich) + Ch_PIPSsurf0(ich)

      ENDIF

      ! --- Forced inflow/ourflow 
      ! inflow: point source
      Ch_LPOPxpnt(ich) = Ch_LPOPpnt(ich,itime) * dt
      Ch_RPOPxpnt(ich) = Ch_RPOPpnt(ich,itime) * dt
      Ch_LDOPxpnt(ich) = Ch_LDOPpnt(ich,itime) * dt
      Ch_RDOPxpnt(ich) = Ch_RDOPpnt(ich,itime) * dt
      Ch_PO4xpnt(ich)  = Ch_PO4pnt(ich,itime)  * dt

      ! inflow: water discharge
      Ch_LPOPxdis(ich) = Ch_LPOPdis(ich,itime) * dt
      Ch_RPOPxdis(ich) = Ch_RPOPdis(ich,itime) * dt
      Ch_LDOPxdis(ich) = Ch_LDOPdis(ich,itime) * dt
      Ch_RDOPxdis(ich) = Ch_RDOPdis(ich,itime) * dt
      Ch_PO4xdis(ich)  = Ch_PO4dis(ich,itime)  * dt
 
      ! update storage (inflows)
      Ch_LPOP(ich) = Ch_LPOP(ich) + Ch_LPOPxpnt(ich) + Ch_LPOPxdis(ich)
      Ch_RPOP(ich) = Ch_RPOP(ich) + Ch_RPOPxpnt(ich) + Ch_RPOPxdis(ich)
      Ch_LDOP(ich) = Ch_LDOP(ich) + Ch_LDOPxpnt(ich) + Ch_LDOPxdis(ich)
      Ch_RDOP(ich) = Ch_RDOP(ich) + Ch_RDOPxpnt(ich) + Ch_RDOPxdis(ich)
      Ch_PO4(ich)  = Ch_PO4(ich)  + Ch_PO4xpnt(ich)  + Ch_PO4xdis(ich)

      ! outflow: water abstraction
      abs_ratio = SedCNP_hydro%Qabs(ich,itime) / (WStorage + Qdsch * dt)
      if(SedCNP_hydro%Qabs(ich,itime) .gt. 0.) then
         Ch_LPOPabs(ich,itime) = Ch_LPOP(ich) * abs_ratio
         Ch_RPOPabs(ich,itime) = Ch_RPOP(ich) * abs_ratio
         Ch_LDOPabs(ich,itime) = Ch_LDOP(ich) * abs_ratio
         Ch_RDOPabs(ich,itime) = Ch_RDOP(ich) * abs_ratio
         Ch_PO4abs(ich,itime)  = Ch_PO4(ich)  * abs_ratio
      endif
      Ch_LPOPxabs(ich) = Ch_LPOPabs(ich,itime) * dt  
      Ch_RPOPxabs(ich) = Ch_RPOPabs(ich,itime) * dt  
      Ch_LDOPxabs(ich) = Ch_LDOPabs(ich,itime) * dt  
      Ch_RDOPxabs(ich) = Ch_RDOPabs(ich,itime) * dt  
      Ch_PO4xabs(ich)  = Ch_PO4abs(ich,itime)  * dt

      ! update storage (outflows)
      Ch_LPOP(ich) = Ch_LPOP(ich) - Ch_LPOPxabs(ich)
      Ch_RPOP(ich) = Ch_RPOP(ich) - Ch_RPOPxabs(ich)
      Ch_LDOP(ich) = Ch_LDOP(ich) - Ch_LDOPxabs(ich)
      Ch_RDOP(ich) = Ch_RDOP(ich) - Ch_RDOPxabs(ich)
      Ch_PO4(ich)  = Ch_PO4(ich)  - Ch_PO4xabs(ich)

      ! Enforce non-negative channel storages after abstraction
      if (Ch_LPOP(ich) < 1.0e-20) Ch_LPOP(ich) = 0.0
      if (Ch_RPOP(ich) < 1.0e-20) Ch_RPOP(ich) = 0.0
      if (Ch_LDOP(ich) < 1.0e-20) Ch_LDOP(ich) = 0.0
      if (Ch_RDOP(ich) < 1.0e-20) Ch_RDOP(ich) = 0.0
      if (Ch_PO4(ich) < 1.0e-20)  Ch_PO4(ich)  = 0.0

      ! ========================================================================
      ! II. CALCULATE-SCALE-UPDATE 
      ! ========================================================================

      ! --- Store initial values for simultaneous flux calculation
      Ch_ALGP_init(:) = Ch_ALGP(ich,:)
      Ch_ZOOP_init    = Ch_ZOOP(ich)
      Ch_LPOP_init    = Ch_LPOP(ich)
      Ch_LDOP_init    = Ch_LDOP(ich)
      Ch_RPOP_init    = Ch_RPOP(ich)
      Ch_RDOP_init    = Ch_RDOP(ich)
      Ch_MBMP_init    = Ch_MBMP(ich)
      Ch_PO4_init     = Ch_PO4(ich)
      Ch_PIPA_init    = Ch_PIPA(ich)
      Ch_PIPS_init    = Ch_PIPS(ich)

      ! === A. CALCULATE ALL FLUXES BASED ON INITIAL STORAGES

      ! --- Pre-calculations for transformation rates
      KTC = EXP(0.115 * TWATER(ich) - 2.88)
      KTHETA = 0.6 ! Assumed constant for in-stream saturated conditions
      dW = DWATER(ich)
      dZ = min(dW, 10.0)  ! take account of algaes in the top 10 m of water
      
      ! ALPHA - light attenuation coeff. (USEPA, 2000)
      if (WStorage > 0.0) then
         dALGCconc = sum(Ch_ALGC0(ich,1:nalg)) / WStorage
         dTSSconc = (sum(channelSed%Cchsd(1:nps,ich)) + &
                     Ch_LPOP_init + Ch_RPOP_init) / WStorage
      else
         dALGCconc = 0.0
         dTSSconc = 0.0
      endif
      ALPHA = ALPHADOC + ALPHACHL * dALGCconc + ALPHATSS * dTSSconc

      ! limiting factors for light, nitrogen, and phosphorus
      monthstr = trim(dateSedCNP%olddate(6:7))
      read(monthstr,*) month

      if (gwid > 0) then
         dXi1 = (1.0 - ALBEDO(month) / 100.0) * &
            SedCNP_hydro%SWDOWNavg(gwid) * exp(-ALPHA * dW) / SOLRADMAX
         dXi2 = (1.0 - ALBEDO(month) / 100.0) * &
            SedCNP_hydro%SWDOWNavg(gwid) * exp(-ALPHA * (dW + dZ)) / SOLRADMAX
      else
         dXi1 = 0.0
         dXi2 = 0.0
      end if

      if ((ALPHA * dZ) > 1e-6) then
         dLambda_L = e / (ALPHA * dZ) * (exp(-dXi2) - exp(-dXi1))
      else
         dLambda_L = 0.0
      end if
      
      if (WStorage > 0.0) then
         dPHI_N = (Ch_NH40(ich) + Ch_NO30(ich)) / WStorage
         dPHI_P = Ch_PO4_init / WStorage
      else
         dPHI_N = 0.0
         dPHI_P = 0.0
      end if

      if ((MONOD_N + dPHI_N) > 0.0) then
         dLambda_N = dPHI_N / (MONOD_N + dPHI_N)
      else
         dLambda_N = 0.0
      end if

      if ((MONOD_P + dPHI_P) > 0.0) then
         dLambda_P = dPHI_P / (MONOD_P + dPHI_P)
      else
         dLambda_P = 0.0
      end if
      
      ! --- Zooplankton 
      KTBZOO = calc_KTB(TWATER(ich), TZOO, MZOO)
      dAZP = sum(Ch_ALGP_init(:)) + Ch_LPOP_init + Ch_ZOOP_init

      ! growth
      if(dAZP + ZHALF * WStorage > 1.0E-9) then
         KZG(ich) = KTBZOO * EZI * KZIMAX * (dAZP - ZLOW * WStorage) / &
                     (dAZP + ZHALF * WStorage)
      else
         KZG(ich) = 0.0
      end if
      RK4SOL = calc_RK4SOL(KZG(ich), Ch_ZOOP_init)
      Ch_ZOOPG(ich) = max(0.0, RK4SOL - Ch_ZOOP_init)
      
      ! being grazed by (other) zooplankton
      if(EZI > 0.0 .and. dAZP > 1.0E-9) then
         Ch_ZOOPZ(ich) = (Ch_ZOOPG(ich) / EZI) * (Ch_ZOOP_init / dAZP)
      else
         Ch_ZOOPZ(ich) = 0.0
      end if
      
      ! respiration
      KZR(ich) = KTBZOO * KZRMAX
      RK4SOL = calc_RK4SOL(-KZR(ich),Ch_ZOOP_init)
      Ch_ZOOPR(ich) = Ch_ZOOP_init - RK4SOL
      
      ! excretion - based on growth to ensure mass balance (BK20251127) 
      if(EZI > 0.0) then
         Ch_ZOOPE(ich) = Ch_ZOOPG(ich) * (1.0 - EZI) / EZI
      else
         Ch_ZOOPE(ich) = 0.0
      end if
      
      ! mortality
      KZM(ich) = KTBZOO * KZMMAX
      RK4SOL = calc_RK4SOL(-KZM(ich),Ch_ZOOP_init)
      Ch_ZOOPM(ich) = Ch_ZOOP_init - RK4SOL
      
      ! settling
      if(dW > 0.0) then
         Ch_ZOOPS(ich) = (1.0 - exp(-OMEGAZOO * dt / dW)) * Ch_ZOOP_init
      else
         Ch_ZOOPS(ich) = 0.0
      end if
      
      ! --- Phytoplanktons
      DO ia = 1, nalg
         do im = 1, ntm
            dTA(im) = TALG(ia,im)
            dMA(im) = MALG(ia,im)
         enddo
         KTBALG = calc_KTB(TWATER(ich), dTA, dMA)

         ! growth
         KAG(ich,ia) = KTBALG * min(dLambda_L, dLambda_N, dLambda_P) * KAGMAX(ia)
         RK4SOL = calc_RK4SOL(KAG(ich,ia), Ch_ALGP_init(ia))
         Ch_ALGPG(ich,ia) = max(0.0, RK4SOL - Ch_ALGP_init(ia))
         
         ! being grazed by zooplankton
         if(EZI > 0.0 .and. dAZP > 1.0E-9) then
            Ch_ALGPZ(ich,ia) = &
               (Ch_ZOOPG(ich) / EZI) * (Ch_ALGP_init(ia) / dAZP)
         else
            Ch_ALGPZ(ich,ia) = 0.0
         end if
         
         ! respiration
         KAR(ich,ia) = KTBALG * KARMAX(ia)
         RK4SOL = calc_RK4SOL(-KAR(ich,ia), Ch_ALGP_init(ia))
         Ch_ALGPR(ich,ia) = Ch_ALGP_init(ia) - RK4SOL
         
         ! excretion
         KAE(ich,ia) = KTBALG *(1.0 - dLambda_L) * KAEMAX(ia)
         RK4SOL = calc_RK4SOL(-KAE(ich,ia), Ch_ALGP_init(ia))
         Ch_ALGPE(ich,ia) = Ch_ALGP_init(ia) - RK4SOL
         
         ! mortality
         KAM(ich,ia) = KTBALG * KAMMAX(ia)
         RK4SOL = calc_RK4SOL(-KAM(ich,ia), Ch_ALGP_init(ia))
         Ch_ALGPM(ich,ia) = Ch_ALGP_init(ia) - RK4SOL
         
         ! settling
         if(dW > 0.0) then
            Ch_ALGPS(ich,ia) = (1.0 - exp(-OMEGAALG * dt / dW)) * &
               Ch_ALGP_init(ia)
         else
            Ch_ALGPS(ich,ia) = 0.0
         end if
      END DO
      
      ! --- Organic Matter 
      KTBMBM = calc_KTB(TWATER(ich), TMBM, MMBM)

      ! LPOP grazed by zooplankton
      if(EZI > 0.0 .and. dAZP > 1.0E-9) then
         dLPOPZ = (Ch_ZOOPG(ich) / EZI) * (Ch_LPOP_init / dAZP)
      else
         dLPOPZ = 0.0
      end if
      
      ! LPOP decomposition
      rate0 = KTBMBM * KTHETA * KLPOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_LPOP_init)
      dLPOPD = Ch_LPOP_init - RK4SOL
      
      ! LPOP settling
      if(dW > 0.0) then
         dLPOPS = (1.0 - exp(-OMEGAPOM * dt / dW)) * Ch_LPOP_init
      else
         dLPOPS = 0.0
      end if
      
      ! RPOP decomposition
      rate0 = KTBMBM * KTHETA * KRPOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_RPOP_init)
      dRPOPD = Ch_RPOP_init - RK4SOL
      
      ! RPOP settling
      if(dW>0.0) then
         dRPOPS=(1.-exp(-OMEGAPOM*dt/dW))*Ch_RPOP_init
      else
         dRPOPS=0.0
      end if

      ! LDOP decomposition
      rate0 = KTBMBM * KTHETA * KLDOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_LDOP_init)
      dLDOPD = Ch_LDOP_init - RK4SOL
      
      ! RDOP decomposition
      rate0 = KTBMBM * KTHETA * KRDOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_RDOP_init)
      dRDOPD = Ch_RDOP_init - RK4SOL
      
      ! MBMP decomposition
      rate0 = KTBMBM * KTHETA * KMBMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_MBMP_init)
      dMBMPD = Ch_MBMP_init - RK4SOL

      ! --- Inorganic P
      ! uptake by phytoplanktons
      Ch_PO4uptk(ich) = sum(Ch_ALGPG(ich,:))
      
      ! PO4 sorption
      if(Ch_PO4_init >= SIGMAPO4 * Ch_PIPA_init) then
         if(Ch_PIPA_init > 0.0) then
            beta=Ch_PO4_init / (SIGMAPO4 * Ch_PIPA_init)
         else
            beta = 100.0
         endif
         Ch_PO4ADS(ich) = KTC * KTHETA * KADSPO4 * (1.0 - exp(-beta)) * &
            (Ch_PO4_init - SIGMAPO4 * Ch_PIPA_init) * dt
      else
         if(Ch_PO4_init > 0.0) then
            beta=(SIGMAPO4 * Ch_PIPA_init) / Ch_PO4_init
         else
            beta = 100.0
         endif
         Ch_PO4ADS(ich) = 0.1 * KTC * KTHETA * KADSPO4 * (1.0 - exp(-beta)) * &
            (Ch_PO4_init - SIGMAPO4 * Ch_PIPA_init) * dt
      endif
      
      ! PIPA sorption
      if(Ch_PIPA_init .ge. SIGMAPIPA * Ch_PIPS_init) then
         if(Ch_PIPS_init > 0.0) then
            beta = Ch_PIPA_init / (SIGMAPIPA * Ch_PIPS_init)
         else
            beta = 100.0
         endif
         Ch_PIPAADS(ich) = KTC * KTHETA * KADSPIPA * (1.0 - exp(-beta)) * &
            (Ch_PIPA_init - SIGMAPIPA * Ch_PIPS_init) * dt
      else
         if(Ch_PIPA_init > 0.0) then
            beta = (SIGMAPIPA * Ch_PIPS_init) / Ch_PIPA_init
         else
            beta = 100.0
         endif
         Ch_PIPAADS(ich) = 0.1 * KTC * KTHETA * KADSPIPA * &
            (1.0 - exp(-beta)) * (Ch_PIPA_init - SIGMAPIPA * Ch_PIPS_init) * dt
      endif
      
      ! --- outflow: downstream discharge
      if(WStorage > 0.0 .and. Qdsch > 0.0) then
         dPI = 1.0 - EXP(-2.3 * (Qdsch * dt) / (WStorage + Qdsch * dt))
      else
         dPI = 0.0
      end if
      Ch_ALGPdsch(ich,:) = Ch_ALGP_init(:) * dPI
      Ch_ZOOPdsch(ich)   = Ch_ZOOP_init * dPI
      Ch_LPOPdsch(ich)   = Ch_LPOP_init * dPI
      Ch_RPOPdsch(ich)   = Ch_RPOP_init * dPI
      Ch_MBMPdsch(ich)   = Ch_MBMP_init * dPI
      Ch_LDOPdsch(ich)   = Ch_LDOP_init * dPI
      Ch_RDOPdsch(ich)   = Ch_RDOP_init * dPI
      Ch_PO4dsch(ich)    = Ch_PO4_init * dPI
      
      ! Particulate Inorganic P transport coupled with sediment transport.   !---BK20260213
      ! This logic apportions total PIPA/PIPS based on the specific surface
      ! area (SSA) of the suspended sediment particles.
      TotalSuspendedSSA = 0.0
      do ips = 1, nps
         TotalSuspendedSSA = TotalSuspendedSSA + (channelSed%Scw(ips, ich) * SSA(ips))
      enddo

      Ch_PIPAdsch(ich) = 0.0
      Ch_PIPSdsch(ich) = 0.0
      if (TotalSuspendedSSA > 1.0e-20) then
         do ips = 1, nps
            ! The flux of P is calculated by assuming P is uniformly distributed
            ! per unit of surface area. The total flux is the concentration of P
            ! per unit area multiplied by the total surface area of sediment
            ! being discharged.
            Ch_PIPAdsch(ich) = Ch_PIPAdsch(ich) + Ch_PIPA_init * &
               (SSA(ips) / TotalSuspendedSSA) * (channelSed%Sdsch(ips,ich) * dt)
            Ch_PIPSdsch(ich) = Ch_PIPSdsch(ich) + Ch_PIPS_init * &
               (SSA(ips) / TotalSuspendedSSA) * (channelSed%Sdsch(ips,ich) * dt)
         enddo
      endif      !---BK20260213

      ! === B. SCALE FLUXES TO ENSURE MASS CONSERVATION

      ! Zooplankton
      TotalFlux_ZOOP = Ch_ZOOPZ(ich) + Ch_ZOOPR(ich) + Ch_ZOOPM(ich) + &
         Ch_ZOOPS(ich) + Ch_ZOOPdsch(ich)
      if(TotalFlux_ZOOP > Ch_ZOOP_init .and. TotalFlux_ZOOP > 0.0) then
         Scale_ZOOP = Ch_ZOOP_init / TotalFlux_ZOOP
      else
         Scale_ZOOP = 1.0
      endif
      Ch_ZOOPZ(ich)    = Ch_ZOOPZ(ich) * Scale_ZOOP
      Ch_ZOOPR(ich)    = Ch_ZOOPR(ich) * Scale_ZOOP
      Ch_ZOOPM(ich)    = Ch_ZOOPM(ich) * Scale_ZOOP
      Ch_ZOOPS(ich)    = Ch_ZOOPS(ich) * Scale_ZOOP
      Ch_ZOOPdsch(ich) = Ch_ZOOPdsch(ich) * Scale_ZOOP
      
      ! LPOP
      TotalFlux_LPOP = dLPOPZ + dLPOPD + dLPOPS + Ch_LPOPdsch(ich)
      if(TotalFlux_LPOP > Ch_LPOP_init .and. TotalFlux_LPOP > 0.0) then
         Scale_LPOP = Ch_LPOP_init / TotalFlux_LPOP
      else
         Scale_LPOP = 1.0
      endif
      dLPOPZ           = dLPOPZ * Scale_LPOP
      dLPOPD           = dLPOPD * Scale_LPOP
      dLPOPS           = dLPOPS * Scale_LPOP
      Ch_LPOPdsch(ich) = Ch_LPOPdsch(ich) * Scale_LPOP

      ! LDOP
      TotalFlux_LDOP = dLDOPD + Ch_LDOPdsch(ich)
      if(TotalFlux_LDOP > Ch_LDOP_init .and. TotalFlux_LDOP > 0.0) then
         Scale_LDOP = Ch_LDOP_init / TotalFlux_LDOP
      else
         Scale_LDOP = 1.0
      endif
      dLDOPD = dLDOPD * Scale_LDOP
      Ch_LDOPdsch(ich) = Ch_LDOPdsch(ich) * Scale_LDOP
      
      ! RPOP
      TotalFlux_RPOP = dRPOPD + dRPOPS + Ch_RPOPdsch(ich)
      if(TotalFlux_RPOP > Ch_RPOP_init .and. TotalFlux_RPOP > 0.0) then
         Scale_RPOP = Ch_RPOP_init / TotalFlux_RPOP
      else
         Scale_RPOP = 1.0
      endif
      dRPOPD = dRPOPD * Scale_RPOP
      dRPOPS = dRPOPS * Scale_RPOP
      Ch_RPOPdsch(ich) = Ch_RPOPdsch(ich) * Scale_RPOP
      
      ! RDOP
      TotalFlux_RDOP = dRDOPD + Ch_RDOPdsch(ich)
      if(TotalFlux_RDOP > Ch_RDOP_init .and. TotalFlux_RDOP > 0.0) then
         Scale_RDOP = Ch_RDOP_init / TotalFlux_RDOP
      else
         Scale_RDOP = 1.0
      endif
      dRDOPD = dRDOPD * Scale_RDOP
      Ch_RDOPdsch(ich) = Ch_RDOPdsch(ich) * Scale_RDOP
      
      ! MBMP
      TotalFlux_MBMP = dMBMPD + Ch_MBMPdsch(ich)
      if(TotalFlux_MBMP > Ch_MBMP_init .and. TotalFlux_MBMP > 0.0) then
         Scale_MBMP = Ch_MBMP_init / TotalFlux_MBMP
      else
         Scale_MBMP = 1.0
      endif
      dMBMPD = dMBMPD * Scale_MBMP
      Ch_MBMPdsch(ich) = Ch_MBMPdsch(ich) * Scale_MBMP

      ! Scale inorganic P fluxes
      if(Ch_PO4ADS(ich) < 0.0) &
         Ch_PO4ADS(ich) = min(abs(Ch_PO4ADS(ich)), Ch_PIPA_init) * (-1.0)
      if(Ch_PIPAADS(ich) < 0.0) &
         Ch_PIPAADS(ich) = min(abs(Ch_PIPAADS(ich)), Ch_PIPS_init) * (-1.0)
      
      ! PO4   
      TotalFlux_PO4 = Ch_PO4uptk(ich) + Ch_PO4dsch(ich)
      if(Ch_PO4ADS(ich) > 0.0) TotalFlux_PO4 = TotalFlux_PO4 + Ch_PO4ADS(ich)
      if(TotalFlux_PO4 > Ch_PO4_init .and. TotalFlux_PO4 > 0.0) then
         Scale_PO4 = Ch_PO4_init / TotalFlux_PO4
      else
         Scale_PO4 = 1.0
      endif
      Ch_PO4uptk(ich) = Ch_PO4uptk(ich) * Scale_PO4
      Ch_PO4dsch(ich) = Ch_PO4dsch(ich) * Scale_PO4
      if(Ch_PO4ADS(ich) > 0.0) Ch_PO4ADS(ich) = Ch_PO4ADS(ich) * Scale_PO4

      ! Scale phytoplankton growth proportionally to nutrient uptake
      Ch_ALGPG(ich,:) = Ch_ALGPG(ich,:) * Scale_PO4
      
      ! Phytoplanktons
      DO ia=1,nalg
         TotalFlux_ALGP(ia) = Ch_ALGPZ(ich,ia) + Ch_ALGPR(ich,ia) + &
            Ch_ALGPE(ich,ia) + Ch_ALGPM(ich,ia) + Ch_ALGPS(ich,ia) + &
            Ch_ALGPdsch(ich,ia)
         if (TotalFlux_ALGP(ia) > (Ch_ALGP_init(ia) + Ch_ALGPG(ich,ia))) then
            if ((Ch_ALGP_init(ia) + Ch_ALGPG(ich,ia)) > 1.0E-20) then
               Scale_ALGP(ia) = (Ch_ALGP_init(ia) + Ch_ALGPG(ich,ia)) / TotalFlux_ALGP(ia)
            else
               Scale_ALGP(ia) = 0.0
            endif
         else
            Scale_ALGP(ia) = 1.0
         endif
         Ch_ALGPZ(ich,ia)    = Ch_ALGPZ(ich,ia) * Scale_ALGP(ia)
         Ch_ALGPR(ich,ia)    = Ch_ALGPR(ich,ia) * Scale_ALGP(ia)
         Ch_ALGPE(ich,ia)    = Ch_ALGPE(ich,ia) * Scale_ALGP(ia)
         Ch_ALGPM(ich,ia)    = Ch_ALGPM(ich,ia) * Scale_ALGP(ia)
         Ch_ALGPS(ich,ia)    = Ch_ALGPS(ich,ia) * Scale_ALGP(ia)
         Ch_ALGPdsch(ich,ia) = Ch_ALGPdsch(ich,ia) * Scale_ALGP(ia)
      ENDDO
      
      ! PIPA
      TotalFlux_PIPA = Ch_PIPAdsch(ich)
      if(Ch_PIPAADS(ich) > 0.0) &
         TotalFlux_PIPA = TotalFlux_PIPA + Ch_PIPAADS(ich)
      if(TotalFlux_PIPA > Ch_PIPA_init .and. TotalFlux_PIPA > 0.0) then
         Scale_PIPA = Ch_PIPA_init / TotalFlux_PIPA
      else
         Scale_PIPA = 1.0
      endif
         Ch_PIPAdsch(ich) = Ch_PIPAdsch(ich) * Scale_PIPA
      if(Ch_PIPAADS(ich) > 0.0) Ch_PIPAADS(ich) = Ch_PIPAADS(ich) * Scale_PIPA

      TotalFlux_PIPS = Ch_PIPSdsch(ich)
      if(TotalFlux_PIPS > Ch_PIPS_init .and. TotalFlux_PIPS > 0.0) then
         Scale_PIPS = Ch_PIPS_init / TotalFlux_PIPS
      else
         Scale_PIPS = 1.0
      endif
      Ch_PIPSdsch(ich) = Ch_PIPSdsch(ich) * Scale_PIPS

      ! === C. UPDATE STORAGES BASED ON SCALED FLUXES

      do ia=1,nalg
      Ch_ALGP(ich,ia) = Ch_ALGP_init(ia) + Ch_ALGPG(ich,ia) - &
                        Ch_ALGPZ(ich,ia) - Ch_ALGPR(ich,ia) - &
                        Ch_ALGPE(ich,ia) - Ch_ALGPM(ich,ia) - &
                        Ch_ALGPS(ich,ia) - Ch_ALGPdsch(ich,ia)
      end do
      
      Ch_ZOOP(ich) = Ch_ZOOP_init + Ch_ZOOPG(ich) - Ch_ZOOPZ(ich) - &
                     Ch_ZOOPR(ich) - Ch_ZOOPM(ich) - &  ! Ch_ZOOPE(ich) excluded!!!
                     Ch_ZOOPS(ich) - Ch_ZOOPdsch(ich)
      Ch_LPOP(ich) = Ch_LPOP_init + sum(Ch_ALGPM(ich,:)) + Ch_ZOOPM(ich) +  &
                     Ch_ZOOPE(ich) - dLPOPZ - dLPOPD - dLPOPS - &
                     Ch_LPOPdsch(ich)
      Ch_LDOP(ich) = Ch_LDOP_init + sum(Ch_ALGPE(ich,:)) + FDEC * &
                     (1.0 - FREF) * dLPOPD + FDEC * dMBMPD - dLDOPD - &
                     Ch_LDOPdsch(ich)
      Ch_RPOP(ich) = Ch_RPOP_init + FDEC * FREF * dLPOPD - dRPOPD - dRPOPS - &
                     Ch_RPOPdsch(ich)
      Ch_RDOP(ich) = Ch_RDOP_init + FDEC * (dLDOPD + dRPOPD + dRDOPD) - &
                     dRDOPD - Ch_RDOPdsch(ich)
      Ch_MBMP(ich) = Ch_MBMP_init + FBIO * (dLPOPD + dRPOPD + dLDOPD + &
                     dRDOPD + dMBMPD) - dMBMPD - Ch_MBMPdsch(ich)
      Ch_PO4(ich)  = Ch_PO4_init + FMET * (dLPOPD + dRPOPD + dLDOPD + &
                     dRDOPD + dMBMPD) + sum(Ch_ALGPR(ich,:)) + & 
                     Ch_ZOOPR(ich) - Ch_PO4uptk(ich) - Ch_PO4ADS(ich) - &
                     Ch_PO4dsch(ich)
      Ch_PIPA(ich) = Ch_PIPA_init + Ch_PO4ADS(ich) - Ch_PIPAADS(ich) - &
                     Ch_PIPAdsch(ich)
      Ch_PIPS(ich) = Ch_PIPS_init + Ch_PIPAADS(ich) - Ch_PIPSdsch(ich)
      Ch_POPDEPOSIT(ich) = Ch_POPDEPOSIT0(ich) + sum(Ch_ALGPS(ich,:)) + &
                     Ch_ZOOPS(ich) + dLPOPS + dRPOPS 

      ! Enforce non-negative channel storages
      do ia=1,nalg
        if (Ch_ALGP(ich,ia) < 1.0e-20) Ch_ALGP(ich,ia) = 0.0
      end do
      if (Ch_ZOOP(ich) < 1.0e-20) Ch_ZOOP(ich) = 0.0
      if (Ch_LPOP(ich) < 1.0e-20) Ch_LPOP(ich) = 0.0
      if (Ch_LDOP(ich) < 1.0e-20) Ch_LDOP(ich) = 0.0
      if (Ch_RPOP(ich) < 1.0e-20) Ch_RPOP(ich) = 0.0
      if (Ch_RDOP(ich) < 1.0e-20) Ch_RDOP(ich) = 0.0
      if (Ch_MBMP(ich) < 1.0e-20) Ch_MBMP(ich) = 0.0
      if (Ch_PO4(ich) < 1.0e-20)  Ch_PO4(ich) = 0.0
      if (Ch_PIPA(ich) < 1.0e-20) Ch_PIPA(ich) = 0.0
      if (Ch_PIPS(ich) < 1.0e-20) Ch_PIPS(ich) = 0.0      
      if (Ch_PONDEPOSIT(ich) < 1.0e-20) Ch_PONDEPOSIT(ich) = 0.0

   end subroutine ChannelPhosphorusCycle
   

   subroutine CheckPMB_So(i, j, lcover)
      
      implicit none

      ! --- Input arguments
      integer, intent(in) :: i, j
      integer, intent(in) :: lcover
      
      ! --- Local variables
      integer          :: isl, imonth
      character(len=2) :: smonth
      real(8)          :: dt

      dt = real(SedCNPmodel%SedCNP_timestep)
      smonth = trim(dateSedCNP%olddate(6:7))
      read(smonth,*) imonth      

      DO isl = 1, domain%nsl
         IF (isl == 1) THEN 
            So_PSC(i,isl,j) = &
                    (So_LPOPLIT(i,j)     - So_LPOPLIT0(i,j))     &
                  + (So_LPOPRES(i,isl,j) - So_LPOPRES0(i,isl,j)) &
                  + (So_LPOPEXC(i,j)     - So_LPOPEXC0(i,j))     &
                  + (So_LPOPMAN(i,j)     - So_LPOPMAN0(i,j))     &
                  + (So_RPOP(i,isl,j)    - So_RPOP0(i,isl,j))    &
                  + (So_LDOP(i,isl,j)    - So_LDOP0(i,isl,j))    &
                  + (So_RDOP(i,isl,j)    - So_RDOP0(i,isl,j))    &
                  + (So_MBMP(i,isl,j)    - So_MBMP0(i,isl,j))    &
                  + (So_PO4(i,isl,j)     - So_PO40(i,isl,j))     &
                  + (So_PIPA(i,isl,j)    - So_PIPA0(i,isl,j))    &
                  + (So_PIPS(i,isl,j)    - So_PIPS0(i,isl,j))

            So_PMBError(i,isl,j) = &
                    So_LPOPLITin(i,j)    + So_LPOPRESin(i,isl,j) &
                  + So_LPOPEXCin(i,j)    + So_LPOPMANin(i,j)     &
                  + So_PO4fxin(i,j)      + So_LDOPgwso0(i,j)     &  !BK20260316
                  + So_RDOPgwso0(i,j)    + So_PO4gwso0(i,j)      &  !BK20260316
                  - So_LPOPsurf(i,j)     - So_RPOPsurf(i,j)      &
                  - So_MBMPsurf(i,j)     - So_LDOPsurf(i,j)      &
                  - So_RDOPsurf(i,j)     - So_PO4surf(i,j)       &
                  - So_PIPAsurf(i,j)     - So_PIPSsurf(i,j)      &
                  - So_LDOPintf(i,isl,j) - So_RDOPintf(i,isl,j)  &
                  - So_PO4intf(i,isl,j)  - So_LDOPperc(i,isl,j)  &
                  - So_RDOPperc(i,isl,j) - So_PO4perc(i,isl,j)   &
                  - So_PO4uptk(i,isl,j)  - So_PSC(i,isl,j)
         ELSE 
            So_PSC(i,isl,j) = &
                    (So_LPOPRES(i,isl,j) - So_LPOPRES0(i,isl,j)) &
                  + (So_RPOP(i,isl,j)    - So_RPOP0(i,isl,j))    &
                  + (So_LDOP(i,isl,j)    - So_LDOP0(i,isl,j))    &
                  + (So_RDOP(i,isl,j)    - So_RDOP0(i,isl,j))    &
                  + (So_MBMP(i,isl,j)    - So_MBMP0(i,isl,j))    &
                  + (So_PO4(i,isl,j)     - So_PO40(i,isl,j))     &
                  + (So_PIPA(i,isl,j)    - So_PIPA0(i,isl,j))    &
                  + (So_PIPS(i,isl,j)    - So_PIPS0(i,isl,j))

            So_PMBError(i,isl,j) = &
                    So_LPOPRESin(i,isl,j)   + So_LDOPperc0(i,isl-1,j) &
                  + So_RDOPperc0(i,isl-1,j) + So_PO4perc0(i,isl-1,j)  &
                  - So_LDOPintf(i,isl,j)    - So_RDOPintf(i,isl,j)    &
                  - So_PO4intf(i,isl,j)     - So_LDOPperc(i,isl,j)    &
                  - So_RDOPperc(i,isl,j)    - So_PO4perc(i,isl,j)     &
                  - So_PO4uptk(i,isl,j)     - So_PSC(i,isl,j)
         END IF
      END DO

   end subroutine CheckPMB_So


   ! subroutine CheckPMB_Aq (i,j)

   !    implicit none

   !    ! --- Input arguments
   !    integer, intent(in) :: i, j

   !    Aq_PSC(i,j) = (Aq_LDOP(i,j) - Aq_LDOP0(i,j)) + &
   !                  (Aq_RDOP(i,j) - Aq_RDOP0(i,j)) + &
   !                  (Aq_PO4(i,j)  - Aq_PO40(i,j))  + &
   !                  (Aq_PIPA(i,j) - Aq_PIPA0(i,j)) + &
   !                  (Aq_PIPS(i,j) - Aq_PIPS0(i,j))
    
   !    Aq_PMBError(i,j) = So_LDOPsoaq0(i,j) + So_RDOPsoaq0(i,j) + &
   !                       So_PO4soaq0(i,j)  - Aq_LDOPaqso(i,j)  - &
   !                       Aq_RDOPaqso(i,j)  - Aq_PO4aqso(i,j)   - &
   !                       Aq_LDOPaqgw(i,j)  - Aq_RDOPaqgw(i,j)  - &
   !                       Aq_PO4aqgw(i,j)   - Aq_PSC(i,j)

   ! end subroutine CheckPMB_Aq


   subroutine CheckPMB_Gw (gwid)

      implicit none

      ! --- Input arguments
      integer, intent(in) :: gwid

      Gw_PSC(gwid) = &
              (Gw_LDOP(gwid) - Gw_LDOP0(gwid)) &
            + (Gw_RDOP(gwid) - Gw_RDOP0(gwid)) &
            + (Gw_PO4(gwid)  - Gw_PO40(gwid))

      Gw_PMBError(gwid) = &
              Gw_LDOPsogw0(gwid) + Gw_RDOPsogw0(gwid) &  !BK20260316
            + Gw_PO4sogw0(gwid)  - Gw_LDOPgwso(gwid)  &  !BK20260316
            - Gw_RDOPgwso(gwid)  - Gw_PO4gwso(gwid)   &
            - Gw_LDOPgwch(gwid)  - Gw_RDOPgwch(gwid)  &
            - Gw_PO4gwch(gwid)   - Gw_PSC(gwid)

   end subroutine CheckPMB_Gw


   subroutine CheckPMB_Ch (ich)

      implicit none

      ! --- Input arguments
      integer, intent(in) :: ich
      
      ! --- Local variables
      integer :: chid, gwid
      real(8) :: dt

      dt = real(SedCNPmodel%SedCNP_timestep)
      chid = SedCNP_hydro%linkID(ich)
      gwid = domain%gwid_ch(chid)

      Ch_PSC(ich) = &
              (Ch_ALGP(ich,1)     - Ch_ALGP0(ich,1))    &
            + (Ch_ALGP(ich,2)     - Ch_ALGP0(ich,2))    &
            + (Ch_ALGP(ich,3)     - Ch_ALGP0(ich,3))    &
            + (Ch_ZOOP(ich)       - Ch_ZOOP0(ich))      &
            + (Ch_MBMP(ich)       - Ch_MBMP0(ich))      &
            + (Ch_LPOP(ich)       - Ch_LPOP0(ich))      &
            + (Ch_RPOP(ich)       - Ch_RPOP0(ich))      &
            + (Ch_LDOP(ich)       - Ch_LDOP0(ich))      &
            + (Ch_RDOP(ich)       - Ch_RDOP0(ich))      &
            + (Ch_PO4(ich)        - Ch_PO40(ich))       &
            + (Ch_PIPA(ich)       - Ch_PIPA0(ich))      &
            + (Ch_PIPS(ich)       - Ch_PIPS0(ich))      &
            + (Ch_POPDEPOSIT(ich) - Ch_POPDEPOSIT0(ich))

      Ch_PMBError(ich) = &
               Ch_ALGPusch0(ich,1) + Ch_ALGPusch0(ich,2) &
            + Ch_ALGPusch0(ich,3)  + Ch_ZOOPusch0(ich)   &
            + Ch_MBMPusch0(ich)    + Ch_LPOPusch0(ich)   &
            + Ch_RPOPusch0(ich)    + Ch_LDOPusch0(ich)   &
            + Ch_RDOPusch0(ich)    + Ch_PO4usch0(ich)    &
            + Ch_PIPAusch0(ich)    + Ch_PIPSusch0(ich)   &
            + Ch_LPOPsurf0(ich)    + Ch_RPOPsurf0(ich)   &
            + Ch_MBMPsurf0(ich)    + Ch_LDOPsurf0(ich)   &
            + Ch_RDOPsurf0(ich)    + Ch_PO4surf0(ich)    &
            + Ch_PIPAsurf0(ich)    + Ch_PIPSsurf0(ich)   &
            + Ch_LDOPintf0(ich)    + Ch_RDOPintf0(ich)   &
            + Ch_PO4intf0(ich)     + Ch_LDOPgwch0(ich)   &
            + Ch_RDOPgwch0(ich)    + Ch_PO4gwch0(ich)    &
            + Ch_LPOPxpnt(ich)     + Ch_RPOPxpnt(ich)    &
            + Ch_LDOPxpnt(ich)     + Ch_RDOPxpnt(ich)    &
            + Ch_PO4xpnt(ich)      + Ch_LPOPxdis(ich)    &
            + Ch_RPOPxdis(ich)     + Ch_LDOPxdis(ich)    &
            + Ch_RDOPxdis(ich)     + Ch_PO4xdis(ich)     &
            - Ch_LPOPxabs(ich)     - Ch_RPOPxabs(ich)    &
            - Ch_LDOPxabs(ich)     - Ch_RDOPxabs(ich)    &
            - Ch_PO4xabs(ich)      - Ch_ALGPdsch(ich,1)  &
            - Ch_ALGPdsch(ich,2)   - Ch_ALGPdsch(ich,3)  &
            - Ch_ZOOPdsch(ich)     - Ch_MBMPdsch(ich)    &
            - Ch_LPOPdsch(ich)     - Ch_RPOPdsch(ich)    &
            - Ch_LDOPdsch(ich)     - Ch_RDOPdsch(ich)    &
            - Ch_PO4dsch(ich)      - Ch_PIPAdsch(ich)    &
            - Ch_PIPSdsch(ich)     - Ch_PSC(ich)

   end subroutine CheckPMB_Ch

end module Phosphorus