module Nitrogen

   use CNPfunctions
   use module_SedCNPvariables
!=====||__WHQ5403q__||=====!
!
   use SedCNP_config, only: SedCNPmodel
!
!=====||__WHQ5403q__||=====!
   use CNPparams

   implicit none

contains

   subroutine SoilNitrogenCycle(i, j, lcover, itime, TSOIL, THETA, &
      THETA_S, THETA_F, So_WStorage, Quptk, Qsurf, Qperc, Qintf)

      implicit none

      ! --- Input arguments
      integer,               intent(in) :: i, j, lcover, itime
      real(8), dimension(4), intent(in) :: TSOIL, THETA, So_WStorage, Quptk
      real(8),               intent(in) :: THETA_S, THETA_F, Qsurf
      real(8), dimension(4), intent(in) :: Qperc
      real(8),               intent(in) :: Qintf

      ! --- Local variables
      integer               :: isl, status, imonth
      character(len=2)      :: smonth
      real(8)               :: dt, rate0, KTBMBM, KTHETA, dPI, RK4SOL
      real(8), allocatable, dimension(:) :: Qintf_layer, Weight  !BK20260612
      real(8)               :: TotalWeight
      real(8)               :: dx    ! grid cell size [m]  !BK20260301
      real(8)               :: q_out, w_strg  !(m3) !BK20260301      
      integer               :: k

      ! Decomposition and transport fluxes
      real(8) :: dLPONLITD, dLPONRESD, dLPONEXCD, dLPONMAND
      real(8) :: dRPOND, dLDOND, dRDOND, dMBMND
      real(8) :: dLPONLITsurf, dLPONRESsurf, dLPONEXCsurf, dLPONMANsurf

      ! Initial State and Flux Scaling Variables
      real(8) :: So_LPONLIT_init, So_LPONRES_init, So_LPONEXC_init
      real(8) :: So_LPONMAN_init, So_RPON_init, So_LDON_init, So_RDON_init
      real(8) :: So_MBMN_init, So_NH4_init, So_NO3_init
      real(8) :: TotalFlux_LPONLIT, TotalFlux_LPONRES, TotalFlux_LPONEXC
      real(8) :: TotalFlux_LPONMAN, TotalFlux_RPON, TotalFlux_LDON
      real(8) :: TotalFlux_RDON, TotalFlux_MBMN, TotalFlux_NH4, TotalFlux_NO3
      real(8) :: Scale_LPONLIT, Scale_LPONRES, Scale_LPONEXC, Scale_LPONMAN
      real(8) :: Scale_RPON, Scale_LDON, Scale_RDON, Scale_MBMN
      real(8) :: Scale_NH4, Scale_NO3


      ! ====================================================================
      ! I. INITIALIZE & UPDATE STORAGES WITH EXTERNAL/INTERNAL INPUTS
      ! ====================================================================

      dt = real(SedCNPmodel%SedCNP_timestep)
      dx = SedCNP_hydro%dx(1)    !BK20260301
      smonth = trim(dateSedCNP%olddate(6:7))
      read(smonth,*) imonth

      allocate(Qintf_layer(domain%nsl), Weight(domain%nsl)) !BK20260612

      ! Apportion total interflow Qintf to layers based on water content above 
      ! field capacity
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

      ! ===================================================================
      ! II. CALCULATE-SCALE-UPDATE LOGIC
      ! ===================================================================

      DO isl = 1, domain%nsl

         ! -----------------------------------------------------------------
         ! A. Store initial values for simultaneous flux calculation
         ! (previous timestep state + current timestep inputs)
         ! -----------------------------------------------------------------
         
         if (isl == 1) then
            So_LPONLITin(i,j) = So_LPONLITinso(lcover,imonth) * dt
            So_LPONEXCin(i,j) = So_LPONEXCinso(lcover,imonth) * dt
            So_LPONMANin(i,j) = So_LPONMANinso(lcover,imonth) * dt
            So_LPONRESin(i,isl,j) = So_LPONRESinso(lcover,imonth) * dt * &
               FROOT(lcover, isl)
            So_NH4fxin(i,j)   = So_NH4fert(lcover,imonth) * dt
            So_NO3fxin(i,j)   = So_NO3fert(lcover,imonth) * dt

            So_LPONLIT_init = So_LPONLIT0(i,j) + So_LPONLITin(i,j)
            So_LPONEXC_init = So_LPONEXC0(i,j) + So_LPONEXCin(i,j)
            So_LPONMAN_init = So_LPONMAN0(i,j) + So_LPONMANin(i,j)
            So_LPONRES_init = So_LPONRES0(i,isl,j) + So_LPONRESin(i,isl,j)
            So_RPON_init    = So_RPON0(i,isl,j)
            So_LDON_init    = So_LDON0(i,isl,j) + So_LDONgwso0(i,j)  !BK20260316
            So_RDON_init    = So_RDON0(i,isl,j) + So_RDONgwso0(i,j)  !BK20260316
            So_MBMN_init    = So_MBMN0(i,isl,j)
            So_NH4_init     = So_NH40(i,isl,j) + So_NH4gwso0(i,j) + So_NH4fxin(i,j)  !BK20260316
            So_NO3_init     = So_NO30(i,isl,j) + So_NO3gwso0(i,j) + So_NO3fxin(i,j)  !BK20260316
         else
            So_LPONRESin(i,isl,j) = So_LPONRESinso(lcover,imonth) * dt * &
               FROOT(lcover, isl)

            So_LPONLIT_init = 0.0  ! Does not exist in lower layers
            So_LPONEXC_init = 0.0  ! Does not exist in lower layers
            So_LPONMAN_init = 0.0  ! Does not exist in lower layers
            So_LPONRES_init = So_LPONRES0(i,isl,j) + So_LPONRESin(i,isl,j)
            So_RPON_init    = So_RPON0(i,isl,j)
            So_LDON_init    = So_LDON0(i,isl,j) + So_LDONperc0(i,isl-1,j)
            So_RDON_init    = So_RDON0(i,isl,j) + So_RDONperc0(i,isl-1,j)
            So_MBMN_init    = So_MBMN0(i,isl,j)
            So_NH4_init     = So_NH40(i,isl,j) + So_NH4perc0(i,isl-1,j)
            So_NO3_init     = So_NO30(i,isl,j) + So_NO3perc0(i,isl-1,j)
         endif

         ! -----------------------------------------------------------------
         ! B. Calculate all transformation and transport fluxes
         ! -----------------------------------------------------------------

         KTBMBM = calc_KTB (TSOIL(isl), TMBM, MMBM)
         KTHETA = calc_KTHETA (THETA(isl), THETA_S, THETA_F)

         ! --- Organic N Decomposition
         if (isl == 1) then
            rate0 = KTBMBM * KTHETA * KLIT(lcover)
            RK4SOL = calc_RK4SOL(-rate0, So_LPONLIT_init)
            dLPONLITD = So_LPONLIT_init - RK4SOL

            rate0 = KTBMBM * KTHETA * KEXC(lcover)
            RK4SOL = calc_RK4SOL(-rate0, So_LPONEXC_init)
            dLPONEXCD = So_LPONEXC_init - RK4SOL

            rate0 = KTBMBM * KTHETA * KMAN(lcover)
            RK4SOL = calc_RK4SOL(-rate0, So_LPONMAN_init)
            dLPONMAND = So_LPONMAN_init - RK4SOL
         endif

         rate0 = KTBMBM * KTHETA * KRES(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_LPONRES_init)
         dLPONRESD = So_LPONRES_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KRPOM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_RPON_init)
         dRPOND = So_RPON_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KLDOM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_LDON_init)
         dLDOND = So_LDON_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KRDOM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_RDON_init)
         dRDOND = So_RDON_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KMBM(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_MBMN_init)
         dMBMND = So_MBMN_init - RK4SOL

         ! --- Inorganic N Transformations & Plant Uptake
         if ((Quptk(isl) > 0.0) .and. (So_WStorage(isl) > 0.0)) then
            So_NH4uptk(i,isl,j) = (1.0 - exp(-Quptk(isl) * domain%areaxy / &
               1000.0 * domain%DT / So_WStorage(isl))) * So_NH4_init
            So_NO3uptk(i,isl,j) = (1.0 - exp(-Quptk(isl) * domain%areaxy / &
               1000.0 * domain%DT / So_WStorage(isl))) * So_NO3_init
         else
            So_NH4uptk(i,isl,j) = 0.0
            So_NO3uptk(i,isl,j) = 0.0
         endif

         rate0 = KTBMBM * KTHETA * KNIT(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_NH4_init)
         So_NITRI(i,isl,j) = So_NH4_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KVOL(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_NH4_init)
         So_NH3VOL(i,isl,j) = So_NH4_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KDEN(lcover)
         RK4SOL = calc_RK4SOL(-rate0, So_NO3_init)
         So_DENIT(i,isl,j) = So_NO3_init - RK4SOL

         ! --- Transport by surface runoff
         if (isl == 1) then
            q_out = Qsurf / 1000. * dx * dx * dt ! (m3)
            w_strg = THETA(1) * Dmix * dx * dx   ! (m3)  !BK20260302
            !POM
            dPI = calc_dPI(etaPOMsurf, q_out, w_strg)  !BK20260302
            dLPONLITsurf = So_LPONLIT_init * (Dmix / DSOIL(isl)) * dPI
            dLPONEXCsurf = So_LPONEXC_init * (Dmix / DSOIL(isl)) * dPI
            dLPONMANsurf = So_LPONMAN_init * (Dmix / DSOIL(isl)) * dPI
            dLPONRESsurf = So_LPONRES_init * (Dmix / DSOIL(isl)) * dPI
            So_LPONsurf(i,j) = dLPONLITsurf + dLPONRESsurf + &
                               dLPONEXCsurf + dLPONMANsurf
            So_RPONsurf(i,j) = So_RPON_init * (Dmix / DSOIL(isl)) * dPI
            So_MBMNsurf(i,j) = So_MBMN_init * (Dmix / DSOIL(isl)) * dPI
            !DOM
            dPI = calc_dPI(etaDOMsurf, q_out, w_strg)  !BK20260302
            So_LDONsurf(i,j) = So_LDON_init * (Dmix / DSOIL(isl)) * dPI
            So_RDONsurf(i,j) = So_RDON_init * (Dmix / DSOIL(isl)) * dPI
            !NH4
            dPI = calc_dPI(etaNH4surf, q_out, w_strg)  !BK20260302
            So_NH4surf(i,j) = So_NH4_init * (Dmix / DSOIL(isl)) * dPI
            !NO3
            dPI = calc_dPI(etaNO3surf, q_out, w_strg)  !BK20260302
            So_NO3surf(i,j) = So_NO3_init * (Dmix / DSOIL(isl)) * dPI
         endif

         ! --- Transport by Interflow
         q_out = Qintf_layer(isl) / 1000. * dx * dx * dt  ! (m3)
         w_strg = THETA(isl) * dx * dx    ! (m3)  !BK20260302   
         !DOM
         dPI = calc_dPI(etaDOMintf, q_out, w_strg)  !BK20260302
         So_LDONintf(i,isl,j) = So_LDON_init * dPI
         So_RDONintf(i,isl,j) = So_RDON_init * dPI
         !NH4
         dPI = calc_dPI(etaNH4intf, q_out, w_strg)  !BK20260302
         So_NH4intf(i,isl,j) = So_NH4_init * dPI
         !NO3
         dPI = calc_dPI(etaNO3intf, q_out, w_strg)  !BK20260302
         So_NO3intf(i,isl,j) = So_NO3_init * dPI

         ! --- Transport by Percolation
         q_out = Qperc(isl) / 1000. * dx * dx * dt  ! (m3)
         w_strg = THETA(isl) * dx * dx    ! (m3)  !BK20260302   
         !DOM
         dPI = calc_dPI(etaDOMperc, q_out, w_strg)  !BK20260302
         So_LDONperc(i,isl,j) = So_LDON_init * dPI
         So_RDONperc(i,isl,j) = So_RDON_init * dPI
         !NH4
         dPI = calc_dPI(etaNH4perc, q_out, w_strg)  !BK20260302
         So_NH4perc(i,isl,j) = So_NH4_init * dPI
         !NO3
         dPI = calc_dPI(etaNO3perc, q_out, w_strg)  !BK20260302
         So_NO3perc(i,isl,j) = So_NO3_init * dPI

         ! -----------------------------------------------------------------
         ! C. Scale fluxes to ensure mass conservation
         ! -----------------------------------------------------------------

         if (isl == 1) then
            TotalFlux_LPONLIT = dLPONLITD + dLPONLITsurf
            if (TotalFlux_LPONLIT > So_LPONLIT_init .and. &
                  TotalFlux_LPONLIT > 0.0) then
               Scale_LPONLIT = So_LPONLIT_init / TotalFlux_LPONLIT
               dLPONLITD = dLPONLITD * Scale_LPONLIT
               dLPONLITsurf = dLPONLITsurf * Scale_LPONLIT
            endif
            
            TotalFlux_LPONEXC = dLPONEXCD + dLPONEXCsurf
            if (TotalFlux_LPONEXC > So_LPONEXC_init .and. &
                  TotalFlux_LPONEXC > 0.0) then
               Scale_LPONEXC = So_LPONEXC_init / TotalFlux_LPONEXC
               dLPONEXCD = dLPONEXCD * Scale_LPONEXC
               dLPONEXCsurf = dLPONEXCsurf * Scale_LPONEXC
            endif
            
            TotalFlux_LPONMAN = dLPONMAND + dLPONMANsurf
            if (TotalFlux_LPONMAN > So_LPONMAN_init .and. &
                  TotalFlux_LPONMAN > 0.0) then
               Scale_LPONMAN = So_LPONMAN_init / TotalFlux_LPONMAN
               dLPONMAND = dLPONMAND * Scale_LPONMAN
               dLPONMANsurf = dLPONMANsurf * Scale_LPONMAN
            endif
         endif

         if (isl == 1) then
            TotalFlux_LPONRES = dLPONRESD + dLPONRESsurf
         else
            TotalFlux_LPONRES = dLPONRESD
         endif
         if (TotalFlux_LPONRES > So_LPONRES_init .and. &
               TotalFlux_LPONRES > 0.0) then
            Scale_LPONRES = So_LPONRES_init / TotalFlux_LPONRES
            dLPONRESD = dLPONRESD * Scale_LPONRES
            if(isl == 1) dLPONRESsurf = dLPONRESsurf * Scale_LPONRES
         endif

         if (isl == 1) then
             TotalFlux_RPON = dRPOND + So_RPONsurf(i,j)
         else
             TotalFlux_RPON = dRPOND
         endif
         if (TotalFlux_RPON > So_RPON_init .and. TotalFlux_RPON > 0.0) then
            Scale_RPON = So_RPON_init / TotalFlux_RPON
            dRPOND = dRPOND * Scale_RPON
            if(isl == 1) So_RPONsurf(i,j) = So_RPONsurf(i,j) * Scale_RPON
         endif

         if (isl == 1) then
            TotalFlux_LDON = dLDOND + So_LDONsurf(i,j) + &
                  So_LDONintf(i,isl,j) + So_LDONperc(i,isl,j)
         else
            TotalFlux_LDON = dLDOND + So_LDONintf(i,isl,j) + &
                  So_LDONperc(i,isl,j)
         endif
         if(TotalFlux_LDON > So_LDON_init .and. TotalFlux_LDON > 0.0) then
            Scale_LDON = So_LDON_init/TotalFlux_LDON
            dLDOND = dLDOND * Scale_LDON
            if(isl == 1) So_LDONsurf(i,j) = So_LDONsurf(i,j) * Scale_LDON
            So_LDONintf(i,isl,j) = So_LDONintf(i,isl,j) * Scale_LDON
            So_LDONperc(i,isl,j) = So_LDONperc(i,isl,j) * Scale_LDON
         endif

         if (isl == 1) then
            TotalFlux_RDON = dRDOND + So_RDONsurf(i,j) + &
                  So_RDONintf(i,isl,j) + So_RDONperc(i,isl,j)
         else
            TotalFlux_RDON = dRDOND + So_RDONintf(i,isl,j) + &
                  So_RDONperc(i,isl,j)
         endif
         if(TotalFlux_RDON > So_RDON_init .and. TotalFlux_RDON > 0.0) then
            Scale_RDON = So_RDON_init/TotalFlux_RDON
            dRDOND = dRDOND * Scale_RDON
            if(isl == 1) So_RDONsurf(i,j) = So_RDONsurf(i,j) * Scale_RDON
            So_RDONintf(i,isl,j) = So_RDONintf(i,isl,j) * Scale_RDON
            So_RDONperc(i,isl,j) = So_RDONperc(i,isl,j) * Scale_RDON
         endif

         if (isl == 1) then
            TotalFlux_MBMN = dMBMND + So_MBMNsurf(i,j)
         else
            TotalFlux_MBMN = dMBMND
         endif
         if (TotalFlux_MBMN > So_MBMN_init .and. TotalFlux_MBMN > 0.0) then
            Scale_MBMN = So_MBMN_init / TotalFlux_MBMN
            dMBMND = dMBMND * Scale_MBMN
            if(isl == 1) So_MBMNsurf(i,j) = So_MBMNsurf(i,j) * Scale_MBMN
         endif

         if (isl == 1) then
            TotalFlux_NH4 = So_NH4uptk(i,isl,j) + So_NITRI(i,isl,j) + &
               So_NH3VOL(i,isl,j) + So_NH4surf(i,j) + So_NH4intf(i,isl,j) + &
               So_NH4perc(i,isl,j)
         else
            TotalFlux_NH4 = So_NH4uptk(i,isl,j) + So_NITRI(i,isl,j) + &
               So_NH3VOL(i,isl,j) + So_NH4intf(i,isl,j) + So_NH4perc(i,isl,j)
         endif
         if (TotalFlux_NH4 > So_NH4_init .and. TotalFlux_NH4 > 0.0) then
            Scale_NH4 = So_NH4_init / TotalFlux_NH4
            So_NH4uptk(i,isl,j) = So_NH4uptk(i,isl,j) * Scale_NH4
            So_NITRI(i,isl,j)   = So_NITRI(i,isl,j)   * Scale_NH4
            So_NH3VOL(i,isl,j)  = So_NH3VOL(i,isl,j)  * Scale_NH4
            if(isl == 1) So_NH4surf(i,j) = So_NH4surf(i,j) * Scale_NH4
            So_NH4intf(i,isl,j) = So_NH4intf(i,isl,j) * Scale_NH4
            So_NH4perc(i,isl,j) = So_NH4perc(i,isl,j) * Scale_NH4
         endif

         if (isl == 1) then
            TotalFlux_NO3 = So_NO3uptk(i,isl,j) + So_DENIT(i,isl,j) + &
               So_NO3surf(i,j) + So_NO3intf(i,isl,j) + So_NO3perc(i,isl,j)
         else
            TotalFlux_NO3 = So_NO3uptk(i,isl,j) + So_DENIT(i,isl,j) + &
               So_NO3intf(i,isl,j) + So_NO3perc(i,isl,j)
         endif
         if (TotalFlux_NO3 > So_NO3_init .and. TotalFlux_NO3 > 0.0) then
            Scale_NO3 = So_NO3_init / TotalFlux_NO3
            So_NO3uptk(i,isl,j) = So_NO3uptk(i,isl,j) * Scale_NO3
            So_DENIT(i,isl,j)   = So_DENIT(i,isl,j)   * Scale_NO3
            if(isl == 1) So_NO3surf(i,j) = So_NO3surf(i,j) * Scale_NO3
            So_NO3intf(i,isl,j) = So_NO3intf(i,isl,j) * Scale_NO3
            So_NO3perc(i,isl,j) = So_NO3perc(i,isl,j) * Scale_NO3
         endif

         ! Recalculate So_LPONsurf with potentially scaled component fluxes
         if (isl == 1) then
            So_LPONsurf(i,j) = dLPONLITsurf + dLPONRESsurf + &
                                 dLPONEXCsurf + dLPONMANsurf
         endif

         ! -----------------------------------------------------------------
         ! D. Update final storages using scaled fluxes
         ! -----------------------------------------------------------------

         if (isl == 1) then
            So_LPONLIT(i,j) = So_LPONLIT_init - dLPONLITD - dLPONLITsurf
            So_LPONEXC(i,j) = So_LPONEXC_init - dLPONEXCD - dLPONEXCsurf
            So_LPONMAN(i,j) = So_LPONMAN_init - dLPONMAND - dLPONMANsurf
            So_LPONRES(i,isl,j) = So_LPONRES_init - dLPONRESD - dLPONRESsurf
            So_RPON(i,isl,j) = So_RPON_init + FDEC * FREF * (dLPONLITD + &
               dLPONRESD + dLPONEXCD + dLPONMAND) - dRPOND - So_RPONsurf(i,j)
            So_LDON(i,isl,j) = So_LDON_init + FDEC * (1. - FREF) * &
               (dLPONLITD + dLPONRESD + dLPONEXCD + dLPONMAND) + FDEC * &
               dMBMND - dLDOND - So_LDONsurf(i,j) - So_LDONintf(i,isl,j) - &
               So_LDONperc(i,isl,j)
            So_RDON(i,isl,j) = So_RDON_init + FDEC * (dLDOND + dRPOND + &
               dRDOND) - dRDOND - So_RDONsurf(i,j) - So_RDONintf(i,isl,j) - &
               So_RDONperc(i,isl,j)
            So_MBMN(i,isl,j) = So_MBMN_init + FBIO * (dLPONLITD + dLPONRESD + &
               dLPONEXCD + dLPONMAND + dRPOND + dLDOND + dRDOND + dMBMND) - &
               dMBMND - So_MBMNsurf(i,j)
            So_NH4(i,isl,j)  = So_NH4_init + FMET * (dLPONLITD + dLPONRESD + &
               dLPONEXCD + dLPONMAND + dRPOND + dLDOND + dRDOND + dMBMND) - &
               So_NH4uptk(i,isl,j) - So_NITRI(i,isl,j) - So_NH3VOL(i,isl,j) - &
               So_NH4surf(i,j) - So_NH4intf(i,isl,j) - So_NH4perc(i,isl,j)
            So_NO3(i,isl,j)  = So_NO3_init + So_NITRI(i,isl,j) - &
               So_NO3uptk(i,isl,j) - So_DENIT(i,isl,j) - So_NO3surf(i,j) - &
               So_NO3intf(i,isl,j) - So_NO3perc(i,isl,j)
         else
            So_LPONRES(i,isl,j) = So_LPONRES_init - dLPONRESD
            So_RPON(i,isl,j) = So_RPON_init + FDEC * FREF * dLPONRESD - dRPOND
            So_LDON(i,isl,j) = So_LDON_init + FDEC * (1. - FREF) * dLPONRESD + &
               FDEC*dMBMND - dLDOND - So_LDONintf(i,isl,j) - So_LDONperc(i,isl,j)
            So_RDON(i,isl,j) = So_RDON_init + FDEC * (dLDOND + dRPOND + &
               dRDOND) - dRDOND - So_RDONintf(i,isl,j) - So_RDONperc(i,isl,j)
            So_MBMN(i,isl,j) = So_MBMN_init + FBIO * (dLPONRESD + dRPOND + &
               dLDOND + dRDOND + dMBMND) - dMBMND
            So_NH4(i,isl,j)  = So_NH4_init + FMET * (dLPONRESD + dRPOND + &
               dLDOND + dRDOND + dMBMND) - So_NH4uptk(i,isl,j) - &
               So_NITRI(i,isl,j) - So_NH3VOL(i,isl,j) - So_NH4intf(i,isl,j) - &
               So_NH4perc(i,isl,j)
            So_NO3(i,isl,j)  = So_NO3_init + So_NITRI(i,isl,j) - &
               So_NO3uptk(i,isl,j) - So_DENIT(i,isl,j) - So_NO3intf(i,isl,j) - &
               So_NO3perc(i,isl,j)
         endif

      END DO

      ! Enforce non-negative soil storages
      if (So_LPONLIT(i,j) < 1.0e-20) So_LPONLIT(i,j) = 0.0
      if (So_LPONEXC(i,j) < 1.0e-20) So_LPONEXC(i,j) = 0.0
      if (So_LPONMAN(i,j) < 1.0e-20) So_LPONMAN(i,j) = 0.0
      DO isl = 1, domain%nsl
        if (So_LPONRES(i,isl,j) < 1.0e-20) So_LPONRES(i,isl,j) = 0.0
        if (So_RPON(i,isl,j) < 1.0e-20) So_RPON(i,isl,j) = 0.0
        if (So_LDON(i,isl,j) < 1.0e-20) So_LDON(i,isl,j) = 0.0
        if (So_RDON(i,isl,j) < 1.0e-20) So_RDON(i,isl,j) = 0.0
        if (So_MBMN(i,isl,j) < 1.0e-20) So_MBMN(i,isl,j) = 0.0
        if (So_NH4(i,isl,j) < 1.0e-20) So_NH4(i,isl,j) = 0.0
        if (So_NO3(i,isl,j) < 1.0e-20) So_NO3(i,isl,j) = 0.0
      END DO

      ! ====================================================================
      ! III. UPDATE OUTPUTS TO OTHER COMPONENTS
      ! ====================================================================

      So_LDONsogw(i,j) = So_LDONperc(i,domain%nsl,j)  !BK20260318
      So_RDONsogw(i,j) = So_RDONperc(i,domain%nsl,j)  !BK20260318
      So_NH4sogw(i,j)  = So_NH4perc(i,domain%nsl,j)  !BK20260318
      So_NO3sogw(i,j)  = So_NO3perc(i,domain%nsl,j)  !BK20260318

   end subroutine SoilNitrogenCycle


   ! subroutine AquiferNitrogenCycle(i, j, Qaqso, Aq_WStorage, Qaqgw)

   !    implicit none

   !    ! --- Input arguments
   !    integer, intent(in) :: i, j
   !    real(8), intent(in) :: Qaqso, Aq_WStorage, Qaqgw

   !    ! --- Local variables
   !    real(8) :: q_out, dx    !BK20260302
   !    real(8) :: dPI, dt

   !    ! --- Variables for calculate-scale-update logic
   !    real(8) :: Aq_LDON_init, Aq_RDON_init, Aq_NH4_init, Aq_NO3_init
   !    real(8) :: TotalFlux_LDON, TotalFlux_RDON, TotalFlux_NH4, TotalFlux_NO3
   !    real(8) :: Scale_LDON, Scale_RDON, Scale_NH4, Scale_NO3

   !    dt = real(SedCNPmodel%SedCNP_timestep)
   !    dx = SedCNP_hydro%dx(1)    !BK20260301

   !    ! ====================================================================
   !    ! I. INITIALIZE STATE AND ADD INFLOWS
   !    ! ====================================================================

   !    Aq_LDON_init = Aq_LDON0(i,j) + So_LDONsoaq0(i,j)
   !    Aq_RDON_init = Aq_RDON0(i,j) + So_RDONsoaq0(i,j)
   !    Aq_NH4_init  = Aq_NH40(i,j)  + So_NH4soaq0(i,j)
   !    Aq_NO3_init  = Aq_NO30(i,j)  + So_NO3soaq0(i,j)

   !    ! ====================================================================
   !    ! II. CALCULATE TRANSPORT FLUXES
   !    ! ====================================================================

   !    ! --- Abstraction for irrigation
   !    if (Aq_WStorage > (Qaqso * domain%areaxy / 1000.0 * dt)) then
   !       Aq_LDONaqso(i,j) = (Qaqso*domain%areaxy/1000.0*dt)*Aq_LDON_init / &
   !                          Aq_WStorage
   !       Aq_RDONaqso(i,j) = (Qaqso*domain%areaxy/1000.0*dt)*Aq_RDON_init / &
   !                          Aq_WStorage
   !       Aq_NH4aqso(i,j)  = (Qaqso*domain%areaxy/1000.0*dt)*Aq_NH4_init / &
   !                          Aq_WStorage
   !       Aq_NO3aqso(i,j)  = (Qaqso*domain%areaxy/1000.0*dt)*Aq_NO3_init / &
   !                          Aq_WStorage
   !    else
   !       Aq_LDONaqso(i,j) = Aq_LDON_init
   !       Aq_RDONaqso(i,j) = Aq_RDON_init
   !       Aq_NH4aqso(i,j)  = Aq_NH4_init
   !       Aq_NO3aqso(i,j)  = Aq_NO3_init
   !    endif

   !    ! --- Aquifer discharge to groundwater
   !    q_out = Qaqgw / 1000. * dx * dx * dt  ! (m3)
   !    !DOM
   !    dPI = calc_dPI(etaDOMaqgw, q_out, Aq_Wstorage)  !BK20260302
   !    Aq_LDONaqgw(i,j) = Aq_LDON_init * dPI
   !    Aq_RDONaqgw(i,j) = Aq_RDON_init * dPI
   !    !NH4
   !    dPI = calc_dPI(etaNH4aqgw, q_out, Aq_Wstorage)  !BK20260302
   !    Aq_NH4aqgw(i,j) = Aq_NH4_init * dPI
   !    !NO3
   !    dPI = calc_dPI(etaNO3aqgw, q_out, Aq_Wstorage)  !BK20260302
   !    Aq_NO3aqgw(i,j) = Aq_NO3_init * dPI

   !    ! ====================================================================
   !    ! III. SCALE FLUXES TO ENSURE MASS CONSERVATION
   !    ! ====================================================================

   !    TotalFlux_LDON = Aq_LDONaqso(i,j) + Aq_LDONaqgw(i,j)
   !    if (TotalFlux_LDON > Aq_LDON_init .and. TotalFlux_LDON > 0.0) then
   !       Scale_LDON = Aq_LDON_init / TotalFlux_LDON
   !       Aq_LDONaqso(i,j) = Aq_LDONaqso(i,j) * Scale_LDON
   !       Aq_LDONaqgw(i,j) = Aq_LDONaqgw(i,j) * Scale_LDON
   !    endif

   !    TotalFlux_RDON = Aq_RDONaqso(i,j) + Aq_RDONaqgw(i,j)
   !    if (TotalFlux_RDON > Aq_RDON_init .and. TotalFlux_RDON > 0.0) then
   !       Scale_RDON = Aq_RDON_init / TotalFlux_RDON
   !       Aq_RDONaqso(i,j) = Aq_RDONaqso(i,j) * Scale_RDON
   !       Aq_RDONaqgw(i,j) = Aq_RDONaqgw(i,j) * Scale_RDON
   !    endif

   !    TotalFlux_NH4 = Aq_NH4aqso(i,j) + Aq_NH4aqgw(i,j)
   !    if (TotalFlux_NH4 > Aq_NH4_init .and. TotalFlux_NH4 > 0.0) then
   !       Scale_NH4 = Aq_NH4_init / TotalFlux_NH4
   !       Aq_NH4aqso(i,j) = Aq_NH4aqso(i,j) * Scale_NH4
   !       Aq_NH4aqgw(i,j) = Aq_NH4aqgw(i,j) * Scale_NH4
   !    endif

   !    TotalFlux_NO3 = Aq_NO3aqso(i,j) + Aq_NO3aqgw(i,j)
   !    if (TotalFlux_NO3 > Aq_NO3_init .and. TotalFlux_NO3 > 0.0) then
   !       Scale_NO3 = Aq_NO3_init / TotalFlux_NO3
   !       Aq_NO3aqso(i,j) = Aq_NO3aqso(i,j) * Scale_NO3
   !       Aq_NO3aqgw(i,j) = Aq_NO3aqgw(i,j) * Scale_NO3
   !    endif

   !    ! ====================================================================
   !    ! IV. UPDATE FINAL STORAGES
   !    ! ====================================================================

   !    Aq_LDON(i,j) = Aq_LDON_init - Aq_LDONaqso(i,j) - Aq_LDONaqgw(i,j)
   !    Aq_RDON(i,j) = Aq_RDON_init - Aq_RDONaqso(i,j) - Aq_RDONaqgw(i,j)
   !    Aq_NH4(i,j)  = Aq_NH4_init  - Aq_NH4aqso(i,j)  - Aq_NH4aqgw(i,j)
   !    Aq_NO3(i,j)  = Aq_NO3_init  - Aq_NO3aqso(i,j)  - Aq_NO3aqgw(i,j)

   !    ! Enforce non-negative aquifer storages
   !    if (Aq_LDON(i,j) < 1.0e-20) Aq_LDON(i,j) = 0.0
   !    if (Aq_RDON(i,j) < 1.0e-20) Aq_RDON(i,j) = 0.0
   !    if (Aq_NH4(i,j) < 1.0e-20) Aq_NH4(i,j) = 0.0
   !    if (Aq_NO3(i,j) < 1.0e-20) Aq_NO3(i,j) = 0.0

   ! end subroutine AquiferNitrogenCycle


   subroutine GroundwaterNitrogenCycle(gwid, Gw_WStorage, Qgwch)  !BK20260318

      implicit none

      ! --- Input arguments
      integer, intent(in) :: gwid
      real(8), intent(in)  :: Gw_WStorage    ![m3]       !BK20260316
      real(8), intent(in)  :: Qgwch          ![mm s-1]   !BK20260316

      ! --- Local variables
      integer :: i, j, istat
      real(8) :: q_out                       ![m3]       !BK20260301  
      real(8) :: dPI, dx, dt                             !BK20260302     
      real(8), allocatable :: Qgwso(:,:)     ![mm s-1]   !BK20260316
      real(8), allocatable :: Qgwso_bas(:)   ![m3]       !BK20260316

      ! --- Variables for calculate-scale-update logic
      real(8) :: Gw_LDON_init, Gw_RDON_init, Gw_NH4_init, Gw_NO3_init
      real(8) :: TotalFlux_LDON, TotalFlux_RDON, TotalFlux_NH4, TotalFlux_NO3
      real(8) :: Scale_LDON, Scale_RDON, Scale_NH4, Scale_NO3
      real(8) :: Scale_gwso  !BK20260320

      dt = real(SedCNPmodel%SedCNP_timestep)
      dx = SedCNP_hydro%dx(1)    !BK20260301

      allocate(Qgwso(domain%ix,domain%jx),stat=istat)  !BK20260316
      allocate(Qgwso_bas(domain%nbasin),stat=istat)    !BK20260316

      ! ====================================================================
      ! I. AGGREGATE INFLOWS FROM SOIL
      ! ====================================================================

      Gw_LDONsogw0(gwid) = 0.0  !BK20260316
      Gw_RDONsogw0(gwid) = 0.0  !BK20260316
      Gw_NH4sogw0(gwid) = 0.0  !BK20260316
      Gw_NO3sogw0(gwid) = 0.0  !BK20260316

      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               Gw_LDONsogw0(gwid) = Gw_LDONsogw0(gwid) + So_LDONsogw0(i,j)  !BK20260316
               Gw_RDONsogw0(gwid) = Gw_RDONsogw0(gwid) + So_RDONsogw0(i,j)  !BK20260316
               Gw_NH4sogw0(gwid)  = Gw_NH4sogw0(gwid) + So_NH4sogw0(i,j)  !BK20260316
               Gw_NO3sogw0(gwid)  = Gw_NO3sogw0(gwid) + So_NO3sogw0(i,j)  !BK20260316
            END IF
         END DO
      END DO

      ! ====================================================================
      ! II. INITIALIZE STATE AND ADD INFLOWS
      ! ====================================================================

      Gw_LDON_init = Gw_LDON0(gwid) + Gw_LDONsogw0(gwid) !BK20260316
      Gw_RDON_init = Gw_RDON0(gwid) + Gw_RDONsogw0(gwid) !BK20260316
      Gw_NH4_init  = Gw_NH40(gwid)  + Gw_NH4sogw0(gwid) !BK20260316
      Gw_NO3_init  = Gw_NO30(gwid)  + Gw_NO3sogw0(gwid) !BK20260316

      ! ====================================================================
      ! III. CALCULATE TRANSPORT FLUXES 
      ! ====================================================================

      ! --- groundwater abstraction for irrigation !BK20260316
      Qgwso_bas(gwid) = 0.0  
      Gw_LDOCgwso(gwid) = 0.0  
      Gw_RDOCgwso(gwid) = 0.0
      Gw_NH4gwso(gwid) = 0.0
      Gw_NO3gwso(gwid) = 0.0

      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               ! Clamp irridepth to >= 0 to exclude fill values (e.g. -9999)
               ! Qgwso(i,j) = SedCNP_hydro%irridepth(i,j) / dt  ! [mm s-1]
               Qgwso(i,j) = max(0.0d0, SedCNP_hydro%irridepth(i,j)) / dt  ! [mm s-1] !BK20260508
               Qgwso_bas(gwid) = Qgwso_bas(gwid) + Qgwso(i,j) / 1000.0 &
                     * domain%areaxy * dt  ! [m3]
            END IF
         END DO
      END DO

      ! if (Gw_WStorage > Qgwso_bas(gwid)) then
      if (Gw_WStorage > 0.0d0 .and. Gw_WStorage > Qgwso_bas(gwid)) then  !BK20260508
         Gw_LDONgwso(gwid) = Gw_LDON_init * Qgwso_bas(gwid) / Gw_WStorage
         Gw_RDONgwso(gwid) = Gw_RDON_init * Qgwso_bas(gwid) / Gw_WStorage
         Gw_NH4gwso(gwid) = Gw_NH4_init * Qgwso_bas(gwid) / Gw_WStorage
         Gw_NO3gwso(gwid) = Gw_NO3_init * Qgwso_bas(gwid) / Gw_WStorage
      ! else
      else if (Gw_WStorage > 0.0d0) then  !BK20260508
         Gw_LDONgwso(gwid) = Gw_LDON_init
         Gw_RDONgwso(gwid) = Gw_RDON_init
         Gw_NH4gwso(gwid) = Gw_NH4_init
         Gw_NO3gwso(gwid) = Gw_NO3_init
      else                         !---BK20260508
         Gw_LDONgwso(gwid) = 0.0d0
         Gw_RDONgwso(gwid) = 0.0d0
         Gw_NH4gwso(gwid) = 0.0d0
         Gw_NO3gwso(gwid) = 0.0d0  !---BK20260508
      endif

      ! --- groundwater discharge to channel
      q_out = Qgwch / 1000. * domain%areaxy * dt  ! (m3)
      !DOM
      dPI = calc_dPI(etaDOMgwch, q_out, Gw_WStorage)  !BK20260302
      Gw_LDONgwch(gwid) = Gw_LDON_init * dPI
      Gw_RDONgwch(gwid) = Gw_RDON_init * dPI
      !NH4
      dPI = calc_dPI(etaNH4gwch, q_out, Gw_WStorage)  !BK20260302
      Gw_NH4gwch(gwid) = Gw_NH4_init * dPI
      !NO3
      dPI = calc_dPI(etaNO3gwch, q_out, Gw_WStorage)  !BK20260302
      Gw_NO3gwch(gwid) = Gw_NO3_init * dPI

      ! ====================================================================
      ! IV. SCALE FLUXES TO ENSURE MASS CONSERVATION
      ! ====================================================================

      ! --- (1) scale flux *gwso(gwid) and *gwch(gwid) at the subbasin-level      
      TotalFlux_LDON = Gw_LDONgwso(gwid) + Gw_LDONgwch(gwid)
      if (TotalFlux_LDON > Gw_LDON_init .and. TotalFlux_LDON > 0.0) then
         Scale_LDON = Gw_LDON_init / TotalFlux_LDON
         Gw_LDONgwso(gwid) = Gw_LDONgwso(gwid) * Scale_LDON
         Gw_LDONgwch(gwid) = Gw_LDONgwch(gwid) * Scale_LDON
      endif

      TotalFlux_RDON = Gw_RDONgwso(gwid) + Gw_RDONgwch(gwid)
      if (TotalFlux_RDON > Gw_RDON_init .and. TotalFlux_RDON > 0.0) then
         Scale_RDON = Gw_RDON_init / TotalFlux_RDON
         Gw_RDONgwso(gwid) = Gw_RDONgwso(gwid) * Scale_RDON
         Gw_RDONgwch(gwid) = Gw_RDONgwch(gwid) * Scale_RDON
      endif

      TotalFlux_NH4 = Gw_NH4gwso(gwid) + Gw_NH4gwch(gwid)
      if (TotalFlux_NH4 > Gw_NH4_init .and. TotalFlux_NH4 > 0.0) then
         Scale_NH4 = Gw_NH4_init / TotalFlux_NH4
         Gw_NH4gwso(gwid) = Gw_NH4gwso(gwid) * Scale_NH4
         Gw_NH4gwch(gwid) = Gw_NH4gwch(gwid) * Scale_NH4
      endif

      TotalFlux_NO3 = Gw_NO3gwso(gwid) + Gw_NO3gwch(gwid)
      if (TotalFlux_NO3 > Gw_NO3_init .and. TotalFlux_NO3 > 0.0) then
         Scale_NO3 = Gw_NO3_init / TotalFlux_NO3
         Gw_NO3gwso(gwid) = Gw_NO3gwso(gwid) * Scale_NO3
         Gw_NO3gwch(gwid) = Gw_NO3gwch(gwid) * Scale_NO3
      endif

      ! --- (2) scale flux *gwso(i,j) at the grid-level  !BK20260316
      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               if (Qgwso_bas(gwid) > 1.0e-9) then
                  Scale_gwso = (Qgwso(i,j) / 1000.0 * domain%areaxy * dt) / Qgwso_bas(gwid)
               else
                  Scale_gwso = 0.0
               endif
               So_LDONgwso(i,j) = Gw_LDONgwso(gwid) * Scale_gwso
               So_RDONgwso(i,j) = Gw_RDONgwso(gwid) * Scale_gwso
               So_NH4gwso(i,j) = Gw_NH4gwso(gwid) * Scale_gwso
               So_NO3gwso(i,j) = Gw_NO3gwso(gwid) * Scale_gwso 
            END IF
         END DO
      END DO

      ! ====================================================================
      ! V. UPDATE FINAL STORAGES
      ! ====================================================================

      Gw_LDON(gwid) = Gw_LDON_init - Gw_LDONgwso(gwid) - Gw_LDONgwch(gwid)
      Gw_RDON(gwid) = Gw_RDON_init - Gw_RDONgwso(gwid) - Gw_RDONgwch(gwid)
      Gw_NH4(gwid)  = Gw_NH4_init  - Gw_NH4gwso(gwid)  - Gw_NH4gwch(gwid)
      Gw_NO3(gwid)  = Gw_NO3_init  - Gw_NO3gwso(gwid)  - Gw_NO3gwch(gwid)

      ! Enforce non-negative groundwater storages
      if (Gw_LDON(gwid) < 1.0e-20) Gw_LDON(gwid) = 0.0
      if (Gw_RDON(gwid) < 1.0e-20) Gw_RDON(gwid) = 0.0
      if (Gw_NH4(gwid) < 1.0e-20) Gw_NH4(gwid) = 0.0
      if (Gw_NO3(gwid) < 1.0e-20) Gw_NO3(gwid) = 0.0

   end subroutine GroundwaterNitrogenCycle


   subroutine ChannelNitrogenCycle(ich, itime, Qdsch, WtopWdth, WStorage)

      implicit none

      ! --- Input arguments
      integer, intent(in) :: ich, itime
      real(8), intent(in)  :: Qdsch, WtopWdth, WStorage

      ! --- Local variables
      integer :: i, j, isl, ips, ia, im, month, stat, chid, gwid
      real(8) :: rate0, dW, dZ, dXi1, dXi2, dLambda_L, dLambda_N, dLambda_P
      real(8) :: dPHI_N, dPHI_P, dLPOND, dLPONZ, dLPONS, dLDOND, dRPOND
      real(8) :: dRPONS, dRDOND, dMBMND, dAZP, dALGCconc, dTSSconc
      real(8) :: dTA(ntm), dMA(ntm), dt, abs_ratio
      real(8) :: RK4SOL, KTBMBM, KTBZOO, KTBALG, KTHETA, ALPHA, dPI
      real(8) :: KZG(domain%nch), KZR(domain%nch), KZM(domain%nch)
      real(8) :: KAG(domain%nch,nalg), KAR(domain%nch,nalg)
      real(8) :: KAE(domain%nch,nalg), KAM(domain%nch,nalg)
      real(8), PARAMETER :: e = 2.718281828459
      character(len=2) :: monthstr

      ! --- Variables for calculate-scale-update logic
      real(8) :: Ch_ALGN_init(nalg), Ch_ZOON_init, Ch_LPON_init
      real(8) :: Ch_LDON_init, Ch_RPON_init, Ch_RDON_init, Ch_MBMN_init
      real(8) :: Ch_NH4_init, Ch_NO3_init
      real(8) :: TotalFlux_ALGN(nalg), Scale_ALGN(nalg)
      real(8) :: TotalFlux_ZOON, Scale_ZOON
      real(8) :: TotalFlux_LPON, Scale_LPON
      real(8) :: TotalFlux_RPON, Scale_RPON
      real(8) :: TotalFlux_LDON, Scale_LDON
      real(8) :: TotalFlux_RDON, Scale_RDON
      real(8) :: TotalFlux_MBMN, Scale_MBMN
      real(8) :: TotalFlux_NH4, Scale_NH4
      real(8) :: TotalFlux_NO3, Scale_NO3

      dt = real(SedCNPmodel%SedCNP_timestep)
      chid = SedCNP_hydro%linkID(ich)
      gwid = domain%gwid_ch(chid)

      ! ================================================================
      ! I. INFLOWS
      ! ================================================================

      ! --- Inflow from upstream channels
      Ch_ALGN(ich,1:nalg) = Ch_ALGN0(ich,1:nalg) + Ch_ALGNusch0(ich,1:nalg)
      Ch_ZOON(ich) = Ch_ZOON0(ich) + Ch_ZOONusch0(ich)
      Ch_MBMN(ich) = Ch_MBMN0(ich) + Ch_MBMNusch0(ich)
      Ch_LPON(ich) = Ch_LPON0(ich) + Ch_LPONusch0(ich)
      Ch_RPON(ich) = Ch_RPON0(ich) + Ch_RPONusch0(ich)
      Ch_LDON(ich) = Ch_LDON0(ich) + Ch_LDONusch0(ich)
      Ch_RDON(ich) = Ch_RDON0(ich) + Ch_RDONusch0(ich)
      Ch_NH4(ich)  = Ch_NH40(ich)  + Ch_NH4usch0(ich)
      Ch_NO3(ich)  = Ch_NO30(ich)  + Ch_NO3usch0(ich)

      ! --- Lateral inflows (surface, interflow, groundwater)
      IF (gwid >= 0) THEN

         ! initialisation
         Ch_LPONsurf0(ich) = 0.0
         Ch_RPONsurf0(ich) = 0.0
         Ch_LDONsurf0(ich) = 0.0
         Ch_RDONsurf0(ich) = 0.0
         Ch_MBMNsurf0(ich) = 0.0
         Ch_NH4surf0(ich)  = 0.0
         Ch_NO3surf0(ich)  = 0.0
         Ch_LDONintf0(ich) = 0.0
         Ch_RDONintf0(ich) = 0.0
         Ch_NH4intf0(ich)  = 0.0
         Ch_NO3intf0(ich)  = 0.0

         DO i = 1, domain%ix
            DO j = 1, domain%jx
               IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN

                  ! surface runoff
                  Ch_LPONsurf0(ich) = Ch_LPONsurf0(ich) + So_LPONsurf0(i,j)
                  Ch_RPONsurf0(ich) = Ch_RPONsurf0(ich) + So_RPONsurf0(i,j)
                  Ch_LDONsurf0(ich) = Ch_LDONsurf0(ich) + So_LDONsurf0(i,j)
                  Ch_RDONsurf0(ich) = Ch_RDONsurf0(ich) + So_RDONsurf0(i,j)
                  Ch_MBMNsurf0(ich) = Ch_MBMNsurf0(ich) + So_MBMNsurf0(i,j)
                  Ch_NH4surf0(ich)  = Ch_NH4surf0(ich)  + So_NH4surf0(i,j)
                  Ch_NO3surf0(ich)  = Ch_NO3surf0(ich)  + So_NO3surf0(i,j)

                  ! interflow
                  DO isl = 1, domain%nsl
                     Ch_LDONintf0(ich) = Ch_LDONintf0(ich) + &
                        So_LDONintf0(i,isl,j)
                     Ch_RDONintf0(ich) = Ch_RDONintf0(ich) + &
                        So_RDONintf0(i,isl,j)
                     Ch_NH4intf0(ich) = Ch_NH4intf0(ich) + &
                        So_NH4intf0(i,isl,j)
                     Ch_NO3intf0(ich) = Ch_NO3intf0(ich) + &
                        So_NO3intf0(i,isl,j)
                  END DO

               END IF
            END DO
         END DO

         ! groundwater inflow
         Ch_LDONgwch0(ich) = Gw_LDONgwch0(gwid)
         Ch_RDONgwch0(ich) = Gw_RDONgwch0(gwid)
         Ch_NH4gwch0(ich)  = Gw_NH4gwch0(gwid)
         Ch_NO3gwch0(ich)  = Gw_NO3gwch0(gwid)

         ! update storage
         Ch_LPON(ich) = Ch_LPON(ich) + Ch_LPONsurf0(ich)
         Ch_RPON(ich) = Ch_RPON(ich) + Ch_RPONsurf0(ich)
         Ch_LDON(ich) = Ch_LDON(ich) + Ch_LDONsurf0(ich) + Ch_LDONintf0(ich) + &
            Ch_LDONgwch0(ich)
         Ch_RDON(ich) = Ch_RDON(ich) + Ch_RDONsurf0(ich) + Ch_RDONintf0(ich) + &
            Ch_RDONgwch0(ich)
         Ch_MBMN(ich) = Ch_MBMN(ich) + Ch_MBMNsurf0(ich)
         Ch_NH4(ich)  = Ch_NH4(ich)  + Ch_NH4surf0(ich) + Ch_NH4intf0(ich) + &
            Ch_NH4gwch0(ich)
         Ch_NO3(ich)  = Ch_NO3(ich)  + Ch_NO3surf0(ich) + Ch_NO3intf0(ich) + &
            Ch_NO3gwch0(ich)

      ENDIF

      ! --- Forced inflow/ourflow 
      ! inflow: point source
      Ch_LPONxpnt(ich) = Ch_LPONpnt(ich,itime) * dt
      Ch_RPONxpnt(ich) = Ch_RPONpnt(ich,itime) * dt
      Ch_LDONxpnt(ich) = Ch_LDONpnt(ich,itime) * dt
      Ch_RDONxpnt(ich) = Ch_RDONpnt(ich,itime) * dt
      Ch_NH4xpnt(ich)  = Ch_NH4pnt(ich,itime)  * dt
      Ch_NO3xpnt(ich)  = Ch_NO3pnt(ich,itime)  * dt

      ! inflow: water discharge
      Ch_LPONxdis(ich) = Ch_LPONdis(ich,itime) * dt
      Ch_RPONxdis(ich) = Ch_RPONdis(ich,itime) * dt
      Ch_LDONxdis(ich) = Ch_LDONdis(ich,itime) * dt
      Ch_RDONxdis(ich) = Ch_RDONdis(ich,itime) * dt
      Ch_NH4xdis(ich)  = Ch_NH4dis(ich,itime)  * dt
      Ch_NO3xdis(ich)  = Ch_NO3dis(ich,itime)  * dt

      ! update storage (inflows)
      Ch_LPON(ich) = Ch_LPON(ich) + Ch_LPONxpnt(ich) + Ch_LPONxdis(ich)
      Ch_RPON(ich) = Ch_RPON(ich) + Ch_RPONxpnt(ich) + Ch_RPONxdis(ich)
      Ch_LDON(ich) = Ch_LDON(ich) + Ch_LDONxpnt(ich) + Ch_LDONxdis(ich)
      Ch_RDON(ich) = Ch_RDON(ich) + Ch_RDONxpnt(ich) + Ch_RDONxdis(ich)
      Ch_NH4(ich)  = Ch_NH4(ich)  + Ch_NH4xpnt(ich)  + Ch_NH4xdis(ich)
      Ch_NO3(ich)  = Ch_NO3(ich)  + Ch_NO3xpnt(ich)  + Ch_NO3xdis(ich)

      ! outflow: water abstraction
      abs_ratio = SedCNP_hydro%Qabs(ich,itime) / (WStorage + Qdsch * dt)
      if(SedCNP_hydro%Qabs(ich,itime) .gt. 0.) then
         Ch_LPONabs(ich,itime) = Ch_LPON(ich) * abs_ratio
         Ch_RPONabs(ich,itime) = Ch_RPON(ich) * abs_ratio
         Ch_LDONabs(ich,itime) = Ch_LDON(ich) * abs_ratio
         Ch_RDONabs(ich,itime) = Ch_RDON(ich) * abs_ratio
         Ch_NH4abs(ich,itime)  = Ch_NH4(ich)  * abs_ratio
         Ch_NO3abs(ich,itime)  = Ch_NO3(ich)  * abs_ratio   
      endif
      Ch_LPONxabs(ich) = Ch_LPONabs(ich,itime) * dt  
      Ch_RPONxabs(ich) = Ch_RPONabs(ich,itime) * dt  
      Ch_LDONxabs(ich) = Ch_LDONabs(ich,itime) * dt  
      Ch_RDONxabs(ich) = Ch_RDONabs(ich,itime) * dt  
      Ch_NH4xabs(ich)  = Ch_NH4abs(ich,itime)  * dt
      Ch_NO3xabs(ich)  = Ch_NO3abs(ich,itime)  * dt

      ! update storage (outflows)
      Ch_LPON(ich) = Ch_LPON(ich) - Ch_LPONxabs(ich)
      Ch_RPON(ich) = Ch_RPON(ich) - Ch_RPONxabs(ich)
      Ch_LDON(ich) = Ch_LDON(ich) - Ch_LDONxabs(ich)
      Ch_RDON(ich) = Ch_RDON(ich) - Ch_RDONxabs(ich)
      Ch_NH4(ich)  = Ch_NH4(ich)  - Ch_NH4xabs(ich)
      Ch_NO3(ich)  = Ch_NO3(ich)  - Ch_NO3xabs(ich)

      ! Enforce non-negative channel storages after abstraction
      if (Ch_LPON(ich) < 1.0e-20) Ch_LPON(ich) = 0.0
      if (Ch_RPON(ich) < 1.0e-20) Ch_RPON(ich) = 0.0
      if (Ch_LDON(ich) < 1.0e-20) Ch_LDON(ich) = 0.0
      if (Ch_RDON(ich) < 1.0e-20) Ch_RDON(ich) = 0.0
      if (Ch_NH4(ich) < 1.0e-20)  Ch_NH4(ich)  = 0.0
      if (Ch_NO3(ich) < 1.0e-20)  Ch_NO3(ich)  = 0.0

      ! ================================================================
      ! II. CALCULATE-SCALE-UPDATE
      ! ================================================================

      ! --- Store initial values for simultaneous flux calculation
      Ch_ALGN_init(:) = Ch_ALGN(ich,:)
      Ch_ZOON_init    = Ch_ZOON(ich)
      Ch_LPON_init    = Ch_LPON(ich)
      Ch_LDON_init    = Ch_LDON(ich)
      Ch_RPON_init    = Ch_RPON(ich)
      Ch_RDON_init    = Ch_RDON(ich)
      Ch_MBMN_init    = Ch_MBMN(ich)
      Ch_NH4_init     = Ch_NH4(ich)
      Ch_NO3_init     = Ch_NO3(ich)

      ! === A. CALCULATE ALL FLUXES BASED ON INITIAL STORAGES

      ! --- Pre-calculations for transformation rates
      KTHETA = 0.6 ! Assumed constant for in-stream saturated conditions
      dW = DWATER(ich)
      dZ = min(dW, 10.0)  ! take account of algaes in the top 10 m of water

      ! ALPHA - light attenuation coeff. (USEPA, 2000)
      if (WStorage > 0.0) then
         dALGCconc = sum(Ch_ALGC0(ich,1:nalg)) / WStorage
         dTSSconc = (sum(channelSed%Cchsd(1:nps,ich)) + Ch_LPOC0(ich) + &
            Ch_RPOC0(ich)) / WStorage
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
         dLambda_L = e / (ALPHA*dZ) * (exp(-dXi2) - exp(-dXi1))
      else
         dLambda_L = 0.0
      endif

      if (WStorage > 0.0) then
         dPHI_N = (Ch_NH40(ich) + Ch_NO30(ich)) / WStorage
         dPHI_P = Ch_PO40(ich) / WStorage
      else
         dPHI_N = 0.0
         dPHI_P = 0.0
      end if

      if ((MONOD_N + dPHI_N) > 0.0) then
         dLambda_N = dPHI_N / (MONOD_N + dPHI_N)
      else
         dLambda_N = 0.0
      endif

      if ((MONOD_P + dPHI_P) > 0.0) then
         dLambda_P = dPHI_P / (MONOD_P + dPHI_P)
      else
         dLambda_P = 0.0
      endif

      ! --- Zooplankton
      KTBZOO = calc_KTB(TWATER(ich), TZOO, MZOO)
      dAZP = sum(Ch_ALGN_init(:)) + Ch_LPON_init + Ch_ZOON_init

      ! growth
      if (dAZP + ZHALF * WStorage > 1.0E-9) then
         KZG(ich) = KTBZOO * EZI * KZIMAX * (dAZP - ZLOW * WStorage) / &
            (dAZP + ZHALF * WStorage)
      else
         KZG(ich) = 0.0
      endif
      RK4SOL = calc_RK4SOL(KZG(ich), Ch_ZOON_init)
      Ch_ZOONG(ich) = max(0.0, RK4SOL - Ch_ZOON_init)

      ! being grazed by (other) zooplankton
      if (EZI > 0.0 .and. dAZP > 1.0E-9) then
         Ch_ZOONZ(ich) = (Ch_ZOONG(ich) / EZI) * (Ch_ZOON_init / dAZP)
      else
         Ch_ZOONZ(ich) = 0.0
      endif

      ! respiration
      KZR(ich) = KTBZOO * KZRMAX
      RK4SOL = calc_RK4SOL(-KZR(ich), Ch_ZOON_init)
      Ch_ZOONR(ich) = Ch_ZOON_init - RK4SOL

      ! excretion - based on growth to ensure mass balance (BK20251127)      
      if (EZI > 0.0) then
         Ch_ZOONE(ich) = Ch_ZOONG(ich) * (1.0 - EZI) / EZI
      else
         Ch_ZOONE(ich) = 0.0
      endif

      ! mortality 
      KZM(ich) = KTBZOO * KZMMAX
      RK4SOL = calc_RK4SOL(-KZM(ich), Ch_ZOON_init)
      Ch_ZOONM(ich) = Ch_ZOON_init - RK4SOL

      ! settling      
      if (dW > 0.0) then
         Ch_ZOONS(ich) = (1.0 - exp(-OMEGAZOO * dt / dW)) * Ch_ZOON_init
      else
         Ch_ZOONS(ich) = 0.0
      endif

      ! --- Phytoplanktons
      DO ia = 1, nalg
         do im = 1, ntm
            dTA(im) = TALG(ia,im)
            dMA(im) = MALG(ia,im)
         enddo         
         KTBALG = calc_KTB(TWATER(ich), dTA, dMA)

         ! growth
         KAG(ich,ia) = KTBALG * min(dLambda_L, dLambda_N, dLambda_P) * KAGMAX(ia)
         RK4SOL = calc_RK4SOL(KAG(ich,ia), Ch_ALGN_init(ia))
         Ch_ALGNG(ich,ia) = max(0.0, RK4SOL - Ch_ALGN_init(ia))

         ! being grazed by zooplankton
         if (EZI > 0.0 .and. dAZP > 1.0E-9) then
            Ch_ALGNZ(ich,ia) = &
               (Ch_ZOONG(ich) / EZI) * (Ch_ALGN_init(ia) / dAZP)
         else
            Ch_ALGNZ(ich,ia) = 0.0
         endif

         ! respiration
         KAR(ich,ia) = KTBALG * KARMAX(ia)
         RK4SOL = calc_RK4SOL(-KAR(ich,ia), Ch_ALGN_init(ia))
         Ch_ALGNR(ich,ia) = Ch_ALGN_init(ia) - RK4SOL

         ! excretion
         KAE(ich,ia) = KTBALG * (1.0 - dLambda_L) * KAEMAX(ia)
         RK4SOL = calc_RK4SOL(-KAE(ich,ia), Ch_ALGN_init(ia))
         Ch_ALGNE(ich,ia) = Ch_ALGN_init(ia) - RK4SOL

         ! mortality
         KAM(ich,ia) = KTBALG * KAMMAX(ia)
         RK4SOL = calc_RK4SOL(-KAM(ich,ia), Ch_ALGN_init(ia))
         Ch_ALGNM(ich,ia) = Ch_ALGN_init(ia) - RK4SOL

         ! settling
         if (dW > 0.0) then
            Ch_ALGNS(ich,ia) = (1.0 - exp(-OMEGAALG * dt / dW)) * &
               Ch_ALGN_init(ia)
         else
            Ch_ALGNS(ich,ia) = 0.0
         endif
      END DO

      ! --- Organic matter
      KTBMBM = calc_KTB(TWATER(ich), TMBM, MMBM)

      ! LPON grazed by zooplankton
      if (EZI > 0.0 .and. dAZP > 1.0E-9) then
         dLPONZ = (Ch_ZOONG(ich) / EZI) * (Ch_LPON_init / dAZP)
      else
         dLPONZ = 0.0
      endif

      ! LPON decomposition
      rate0 = KTBMBM * KTHETA * KLPOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_LPON_init)
      dLPOND = Ch_LPON_init - RK4SOL

      ! LPON settling
      if (dW > 0.0) then
         dLPONS = (1.0 - exp(-OMEGAPOM * dt / dW)) * Ch_LPON_init
      else
         dLPONS = 0.0
      endif

      ! RPON decomposition
      rate0 = KTBMBM * KTHETA * KRPOMW
      RK4SOL = calc_RK4SOL(-rate0,Ch_RPON_init)
      dRPOND = Ch_RPON_init - RK4SOL

      ! RPON settling
      if (dW > 0.0) then
         dRPONS = (1.0 - exp(-OMEGAPOM * dt / dW)) * Ch_RPON_init
      else
         dRPONS = 0.0
      endif

      ! LDON decomposition
      rate0 = KTBMBM * KTHETA * KLDOMW
      RK4SOL = calc_RK4SOL(-rate0,Ch_LDON_init)
      dLDOND = Ch_LDON_init - RK4SOL

      ! RDON decomposition
      rate0 = KTBMBM * KTHETA * KRDOMW
      RK4SOL = calc_RK4SOL(-rate0,Ch_RDON_init)
      dRDOND = Ch_RDON_init - RK4SOL

      ! MBMN decompositio
      rate0 = KTBMBM * KTHETA * KMBMW
      RK4SOL = calc_RK4SOL(-rate0,Ch_MBMN_init)
      dMBMND = Ch_MBMN_init - RK4SOL

      ! --- Inorganic N 
      ! uptake by phytoplanktons
      Ch_NH4uptk(ich) = 0.0
      Ch_NO3uptk(ich) = 0.0
      if ((Ch_NH4_init + Ch_NO3_init) > 0.0) then
         DO ia = 1, nalg
            Ch_NH4uptk(ich) = Ch_NH4uptk(ich) + Ch_ALGNG(ich,ia) * &
               Ch_NH4_init / (Ch_NH4_init + Ch_NO3_init)
            Ch_NO3uptk(ich) = Ch_NO3uptk(ich) + Ch_ALGNG(ich,ia) * &
               Ch_NO3_init / (Ch_NH4_init + Ch_NO3_init)
         END DO
      endif

      ! nitrification
      rate0 = KTBMBM * KTHETA * KNITW
      RK4SOL = calc_RK4SOL(-rate0, Ch_NH4_init)
      Ch_NITRI(ich) = Ch_NH4_init - RK4SOL

      ! NH3 volatilization
      rate0 = KTBMBM * KTHETA * KVOLW
      RK4SOL = calc_RK4SOL(-rate0, Ch_NH4_init)
      Ch_NH3VOL(ich) = Ch_NH4_init - RK4SOL

      ! denitrification
      rate0 = KTBMBM * KTHETA * KDENW
      RK4SOL = calc_RK4SOL(-rate0, Ch_NO3_init)
      Ch_DENIT(ich) = Ch_NO3_init - RK4SOL

      ! --- outflow: downstream discharge
      if (WStorage > 0.0 .and. Qdsch > 0.0) then
         dPI = 1.0 - EXP(-2.3 * (Qdsch * dt) / (WStorage + Qdsch * dt))
      else
         dPI = 0.0
      endif
      Ch_ALGNdsch(ich,:) = Ch_ALGN_init(:) * dPI
      Ch_ZOONdsch(ich)   = Ch_ZOON_init * dPI
      Ch_LPONdsch(ich)   = Ch_LPON_init * dPI
      Ch_RPONdsch(ich)   = Ch_RPON_init * dPI
      Ch_MBMNdsch(ich)   = Ch_MBMN_init * dPI
      Ch_LDONdsch(ich)   = Ch_LDON_init * dPI
      Ch_RDONdsch(ich)   = Ch_RDON_init * dPI
      Ch_NH4dsch(ich)    = Ch_NH4_init * dPI
      Ch_NO3dsch(ich)    = Ch_NO3_init * dPI

      ! === B. SCALE FLUXES TO ENSURE MASS CONSERVATION

      ! Zooplankton
      TotalFlux_ZOON = Ch_ZOONZ(ich) + Ch_ZOONR(ich) + Ch_ZOONM(ich) + &
         Ch_ZOONE(ich) + Ch_ZOONS(ich) + Ch_ZOONdsch(ich)
      if (TotalFlux_ZOON > Ch_ZOON_init .and. TotalFlux_ZOON > 0.0) then
         Scale_ZOON = Ch_ZOON_init / TotalFlux_ZOON
      else
         Scale_ZOON = 1.0
      endif
      Ch_ZOONZ(ich)    = Ch_ZOONZ(ich) * Scale_ZOON
      Ch_ZOONR(ich)    = Ch_ZOONR(ich) * Scale_ZOON
      Ch_ZOONE(ich)    = Ch_ZOONE(ich) * Scale_ZOON
      Ch_ZOONM(ich)    = Ch_ZOONM(ich) * Scale_ZOON
      Ch_ZOONS(ich)    = Ch_ZOONS(ich) * Scale_ZOON
      Ch_ZOONdsch(ich) = Ch_ZOONdsch(ich) * Scale_ZOON

      ! LPON
      TotalFlux_LPON = dLPONZ + dLPOND + dLPONS + Ch_LPONdsch(ich)
      if (TotalFlux_LPON > Ch_LPON_init .and. TotalFlux_LPON > 0.0) then
         Scale_LPON = Ch_LPON_init / TotalFlux_LPON
      else
         Scale_LPON = 1.0
      endif
      dLPONZ           = dLPONZ * Scale_LPON
      dLPOND           = dLPOND * Scale_LPON
      dLPONS           = dLPONS * Scale_LPON
      Ch_LPONdsch(ich) = Ch_LPONdsch(ich) * Scale_LPON

      ! LDON
      TotalFlux_LDON = dLDOND + Ch_LDONdsch(ich)
      if (TotalFlux_LDON > Ch_LDON_init .and. TotalFlux_LDON > 0.0) then
         Scale_LDON = Ch_LDON_init / TotalFlux_LDON
      else
         Scale_LDON = 1.0
      endif
      dLDOND           = dLDOND * Scale_LDON
      Ch_LDONdsch(ich) = Ch_LDONdsch(ich) * Scale_LDON

      ! RPON
      TotalFlux_RPON = dRPOND + dRPONS + Ch_RPONdsch(ich)
      if (TotalFlux_RPON > Ch_RPON_init .and. TotalFlux_RPON > 0.0) then
         Scale_RPON = Ch_RPON_init / TotalFlux_RPON
      else
         Scale_RPON = 1.0
      endif
      dRPOND           = dRPOND * Scale_RPON
      dRPONS           = dRPONS * Scale_RPON
      Ch_RPONdsch(ich) = Ch_RPONdsch(ich) * Scale_RPON

      ! RDON
      TotalFlux_RDON = dRDOND + Ch_RDONdsch(ich)
      if (TotalFlux_RDON > Ch_RDON_init .and. TotalFlux_RDON > 0.0) then
         Scale_RDON = Ch_RDON_init / TotalFlux_RDON
      else
         Scale_RDON = 1.0
      endif
      dRDOND           = dRDOND * Scale_RDON
      Ch_RDONdsch(ich) = Ch_RDONdsch(ich) * Scale_RDON

      ! MBMN
      TotalFlux_MBMN = dMBMND + Ch_MBMNdsch(ich)
      if (TotalFlux_MBMN > Ch_MBMN_init .and. TotalFlux_MBMN > 0.0) then
         Scale_MBMN = Ch_MBMN_init / TotalFlux_MBMN
      else
         Scale_MBMN = 1.0
      endif
      dMBMND           = dMBMND * Scale_MBMN
      Ch_MBMNdsch(ich) = Ch_MBMNdsch(ich) * Scale_MBMN

      ! NH4
      TotalFlux_NH4 = Ch_NH4uptk(ich) + Ch_NITRI(ich) + Ch_NH3VOL(ich) + &
         Ch_NH4dsch(ich)
      if (TotalFlux_NH4 > Ch_NH4_init .and. TotalFlux_NH4 > 0.0) then
         Scale_NH4 = Ch_NH4_init / TotalFlux_NH4
         Ch_NH4uptk(ich) = Ch_NH4uptk(ich) * Scale_NH4
         Ch_NITRI(ich)   = Ch_NITRI(ich)   * Scale_NH4
         Ch_NH3VOL(ich)  = Ch_NH3VOL(ich)  * Scale_NH4
         Ch_NH4dsch(ich) = Ch_NH4dsch(ich) * Scale_NH4
      else
         Scale_NH4 = 1.0
      endif

      ! NO3
      TotalFlux_NO3 = Ch_NO3uptk(ich) + Ch_DENIT(ich) + Ch_NO3dsch(ich)
      if (TotalFlux_NO3 > Ch_NO3_init .and. TotalFlux_NO3 > 0.0) then
         Scale_NO3 = Ch_NO3_init / TotalFlux_NO3
         Ch_NO3uptk(ich) = Ch_NO3uptk(ich) * Scale_NO3
         Ch_DENIT(ich)   = Ch_DENIT(ich)   * Scale_NO3
         Ch_NO3dsch(ich) = Ch_NO3dsch(ich) * Scale_NO3
      else
         Scale_NO3 = 1.0
      endif

      ! Scale phytoplankton growth proportionally to nutrient uptake
      if ((Ch_NH4_init + Ch_NO3_init) > 0.0) then
         Ch_ALGNG(ich,:) = Ch_ALGNG(ich,:) * (Ch_NH4_init / (Ch_NH4_init + Ch_NO3_init) * Scale_NH4 + &
                                             Ch_NO3_init / (Ch_NH4_init + Ch_NO3_init) * Scale_NO3)
      else
         Ch_ALGNG(ich,:) = 0.0
      endif

      ! Phytoplanktons
      DO ia=1,nalg
         TotalFlux_ALGN(ia) = Ch_ALGNZ(ich,ia) + Ch_ALGNR(ich,ia) + &
            Ch_ALGNE(ich,ia) + Ch_ALGNM(ich,ia) + Ch_ALGNS(ich,ia) + &
            Ch_ALGNdsch(ich,ia)
        if (TotalFlux_ALGN(ia) > (Ch_ALGN_init(ia) + Ch_ALGNG(ich,ia))) then
            if ((Ch_ALGN_init(ia) + Ch_ALGNG(ich,ia)) > 1.0E-20) then
               Scale_ALGN(ia) = (Ch_ALGN_init(ia) + Ch_ALGNG(ich,ia)) / TotalFlux_ALGN(ia)
            else
               Scale_ALGN(ia) = 0.0
            endif
         else
            Scale_ALGN(ia) = 1.0
         endif
         Ch_ALGNZ(ich,ia)    = Ch_ALGNZ(ich,ia) * Scale_ALGN(ia)
         Ch_ALGNR(ich,ia)    = Ch_ALGNR(ich,ia) * Scale_ALGN(ia)
         Ch_ALGNE(ich,ia)    = Ch_ALGNE(ich,ia) * Scale_ALGN(ia)
         Ch_ALGNM(ich,ia)    = Ch_ALGNM(ich,ia) * Scale_ALGN(ia)
         Ch_ALGNS(ich,ia)    = Ch_ALGNS(ich,ia) * Scale_ALGN(ia)
         Ch_ALGNdsch(ich,ia) = Ch_ALGNdsch(ich,ia) * Scale_ALGN(ia)
      ENDDO

      ! === C. UPDATE STORAGES BASED ON SCALED FLUXES

      do ia=1,nalg
         Ch_ALGN(ich,ia) = Ch_ALGN_init(ia) + Ch_ALGNG(ich,ia) - &
                           Ch_ALGNZ(ich,ia) - Ch_ALGNR(ich,ia) - &
                           Ch_ALGNE(ich,ia) - Ch_ALGNM(ich,ia) - &
                           Ch_ALGNS(ich,ia) - Ch_ALGNdsch(ich,ia)
      end do
      Ch_ZOON(ich) = Ch_ZOON_init + Ch_ZOONG(ich) - Ch_ZOONZ(ich) - &
                     Ch_ZOONR(ich) - Ch_ZOONM(ich) - &  ! Ch_ZOONE(ich) excluded!!!
                     Ch_ZOONS(ich) - Ch_ZOONdsch(ich)
      Ch_LPON(ich) = Ch_LPON_init + sum(Ch_ALGNM(ich,:)) + Ch_ZOONM(ich) + &
                     Ch_ZOONE(ich) - dLPONZ - dLPOND - dLPONS - &
                     Ch_LPONdsch(ich)
      Ch_LDON(ich) = Ch_LDON_init + sum(Ch_ALGNE(ich,:)) + FDEC * &
                     (1.0 - FREF) * dLPOND + FDEC * dMBMND - dLDOND - &
                     Ch_LDONdsch(ich)
      Ch_RPON(ich) = Ch_RPON_init + FDEC * FREF * dLPOND - dRPOND - dRPONS - &
                     Ch_RPONdsch(ich)
      Ch_RDON(ich) = Ch_RDON_init + FDEC * (dLDOND + dRPOND + dRDOND) - &
                     dRDOND - Ch_RDONdsch(ich)
      Ch_MBMN(ich) = Ch_MBMN_init + FBIO * (dLPOND + dRPOND + dLDOND + &
                     dRDOND + dMBMND) - dMBMND - Ch_MBMNdsch(ich)
      Ch_NH4(ich)  = Ch_NH4_init + FMET * (dLPOND + dRPOND + dLDOND + &
                     dRDOND + dMBMND) + sum(Ch_ALGNR(ich,:)) + &
                     Ch_ZOONR(ich) - (Ch_NH4uptk(ich) + Ch_NITRI(ich) + &
                     Ch_NH3VOL(ich) + Ch_NH4dsch(ich))
      Ch_NO3(ich)  = Ch_NO3_init + Ch_NITRI(ich) - (Ch_NO3uptk(ich) + &
                     Ch_DENIT(ich) + Ch_NO3dsch(ich))
      Ch_PONDEPOSIT(ich) = Ch_PONDEPOSIT0(ich) + sum(Ch_ALGNS(ich,:)) + &
                     Ch_ZOONS(ich) + dLPONS + dRPONS

      ! Enforce non-negative channel storages
      do ia=1,nalg
        if (Ch_ALGN(ich,ia) < 1.0e-20) Ch_ALGN(ich,ia) = 0.0
      end do
      if (Ch_ZOON(ich) < 1.0e-20) Ch_ZOON(ich) = 0.0
      if (Ch_LPON(ich) < 1.0e-20) Ch_LPON(ich) = 0.0
      if (Ch_LDON(ich) < 1.0e-20) Ch_LDON(ich) = 0.0
      if (Ch_RPON(ich) < 1.0e-20) Ch_RPON(ich) = 0.0
      if (Ch_RDON(ich) < 1.0e-20) Ch_RDON(ich) = 0.0
      if (Ch_MBMN(ich) < 1.0e-20) Ch_MBMN(ich) = 0.0
      if (Ch_NH4(ich) < 1.0e-20)  Ch_NH4(ich) = 0.0
      if (Ch_NO3(ich) < 1.0e-20)  Ch_NO3(ich) = 0.0
      if (Ch_PONDEPOSIT(ich) < 1.0e-20) Ch_PONDEPOSIT(ich) = 0.0

   end subroutine ChannelNitrogenCycle


   subroutine CheckNMB_So(i, j, lcover)

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
            So_NSC(i,isl,j) = &
                    (So_LPONLIT(i,j)     - So_LPONLIT0(i,j))     &
                  + (So_LPONRES(i,isl,j) - So_LPONRES0(i,isl,j)) &
                  + (So_LPONEXC(i,j)     - So_LPONEXC0(i,j))     &
                  + (So_LPONMAN(i,j)     - So_LPONMAN0(i,j))     &
                  + (So_RPON(i,isl,j)    - So_RPON0(i,isl,j))    &
                  + (So_LDON(i,isl,j)    - So_LDON0(i,isl,j))    &
                  + (So_RDON(i,isl,j)    - So_RDON0(i,isl,j))    &
                  + (So_MBMN(i,isl,j)    - So_MBMN0(i,isl,j))    &
                  + (So_NH4(i,isl,j)     - So_NH40(i,isl,j))     &
                  + (So_NO3(i,isl,j)     - So_NO30(i,isl,j))

            So_NMBError(i,isl,j) = &
                    So_LPONLITin(i,j)    + So_LPONRESin(i,isl,j) &
                  + So_LPONEXCin(i,j)    + So_LPONMANin(i,j)     &
                  + So_NH4fxin(i,j)      + So_NO3fxin(i,j)       &
                  + So_LDONgwso0(i,j)    + So_RDONgwso0(i,j)     &  !BK20260316
                  + So_NH4gwso0(i,j)     + So_NO3gwso0(i,j)      &  !BK20260316
                  - So_LPONsurf(i,j)     - So_RPONsurf(i,j)      &
                  - So_MBMNsurf(i,j)     - So_LDONsurf(i,j)      &
                  - So_RDONsurf(i,j)     - So_NH4surf(i,j)       &
                  - So_NO3surf(i,j)      - So_LDONintf(i,isl,j)  &
                  - So_RDONintf(i,isl,j) - So_NH4intf(i,isl,j)   &
                  - So_NO3intf(i,isl,j)  - So_LDONperc(i,isl,j)  &
                  - So_RDONperc(i,isl,j) - So_NH4perc(i,isl,j)   &
                  - So_NO3perc(i,isl,j)  - So_NH4uptk(i,isl,j)   &
                  - So_NO3uptk(i,isl,j)  - So_NH3VOL(i,isl,j)    &
                  - So_DENIT(i,isl,j)    - So_NSC(i,isl,j)
         ELSE 
            So_NSC(i,isl,j) = &
                    (So_LPONRES(i,isl,j) - So_LPONRES0(i,isl,j)) &
                  + (So_RPON(i,isl,j)    - So_RPON0(i,isl,j))    &
                  + (So_LDON(i,isl,j)    - So_LDON0(i,isl,j))    &
                  + (So_RDON(i,isl,j)    - So_RDON0(i,isl,j))    &
                  + (So_MBMN(i,isl,j)    - So_MBMN0(i,isl,j))    &
                  + (So_NH4(i,isl,j)     - So_NH40(i,isl,j))     &
                  + (So_NO3(i,isl,j)     - So_NO30(i,isl,j))

            So_NMBError(i,isl,j) = &
                    So_LPONRESin(i,isl,j)   + So_LDONperc0(i,isl-1,j) &
                  + So_RDONperc0(i,isl-1,j) + So_NH4perc0(i,isl-1,j)  &
                  + So_NO3perc0(i,isl-1,j)  - So_LDONintf(i,isl,j)    &
                  - So_RDONintf(i,isl,j)    - So_NH4intf(i,isl,j)     &
                  - So_NO3intf(i,isl,j)     - So_LDONperc(i,isl,j)    &
                  - So_RDONperc(i,isl,j)    - So_NH4perc(i,isl,j)     &
                  - So_NO3perc(i,isl,j)     - So_NH4uptk(i,isl,j)     &
                  - So_NO3uptk(i,isl,j)     - So_NH3VOL(i,isl,j)      &
                  - So_DENIT(i,isl,j)       - So_NSC(i,isl,j)
         END IF
      END DO

   end subroutine CheckNMB_So


   subroutine CheckNMB_Gw (gwid)

      implicit none
      
      ! --- Input arguments
      integer, intent(in) :: gwid

      Gw_NSC(gwid) = &
              (Gw_LDON(gwid) - Gw_LDON0(gwid)) &
            + (Gw_RDON(gwid) - Gw_RDON0(gwid)) &
            + (Gw_NH4(gwid)  - Gw_NH40(gwid))  &
            + (Gw_NO3(gwid)  - Gw_NO30(gwid))

      Gw_NMBError(gwid) = &
              Gw_LDONsogw0(gwid) + Gw_RDONsogw0(gwid) &  !BK20260316
            + Gw_NH4sogw0(gwid)  + Gw_NO3sogw0(gwid)  &  !BK20260316
            - Gw_LDONgwso(gwid)  - Gw_RDONgwso(gwid)  &  !BK20260320
            - Gw_NH4gwso(gwid)   - Gw_NO3gwso(gwid)   &  !BK20260320
            - Gw_LDONgwch(gwid)  - Gw_RDONgwch(gwid)  &
            - Gw_NH4gwch(gwid)   - Gw_NO3gwch(gwid)   &
            - Gw_NSC(gwid)

   end subroutine CheckNMB_Gw


   subroutine CheckNMB_Ch (ich)

      implicit none

      ! --- Input arguments
      integer, intent(in) :: ich

      ! --- Local variables
      integer :: chid, gwid
      real(8) :: dt

      dt = real(SedCNPmodel%SedCNP_timestep)
      chid = SedCNP_hydro%linkID(ich)
      gwid = domain%gwid_ch(chid)

      Ch_NSC(ich)  = &
              (Ch_ALGN(ich,1)     - Ch_ALGN0(ich,1))     &
            + (Ch_ALGN(ich,2)     - Ch_ALGN0(ich,2))     &
            + (Ch_ALGN(ich,3)     - Ch_ALGN0(ich,3))     &
            + (Ch_ZOON(ich)       - Ch_ZOON0(ich))       &
            + (Ch_MBMN(ich)       - Ch_MBMN0(ich))       &
            + (Ch_LPON(ich)       - Ch_LPON0(ich))       &
            + (Ch_RPON(ich)       - Ch_RPON0(ich))       &
            + (Ch_LDON(ich)       - Ch_LDON0(ich))       &
            + (Ch_RDON(ich)       - Ch_RDON0(ich))       &
            + (Ch_NH4(ich)        - Ch_NH40(ich))        &
            + (Ch_NO3(ich)        - Ch_NO30(ich))        &
            + (Ch_PONDEPOSIT(ich) - Ch_PONDEPOSIT0(ich))

      Ch_NMBError(ich) = &
              Ch_ALGNusch0(ich,1) + Ch_ALGNusch0(ich,2) &
            + Ch_ALGNusch0(ich,3) + Ch_ZOONusch0(ich)   &
            + Ch_MBMNusch0(ich)   + Ch_LPONusch0(ich)   &
            + Ch_RPONusch0(ich)   + Ch_LDONusch0(ich)   &
            + Ch_RDONusch0(ich)   + Ch_NH4usch0(ich)    &
            + Ch_NO3usch0(ich)    + Ch_LPONsurf0(ich)   &
            + Ch_RPONsurf0(ich)   + Ch_MBMNsurf0(ich)   &
            + Ch_LDONsurf0(ich)   + Ch_RDONsurf0(ich)   &
            + Ch_NH4surf0(ich)    + Ch_NO3surf0(ich)    &
            + Ch_LDONintf0(ich)   + Ch_RDONintf0(ich)   &
            + Ch_NH4intf0(ich)    + Ch_NO3intf0(ich)    &
            + Ch_LDONgwch0(ich)   + Ch_RDONgwch0(ich)   &
            + Ch_NH4gwch0(ich)    + Ch_NO3gwch0(ich)    &
            + Ch_LPONxpnt(ich)    + Ch_RPONxpnt(ich)    &
            + Ch_LDONxpnt(ich)    + Ch_RDONxpnt(ich)    &
            + Ch_NH4xpnt(ich)     + Ch_NO3xpnt(ich)     &
            + Ch_LPONxdis(ich)    + Ch_RPONxdis(ich)    &
            + Ch_LDONxdis(ich)    + Ch_RDONxdis(ich)    &
            + Ch_NH4xdis(ich)     + Ch_NO3xdis(ich)     &
            - Ch_LPONxabs(ich)    - Ch_RPONxabs(ich)    &
            - Ch_LDONxabs(ich)    - Ch_RDONxabs(ich)    &
            - Ch_NH4xabs(ich)     - Ch_NO3xabs(ich)     &
            - Ch_ALGNdsch(ich,1)  - Ch_ALGNdsch(ich,2)  &
            - Ch_ALGNdsch(ich,3)  - Ch_ZOONdsch(ich)    &
            - Ch_MBMNdsch(ich)    - Ch_LPONdsch(ich)    &
            - Ch_RPONdsch(ich)    - Ch_LDONdsch(ich)    &
            - Ch_RDONdsch(ich)    - Ch_NH4dsch(ich)     &
            - Ch_NO3dsch(ich)     - Ch_NH3VOL(ich)      &
            - Ch_DENIT(ich)       - Ch_NSC(ich)

   end subroutine CheckNMB_Ch

end module Nitrogen
