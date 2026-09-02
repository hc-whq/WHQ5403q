module Carbon

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

   subroutine SoilCarbonCycle(i, j, lcover, itime, TSOIL, THETA, THETA_S, &
        THETA_F, So_WStorage, Quptk, Qsurf, Qperc, Qintf)  

      implicit none

      ! --- Input arguments
      integer,               intent(in) :: i, j, lcover, itime
      real(8), dimension(4), intent(in) :: TSOIL, THETA, So_WStorage, Quptk
      real(8),               intent(in) :: THETA_S, THETA_F, Qsurf
      real(8), dimension(4), intent(in) :: Qperc
      real(8),               intent(in) :: Qintf

      ! --- Local Variables
      character(len=2)      :: smonth
      integer               :: imonth, isl, status
      real(8)               :: dt, rate0, KTBMBM, KTHETA, dPI, RK4SOL
      real(8), allocatable, dimension(:) :: Qintf_layer, Weight  !BK20260612 
      real(8)               :: TotalWeight
      real(8)               :: dx    ! grid cell size [m]  !BK20260301
      real(8)               :: q_out, w_strg  !(m3) !BK20260301  
      integer               :: k

      ! Decomposition and transport fluxes
      real(8) :: dLPOCLITD, dLPOCRESD, dLPOCEXCD, dLPOCMAND
      real(8) :: dLDOCD, dRPOCD, dRDOCD, dMBMCD
      real(8) :: dLPOCLITsurf, dLPOCRESsurf, dLPOCEXCsurf, dLPOCMANsurf

      ! Initial State and Flux Scaling Variables
      real(8) :: So_LPOCLIT_init, So_LPOCRES_init, So_LPOCEXC_init
      real(8) :: So_LPOCMAN_init, So_RPOC_init, So_LDOC_init
      real(8) :: So_RDOC_init, So_MBMC_init
      real(8) :: TotalFlux_LPOCLIT, TotalFlux_LPOCRES, TotalFlux_LPOCEXC
      real(8) :: TotalFlux_LPOCMAN, TotalFlux_RPOC, TotalFlux_LDOC
      real(8) :: TotalFlux_RDOC, TotalFlux_MBMC
      real(8) :: Scale_LPOCLIT, Scale_LPOCRES, Scale_LPOCEXC, Scale_LPOCMAN
      real(8) :: Scale_RPOC, Scale_LDOC, Scale_RDOC, Scale_MBMC

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

      ! =================================================================
      ! II. CALCULATE-SCALE-UPDATE LOGIC
      ! =================================================================

      DO isl = 1, domain%nsl

         ! -----------------------------------------------------------------
         ! A. Store initial values for simultaneous flux calculation
         ! (previous timestep state + current timestep inputs)
         ! -----------------------------------------------------------------
         
         if (isl == 1) then
            So_LPOCLITin(i,j) = So_LPOCLITinso(lcover,imonth) * dt
            So_LPOCEXCin(i,j) = So_LPOCEXCinso(lcover,imonth) * dt
            So_LPOCMANin(i,j) = So_LPOCMANinso(lcover,imonth) * dt
            So_LPOCRESin(i,isl,j) = So_LPOCRESinso(lcover,imonth) * dt * &
               FROOT(lcover, isl)

            So_LPOCLIT_init = So_LPOCLIT0(i,j) + So_LPOCLITin(i,j)
            So_LPOCEXC_init = So_LPOCEXC0(i,j) + So_LPOCEXCin(i,j)
            So_LPOCMAN_init = So_LPOCMAN0(i,j) + So_LPOCMANin(i,j)
            So_LPOCRES_init = So_LPOCRES0(i,isl,j) + So_LPOCRESin(i,isl,j)
            So_RPOC_init    = So_RPOC0(i,isl,j)
            So_LDOC_init    = So_LDOC0(i,isl,j) + So_LDOCgwso0(i,j)  !BK20260316
            So_RDOC_init    = So_RDOC0(i,isl,j) + So_RDOCgwso0(i,j)  !BK20260316
            So_MBMC_init    = So_MBMC0(i,isl,j)
         else
            So_LPOCRESin(i,isl,j) = So_LPOCRESinso(lcover,imonth) * dt * &
               FROOT(lcover, isl)

            So_LPOCLIT_init = 0.0  ! Does not exist in lower layers
            So_LPOCEXC_init = 0.0  ! Does not exist in lower layers
            So_LPOCMAN_init = 0.0  ! Does not exist in lower layers
            So_LPOCRES_init = So_LPOCRES0(i,isl,j) + So_LPOCRESin(i,isl,j)
            So_RPOC_init    = So_RPOC0(i,isl,j)
            So_LDOC_init    = So_LDOC0(i,isl,j) + So_LDOCperc0(i,isl-1,j)
            So_RDOC_init    = So_RDOC0(i,isl,j) + So_RDOCperc0(i,isl-1,j)
            So_MBMC_init    = So_MBMC0(i,isl,j)
         endif

         ! -----------------------------------------------------------------
         ! B. Calculate all transformation and transport fluxes
         ! -----------------------------------------------------------------

         KTBMBM = calc_KTB (TSOIL(isl), TMBM, MMBM)
         KTHETA = calc_KTHETA (THETA(isl), THETA_S, THETA_F)

         ! --- Organic C Decomposition ---
         if (isl == 1) then
            rate0 = KTBMBM * KTHETA * KLIT(lcover)
            RK4SOL = calc_RK4SOL (-rate0, So_LPOCLIT_init)
            dLPOCLITD = So_LPOCLIT_init - RK4SOL

            rate0 = KTBMBM * KTHETA * KEXC(lcover)
            RK4SOL = calc_RK4SOL (-rate0, So_LPOCEXC_init)
            dLPOCEXCD = So_LPOCEXC_init - RK4SOL

            rate0 = KTBMBM * KTHETA * KMAN(lcover)
            RK4SOL = calc_RK4SOL (-rate0, So_LPOCMAN_init)
            dLPOCMAND = So_LPOCMAN_init - RK4SOL
         endif

         rate0 = KTBMBM * KTHETA * KRES(lcover)
         RK4SOL = calc_RK4SOL (-rate0, So_LPOCRES_init)
         dLPOCRESD = So_LPOCRES_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KRPOM(lcover)
         RK4SOL = calc_RK4SOL (-rate0, So_RPOC_init)
         dRPOCD = So_RPOC_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KLDOM(lcover)
         RK4SOL = calc_RK4SOL (-rate0, So_LDOC_init)
         dLDOCD = So_LDOC_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KRDOM(lcover)
         RK4SOL = calc_RK4SOL (-rate0, So_RDOC_init)
         dRDOCD = So_RDOC_init - RK4SOL

         rate0 = KTBMBM * KTHETA * KMBM(lcover)
         RK4SOL = calc_RK4SOL (-rate0, So_MBMC_init)
         dMBMCD = So_MBMC_init - RK4SOL

         ! --- Transport by surface runoff
         if (isl == 1) then
            q_out = Qsurf / 1000. * dx * dx * dt ! (m3)
            w_strg = THETA(1) * Dmix * dx * dx   ! (m3)  !BK20260301
            !POM
            dPI = calc_dPI(etaPOMsurf, q_out, w_strg)  !BK20260301
            dLPOCLITsurf = So_LPOCLIT_init * (Dmix / DSOIL(isl)) * dPI
            dLPOCRESsurf = So_LPOCRES_init * (Dmix / DSOIL(isl)) * dPI
            dLPOCEXCsurf = So_LPOCEXC_init * (Dmix / DSOIL(isl)) * dPI
            dLPOCMANsurf = So_LPOCMAN_init * (Dmix / DSOIL(isl)) * dPI
            So_LPOCsurf(i,j) = dLPOCLITsurf + dLPOCRESsurf + &
                               dLPOCEXCsurf + dLPOCMANsurf
            So_RPOCsurf(i,j) = So_RPOC_init * (Dmix / DSOIL(isl)) * dPI
            So_MBMCsurf(i,j) = So_MBMC_init * (Dmix / DSOIL(isl)) * dPI
            !DOM
            dPI = calc_dPI(etaDOMsurf, q_out, w_strg)  !BK20260302
            So_LDOCsurf(i,j) = So_LDOC_init * (Dmix / DSOIL(isl)) * dPI
            So_RDOCsurf(i,j) = So_RDOC_init * (Dmix / DSOIL(isl)) * dPI
         endif

         ! --- Transport by Interflow ---
         q_out = Qintf_layer(isl) / 1000. * dx * dx * dt  ! (m3)
         w_strg = THETA(isl) * dx * dx    ! (m3)  !BK20260302
         !DOM
         dPI = calc_dPI(etaDOMintf, q_out, w_strg)  !BK20260302
         So_LDOCintf(i,isl,j) = So_LDOC_init * dPI
         So_RDOCintf(i,isl,j) = So_RDOC_init * dPI

         ! --- Transport by Percolation ---         
         q_out = Qperc(isl) / 1000. * dx * dx * dt  ! (m3)
         w_strg = THETA(isl) * dx * dx   ! (m3)  !BK20260302
         !DOM
         dPI = calc_dPI(etaDOMperc, q_out, w_strg)  !BK20260302
         So_LDOCperc(i,isl,j) = So_LDOC_init * dPI
         So_RDOCperc(i,isl,j) = So_RDOC_init * dPI

         ! -----------------------------------------------------------------
         ! C. Scale fluxes to ensure mass conservation
         ! -----------------------------------------------------------------

         if (isl == 1) then
            TotalFlux_LPOCLIT = dLPOCLITD + dLPOCLITsurf
            if (TotalFlux_LPOCLIT > So_LPOCLIT_init .and. &
                  TotalFlux_LPOCLIT > 0.0) then
               Scale_LPOCLIT = So_LPOCLIT_init / TotalFlux_LPOCLIT
               dLPOCLITD = dLPOCLITD * Scale_LPOCLIT
               dLPOCLITsurf = dLPOCLITsurf * Scale_LPOCLIT
            endif
            
            TotalFlux_LPOCEXC = dLPOCEXCD + dLPOCEXCsurf
            if (TotalFlux_LPOCEXC > So_LPOCEXC_init .and. &
                  TotalFlux_LPOCEXC > 0.0) then
               Scale_LPOCEXC = So_LPOCEXC_init / TotalFlux_LPOCEXC
               dLPOCEXCD = dLPOCEXCD * Scale_LPOCEXC
               dLPOCEXCsurf = dLPOCEXCsurf * Scale_LPOCEXC
            endif
            
            TotalFlux_LPOCMAN = dLPOCMAND + dLPOCMANsurf
            if (TotalFlux_LPOCMAN > So_LPOCMAN_init .and. &
                  TotalFlux_LPOCMAN > 0.0) then
               Scale_LPOCMAN = So_LPOCMAN_init / TotalFlux_LPOCMAN
               dLPOCMAND = dLPOCMAND * Scale_LPOCMAN
               dLPOCMANsurf = dLPOCMANsurf * Scale_LPOCMAN
            endif
         endif
         
         if (isl == 1) then
            TotalFlux_LPOCRES = dLPOCRESD + dLPOCRESsurf
         else
            TotalFlux_LPOCRES = dLPOCRESD
         endif
         if (TotalFlux_LPOCRES > So_LPOCRES_init .and. &
               TotalFlux_LPOCRES > 0.0) then
            Scale_LPOCRES = So_LPOCRES_init / TotalFlux_LPOCRES
            dLPOCRESD = dLPOCRESD * Scale_LPOCRES
            if(isl == 1) dLPOCRESsurf = dLPOCRESsurf * Scale_LPOCRES
         endif
         
         if (isl == 1) then
             TotalFlux_RPOC = dRPOCD + So_RPOCsurf(i,j)
         else
             TotalFlux_RPOC = dRPOCD
         endif
         if (TotalFlux_RPOC > So_RPOC_init .and. TotalFlux_RPOC > 0.0) then
            Scale_RPOC = So_RPOC_init / TotalFlux_RPOC
            dRPOCD = dRPOCD * Scale_RPOC
            if(isl == 1) So_RPOCsurf(i,j) = So_RPOCsurf(i,j) * Scale_RPOC
         endif

         if (isl == 1) then
            TotalFlux_LDOC = dLDOCD + So_LDOCsurf(i,j) + &
                  So_LDOCintf(i,isl,j) + So_LDOCperc(i,isl,j)
         else
            TotalFlux_LDOC = dLDOCD + So_LDOCintf(i,isl,j) + &
                  So_LDOCperc(i,isl,j)
         endif
         if(TotalFlux_LDOC > So_LDOC_init .and. TotalFlux_LDOC > 0.0) then
            Scale_LDOC = So_LDOC_init/TotalFlux_LDOC
            dLDOCD = dLDOCD * Scale_LDOC
            if(isl == 1) So_LDOCsurf(i,j) = So_LDOCsurf(i,j) * Scale_LDOC
            So_LDOCintf(i,isl,j) = So_LDOCintf(i,isl,j) * Scale_LDOC
            So_LDOCperc(i,isl,j) = So_LDOCperc(i,isl,j) * Scale_LDOC
         endif
         
         if (isl == 1) then
            TotalFlux_RDOC = dRDOCD + So_RDOCsurf(i,j) + &
                  So_RDOCintf(i,isl,j) + So_RDOCperc(i,isl,j)
         else
            TotalFlux_RDOC = dRDOCD + So_RDOCintf(i,isl,j) + &
                  So_RDOCperc(i,isl,j)
         endif
         if(TotalFlux_RDOC > So_RDOC_init .and. TotalFlux_RDOC > 0.0) then
            Scale_RDOC = So_RDOC_init/TotalFlux_RDOC
            dRDOCD = dRDOCD * Scale_RDOC
            if(isl == 1) So_RDOCsurf(i,j) = So_RDOCsurf(i,j) * Scale_RDOC
            So_RDOCintf(i,isl,j) = So_RDOCintf(i,isl,j) * Scale_RDOC
            So_RDOCperc(i,isl,j) = So_RDOCperc(i,isl,j) * Scale_RDOC
         endif

         if (isl == 1) then
            TotalFlux_MBMC = dMBMCD + So_MBMCsurf(i,j)
         else
            TotalFlux_MBMC = dMBMCD
         endif
         if (TotalFlux_MBMC > So_MBMC_init .and. TotalFlux_MBMC > 0.0) then
            Scale_MBMC = So_MBMC_init / TotalFlux_MBMC
            dMBMCD = dMBMCD * Scale_MBMC
            if(isl == 1) So_MBMCsurf(i,j) = So_MBMCsurf(i,j) * Scale_MBMC
         endif

         ! Recalculate So_LPOCsurf with potentially scaled component fluxes
         if (isl == 1) then
            So_LPOCsurf(i,j) = dLPOCLITsurf + dLPOCRESsurf + &
                                 dLPOCEXCsurf + dLPOCMANsurf
         endif

         ! -----------------------------------------------------------------
         ! D. Update final storages using scaled fluxes
         ! -----------------------------------------------------------------

         if (isl == 1) then
            So_LPOCLIT(i,j) = So_LPOCLIT_init - dLPOCLITD - dLPOCLITsurf
            So_LPOCEXC(i,j) = So_LPOCEXC_init - dLPOCEXCD - dLPOCEXCsurf
            So_LPOCMAN(i,j) = So_LPOCMAN_init - dLPOCMAND - dLPOCMANsurf
            So_LPOCRES(i,isl,j) = So_LPOCRES_init - dLPOCRESD - dLPOCRESsurf
            So_RPOC(i,isl,j) = So_RPOC_init + FDEC * FREF * (dLPOCLITD + &
               dLPOCRESD + dLPOCEXCD + dLPOCMAND) - dRPOCD - So_RPOCsurf(i,j)
            So_LDOC(i,isl,j) = So_LDOC_init + FDEC * (1. - FREF) * &
               (dLPOCLITD + dLPOCRESD + dLPOCEXCD + dLPOCMAND) + &
               FDEC * dMBMCD - dLDOCD - So_LDOCsurf(i,j) - &
               So_LDOCintf(i,isl,j) - So_LDOCperc(i,isl,j)
            So_RDOC(i,isl,j) = So_RDOC_init + FDEC * (dLDOCD + dRPOCD + &
               dRDOCD) - dRDOCD - So_RDOCsurf(i,j) - So_RDOCintf(i,isl,j) - &
               So_RDOCperc(i,isl,j)
            So_MBMC(i,isl,j) = So_MBMC_init + FBIO * (dLPOCLITD + dLPOCRESD + &
               dLPOCEXCD + dLPOCMAND + dRPOCD + dLDOCD + dRDOCD + dMBMCD) - &
               dMBMCD - So_MBMCsurf(i,j)
            So_CO2(i,isl,j) = FMET * (dLPOCLITD + dLPOCRESD + dLPOCEXCD + &
               dLPOCMAND + dRPOCD + dLDOCD + dRDOCD + dMBMCD)
         else
            So_LPOCRES(i,isl,j) = So_LPOCRES_init - dLPOCRESD
            So_RPOC(i,isl,j) = So_RPOC_init + FDEC * FREF * dLPOCRESD - dRPOCD
            So_LDOC(i,isl,j) = So_LDOC_init + FDEC * (1. - FREF) * dLPOCRESD + &
               FDEC * dMBMCD - dLDOCD - So_LDOCintf(i,isl,j) - &
               So_LDOCperc(i,isl,j)
            So_RDOC(i,isl,j) = So_RDOC_init + FDEC * (dLDOCD + dRPOCD + &
               dRDOCD) - dRDOCD - So_RDOCintf(i,isl,j) - So_RDOCperc(i,isl,j)
            So_MBMC(i,isl,j) = So_MBMC_init + FBIO * (dLPOCRESD + dRPOCD + &
               dLDOCD + dRDOCD + dMBMCD) - dMBMCD
            So_CO2(i,isl,j) = FMET * (dLPOCRESD + dRPOCD + dLDOCD + dRDOCD + &
               dMBMCD)
         endif
      END DO

      ! Enforce non-negative soil storages
      if (So_LPOCLIT(i,j) < 1.0e-20) So_LPOCLIT(i,j) = 0.0
      if (So_LPOCEXC(i,j) < 1.0e-20) So_LPOCEXC(i,j) = 0.0
      if (So_LPOCMAN(i,j) < 1.0e-20) So_LPOCMAN(i,j) = 0.0
      DO isl = 1, domain%nsl
        if (So_LPOCRES(i,isl,j) < 1.0e-20) So_LPOCRES(i,isl,j) = 0.0
        if (So_RPOC(i,isl,j) < 1.0e-20) So_RPOC(i,isl,j) = 0.0
        if (So_LDOC(i,isl,j) < 1.0e-20) So_LDOC(i,isl,j) = 0.0
        if (So_RDOC(i,isl,j) < 1.0e-20) So_RDOC(i,isl,j) = 0.0
        if (So_MBMC(i,isl,j) < 1.0e-20) So_MBMC(i,isl,j) = 0.0
      END DO

      ! ====================================================================
      ! III. UPDATE OUTPUTS TO OTHER COMPONENTS
      ! ====================================================================

      So_LDOCsogw(i,j) = So_LDOCperc(i,domain%nsl,j)  !BK20260318
      So_RDOCsogw(i,j) = So_RDOCperc(i,domain%nsl,j)  !BK20260318

   end subroutine SoilCarbonCycle


   ! subroutine AquiferCarbonCycle(i, j, Qaqso, Aq_WStorage, Qaqgw)

   !    implicit none

   !    ! --- Input arguments
   !    integer, intent(in) :: i, j
   !    real(8), intent(in) :: Qaqso
   !    real(8), intent(in) :: Aq_WStorage
   !    real(8), intent(in) :: Qaqgw

   !    ! --- Local variables
   !    real(8) :: q_out, dx    !BK20260302      
   !    real(8) :: dPI, dt

   !    ! --- Variables for calculate-scale-update logic
   !    real(8) :: Aq_LDOC_init, Aq_RDOC_init
   !    real(8) :: TotalFlux_LDOC, TotalFlux_RDOC
   !    real(8) :: Scale_LDOC, Scale_RDOC

   !    dt = real(SedCNPmodel%SedCNP_timestep)
   !    dx = SedCNP_hydro%dx(1)    !BK20260301

   !    ! ====================================================================
   !    ! I. INITIALIZE STATE AND ADD INFLOWS
   !    ! ====================================================================

   !    Aq_LDOC_init = Aq_LDOC0(i,j) + So_LDOCsoaq0(i,j)
   !    Aq_RDOC_init = Aq_RDOC0(i,j) + So_RDOCsoaq0(i,j)

   !    ! ====================================================================
   !    ! II. CALCULATE TRANSPORT FLUXES
   !    ! ====================================================================

   !    ! --- Abstraction for irrigation
   !    if (Aq_WStorage > (Qaqso * domain%areaxy / 1000.0 * dt)) then
   !       Aq_LDOCaqso(i,j) = (Qaqso * domain%areaxy / 1000.0 * dt) * &
   !                          Aq_LDOC_init / Aq_WStorage
   !       Aq_RDOCaqso(i,j) = (Qaqso * domain%areaxy / 1000.0 * dt) * &
   !                          Aq_RDOC_init / Aq_WStorage
   !    else
   !       Aq_LDOCaqso(i,j) = Aq_LDOC_init
   !       Aq_RDOCaqso(i,j) = Aq_RDOC_init
   !    endif

   !    ! --- Aquifer discharge to groundwater
   !    q_out = Qaqgw / 1000. * dx * dx * dt  ! (m3)    !BK20260302
   !    dPI = calc_dPI(etaDOMaqgw, q_out, Aq_WStorage)  !BK20260302
   !    Aq_LDOCaqgw(i,j) = Aq_LDOC_init * dPI
   !    Aq_RDOCaqgw(i,j) = Aq_RDOC_init * dPI

   !    ! ====================================================================
   !    ! III. SCALE FLUXES TO ENSURE MASS CONSERVATION
   !    ! ====================================================================

   !    TotalFlux_LDOC = Aq_LDOCaqso(i,j) + Aq_LDOCaqgw(i,j)
   !    if (TotalFlux_LDOC > Aq_LDOC_init .and. TotalFlux_LDOC > 0.0) then
   !       Scale_LDOC = Aq_LDOC_init / TotalFlux_LDOC
   !       Aq_LDOCaqso(i,j) = Aq_LDOCaqso(i,j) * Scale_LDOC
   !       Aq_LDOCaqgw(i,j) = Aq_LDOCaqgw(i,j) * Scale_LDOC
   !    endif

   !    TotalFlux_RDOC = Aq_RDOCaqso(i,j) + Aq_RDOCaqgw(i,j)
   !    if (TotalFlux_RDOC > Aq_RDOC_init .and. TotalFlux_RDOC > 0.0) then
   !       Scale_RDOC = Aq_RDOC_init / TotalFlux_RDOC
   !       Aq_RDOCaqso(i,j) = Aq_RDOCaqso(i,j) * Scale_RDOC
   !       Aq_RDOCaqgw(i,j) = Aq_RDOCaqgw(i,j) * Scale_RDOC
   !    endif

   !    ! ====================================================================
   !    ! IV. UPDATE FINAL STORAGES
   !    ! ====================================================================

   !    Aq_LDOC(i,j) = Aq_LDOC_init - Aq_LDOCaqso(i,j) - Aq_LDOCaqgw(i,j)
   !    Aq_RDOC(i,j) = Aq_RDOC_init - Aq_RDOCaqso(i,j) - Aq_RDOCaqgw(i,j)

   !    ! Enforce non-negative aquifer storages
   !    if (Aq_LDOC(i,j) < 1.0e-20) Aq_LDOC(i,j) = 0.0
   !    if (Aq_RDOC(i,j) < 1.0e-20) Aq_RDOC(i,j) = 0.0

   ! end subroutine AquiferCarbonCycle


   subroutine GroundwaterCarbonCycle(gwid, Gw_WStorage, Qgwch)  !BK20260316

      implicit none

      ! --- Input arguments
      integer, intent(in)  :: gwid
      real(8), intent(in)  :: Gw_WStorage    ![m3]       !BK20260316
      real(8), intent(in)  :: Qgwch          ![mm s-1]   !BK20260316

      ! --- Local variables
      integer :: i, j, istat
      real(8) :: q_out                       ![m3]       !BK20260301  
      real(8) :: dPI, dx, dt                             !BK20260302     
      real(8), allocatable :: Qgwso(:,:)     ![mm s-1]   !BK20260316
      real(8), allocatable :: Qgwso_bas(:)   ![m3]       !BK20260316

      ! --- Variables for calculate-scale-update logic
      real(8) :: Gw_LDOC_init, Gw_RDOC_init
      real(8) :: TotalFlux_LDOC, TotalFlux_RDOC
      real(8) :: Scale_LDOC, Scale_RDOC, Scale_gwso

      dt = real(SedCNPmodel%SedCNP_timestep)
      dx = SedCNP_hydro%dx(1)    !BK20260301

      allocate(Qgwso(domain%ix,domain%jx),stat=istat)  !BK20260316
      allocate(Qgwso_bas(domain%nbasin),stat=istat)    !BK20260316

      ! ====================================================================
      ! I. AGGREGATE INFLOWS FROM SOIL
      ! ====================================================================

      Gw_LDOCsogw0(gwid) = 0.0  !BK20260316
      Gw_RDOCsogw0(gwid) = 0.0  !BK20260316

      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               Gw_LDOCsogw0(gwid) = Gw_LDOCsogw0(gwid) + So_LDOCsogw0(i,j)  !BK20260316
               Gw_RDOCsogw0(gwid) = Gw_RDOCsogw0(gwid) + So_RDOCsogw0(i,j)  !BK20260316
            END IF
         END DO
      END DO

      ! ====================================================================
      ! II. INITIALIZE STATE AND ADD INFLOWS
      ! ====================================================================

      Gw_LDOC_init = Gw_LDOC0(gwid) + Gw_LDOCsogw0(gwid) 
      Gw_RDOC_init = Gw_RDOC0(gwid) + Gw_RDOCsogw0(gwid) 

      ! ====================================================================
      ! III. CALCULATE TRANSPORT FLUXES 
      ! ====================================================================

      ! --- groundwater abstraction for irrigation !BK20260316
      Qgwso_bas(gwid) = 0.0  
      Gw_LDOCgwso(gwid) = 0.0  
      Gw_RDOCgwso(gwid) = 0.0

      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               ! Clamp irridepth to >= 0 to exclude fill values (e.g. -9999)
               ! Qgwso(i,j) = SedCNP_hydro%irridepth(i,j) / dt        ! [mm s-1]
               Qgwso(i,j) = max(0.0d0, SedCNP_hydro%irridepth(i,j)) / dt  ! [mm s-1] !BK20260508
               Qgwso_bas(gwid) = Qgwso_bas(gwid) + Qgwso(i,j) / 1000.0 &
                     * domain%areaxy * dt  ! [m3]
            END IF
         END DO
      END DO

      ! if (Gw_WStorage > Qgwso_bas(gwid)) then
      if (Gw_WStorage > 0.0d0 .and. Gw_WStorage > Qgwso_bas(gwid)) then !BK20260508
         Gw_LDOCgwso(gwid) = Gw_LDOC_init * Qgwso_bas(gwid) / Gw_WStorage
         Gw_RDOCgwso(gwid) = Gw_RDOC_init * Qgwso_bas(gwid) / Gw_WStorage
      else if (Gw_WStorage > 0.0d0) then  !---BK20260508
         Gw_LDOCgwso(gwid) = Gw_LDOC_init 
         Gw_RDOCgwso(gwid) = Gw_RDOC_init !---BK20260508
      else
         !Gw_LDOCgwso(gwid) = Gw_LDOC_init 
         !Gw_RDOCgwso(gwid) = Gw_RDOC_init
         Gw_LDOCgwso(gwid) = 0.0d0  !BK20260508
         Gw_RDOCgwso(gwid) = 0.0d0  !BK20260508
      endif

      ! --- groundwater discharge to channel
      q_out = Qgwch / 1000. * domain%areaxy * dt  ! [m3]    !BK20260302
      !DOM
      dPI = calc_dPI(etaDOMgwch, q_out, Gw_WStorage)        !BK20260302
      Gw_LDOCgwch(gwid) = Gw_LDOC_init * dPI
      Gw_RDOCgwch(gwid) = Gw_RDOC_init * dPI

      ! ====================================================================
      ! IV. SCALE FLUXES TO ENSURE MASS CONSERVATION
      ! ====================================================================
 
      ! --- (1) scale flux *gwso(gwid) and *gwch(gwid) at the subbasin-level
      TotalFlux_LDOC = Gw_LDOCgwso(gwid) + Gw_LDOCgwch(gwid)
      if (TotalFlux_LDOC > Gw_LDOC_init .and. TotalFlux_LDOC > 0.0) then
         Scale_LDOC = Gw_LDOC_init / TotalFlux_LDOC
         Gw_LDOCgwso(gwid) = Gw_LDOCgwso(gwid) * Scale_LDOC
         Gw_LDOCgwch(gwid) = Gw_LDOCgwch(gwid) * Scale_LDOC
      endif

      TotalFlux_RDOC = Gw_RDOCgwso(gwid) + Gw_RDOCgwch(gwid)
      if (TotalFlux_RDOC > Gw_RDOC_init .and. TotalFlux_RDOC > 0.0) then
         Scale_RDOC = Gw_RDOC_init / TotalFlux_RDOC
         Gw_RDOCgwso(gwid) = Gw_RDOCgwso(gwid) * Scale_RDOC
         Gw_RDOCgwch(gwid) = Gw_RDOCgwch(gwid) * Scale_RDOC
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
               So_LDOCgwso(i,j) = Gw_LDOCgwso(gwid) * Scale_gwso
               So_RDOCgwso(i,j) = Gw_RDOCgwso(gwid) * Scale_gwso                    
            END IF
         END DO
      END DO

      ! ====================================================================
      ! V. UPDATE FINAL STORAGES
      ! ====================================================================

      Gw_LDOC(gwid) = Gw_LDOC_init - Gw_LDOCgwso(gwid) - Gw_LDOCgwch(gwid)
      Gw_RDOC(gwid) = Gw_RDOC_init - Gw_RDOCgwso(gwid) - Gw_RDOCgwch(gwid)

      ! Enforce non-negative groundwater storages
      if (Gw_LDOC(gwid) < 1.0e-20) Gw_LDOC(gwid) = 0.0
      if (Gw_RDOC(gwid) < 1.0e-20) Gw_RDOC(gwid) = 0.0

   end subroutine GroundwaterCarbonCycle


   subroutine ChannelCarbonCycle(ich, itime, Qdsch, WtopWdth, WStorage)

      implicit none

      ! --- Input arguments
      integer, intent(in) :: ich, itime
      real(8), intent(in) :: Qdsch, WtopWdth, WStorage

      ! --- Local variables
      integer :: i, j, isl, ips, ia, im, month, stat, chid, gwid
      real(8) :: rate0, dW, dZ, dXi1, dXi2, dLambda_L, dLambda_N, dLambda_P
      real(8) :: dPHI_N, dPHI_P, dLPOCD, dLPOCZ, dLPOCS, dLDOCD, dRPOCD
      real(8) :: dRPOCS, dRDOCD, dMBMCD, dAZP, dALGCconc, dTSSconc
      real(8) :: dTA(ntm), dMA(ntm), dt, abs_ratio
      real(8) :: SatCO2, KCO2, pKH, PCO2, KO2, U, H, KCO2_20, KCO2_T
      real(8) :: RK4SOL, KTBMBM, KTBZOO, KTBALG, KTHETA, ALPHA, dPI
      real(8) :: KZG(domain%nch), KZR(domain%nch), KZM(domain%nch)
      real(8) :: KAG(domain%nch,nalg), KAR(domain%nch,nalg)
      real(8) :: KAE(domain%nch,nalg), KAM(domain%nch,nalg)
      real(8), parameter :: e = 2.718281828459
      character(len=2) :: monthstr

      ! --- Variables for calculate-scale-update logic
      real(8) :: Ch_ALGC_init(nalg), Ch_ZOOC_init, Ch_LPOC_init, Ch_LDOC_init
      real(8) :: Ch_RPOC_init, Ch_RDOC_init, Ch_MBMC_init, Ch_DIC_init
      real(8) :: TotalFlux_ALGC(nalg), Scale_ALGC(nalg)
      real(8) :: TotalFlux_ZOOC, Scale_ZOOC
      real(8) :: TotalFlux_LPOC, Scale_LPOC
      real(8) :: TotalFlux_RPOC, Scale_RPOC
      real(8) :: TotalFlux_LDOC, Scale_LDOC
      real(8) :: TotalFlux_RDOC, Scale_RDOC
      real(8) :: TotalFlux_MBMC, Scale_MBMC
      real(8) :: TotalFlux_DIC, Scale_DIC

      dt = real(SedCNPmodel%SedCNP_timestep)
      chid = SedCNP_hydro%linkID(ich)
      gwid = domain%gwid_ch(chid)

      !================================================================
      ! I. INFLOWS
      !================================================================

      ! --- Inflow from upstream channels
      Ch_ALGC(ich,1:nalg) = Ch_ALGC0(ich,1:nalg) + Ch_ALGCusch0(ich,1:nalg)
      Ch_ZOOC(ich) = Ch_ZOOC0(ich) + Ch_ZOOCusch0(ich)
      Ch_MBMC(ich) = Ch_MBMC0(ich) + Ch_MBMCusch0(ich)
      Ch_DIC(ich)  = Ch_DIC0(ich)  + Ch_DICusch0(ich)
      Ch_LPOC(ich) = Ch_LPOC0(ich) + Ch_LPOCusch0(ich)
      Ch_RPOC(ich) = Ch_RPOC0(ich) + Ch_RPOCusch0(ich)
      Ch_LDOC(ich) = Ch_LDOC0(ich) + Ch_LDOCusch0(ich)
      Ch_RDOC(ich) = Ch_RDOC0(ich) + Ch_RDOCusch0(ich)

      ! --- Lateral inflows (surface, interflow, groundwater)
      IF (gwid >= 0) THEN  ! if 'genuine' channel to which gwbasin drains

         ! initialisation
         Ch_LPOCsurf0(ich) = 0.0
         Ch_RPOCsurf0(ich) = 0.0
         Ch_LDOCsurf0(ich) = 0.0
         Ch_RDOCsurf0(ich) = 0.0
         Ch_MBMCsurf0(ich) = 0.0
         Ch_LDOCintf0(ich) = 0.0
         Ch_RDOCintf0(ich) = 0.0

         DO i = 1, domain%ix
            DO j = 1, domain%jx
               IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN

                  ! surface runoff
                  Ch_LPOCsurf0(ich) = Ch_LPOCsurf0(ich) + So_LPOCsurf0(i,j)
                  Ch_RPOCsurf0(ich) = Ch_RPOCsurf0(ich) + So_RPOCsurf0(i,j)
                  Ch_LDOCsurf0(ich) = Ch_LDOCsurf0(ich) + So_LDOCsurf0(i,j)
                  Ch_RDOCsurf0(ich) = Ch_RDOCsurf0(ich) + So_RDOCsurf0(i,j)
                  Ch_MBMCsurf0(ich) = Ch_MBMCsurf0(ich) + So_MBMCsurf0(i,j)
                  
                  ! interflow
                  DO isl = 1, domain%nsl
                     Ch_LDOCintf0(ich) = Ch_LDOCintf0(ich) + So_LDOCintf0(i,isl,j)
                     Ch_RDOCintf0(ich) = Ch_RDOCintf0(ich) + So_RDOCintf0(i,isl,j)
                  END DO

               END IF
            END DO
         END DO

         ! groundwater inflow
         Ch_LDOCgwch0(ich) = Gw_LDOCgwch0(gwid)
         Ch_RDOCgwch0(ich) = Gw_RDOCgwch0(gwid)

         ! update storage
         Ch_LPOC(ich) = Ch_LPOC(ich) + Ch_LPOCsurf0(ich)
         Ch_RPOC(ich) = Ch_RPOC(ich) + Ch_RPOCsurf0(ich)
         Ch_LDOC(ich) = Ch_LDOC(ich) + Ch_LDOCsurf0(ich) + Ch_LDOCintf0(ich) + &
            Ch_LDOCgwch0(ich)
         Ch_RDOC(ich) = Ch_RDOC(ich) + Ch_RDOCsurf0(ich) + Ch_RDOCintf0(ich) + &
            Ch_RDOCgwch0(ich)
         Ch_MBMC(ich) = Ch_MBMC(ich) + Ch_MBMCsurf0(ich)

      ENDIF

      ! --- Forced inflow/ourflow 
      ! inflow: point source
      Ch_LPOCxpnt(ich) = Ch_LPOCpnt(ich,itime) * dt
      Ch_RPOCxpnt(ich) = Ch_RPOCpnt(ich,itime) * dt
      Ch_LDOCxpnt(ich) = Ch_LDOCpnt(ich,itime) * dt
      Ch_RDOCxpnt(ich) = Ch_RDOCpnt(ich,itime) * dt

      ! inflow: water discharge
      Ch_LPOCxdis(ich) = Ch_LPOCdis(ich,itime) * dt  
      Ch_RPOCxdis(ich) = Ch_RPOCdis(ich,itime) * dt  
      Ch_LDOCxdis(ich) = Ch_LDOCdis(ich,itime) * dt  
      Ch_RDOCxdis(ich) = Ch_RDOCdis(ich,itime) * dt  

      ! update storage (inflows)
      Ch_LPOC(ich) = Ch_LPOC(ich) + Ch_LPOCxpnt(ich) + Ch_LPOCxdis(ich)
      Ch_RPOC(ich) = Ch_RPOC(ich) + Ch_RPOCxpnt(ich) + Ch_RPOCxdis(ich)
      Ch_LDOC(ich) = Ch_LDOC(ich) + Ch_LDOCxpnt(ich) + Ch_LDOCxdis(ich)
      Ch_RDOC(ich) = Ch_RDOC(ich) + Ch_RDOCxpnt(ich) + Ch_RDOCxdis(ich)

      ! outflow: water abstraction
      abs_ratio = SedCNP_hydro%Qabs(ich,itime) / (WStorage + Qdsch * dt)
      if(SedCNP_hydro%Qabs(ich,itime) .gt. 0.) then
         Ch_LPOCabs(ich,itime) = Ch_LPOC(ich) * abs_ratio
         Ch_RPOCabs(ich,itime) = Ch_RPOC(ich) * abs_ratio
         Ch_LDOCabs(ich,itime) = Ch_LDOC(ich) * abs_ratio
         Ch_RDOCabs(ich,itime) = Ch_RDOC(ich) * abs_ratio
      endif
      Ch_LPOCxabs(ich) = Ch_LPOCabs(ich,itime) * dt  
      Ch_RPOCxabs(ich) = Ch_RPOCabs(ich,itime) * dt  
      Ch_LDOCxabs(ich) = Ch_LDOCabs(ich,itime) * dt  
      Ch_RDOCxabs(ich) = Ch_RDOCabs(ich,itime) * dt  

      ! update storage (outflows)
      Ch_LPOC(ich) = Ch_LPOC(ich) - Ch_LPOCxabs(ich)
      Ch_RPOC(ich) = Ch_RPOC(ich) - Ch_RPOCxabs(ich)
      Ch_LDOC(ich) = Ch_LDOC(ich) - Ch_LDOCxabs(ich)
      Ch_RDOC(ich) = Ch_RDOC(ich) - Ch_RDOCxabs(ich)

      ! Enforce non-negative channel storages after abstraction
      if (Ch_LPOC(ich) < 1.0e-20) Ch_LPOC(ich) = 0.0
      if (Ch_RPOC(ich) < 1.0e-20) Ch_RPOC(ich) = 0.0
      if (Ch_LDOC(ich) < 1.0e-20) Ch_LDOC(ich) = 0.0
      if (Ch_RDOC(ich) < 1.0e-20) Ch_RDOC(ich) = 0.0

      !================================================================
      ! II. CALCULATE-SCALE-UPDATE
      !================================================================

      ! --- Store initial values for simultaneous flux calculation
      Ch_ALGC_init(:) = Ch_ALGC(ich,:)
      Ch_ZOOC_init    = Ch_ZOOC(ich)
      Ch_LPOC_init    = Ch_LPOC(ich)
      Ch_LDOC_init    = Ch_LDOC(ich)
      Ch_RPOC_init    = Ch_RPOC(ich)
      Ch_RDOC_init    = Ch_RDOC(ich)
      Ch_MBMC_init    = Ch_MBMC(ich)
      Ch_DIC_init     = Ch_DIC(ich)

      ! === A. CALCULATE ALL FLUXES BASED ON INITIAL STORAGES
      
      ! --- Pre-calculations for transformation rates
      KTHETA = 0.6 ! Assumed constant for in-stream saturated conditions
      dW = DWATER(ich)
      dZ = min(dW, 10.0)  ! take account of algaes in the top 10 m of water

      ! ALPHA - light attenuation coeff. (Chesapeake Bay Program, 2000)
      if (WStorage > 0.0) then
         dALGCconc = sum(Ch_ALGC_init(:)) / WStorage
         dTSSconc = (sum(channelSed%Cchsd(1:nps,ich)) + Ch_LPOC_init + &
                     Ch_RPOC_init) / WStorage
      else
         dALGCconc = 0.0
         dTSSconc = 0.0
      endif
      ALPHA = ALPHADOC + ALPHACHL * dALGCconc + ALPHATSS * dTSSconc

      ! limiting factors for light, nitrogen, and phosphorus
      monthstr = trim(dateSedCNP%olddate(6:7))
      read(monthstr,*) month

      if (gwid > 0) then
         dXi1 = (1.0 - ALBEDO(month)/100.0) * &
            SedCNP_hydro%SWDOWNavg(gwid) * exp(-ALPHA * dW) / SOLRADMAX
         dXi2 = (1.0 - ALBEDO(month)/100.0) * &
            SedCNP_hydro%SWDOWNavg(gwid) * exp(-ALPHA * (dW + dZ)) / SOLRADMAX
      else
         dXi1 = 0.0
         dXi2 = 0.0
      end if

      if ((ALPHA*dZ) > 1e-6) then
         dLambda_L = e / (ALPHA*dZ) * (exp(-dXi2) - exp(-dXi1))
      else
         dLambda_L = 0.0
      end if

      if (WStorage > 0.0) then
         dPHI_N = (Ch_NH40(ich) + Ch_NO30(ich)) / WStorage
         dPHI_P = Ch_PO40(ich) / WStorage
      else
         dPHI_N = 0.0
         dPHI_P = 0.0
      end if

      if ((MONOD_N + dPHI_N) > 0.) then  !Y.Kwon
         dLambda_N = dPHI_N / (MONOD_N + dPHI_N)  !Y.Kwon
      else
         dLambda_N = 0.0
      endif

      if ((MONOD_P + dPHI_P) > 0.0) then  !Y.Kwon
         dLambda_P = dPHI_P / (MONOD_P + dPHI_P)  !Y.Kwon
      else
         dLambda_P = 0.0
      endif

      ! --- Zooplankton
      KTBZOO = calc_KTB(TWATER(ich), TZOO, MZOO)
      dAZP = sum(Ch_ALGC_init(:)) + Ch_LPOC_init + Ch_ZOOC_init

      ! growth 
      if (dAZP + ZHALF * WStorage > 1.0E-9) then
         KZG(ich) = KTBZOO * EZI * KZIMAX * (dAZP - ZLOW * WStorage) / &
                     (dAZP + ZHALF * WStorage)
      else
         KZG(ich) = 0.0
      end if
      RK4SOL = calc_RK4SOL(KZG(ich), Ch_ZOOC_init)
      Ch_ZOOCG(ich) = max(0.0, RK4SOL - Ch_ZOOC_init)
      
      ! being grazed by (other) zooplankton
      if (EZI > 0.0 .and. dAZP > 0.0) then
         Ch_ZOOCZ(ich) = (Ch_ZOOCG(ich) / EZI) * (Ch_ZOOC_init / dAZP)
      else
         Ch_ZOOCZ(ich) = 0.0
      end if

      ! respiration 
      KZR(ich) = KTBZOO * KZRMAX
      RK4SOL = calc_RK4SOL(-KZR(ich), Ch_ZOOC_init)
      Ch_ZOOCR(ich) = Ch_ZOOC_init - RK4SOL

      ! excretion - based on growth to ensure mass balance (BK20251127)
      if (EZI > 0.0) then
         Ch_ZOOCE(ich) = (Ch_ZOOCG(ich) / EZI) * (1.0 - EZI) 
      else
         Ch_ZOOCE(ich) = 0.0
      end if

      ! mortality 
      KZM(ich) = KTBZOO * KZMMAX
      RK4SOL = calc_RK4SOL(-KZM(ich), Ch_ZOOC_init)
      Ch_ZOOCM(ich) = Ch_ZOOC_init - RK4SOL

      ! settling
      if (dW > 0.0) then
         Ch_ZOOCS(ich) = (1.0 - exp(-OMEGAZOO * dt / dW)) * Ch_ZOOC_init
      else
         Ch_ZOOCS(ich) = 0.0
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
         RK4SOL = calc_RK4SOL(KAG(ich,ia), Ch_ALGC_init(ia))
         Ch_ALGCG(ich,ia) = max(0.0, RK4SOL - Ch_ALGC_init(ia))

         ! grazing by zooplankton
         if (EZI > 0.0 .and. dAZP > 0.0) then
         Ch_ALGCZ(ich,ia) = &
            (Ch_ZOOCG(ich) / EZI) * (Ch_ALGC_init(ia) / dAZP)
         else
         Ch_ALGCZ(ich,ia) = 0.0
         end if

         ! respiration
         KAR(ich,ia) = KTBALG * KARMAX(ia)
         RK4SOL = calc_RK4SOL(-KAR(ich,ia), Ch_ALGC_init(ia))
         Ch_ALGCR(ich,ia) = Ch_ALGC_init(ia) - RK4SOL

         ! excretion
         KAE(ich,ia) = KTBALG * (1.0 - dLambda_L) * KAEMAX(ia)
         RK4SOL = calc_RK4SOL(-KAE(ich,ia), Ch_ALGC_init(ia))
         Ch_ALGCE(ich,ia) = Ch_ALGC_init(ia) - RK4SOL

         ! mortality
         KAM(ich,ia) = KTBALG * KAMMAX(ia)
         RK4SOL = calc_RK4SOL(-KAM(ich,ia), Ch_ALGC_init(ia))
         Ch_ALGCM(ich,ia) = Ch_ALGC_init(ia) - RK4SOL

         ! settling
         if (dW > 0.0) then
         Ch_ALGCS(ich,ia) = (1.0 - exp(-OMEGAALG*dt/dW)) * Ch_ALGC_init(ia)
         else
         Ch_ALGCS(ich,ia) = 0.0
         end if
      END DO

      ! --- Organic matter
      KTBMBM = calc_KTB(TWATER(ich), TMBM, MMBM)

      ! LPOC grazed by zooplankton
      if (EZI > 0.0 .and. dAZP > 0.0) then
         dLPOCZ = (Ch_ZOOCG(ich) / EZI) * (Ch_LPOC_init / dAZP)
      else
         dLPOCZ = 0.0
      end if

      ! LPOC decomposition
      rate0 = KTBMBM * KTHETA * KLPOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_LPOC_init)
      dLPOCD = Ch_LPOC_init - RK4SOL

      ! LPOC settling
      if (dW > 0.0) then
         dLPOCS = (1.0 - exp(-OMEGAPOM * dt / dW)) * Ch_LPOC_init
      else
         dLPOCS = 0.0
      end if

      ! RPOC decomposition
      rate0 = KTBMBM * KTHETA * KRPOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_RPOC_init)
      dRPOCD = Ch_RPOC_init - RK4SOL

      ! RPOC settling
      if (dW > 0.0) then
         dRPOCS = (1.0 - exp(-OMEGAPOM * dt / dW)) * Ch_RPOC_init
      else
         dRPOCS = 0.0
      end if

      ! LDOC decomposition
      rate0 = KTBMBM * KTHETA * KLDOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_LDOC_init)
      dLDOCD = Ch_LDOC_init - RK4SOL

      ! RDOC decomposition
      rate0 = KTBMBM * KTHETA * KRDOMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_RDOC_init)
      dRDOCD = Ch_RDOC_init - RK4SOL

      ! MBMC decomposition
      rate0 = KTBMBM * KTHETA * KMBMW
      RK4SOL = calc_RK4SOL(-rate0, Ch_MBMC_init)
      dMBMCD = Ch_MBMC_init - RK4SOL

      ! --- outflow: downstream discharge
      ! if (WStorage > 0.0 .and. Qdsch > 0.0) then
      !    dPI = 1.0 - EXP(-2.3 * (Qdsch * dt) / (WStorage + Qdsch * dt))
      ! else
      !    dPI = 0.0
      ! end if
      
      !POM
      dPI = calc_dPI(etaPOMdsch, Qdsch, WStorage)  !BK20260302
      Ch_ALGCdsch(ich,:) = Ch_ALGC_init(:) * dPI
      Ch_ZOOCdsch(ich)   = Ch_ZOOC_init * dPI
      Ch_LPOCdsch(ich)   = Ch_LPOC_init * dPI
      Ch_RPOCdsch(ich)   = Ch_RPOC_init * dPI
      Ch_MBMCdsch(ich)   = Ch_MBMC_init * dPI
      !DOM, DIC
      dPI = calc_dPI(etaDOMdsch, Qdsch, WStorage)  !BK20260302
      Ch_LDOCdsch(ich)   = Ch_LDOC_init * dPI
      Ch_RDOCdsch(ich)   = Ch_RDOC_init * dPI
      Ch_DICdsch(ich)    = Ch_DIC_init * dPI

      ! --- CO2 Evasion (Corrected 20251202)
      if (WStorage > 0.0) then
            ! Henry's Law Constant for CO2
            ! Based on simplified temperature relationship valid for 0-30°C
            ! pKH = -log10(KH) where KH is Henry's constant in mol/(L·atm)
            ! KH = 10^(-pKH) [mol/(L·atm)]
            ! Reference: Simplified from Weiss (1974) for freshwater systems
         pKH = 1.47 - 0.0129 * TWATER(ich)

            ! Saturation Concentration, SatCO2 [kg C/m3]
            ! PCO2 = atmospheric CO2 partial pressure = 4.2e-4 atm (420 ppm)
            ! SatCO2 = KH * PCO2 * MW_C [kg C/m3]
            ! where MW_C = 12 g/mol and unit conversion: g/L -> kg/m3 (factor of 1)
         PCO2 = 4.2e-4  !global average CO2 partial pressure hard-coded
         SatCO2 = (10.0**(-pKH)) * PCO2 * 12.0

         ! Reaeration Coefficient, KO2 [s-1]
         U = SedCNP_hydro%chVa(ich)
         H = dW
         if (H > 0.01) then
         KO2 = 3.93 * (U**0.5) / (H**1.5)  ! O’Connor and Dobbins (1958) 
         else
         KO2 = 0.0
         end if

         ! Reaeration Coefficient - adjust for CO2 and Temperature, KCO2 [s-1]
         KCO2_20 = KO2 * 0.923  ! Thibedoux (1996)
         KCO2_T = KCO2_20 * (1.024**(TWATER(ich) - 20.0))
         KCO2 = KCO2_T / 86400.0

         ! Evasion Flux, Ch_CO2EVAS [kg]
         ! Use analytical solution to prevent instability when KCO2*dt > 1
         ! SIGN CONVENTION (Corrected 20251202):
         !   - Ch_CO2EVAS > 0: Evasion (CO2 leaves water to atmosphere) - TYPICAL
         !   - Ch_CO2EVAS < 0: Invasion (CO2 enters water from atmosphere) - RARE
         ! Most streams are supersaturated with CO2, so Ch_CO2EVAS should be POSITIVE
         ! Typical values: 0.1-10 g C/m2/day or ~0.001-0.1 kg C per timestep per channel
         Ch_CO2EVAS(ich) = (Ch_DIC_init - SatCO2 * WStorage) * &
                           (1.0 - exp(-KCO2 * dt))
      else
         Ch_CO2EVAS(ich) = 0.0
      end if

      ! === B. SCALE FLUXES TO ENSURE MASS CONSERVATION

      ! Phytoplanktons
      DO ia = 1, nalg
         TotalFlux_ALGC(ia) = Ch_ALGCZ(ich,ia) + Ch_ALGCR(ich,ia) + &
            Ch_ALGCE(ich,ia) + Ch_ALGCM(ich,ia) + Ch_ALGCS(ich,ia) + &
            Ch_ALGCdsch(ich,ia)
         if (TotalFlux_ALGC(ia) > Ch_ALGC_init(ia) .and. &
            TotalFlux_ALGC(ia) > 0.0) then
            Scale_ALGC(ia) = Ch_ALGC_init(ia) / TotalFlux_ALGC(ia)
         else
            Scale_ALGC(ia) = 1.0
         end if
         Ch_ALGCZ(ich,ia)    = Ch_ALGCZ(ich,ia) * Scale_ALGC(ia)
         Ch_ALGCR(ich,ia)    = Ch_ALGCR(ich,ia) * Scale_ALGC(ia)
         Ch_ALGCE(ich,ia)    = Ch_ALGCE(ich,ia) * Scale_ALGC(ia)
         Ch_ALGCM(ich,ia)    = Ch_ALGCM(ich,ia) * Scale_ALGC(ia)
         Ch_ALGCS(ich,ia)    = Ch_ALGCS(ich,ia) * Scale_ALGC(ia)
         Ch_ALGCdsch(ich,ia) = Ch_ALGCdsch(ich,ia) * Scale_ALGC(ia)
      ENDDO

      ! Zooplankton
      TotalFlux_ZOOC = Ch_ZOOCZ(ich) + Ch_ZOOCR(ich) + Ch_ZOOCM(ich) + &
         Ch_ZOOCE(ich) + Ch_ZOOCS(ich) + Ch_ZOOCdsch(ich)
      if (TotalFlux_ZOOC > Ch_ZOOC_init .and. TotalFlux_ZOOC > 0.0) then
         Scale_ZOOC = Ch_ZOOC_init / TotalFlux_ZOOC
      else
         Scale_ZOOC = 1.0
      end if
      Ch_ZOOCZ(ich)     = Ch_ZOOCZ(ich) * Scale_ZOOC
      Ch_ZOOCR(ich)     = Ch_ZOOCR(ich) * Scale_ZOOC
      Ch_ZOOCE(ich)     = Ch_ZOOCE(ich) * Scale_ZOOC    
      Ch_ZOOCM(ich)     = Ch_ZOOCM(ich) * Scale_ZOOC
      Ch_ZOOCS(ich)     = Ch_ZOOCS(ich) * Scale_ZOOC
      Ch_ZOOCdsch(ich)  = Ch_ZOOCdsch(ich) * Scale_ZOOC

      ! LPOC
      TotalFlux_LPOC = dLPOCZ + dLPOCD + dLPOCS + Ch_LPOCdsch(ich)
      if (TotalFlux_LPOC > Ch_LPOC_init .and. TotalFlux_LPOC > 0.0) then
         Scale_LPOC = Ch_LPOC_init / TotalFlux_LPOC
      else
         Scale_LPOC = 1.0
      end if
      dLPOCZ           = dLPOCZ * Scale_LPOC
      dLPOCD           = dLPOCD * Scale_LPOC
      dLPOCS           = dLPOCS * Scale_LPOC
      Ch_LPOCdsch(ich) = Ch_LPOCdsch(ich) * Scale_LPOC

      ! LDOC
      TotalFlux_LDOC = dLDOCD + Ch_LDOCdsch(ich)
      if (TotalFlux_LDOC > Ch_LDOC_init .and. TotalFlux_LDOC > 0.0) then
         Scale_LDOC = Ch_LDOC_init / TotalFlux_LDOC
      else
         Scale_LDOC = 1.0
      end if
      dLDOCD           = dLDOCD * Scale_LDOC
      Ch_LDOCdsch(ich) = Ch_LDOCdsch(ich) * Scale_LDOC

      ! RPOC
      TotalFlux_RPOC = dRPOCD + dRPOCS + Ch_RPOCdsch(ich)
      if (TotalFlux_RPOC > Ch_RPOC_init .and. TotalFlux_RPOC > 0.0) then
         Scale_RPOC = Ch_RPOC_init / TotalFlux_RPOC
      else
         Scale_RPOC = 1.0
      end if
      dRPOCD           = dRPOCD * Scale_RPOC
      dRPOCS           = dRPOCS * Scale_RPOC
      Ch_RPOCdsch(ich) = Ch_RPOCdsch(ich) * Scale_RPOC

      ! RDOC
      TotalFlux_RDOC = dRDOCD + Ch_RDOCdsch(ich)
      if (TotalFlux_RDOC > Ch_RDOC_init .and. TotalFlux_RDOC > 0.0) then
         Scale_RDOC = Ch_RDOC_init / TotalFlux_RDOC
      else
         Scale_RDOC = 1.0
      end if
      dRDOCD = dRDOCD * Scale_RDOC
      Ch_RDOCdsch(ich) = Ch_RDOCdsch(ich) * Scale_RDOC

      ! MBMC
      TotalFlux_MBMC = dMBMCD + Ch_MBMCdsch(ich)
      if (TotalFlux_MBMC > Ch_MBMC_init .and. TotalFlux_MBMC > 0.0) then
         Scale_MBMC = Ch_MBMC_init / TotalFlux_MBMC
      else
         Scale_MBMC = 1.0
      end if
      dMBMCD           = dMBMCD * Scale_MBMC
      Ch_MBMCdsch(ich) = Ch_MBMCdsch(ich) * Scale_MBMC

      ! DIC
      TotalFlux_DIC = sum(Ch_ALGCG(ich,:)) + Ch_DICdsch(ich)
      if (Ch_CO2EVAS(ich) > 0.0) TotalFlux_DIC = TotalFlux_DIC + Ch_CO2EVAS(ich)
      if (TotalFlux_DIC > Ch_DIC_init .and. TotalFlux_DIC > 0.0) then
         Scale_DIC = Ch_DIC_init / TotalFlux_DIC
      else
         Scale_DIC = 1.0
      end if
      Ch_ALGCG(ich,:) = Ch_ALGCG(ich,:) * Scale_DIC
      Ch_DICdsch(ich) = Ch_DICdsch(ich) * Scale_DIC
      if (Ch_CO2EVAS(ich)>0.) Ch_CO2EVAS(ich) = Ch_CO2EVAS(ich) * Scale_DIC

      ! === C. UPDATE STORAGES BASED ON SCALED FLUXES
      
      do ia=1,nalg
      Ch_ALGC(ich,ia) = Ch_ALGC_init(ia) + Ch_ALGCG(ich,ia) - &
                        Ch_ALGCZ(ich,ia) - Ch_ALGCR(ich,ia) - &
                        Ch_ALGCE(ich,ia) - Ch_ALGCM(ich,ia) - &
                        Ch_ALGCS(ich,ia) - Ch_ALGCdsch(ich,ia)
      end do
      Ch_ZOOC(ich) = Ch_ZOOC_init + Ch_ZOOCG(ich) - Ch_ZOOCZ(ich) - &
                     Ch_ZOOCR(ich) - Ch_ZOOCM(ich) - &  ! Ch_ZOOCE(ich) excluded!!!
                     Ch_ZOOCS(ich) - Ch_ZOOCdsch(ich)
      Ch_LPOC(ich) = Ch_LPOC_init + sum(Ch_ALGCM(ich,:)) + Ch_ZOOCM(ich) + &
                     Ch_ZOOCE(ich) - dLPOCZ - dLPOCD - dLPOCS - &
                     Ch_LPOCdsch(ich)
      Ch_LDOC(ich) = Ch_LDOC_init + sum(Ch_ALGCE(ich,:)) + FDEC * &
                     (1.0 - FREF) * dLPOCD + FDEC * dMBMCD - dLDOCD - &
                     Ch_LDOCdsch(ich)
      Ch_RPOC(ich) = Ch_RPOC_init + FDEC*FREF * dLPOCD - dRPOCD - dRPOCS - &
                     Ch_RPOCdsch(ich)
      Ch_RDOC(ich) = Ch_RDOC_init + FDEC * (dLDOCD + dRPOCD + dRDOCD) - &
                     dRDOCD - Ch_RDOCdsch(ich)
      Ch_MBMC(ich) = Ch_MBMC_init + FBIO * (dLPOCD + dRPOCD + dLDOCD + &
                     dRDOCD + dMBMCD) - dMBMCD - Ch_MBMCdsch(ich)
      Ch_DIC(ich)  = Ch_DIC_init + FMET * (dLPOCD + dRPOCD + dLDOCD + &
                     dRDOCD + dMBMCD) + Ch_ZOOCR(ich) + &
                     sum(Ch_ALGCR(ich,:)) - sum(Ch_ALGCG(ich,:)) - &
                     Ch_CO2EVAS(ich) - Ch_DICdsch(ich)
      Ch_POCDEPOSIT(ich) = Ch_POCDEPOSIT0(ich) + sum(Ch_ALGCS(ich,:)) + &
                     Ch_ZOOCS(ich) + dLPOCS + dRPOCS

      ! Enforce non-negative channel storages
      do ia=1,nalg
         if (Ch_ALGC(ich,ia) < 1.0e-20) Ch_ALGC(ich,ia) = 0.0
      end do
      if (Ch_ZOOC(ich) < 1.0e-20) Ch_ZOOC(ich) = 0.0
      if (Ch_LPOC(ich) < 1.0e-20) Ch_LPOC(ich) = 0.0
      if (Ch_LDOC(ich) < 1.0e-20) Ch_LDOC(ich) = 0.0
      if (Ch_RPOC(ich) < 1.0e-20) Ch_RPOC(ich) = 0.0
      if (Ch_RDOC(ich) < 1.0e-20) Ch_RDOC(ich) = 0.0
      if (Ch_MBMC(ich) < 1.0e-20) Ch_MBMC(ich) = 0.0
      if (Ch_DIC(ich) < 1.0e-20) Ch_DIC(ich) = 0.0
      if (Ch_POCDEPOSIT(ich) < 1.0e-20) Ch_POCDEPOSIT(ich) = 0.0
    
   end subroutine ChannelCarbonCycle

   
   subroutine CheckCMB_So(i, j, lcover)

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
            So_CSC(i,isl,j) = &
                    (So_LPOCLIT(i,j)     - So_LPOCLIT0(i,j))     &
                  + (So_LPOCRES(i,isl,j) - So_LPOCRES0(i,isl,j)) &
                  + (So_LPOCEXC(i,j)     - So_LPOCEXC0(i,j))     &
                  + (So_LPOCMAN(i,j)     - So_LPOCMAN0(i,j))     &
                  + (So_RPOC(i,isl,j)    - So_RPOC0(i,isl,j))    &
                  + (So_LDOC(i,isl,j)    - So_LDOC0(i,isl,j))    &
                  + (So_RDOC(i,isl,j)    - So_RDOC0(i,isl,j))    &
                  + (So_MBMC(i,isl,j)    - So_MBMC0(i,isl,j))

            So_CMBError(i,isl,j) = &
                    So_LPOCLITin(i,j)    + So_LPOCRESin(i,isl,j) &
                  + So_LPOCEXCin(i,j)    + So_LPOCMANin(i,j)     &
                  + So_LDOCgwso0(i,j)    + So_RDOCgwso0(i,j)     &  !BK20260316
                  - So_LPOCsurf(i,j)     - So_RPOCsurf(i,j)      & 
                  - So_MBMCsurf(i,j)     - So_LDOCsurf(i,j)      &
                  - So_RDOCsurf(i,j)     - So_LDOCintf(i,isl,j)  &
                  - So_RDOCintf(i,isl,j) - So_LDOCperc(i,isl,j)  &
                  - So_RDOCperc(i,isl,j) - So_CO2(i,isl,j)       &
                  - So_CSC(i,isl,j)
        ELSE  
            So_CSC(i,isl,j) = &
                    (So_LPOCRES(i,isl,j) - So_LPOCRES0(i,isl,j)) &
                  + (So_RPOC(i,isl,j)    - So_RPOC0(i,isl,j))    &
                  + (So_LDOC(i,isl,j)    - So_LDOC0(i,isl,j))    &
                  + (So_RDOC(i,isl,j)    - So_RDOC0(i,isl,j))    &
                  + (So_MBMC(i,isl,j)    - So_MBMC0(i,isl,j))

            So_CMBError(i,isl,j) = &
                    So_LPOCRESin(i,isl,j)   + So_LDOCperc0(i,isl-1,j)  &
                  + So_RDOCperc0(i,isl-1,j) - So_LDOCintf(i,isl,j)     &
                  - So_RDOCintf(i,isl,j)    - So_LDOCperc(i,isl,j)     &
                  - So_RDOCperc(i,isl,j)    - So_CO2(i,isl,j)          &
                  - So_CSC(i,isl,j)
        END IF
      END DO

   end subroutine CheckCMB_So


   subroutine CheckCMB_Gw (gwid)

      implicit none
      
      ! --- Input arguments
      integer, intent(in) :: gwid

      Gw_CSC(gwid) = &
              (Gw_LDOC(gwid) - Gw_LDOC0(gwid)) &
            + (Gw_RDOC(gwid) - Gw_RDOC0(gwid))

      Gw_CMBError(gwid) = &
              Gw_LDOCsogw0(gwid) + Gw_RDOCsogw0(gwid) &  !BK20260316
            - Gw_LDOCgwso(gwid)  - Gw_RDOCgwso(gwid)  &  !BK20260316
            - Gw_LDOCgwch(gwid)  - Gw_RDOCgwch(gwid)  &
            - Gw_CSC(gwid)

   end subroutine CheckCMB_Gw

   
   subroutine CheckCMB_Ch (ich)

      implicit none

      ! --- Input arguments
      integer, intent(in) :: ich

      ! --- Local variables
      integer :: chid, gwid
      real(8) :: dt

      dt = real(SedCNPmodel%SedCNP_timestep)
      chid = SedCNP_hydro%linkID(ich)
      gwid = domain%gwid_ch(chid)

      Ch_CSC(ich) = &
              (Ch_ALGC(ich,1)     - Ch_ALGC0(ich,1))     &
            + (Ch_ALGC(ich,2)     - Ch_ALGC0(ich,2))     &
            + (Ch_ALGC(ich,3)     - Ch_ALGC0(ich,3))     &
            + (Ch_ZOOC(ich)       - Ch_ZOOC0(ich))       &
            + (Ch_MBMC(ich)       - Ch_MBMC0(ich))       &
            + (Ch_DIC(ich)        - Ch_DIC0(ich))        &
            + (Ch_LPOC(ich)       - Ch_LPOC0(ich))       &
            + (Ch_RPOC(ich)       - Ch_RPOC0(ich))       &
            + (Ch_LDOC(ich)       - Ch_LDOC0(ich))       &
            + (Ch_RDOC(ich)       - Ch_RDOC0(ich))       &
            + (Ch_POCDEPOSIT(ich) - Ch_POCDEPOSIT0(ich))

      Ch_CMBError(ich) = &
              Ch_ALGCusch0(ich,1) + Ch_ALGCusch0(ich,2) &
            + Ch_ALGCusch0(ich,3) + Ch_ZOOCusch0(ich)   &
            + Ch_MBMCusch0(ich)   + Ch_DICusch0(ich)    &
            + Ch_LPOCusch0(ich)   + Ch_RPOCusch0(ich)   &
            + Ch_LDOCusch0(ich)   + Ch_RDOCusch0(ich)   &
            + Ch_LPOCsurf0(ich)   &
            + Ch_RPOCsurf0(ich)   + Ch_MBMCsurf0(ich)   &
            + Ch_LDOCsurf0(ich)   + Ch_RDOCsurf0(ich)   &
            + Ch_LDOCintf0(ich)   + Ch_RDOCintf0(ich)   &
            + Ch_LDOCgwch0(ich)   + Ch_RDOCgwch0(ich)   &
            + Ch_LPOCxpnt(ich)    + Ch_RPOCxpnt(ich)    &
            + Ch_LDOCxpnt(ich)    + Ch_RDOCxpnt(ich)    &
            + Ch_LPOCxdis(ich)    + Ch_RPOCxdis(ich)    &
            + Ch_LDOCxdis(ich)    + Ch_RDOCxdis(ich)    &
            - Ch_LPOCxabs(ich)    - Ch_RPOCxabs(ich)    &
            - Ch_LDOCxabs(ich)    - Ch_RDOCxabs(ich)    &
            - Ch_ALGCdsch(ich,1)  - Ch_ALGCdsch(ich,2)  &
            - Ch_ALGCdsch(ich,3)  - Ch_ZOOCdsch(ich)    &
            - Ch_MBMCdsch(ich)    - Ch_DICdsch(ich)     &
            - Ch_LPOCdsch(ich)    - Ch_RPOCdsch(ich)    &
            - Ch_LDOCdsch(ich)    - Ch_RDOCdsch(ich)    &
            - Ch_CO2EVAS(ich)     - Ch_CSC(ich)

   end subroutine CheckCMB_Ch

end module Carbon