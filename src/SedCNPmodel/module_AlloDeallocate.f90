module module_AlloDeallocate

  use module_SedCNPvariables

  contains

  !> @brief Allocates temporary variables for the SedCNP model.
  !! @details This subroutine allocates memory for arrays that are used within
  !! a single time step. Their values are not preserved between time steps.
  subroutine Allocate_Temp()

     implicit none

     integer              :: istat

     ! Allocate variables whose previous time-step values are not needed for
     ! the next time step

     ! --------------------------------------------------------------------
     ! Hydrology Output Variables
     ! --------------------------------------------------------------------
     allocate(SedCNP_hydro%rainfall(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%runsrf(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%olwdepth(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%irridepth(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%infltrate(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%evapwater(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%FVEG(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%CanHt(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%vgtyp(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%sltyp(domain%ix,domain%jx), stat=istat)   !BK20240205
     allocate(SedCNP_hydro%slope(domain%ixrt,domain%jxrt), stat=istat) !slope (from RT output)
     allocate(SedCNP_hydro%runsrfRT(domain%ixrt,domain%jxrt), stat=istat) !from RT output  !Y.Kwon 20230427
     allocate(SedCNP_hydro%chVa(domain%nch), stat=istat)
     allocate(SedCNP_hydro%t_soil(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(SedCNP_hydro%theta(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(SedCNP_hydro%percolate(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(SedCNP_hydro%ETRANLayer(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(SedCNP_hydro%WA(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%GWDischarge(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%EvapCanopy(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%EvapSoil(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%Etran(domain%ix,domain%jx), stat=istat)
     allocate(SedCNP_hydro%swdown(domain%ix,domain%jx), stat=istat)  !BK20240129
     allocate(SedCNP_hydro%T2m(domain%ix,domain%jx), stat=istat)   !Y.Kwon 20240515
     !allocate(SedCNP_hydro%q_aqgw(domain%ixrt,domain%jxrt), stat=istat) ! (from RT output)  !Y.Kwon 20230427
     allocate(SedCNP_hydro%q_sogw(domain%ixrt,domain%jxrt), stat=istat) ! (from RT output)  !BK20260313
     allocate(SedCNP_hydro%q_intf(domain%ixrt,domain%jxrt), stat=istat) ! (from RT output)  !Y.Kwon 20230501
     allocate(SedCNP_hydro%q_lateral(domain%nch), stat=istat) ! (from CHRT output)  !BK20240726
     allocate(SedCNP_hydro%q_surf(domain%nch), stat=istat) ! (from CHRT output)  !BK20240620
     allocate(SedCNP_hydro%head(domain%nch), stat=istat) ! (from CHRT output) !Y.Kwon 20230427
     !allocate(SedCNP_hydro%q_base(domain%nch))  ! (from CHRT output) !Y.Kwon 20230427 !BK20240725
     allocate(SedCNP_hydro%q_dsch(domain%nch), stat=istat)   ! (from CHRT output) Y.Kwon 20230427
     allocate(SedCNP_hydro%z_gwbas(domain%nbasin), stat=istat) ! (from GWOUT output) !BK20260305
     allocate(SedCNP_hydro%q_gwch(domain%nbasin), stat=istat) ! (from GWOUT output) !BK20240620 !BK20240726
     allocate(SedCNP_hydro%AREAbasin(domain%nbasin), stat=istat)  ! (from GWBUCKPARM.nc) !Y.Kwon 20230604 !BK20231011
     !allocate(SedCNP_hydro%AREAbasin(domain%nch))  ! (from GWBUCKPARM.nc) !Y.Kwon 20230604
     allocate(SedCNP_hydro%Bucket_max(domain%nbasin), stat=istat)  ! (from GWBUCKPARM.nc) !BK20240510
     allocate(SedCNP_hydro%Bucket_ini(domain%nbasin), stat=istat)  ! (from GWBUCKPARM.nc) !BK20240510
     allocate(SedCNP_hydro%linkID(domain%nch), stat=istat)  ! (from Route_Link.nc) !Y.Kwon 20230604
     allocate(SedCNP_hydro%dwnstrm_linkID(domain%nch), stat=istat)  ! (from Route_Link.nc) !Y.Kwon 20230604
     !allocate(SedCNP_hydro%linkID_grid(domain%ixrt,domain%jxrt))  !BK20240619
     !allocate(SedCNP_hydro%gwbasin(domain%ix,domain%jx)) !!BK20231011  !BK20240619
     !allocate(SedCNP_hydro%gwbasin_ch(domain%ixrt,domain%jxrt))  !BK20231011  !BK20240619
     allocate(SedCNP_hydro%q_usch0(domain%nch), stat=istat)
     allocate(SedCNP_hydro%Wstg_st(domain%nch), stat=istat)
     !allocate(SedCNP_hydro%Wstg_st0(domain%nch))
     allocate(SedCNP_hydro%SFCTMPavg(domain%nbasin), stat=istat)  !Y.Kwon 20240510
     allocate(SedCNP_hydro%SWDOWNavg(domain%nbasin), stat=istat)  !Y.Kwon 20240830
     allocate(SedCNP_hydro%Qpnt(domain%nch,domain%ntime_sedcnp), stat=istat)     !BK20231016
     allocate(SedCNP_hydro%Qabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(SedCNP_hydro%Qdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     !allocate(THETA_S(nst), THETA_F(nst))  !BK20231102
     allocate(DWATER(domain%nch), stat=istat)  !BK20240129
     allocate(TWATER(domain%nch), stat=istat)  !BK20240129
     allocate(DSOIL(domain%nsl), stat=istat)
     allocate(Qsogw(domain%ix,domain%jx), stat=istat)  !BK20260313
     !allocate(st(domain%ix,domain%jx))            !BK20240205
     !allocate(domain%sbid_cell(domain%ix,domain%jx))   !BK20231011  !BK20240610
     !allocate(domain%gwid_ch(domain%nch))              !BK20231011  !BK20240610 !BK20240726 !BK20250701
     !allocate(domain%x_ch(domain%nch))              !BK20231011  !BK20240509
     !allocate(domain%y_ch(domain%nch))              !BK20231011  !BK20240509

     ! --------------------------------------------------------------------
     ! Soil Carbon
     ! --------------------------------------------------------------------
     allocate(So_LPOCLITinso(nlc,12), stat=istat)
     allocate(So_LPOCRESinso(nlc,12), stat=istat)
     allocate(So_LPOCEXCinso(nlc,12), stat=istat)
     allocate(So_LPOCMANinso(nlc,12), stat=istat)  !BK20240117
     allocate(So_LPOCLITin(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOCRESin(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20240630
     allocate(So_LPOCEXCin(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOCMANin(domain%ix,domain%jx), stat=istat)  !BK20240630
     allocate(So_LPOCLIT(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOCRES(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LPOCEXC(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOCMAN(domain%ix,domain%jx), stat=istat)
     allocate(So_RPOC(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LDOC(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_RDOC(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_MBMC(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_CO2(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LPOCsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_RPOCsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_LDOCsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_RDOCsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_MBMCsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_LDOCintf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_RDOCintf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LDOCperc(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20250930 !BK20251018
     allocate(So_RDOCperc(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20240718 !BK20250930 !BK20251018
     allocate(So_LDOCsogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     allocate(So_RDOCsogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     So_LDOCsogw = 0.0; So_RDOCsogw = 0.0  !BK20260508
     allocate(So_LDOCgwso(domain%ix,domain%jx), stat=istat)  !BK20260322
     allocate(So_RDOCgwso(domain%ix,domain%jx), stat=istat)  !BK20260322
     So_LDOCgwso = 0.0; So_RDOCgwso = 0.0  !BK20260508
     allocate(So_CSC(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_CMBError(domain%ix,domain%nsl,domain%jx), stat=istat)

     ! --------------------------------------------------------------------
     ! Soil Nitrogen
     ! --------------------------------------------------------------------
     allocate(So_LPONLITinso(nlc,12), stat=istat)
     allocate(So_LPONRESinso(nlc,12), stat=istat)
     allocate(So_LPONEXCinso(nlc,12), stat=istat)
     allocate(So_LPONMANinso(nlc,12), stat=istat)  !BK20240117
     allocate(So_NH4fert(nlc,12), stat=istat)
     allocate(So_NO3fert(nlc,12), stat=istat)  !BK20240117
     allocate(So_LPONLITin(domain%ix,domain%jx), stat=istat)
     allocate(So_LPONRESin(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20240630
     allocate(So_LPONEXCin(domain%ix,domain%jx), stat=istat)
     allocate(So_LPONMANin(domain%ix,domain%jx), stat=istat)  !BK20240630
     allocate(So_NH4fxin(domain%ix,domain%jx), stat=istat)
     allocate(So_NO3fxin(domain%ix,domain%jx), stat=istat)  !BK20240630
     allocate(So_LPONLIT(domain%ix,domain%jx), stat=istat)
     allocate(So_LPONRES(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LPONEXC(domain%ix,domain%jx), stat=istat)
     allocate(So_LPONMAN(domain%ix,domain%jx), stat=istat)
     allocate(So_RPON(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LDON(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_RDON(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_MBMN(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NH4(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NO3(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NITRI(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_DENIT(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NH3VOL(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NH4uptk(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NO3uptk(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LPONsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_RPONsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_LDONsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_RDONsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_MBMNsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_NH4surf(domain%ix,domain%jx), stat=istat)
     allocate(So_NO3surf(domain%ix,domain%jx), stat=istat)
     allocate(So_LDONintf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_RDONintf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NH4intf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NO3intf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LDONperc(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20250930
     allocate(So_RDONperc(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20250930
     allocate(So_NH4perc(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20250930
     allocate(So_NO3perc(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20240718 !BK20250930
     allocate(So_LDONsogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     allocate(So_RDONsogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     allocate(So_NH4sogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     allocate(So_NO3sogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     So_LDONsogw = 0.0; So_RDONsogw = 0.0  !BK20260508
     So_NH4sogw  = 0.0; So_NO3sogw  = 0.0  !BK20260508
     allocate(So_LDONgwso(domain%ix,domain%jx), stat=istat) !BK20260322
     allocate(So_RDONgwso(domain%ix,domain%jx), stat=istat) !BK20260322
     allocate(So_NH4gwso(domain%ix,domain%jx), stat=istat)  !BK20260322
     allocate(So_NO3gwso(domain%ix,domain%jx), stat=istat)  !BK20260322
     So_LDONgwso = 0.0; So_RDONgwso = 0.0  !BK20260508
     So_NH4gwso  = 0.0; So_NO3gwso  = 0.0  !BK20260508
     allocate(So_NSC(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_NMBError(domain%ix,domain%nsl,domain%jx), stat=istat)

     ! --------------------------------------------------------------------
     ! Soil Phosphorus
     ! --------------------------------------------------------------------
     allocate(So_LPOPLITinso(nlc,12), stat=istat)
     allocate(So_LPOPRESinso(nlc,12), stat=istat)
     allocate(So_LPOPEXCinso(nlc,12), stat=istat)
     allocate(So_LPOPMANinso(nlc,12), stat=istat)  !BK20240117
     allocate(So_PO4fert(nlc,12),stat=istat)
     allocate(So_LPOPLITin(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOPRESin(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20240630
     allocate(So_LPOPEXCin(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOPMANin(domain%ix,domain%jx), stat=istat)  !BK20240630
     allocate(So_PO4fxin(domain%ix,domain%jx), stat=istat)  !BK20240630
     allocate(So_LPOPLIT(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOPRES(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LPOPEXC(domain%ix,domain%jx), stat=istat)
     allocate(So_LPOPMAN(domain%ix,domain%jx), stat=istat)
     allocate(So_RPOP(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LDOP(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_RDOP(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_MBMP(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PO4(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PIPA(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PIPS(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PO4uptk(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PO4ADS(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PIPAADS(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LPOPsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_RPOPsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_LDOPsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_RDOPsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_MBMPsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_PO4surf(domain%ix,domain%jx), stat=istat)
     allocate(So_PIPAsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_PIPSsurf(domain%ix,domain%jx), stat=istat)
     allocate(So_LDOPintf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_RDOPintf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PO4intf(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_LDOPperc(domain%ix,domain%nsl,domain%jx), stat=istat) !BK20250930
     allocate(So_RDOPperc(domain%ix,domain%nsl,domain%jx), stat=istat) !BK20250930
     allocate(So_PO4perc(domain%ix,domain%nsl,domain%jx), stat=istat)  !BK20240718 !BK20250930
     allocate(So_LDOPsogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     allocate(So_RDOPsogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     allocate(So_PO4sogw(domain%ix,domain%jx), stat=istat)  !BK20260318
     So_LDOPsogw = 0.0; So_RDOPsogw = 0.0  !BK20260508
     So_PO4sogw  = 0.0                      !BK20260508
     allocate(So_LDOPgwso(domain%ix,domain%jx), stat=istat) !BK20260322
     allocate(So_RDOPgwso(domain%ix,domain%jx), stat=istat) !BK20260322
     allocate(So_PO4gwso(domain%ix,domain%jx), stat=istat)  !BK20260322
     So_LDOPgwso = 0.0; So_RDOPgwso = 0.0  !BK20260508
     So_PO4gwso  = 0.0                      !BK20260508
     allocate(So_PSC(domain%ix,domain%nsl,domain%jx), stat=istat)
     allocate(So_PMBError(domain%ix,domain%nsl,domain%jx), stat=istat)

    !  ! --------------------------------------------------------------------
    !  ! Aquifer Carbon
    !  ! --------------------------------------------------------------------
    !  allocate(Aq_LDOC(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDOC(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_LDOCaqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDOCaqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_LDOCaqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDOCaqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_CSC(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_CMBError(domain%ix,domain%jx), stat=istat)

    !  ! --------------------------------------------------------------------
    !  ! Aquifer Nitrogen
    !  ! --------------------------------------------------------------------
    !  allocate(Aq_LDON(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDON(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NH4(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NO3(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_LDONaqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDONaqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NH4aqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NO3aqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_LDONaqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDONaqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NH4aqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NO3aqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NSC(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_NMBError(domain%ix,domain%jx), stat=istat)

    !  ! --------------------------------------------------------------------
    !  ! Aquifer Phosphorus
    !  ! --------------------------------------------------------------------
    !  allocate(Aq_LDOP(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDOP(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PO4(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PIPA(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PIPS(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PO4ADS(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PIPAADS(domain%ix,domain%jx), stat=istat)  !BK20241216
    !  allocate(Aq_LDOPaqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDOPaqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PO4aqso(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_LDOPaqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_RDOPaqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PO4aqgw(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PSC(domain%ix,domain%jx), stat=istat)
    !  allocate(Aq_PMBError(domain%ix,domain%jx), stat=istat)

     ! --------------------------------------------------------------------
     ! Groundwater Carbon
     ! --------------------------------------------------------------------
     allocate(Gw_LDOC(domain%nbasin), stat=istat)
     allocate(Gw_RDOC(domain%nbasin), stat=istat)
     allocate(Gw_LDOCgwso(domain%nbasin), stat=istat) !BK20260322
     allocate(Gw_RDOCgwso(domain%nbasin), stat=istat) !BK20260322
     allocate(Gw_LDOCgwch(domain%nbasin), stat=istat)
     allocate(Gw_RDOCgwch(domain%nbasin), stat=istat)
     allocate(Gw_CSC(domain%nbasin), stat=istat)
     allocate(Gw_CMBError(domain%nbasin), stat=istat)

     ! --------------------------------------------------------------------
     ! Groundwater Nitrogen
     ! --------------------------------------------------------------------
     allocate(Gw_LDON(domain%nbasin), stat=istat)
     allocate(Gw_RDON(domain%nbasin), stat=istat)
     allocate(Gw_NH4(domain%nbasin), stat=istat)
     allocate(Gw_NO3(domain%nbasin), stat=istat)
     allocate(Gw_LDONgwso(domain%nbasin), stat=istat) !BK20260322     
     allocate(Gw_RDONgwso(domain%nbasin), stat=istat) !BK20260322
     allocate(Gw_NH4gwso(domain%nbasin), stat=istat)  !BK20260322
     allocate(Gw_NO3gwso(domain%nbasin), stat=istat)  !BK20260322  
     allocate(Gw_LDONgwch(domain%nbasin), stat=istat) 
     allocate(Gw_RDONgwch(domain%nbasin), stat=istat)
     allocate(Gw_NH4gwch(domain%nbasin), stat=istat)
     allocate(Gw_NO3gwch(domain%nbasin), stat=istat)
     allocate(Gw_NSC(domain%nbasin), stat=istat)
     allocate(Gw_NMBError(domain%nbasin), stat=istat)

     ! --------------------------------------------------------------------
     ! Groundwater Phosphorus
     ! --------------------------------------------------------------------
     allocate(Gw_LDOP(domain%nbasin), stat=istat)
     allocate(Gw_RDOP(domain%nbasin), stat=istat)
     allocate(Gw_PO4(domain%nbasin), stat=istat)
     allocate(Gw_LDOPgwso(domain%nbasin), stat=istat) !BK20260322
     allocate(Gw_RDOPgwso(domain%nbasin), stat=istat) !BK20260322
     allocate(Gw_PO4gwso(domain%nbasin), stat=istat)  !BK20260322
     allocate(Gw_LDOPgwch(domain%nbasin), stat=istat)
     allocate(Gw_RDOPgwch(domain%nbasin), stat=istat)
     allocate(Gw_PO4gwch(domain%nbasin), stat=istat)
     allocate(Gw_PSC(domain%nbasin), stat=istat)
     allocate(Gw_PMBError(domain%nbasin), stat=istat)

     ! --------------------------------------------------------------------
     ! Channel Carbon
     ! --------------------------------------------------------------------
     allocate(Ch_ALGC(domain%nch,nalg), stat=istat)
     allocate(Ch_ZOOC(domain%nch), stat=istat)
     allocate(Ch_MBMC(domain%nch), stat=istat)
     allocate(Ch_DIC(domain%nch), stat=istat)
     allocate(Ch_LPOC(domain%nch), stat=istat)
     allocate(Ch_RPOC(domain%nch), stat=istat)
     allocate(Ch_LDOC(domain%nch), stat=istat)
     allocate(Ch_RDOC(domain%nch), stat=istat)
     allocate(Ch_POCDEPOSIT(domain%nch), stat=istat)
     allocate(Ch_ALGCG(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGCR(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGCE(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGCM(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGCZ(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGCS(domain%nch,nalg), stat=istat)
     allocate(Ch_ZOOCG(domain%nch), stat=istat)
     allocate(Ch_ZOOCR(domain%nch), stat=istat)
     allocate(Ch_ZOOCE(domain%nch), stat=istat)
     allocate(Ch_ZOOCM(domain%nch), stat=istat)
     allocate(Ch_ZOOCZ(domain%nch), stat=istat)
     allocate(Ch_ZOOCS(domain%nch), stat=istat)
   !   allocate(Ch_ALGCdsch(domain%nch,nalg), stat=istat) !BK20250917 moved to SedCNP_driver_ini
   !   allocate(Ch_ZOOCdsch(domain%nch), stat=istat)
   !   allocate(Ch_MBMCdsch(domain%nch), stat=istat)
   !   allocate(Ch_LPOCdsch(domain%nch), stat=istat)
   !   allocate(Ch_RPOCdsch(domain%nch), stat=istat)
   !   allocate(Ch_LDOCdsch(domain%nch), stat=istat)
   !   allocate(Ch_RDOCdsch(domain%nch), stat=istat)
     allocate(Ch_LPOCpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPOCpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDOCpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDOCpnt(domain%nch,domain%ntime_sedcnp), stat=istat)  !BK20231015
     allocate(Ch_LPOCxpnt(domain%nch), stat=istat)
     allocate(Ch_RPOCxpnt(domain%nch), stat=istat)
     allocate(Ch_LDOCxpnt(domain%nch), stat=istat)
     allocate(Ch_RDOCxpnt(domain%nch), stat=istat)     !BK20240701
     allocate(Ch_LPOCabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPOCabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDOCabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDOCabs(domain%nch,domain%ntime_sedcnp), stat=istat)     !BK20240707
     allocate(Ch_LPOCxabs(domain%nch), stat=istat)
     allocate(Ch_RPOCxabs(domain%nch), stat=istat)
     allocate(Ch_LDOCxabs(domain%nch), stat=istat)
     allocate(Ch_RDOCxabs(domain%nch), stat=istat)     !BK20240729
     allocate(Ch_LPOCdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPOCdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDOCdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDOCdis(domain%nch,domain%ntime_sedcnp), stat=istat)     !BK20240710
     allocate(Ch_LPOCxdis(domain%nch), stat=istat)
     allocate(Ch_RPOCxdis(domain%nch), stat=istat)
     allocate(Ch_LDOCxdis(domain%nch), stat=istat)
     allocate(Ch_RDOCxdis(domain%nch), stat=istat)     !BK20240729
     allocate(Ch_CSC(domain%nch), stat=istat)
     allocate(Ch_CMBError(domain%nch), stat=istat)

     ! --------------------------------------------------------------------
     ! Channel Common
     ! --------------------------------------------------------------------
     allocate(KAG(domain%nch,nalg), stat=istat)
     allocate(KAR(domain%nch,nalg), stat=istat)
     allocate(KAE(domain%nch,nalg), stat=istat)
     allocate(KAM(domain%nch,nalg), stat=istat)  !BK20231120  !BK20240115
     allocate(KZG(domain%nch), stat=istat)
     allocate(KZR(domain%nch), stat=istat)
     !allocate(KZE(domain%nch), stat=istat)  !BK20251127
     allocate(KZM(domain%nch), stat=istat)  !BK20231120

     ! --------------------------------------------------------------------
     ! Channel Nitrogen
     ! --------------------------------------------------------------------
     allocate(Ch_ALGN(domain%nch,nalg), stat=istat)
     allocate(Ch_ZOON(domain%nch), stat=istat)
     allocate(Ch_MBMN(domain%nch), stat=istat)
     allocate(Ch_NH4(domain%nch), stat=istat)
     allocate(Ch_NO3(domain%nch), stat=istat)
     allocate(Ch_LPON(domain%nch), stat=istat)
     allocate(Ch_RPON(domain%nch), stat=istat)
     allocate(Ch_LDON(domain%nch), stat=istat)
     allocate(Ch_RDON(domain%nch), stat=istat)
     allocate(Ch_PONDEPOSIT(domain%nch), stat=istat)
     allocate(Ch_NITRI(domain%nch), stat=istat)
     allocate(Ch_DENIT(domain%nch), stat=istat)
     allocate(Ch_NH3VOL(domain%nch), stat=istat)
     allocate(Ch_NH4uptk(domain%nch), stat=istat)
     allocate(Ch_NO3uptk(domain%nch), stat=istat)
     allocate(Ch_ALGNG(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGNR(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGNE(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGNM(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGNZ(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGNS(domain%nch,nalg), stat=istat)
     allocate(Ch_ZOONG(domain%nch), stat=istat)
     allocate(Ch_ZOONR(domain%nch), stat=istat)
     allocate(Ch_ZOONE(domain%nch), stat=istat)
     allocate(Ch_ZOONM(domain%nch), stat=istat)
     allocate(Ch_ZOONZ(domain%nch), stat=istat)
     allocate(Ch_ZOONS(domain%nch), stat=istat)
   !   allocate(Ch_ALGNdsch(domain%nch,nalg), stat=istat)  !BK20250917 moved to SedCNP_driver_ini
   !   allocate(Ch_ZOONdsch(domain%nch), stat=istat)
   !   allocate(Ch_MBMNdsch(domain%nch), stat=istat)
   !   allocate(Ch_NH4dsch(domain%nch), stat=istat)
   !   allocate(Ch_NO3dsch(domain%nch), stat=istat)
   !   allocate(Ch_LPONdsch(domain%nch), stat=istat)
   !   allocate(Ch_RPONdsch(domain%nch), stat=istat)
   !   allocate(Ch_LDONdsch(domain%nch), stat=istat)
   !   allocate(Ch_RDONdsch(domain%nch), stat=istat)
     allocate(Ch_LPONpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPONpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDONpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDONpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_NH4pnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_NO3pnt(domain%nch,domain%ntime_sedcnp), stat=istat)    !BK20231015
     allocate(Ch_LPONxpnt(domain%nch), stat=istat)
     allocate(Ch_RPONxpnt(domain%nch), stat=istat)
     allocate(Ch_LDONxpnt(domain%nch), stat=istat)
     allocate(Ch_RDONxpnt(domain%nch), stat=istat)
     allocate(Ch_NH4xpnt(domain%nch), stat=istat)
     allocate(Ch_NO3xpnt(domain%nch), stat=istat)     !BK20240701
     allocate(Ch_LPONabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPONabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDONabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDONabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_NH4abs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_NO3abs(domain%nch,domain%ntime_sedcnp), stat=istat)     !BK20240707
     allocate(Ch_LPONxabs(domain%nch), stat=istat)
     allocate(Ch_RPONxabs(domain%nch), stat=istat)
     allocate(Ch_LDONxabs(domain%nch), stat=istat)
     allocate(Ch_RDONxabs(domain%nch), stat=istat)
     allocate(Ch_NH4xabs(domain%nch), stat=istat)
     allocate(Ch_NO3xabs(domain%nch), stat=istat)     !BK20240729
     allocate(Ch_LPONdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPONdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDONdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDONdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_NH4dis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_NO3dis(domain%nch,domain%ntime_sedcnp), stat=istat)     !BK20240710
     allocate(Ch_LPONxdis(domain%nch), stat=istat)
     allocate(Ch_RPONxdis(domain%nch), stat=istat)
     allocate(Ch_LDONxdis(domain%nch), stat=istat)
     allocate(Ch_RDONxdis(domain%nch), stat=istat)
     allocate(Ch_NH4xdis(domain%nch), stat=istat)
     allocate(Ch_NO3xdis(domain%nch), stat=istat)     !BK20240729
     allocate(Ch_NSC(domain%nch), stat=istat)
     allocate(Ch_NMBError(domain%nch), stat=istat)

     ! --------------------------------------------------------------------
     ! Channel Phosphorus
     ! --------------------------------------------------------------------
     allocate(Ch_ALGP(domain%nch,nalg), stat=istat)
     allocate(Ch_ZOOP(domain%nch), stat=istat)
     allocate(Ch_MBMP(domain%nch), stat=istat)
     allocate(Ch_PO4(domain%nch), stat=istat)
     allocate(Ch_PIPA(domain%nch), stat=istat)
     allocate(Ch_PIPS(domain%nch), stat=istat)
     allocate(Ch_LPOP(domain%nch), stat=istat)
     allocate(Ch_RPOP(domain%nch), stat=istat)
     allocate(Ch_LDOP(domain%nch), stat=istat)
     allocate(Ch_RDOP(domain%nch), stat=istat)
     allocate(Ch_POPDEPOSIT(domain%nch), stat=istat)
     allocate(Ch_PO4uptk(domain%nch), stat=istat)
     allocate(Ch_PO4ADS(domain%nch), stat=istat)
     allocate(Ch_PIPAADS(domain%nch), stat=istat)
     allocate(Ch_ALGPG(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGPR(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGPE(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGPM(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGPZ(domain%nch,nalg), stat=istat)
     allocate(Ch_ALGPS(domain%nch,nalg), stat=istat)
     allocate(Ch_ZOOPG(domain%nch), stat=istat)
     allocate(Ch_ZOOPR(domain%nch), stat=istat)
     allocate(Ch_ZOOPE(domain%nch), stat=istat)
     allocate(Ch_ZOOPM(domain%nch), stat=istat)
     allocate(Ch_ZOOPZ(domain%nch), stat=istat)
     allocate(Ch_ZOOPS(domain%nch), stat=istat)
   !   allocate(Ch_ALGPdsch(domain%nch,nalg), stat=istat)  !BK20250917 moved to SedCNP_driver_ini
   !   allocate(Ch_ZOOPdsch(domain%nch), stat=istat)
   !   allocate(Ch_MBMPdsch(domain%nch), stat=istat)
   !   allocate(Ch_PO4dsch(domain%nch), stat=istat)
   !   allocate(Ch_PIPAdsch(domain%nch), stat=istat)
   !   allocate(Ch_PIPSdsch(domain%nch), stat=istat)
   !   allocate(Ch_LPOPdsch(domain%nch), stat=istat)
   !   allocate(Ch_RPOPdsch(domain%nch), stat=istat)
   !   allocate(Ch_LDOPdsch(domain%nch), stat=istat)
   !   allocate(Ch_RDOPdsch(domain%nch), stat=istat)
     allocate(Ch_LPOPpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPOPpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDOPpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDOPpnt(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_PO4pnt(domain%nch,domain%ntime_sedcnp), stat=istat)    !BK20231015
     allocate(Ch_LPOPxpnt(domain%nch), stat=istat)
     allocate(Ch_RPOPxpnt(domain%nch), stat=istat)
     allocate(Ch_LDOPxpnt(domain%nch), stat=istat)
     allocate(Ch_RDOPxpnt(domain%nch), stat=istat)
     allocate(Ch_PO4xpnt(domain%nch), stat=istat)       !BK20240701
     allocate(Ch_LPOPabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPOPabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDOPabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDOPabs(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_PO4abs(domain%nch,domain%ntime_sedcnp), stat=istat)     !BK20240707
     allocate(Ch_LPOPxabs(domain%nch), stat=istat)
     allocate(Ch_RPOPxabs(domain%nch), stat=istat)
     allocate(Ch_LDOPxabs(domain%nch), stat=istat)
     allocate(Ch_RDOPxabs(domain%nch), stat=istat)
     allocate(Ch_PO4xabs(domain%nch), stat=istat)       !BK20240729
     allocate(Ch_LPOPdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RPOPdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_LDOPdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_RDOPdis(domain%nch,domain%ntime_sedcnp), stat=istat)
     allocate(Ch_PO4dis(domain%nch,domain%ntime_sedcnp), stat=istat)     !BK20240707
     allocate(Ch_LPOPxdis(domain%nch), stat=istat)
     allocate(Ch_RPOPxdis(domain%nch), stat=istat)
     allocate(Ch_LDOPxdis(domain%nch), stat=istat)
     allocate(Ch_RDOPxdis(domain%nch), stat=istat)
     allocate(Ch_PO4xdis(domain%nch), stat=istat)       !BK20240729
     allocate(Ch_PSC(domain%nch), stat=istat)
     allocate(Ch_PMBError(domain%nch), stat=istat)

     allocate(SSA(nps),stat=istat)
     allocate(SedCNP_hydro%Gw_WStorage(domain%nbasin))   !BK20231011
     !allocate(SedCNP_hydro%Gw_WStorage(domain%nch))

     ! --------------------------------------------------------------------
     ! Initialize allocated variables
     ! --------------------------------------------------------------------
     SedCNP_hydro%rainfall   = 0.0
     SedCNP_hydro%runsrf     = 0.0
     SedCNP_hydro%runsrfRT   = 0.0  !Y.Kwon 20230427
     SedCNP_hydro%olwdepth   = 0.0
     SedCNP_hydro%irridepth  = 0.0
     SedCNP_hydro%infltrate  = 0.0
     SedCNP_hydro%evapwater  = 0.0
     SedCNP_hydro%FVEG       = 0.0
     SedCNP_hydro%CanHt      = 0.0
     SedCNP_hydro%t_soil     = 0.0
     SedCNP_hydro%theta      = 0.0
     SedCNP_hydro%percolate  = 0.0
     SedCNP_hydro%ETRANLayer = 0.0
     SedCNP_hydro%WA         = 0.0
     SedCNP_hydro%GWDischarge= 0.0
     SedCNP_hydro%EvapCanopy = 0.0
     SedCNP_hydro%EvapSoil   = 0.0
     SedCNP_hydro%Etran      = 0.0
     SedCNP_hydro%swdown     = 0.0  !BK20240129
     SedCNP_hydro%chVa       = 0.0
     SedCNP_hydro%T2m        = 0.0  !Y.Kwon 20240510
     SedCNP_hydro%q_lateral  = 0.0  !BK20240726
     SedCNP_hydro%q_surf     = 0.0  !BK20240620
     !SedCNP_hydro%q_aqgw     = 0.0 !Y.Kwon 20230427
     SedCNP_hydro%q_sogw     = 0.0  !BK20260313
     SedCNP_hydro%q_intf     = 0.0 !Y.Kwon 20230501
     SedCNP_hydro%head       = 0.0 !Y.Kwon 20230427
     !SedCNP_hydro%q_base     = 0.0 !Y.Kwon 20230427 !BK20240620 !BK20240725
     SedCNP_hydro%q_gwch     = 0.0   !BK20240620 !BK20240726
     SedCNP_hydro%z_gwbas    = 0.0   !BK20260305
     SedCNP_hydro%q_dsch     = 0.0 !Y.Kwon 20230604
     SedCNP_hydro%AREAbasin  = 0.0 !Y.Kwon 20230604
     SedCNP_hydro%q_usch0    = 0.0
     SedCNP_hydro%Wstg_st    = 0.0
     !SedCNP_hydro%Wstg_st0   = 0.0
     SedCNP_hydro%linkID     = -999   !Y.Kwon 20230604
     SedCNP_hydro%dwnstrm_linkID = -999 !Y.Kwon 20230604
     SedCNP_hydro%SFCTMPavg  = -999.0  !Y.Kwon 20240510
     SedCNP_hydro%SWDOWNavg  = 0.0  !Y.Kwon 20240830

  end subroutine Allocate_Temp

  !> @brief Deallocates temporary variables for the SedCNP model.
  !! @details This subroutine deallocates memory for arrays that were used
  !! within a single time step.
  subroutine Deallocate_Temp()

     implicit none

     ! Deallocate variables whose previous time-step values are not needed
     ! for the next time step

     ! --------------------------------------------------------------------
     ! Hydrology Output Variables
     ! --------------------------------------------------------------------
     deallocate(SedCNP_hydro%rainfall)
     deallocate(SedCNP_hydro%runsrf)
     deallocate(SedCNP_hydro%olwdepth)
     deallocate(SedCNP_hydro%irridepth)
     deallocate(SedCNP_hydro%infltrate)
     deallocate(SedCNP_hydro%evapwater)
     deallocate(SedCNP_hydro%FVEG)
     deallocate(SedCNP_hydro%CanHt)
     deallocate(SedCNP_hydro%vgtyp)
     deallocate(SedCNP_hydro%sltyp)   !BK20240205
     deallocate(SedCNP_hydro%slope) !slope (from RT output)
     deallocate(SedCNP_hydro%runsrfRT) !from RT output  !Y.Kwon 20230427
     deallocate(SedCNP_hydro%chVa)
     deallocate(SedCNP_hydro%t_soil)
     deallocate(SedCNP_hydro%theta)
     deallocate(SedCNP_hydro%percolate)
     deallocate(SedCNP_hydro%ETRANLayer)
     deallocate(SedCNP_hydro%WA)
     deallocate(SedCNP_hydro%GWDischarge)
     deallocate(SedCNP_hydro%EvapCanopy)
     deallocate(SedCNP_hydro%EvapSoil)
     deallocate(SedCNP_hydro%Etran)
     deallocate(SedCNP_hydro%swdown)  !BK20240129
     deallocate(SedCNP_hydro%T2m)  !Y.Kwon 20240510
     !deallocate(SedCNP_hydro%q_aqgw)
     deallocate(SedCNP_hydro%q_sogw) !WHQ 20260313
     deallocate(SedCNP_hydro%q_lateral)  !BK20240726
     deallocate(SedCNP_hydro%q_surf)  !BK20240620
     deallocate(SedCNP_hydro%q_intf)
     deallocate(SedCNP_hydro%head)
     !deallocate(SedCNP_hydro%q_base)  !BK20240620 !BK20240725
     deallocate(SedCNP_hydro%q_gwch)  !BK20240620 !BK20240726
     deallocate(SedCNP_hydro%z_gwbas)  !BK20260305
     deallocate(SedCNP_hydro%q_dsch)
     deallocate(SedCNP_hydro%AREAbasin)
     deallocate(SedCNP_hydro%Bucket_max)  !BK20240510
     deallocate(SedCNP_hydro%Bucket_ini)  !BK20240510
     deallocate(SedCNP_hydro%linkID)
     deallocate(SedCNP_hydro%dwnstrm_linkID)
     deallocate(SedCNP_hydro%q_usch0)  !BK20250702
     deallocate(SedCNP_hydro%Wstg_st)
     deallocate(SedCNP_hydro%SFCTMPavg)  !Y.Kwon 20240510
     deallocate(SedCNP_hydro%SWDOWNavg)  !Y.Kwon 20240830
     deallocate(SedCNP_hydro%Qpnt)  !BK20250702
     deallocate(SedCNP_hydro%Qabs)  !BK20250702
     deallocate(SedCNP_hydro%Qdis)  !BK20250702

     ! --------------------------------------------------------------------
     ! Sediment and CNP common variables
     ! --------------------------------------------------------------------
     !deallocate(channelSed%Spnt)   !BK20231015
     !deallocate(channelSed%Sxpnt)    !BK20240701
     deallocate(DWATER, TWATER)  !BK20240129
     deallocate(DSOIL)
     deallocate(Qsogw)  !BK20260313
     !deallocate(st)  !BK20240205

     ! --------------------------------------------------------------------
     ! Soil Carbon
     ! --------------------------------------------------------------------
     deallocate(So_LPOCLIT, So_LPOCRES, So_LPOCEXC, So_LPOCMAN)
     deallocate(So_RPOC, So_LDOC, So_RDOC, So_MBMC, So_CO2)
     deallocate(So_LPOCLITinso, So_LPOCRESinso, So_LPOCEXCinso, &
                So_LPOCMANinso)
     deallocate(So_LPOCLITin, So_LPOCRESin, So_LPOCEXCin, So_LPOCMANin)  !BK20240630
     deallocate(So_LDOCperc, So_RDOCperc)
     deallocate(So_LPOCsurf, So_RPOCsurf, So_LDOCsurf, So_RDOCsurf, &
                So_MBMCsurf)
     deallocate(So_LDOCintf, So_RDOCintf)
     deallocate(So_LDOCsogw, So_RDOCsogw)    !BK20260318
     deallocate(So_LDOCgwso, So_RDOCgwso) !BK20260322
     deallocate(So_CSC, So_CMBError)

     ! --------------------------------------------------------------------
     ! Soil Nitrogen
     ! --------------------------------------------------------------------
     deallocate(So_LPONLIT, So_LPONRES, So_LPONEXC, So_LPONMAN)
     deallocate(So_RPON, So_LDON, So_RDON, So_MBMN)
     deallocate(So_NH4, So_NO3)
     deallocate(So_LPONLITinso, So_LPONRESinso, So_LPONEXCinso, &
                So_LPONMANinso)
     deallocate(So_NH4fert, So_NO3fert)
     deallocate(So_LPONLITin, So_LPONRESin, So_LPONEXCin, So_LPONMANin)  !BK20240630
     deallocate(So_NH4fxin, So_NO3fxin)  !BK20240630
     deallocate(So_NITRI, So_DENIT, So_NH3VOL, So_NH4uptk, So_NO3uptk)
     deallocate(So_LDONperc, So_RDONperc, So_NH4perc, So_NO3perc)
     deallocate(So_LPONsurf, So_RPONsurf, So_LDONsurf, So_RDONsurf, &
                So_MBMNsurf)
     deallocate(So_NH4surf, So_NO3surf)
     deallocate(So_LDONintf, So_RDONintf, So_NH4intf, So_NO3intf)
     deallocate(So_LDONsogw, So_RDONsogw, So_NH4sogw, So_NO3sogw)    !BK20260318
     deallocate(So_LDONgwso, So_RDONgwso, So_NH4gwso, So_NO3gwso) !BK20260322
     deallocate(So_NSC, So_NMBError)

     ! --------------------------------------------------------------------
     ! Soil Phosphorus
     ! --------------------------------------------------------------------
     deallocate(So_LPOPLIT, So_LPOPRES, So_LPOPEXC, So_LPOPMAN)
     deallocate(So_RPOP, So_LDOP, So_RDOP, So_MBMP)
     deallocate(So_PO4, So_PIPA, So_PIPS)
     deallocate(So_LPOPLITinso, So_LPOPRESinso, So_LPOPEXCinso, &
                So_LPOPMANinso)
     deallocate(So_PO4fert)
     deallocate(So_LPOPLITin, So_LPOPRESin, So_LPOPEXCin, So_LPOPMANin)  !BK20240630
     deallocate(So_PO4fxin)  !BK20240630
     deallocate(So_PO4uptk, So_PO4ADS, So_PIPAADS)
     deallocate(So_LDOPperc, So_RDOPperc, So_PO4perc)
     deallocate(So_LPOPsurf, So_RPOPsurf, So_LDOPsurf, So_RDOPsurf, &
                So_MBMPsurf)
     deallocate(So_PO4surf, So_PIPAsurf, So_PIPSsurf)
     deallocate(So_LDOPintf, So_RDOPintf, So_PO4intf)
     deallocate(So_LDOPsogw, So_RDOPsogw, So_PO4sogw)  !BK20260318
     deallocate(So_LDOPgwso, So_RDOPgwso, So_PO4gwso) !BK20260322
     deallocate(So_PSC, So_PMBError)

    !  ! --------------------------------------------------------------------
    !  ! Aquifer Carbon
    !  ! --------------------------------------------------------------------
    !  deallocate(Aq_LDOC, Aq_RDOC)
    !  deallocate(Aq_LDOCaqso, Aq_RDOCaqso)
    !  deallocate(Aq_LDOCaqgw, Aq_RDOCaqgw)
    !  deallocate(Aq_CSC, Aq_CMBError)

    !  ! --------------------------------------------------------------------
    !  ! Aquifer Nitrogen
    !  ! --------------------------------------------------------------------
    !  deallocate(Aq_LDON, Aq_RDON, Aq_NH4, Aq_NO3)
    !  deallocate(Aq_LDONaqso, Aq_RDONaqso, Aq_NH4aqso, Aq_NO3aqso)
    !  deallocate(Aq_LDONaqgw, Aq_RDONaqgw, Aq_NH4aqgw, Aq_NO3aqgw)
    !  deallocate(Aq_NSC, Aq_NMBError)

    !  ! --------------------------------------------------------------------
    !  ! Aquifer Phosphorus
    !  ! --------------------------------------------------------------------
    !  deallocate(Aq_LDOP, Aq_RDOP, Aq_PO4, Aq_PIPA, Aq_PIPS)
    !  deallocate(Aq_PO4ADS, Aq_PIPAADS)  !BK20241216
    !  deallocate(Aq_LDOPaqso, Aq_RDOPaqso, Aq_PO4aqso)
    !  deallocate(Aq_LDOPaqgw, Aq_RDOPaqgw, Aq_PO4aqgw)
    !  deallocate(Aq_PSC, Aq_PMBError)

     ! --------------------------------------------------------------------
     ! Groundwater Carbon
     ! --------------------------------------------------------------------
     deallocate(Gw_LDOC, Gw_RDOC)
     deallocate(Gw_LDOCgwch, Gw_RDOCgwch)
     deallocate(Gw_LDOCgwso, Gw_RDOCgwso) !BK20260322
     deallocate(Gw_CSC, Gw_CMBError)

     ! --------------------------------------------------------------------
     ! Groundwater Nitrogen
     ! --------------------------------------------------------------------
     deallocate(Gw_LDON, Gw_RDON, Gw_NH4, Gw_NO3)
     deallocate(Gw_LDONgwso, Gw_RDONgwso, Gw_NH4gwso, Gw_NO3gwso) !BK20260322
     deallocate(Gw_LDONgwch, Gw_RDONgwch, Gw_NH4gwch, Gw_NO3gwch)
     deallocate(Gw_NSC, Gw_NMBError)

     ! --------------------------------------------------------------------
     ! Groundwater Phosphorus
     ! --------------------------------------------------------------------
     deallocate(Gw_LDOP, Gw_RDOP, Gw_PO4)
     deallocate(Gw_LDOPgwso, Gw_RDOPgwso, Gw_PO4gwso) !BK20260322
     deallocate(Gw_LDOPgwch, Gw_RDOPgwch, Gw_PO4gwch)
     deallocate(Gw_PSC, Gw_PMBError)

     ! --------------------------------------------------------------------
     ! Channel Common
     ! --------------------------------------------------------------------
     !deallocate(abs2dis_ich)
     !deallocate(abs2dis_lag)

     ! --------------------------------------------------------------------
     ! Channel Carbon
     ! --------------------------------------------------------------------
     deallocate(Ch_ALGC, Ch_ZOOC, Ch_MBMC, Ch_DIC)
     deallocate(Ch_LPOC, Ch_RPOC, Ch_LDOC, Ch_RDOC, Ch_POCDEPOSIT)
     deallocate(Ch_ALGCG, Ch_ALGCR, Ch_ALGCE, Ch_ALGCM, Ch_ALGCZ)
     deallocate(Ch_ALGCS)
     deallocate(Ch_ZOOCG, Ch_ZOOCR, Ch_ZOOCE, Ch_ZOOCM, Ch_ZOOCZ, Ch_ZOOCS)
     deallocate(Ch_ALGCdsch, Ch_ZOOCdsch, Ch_MBMCdsch, Ch_DICdsch)
     deallocate(Ch_LPOCdsch, Ch_RPOCdsch, Ch_LDOCdsch, Ch_RDOCdsch)
     deallocate(Ch_LPOCpnt, Ch_RPOCpnt, Ch_LDOCpnt, Ch_RDOCpnt)    !BK20250702
     deallocate(Ch_LPOCabs, Ch_RPOCabs, Ch_LDOCabs, Ch_RDOCabs)    !BK20250702
     deallocate(Ch_LPOCdis, Ch_RPOCdis, Ch_LDOCdis, Ch_RDOCdis)    !BK20250702
     deallocate(Ch_LPOCxpnt, Ch_RPOCxpnt, Ch_LDOCxpnt, Ch_RDOCxpnt)    !BK20240701
     deallocate(Ch_LPOCxabs, Ch_RPOCxabs, Ch_LDOCxabs, Ch_RDOCxabs)    !BK20240729
     deallocate(Ch_LPOCxdis, Ch_RPOCxdis, Ch_LDOCxdis, Ch_RDOCxdis)    !BK20240729
     deallocate(Ch_CSC, Ch_CMBError)
     deallocate(KAG, KAR, KAE, KAM)  !BK20231120
     !deallocate(KZG, KZR, KZE, KZM)  !BK20231120
     deallocate(KZG, KZR, KZM)  !BK20231120 !BK20251127

     ! --------------------------------------------------------------------
     ! Channel Nitrogen
     ! --------------------------------------------------------------------
     deallocate(Ch_ALGN, Ch_ZOON, Ch_MBMN, Ch_NH4, Ch_NO3)
     deallocate(Ch_LPON, Ch_RPON, Ch_LDON, Ch_RDON, Ch_PONDEPOSIT)
     deallocate(Ch_NITRI, Ch_DENIT, Ch_NH3VOL, Ch_NH4uptk, Ch_NO3uptk)
     deallocate(Ch_ALGNusch0, Ch_ZOONusch0, Ch_MBMNusch0, Ch_NH4usch0, &
                Ch_NO3usch0)
     deallocate(Ch_ALGNG, Ch_ALGNR, Ch_ALGNE, Ch_ALGNM, Ch_ALGNZ)
     deallocate(Ch_ALGNS)
     deallocate(Ch_ZOONG, Ch_ZOONR, Ch_ZOONE, Ch_ZOONM, Ch_ZOONZ, Ch_ZOONS)
     deallocate(Ch_ALGNdsch, Ch_ZOONdsch, Ch_MBMNdsch, Ch_NH4dsch, &
                Ch_NO3dsch)
     deallocate(Ch_LPONdsch, Ch_RPONdsch, Ch_LDONdsch, Ch_RDONdsch)
     deallocate(Ch_LPONpnt, Ch_RPONpnt, Ch_LDONpnt, Ch_RDONpnt, Ch_NH4pnt, &
                Ch_NO3pnt)    !BK20250702
     deallocate(Ch_LPONabs, Ch_RPONabs, Ch_LDONabs, Ch_RDONabs, Ch_NH4abs, &
                Ch_NO3abs)    !BK20250702
     deallocate(Ch_LPONdis, Ch_RPONdis, Ch_LDONdis, Ch_RDONdis, Ch_NH4dis, &
                Ch_NO3dis)    !BK20250702
     deallocate(Ch_LPONxpnt, Ch_RPONxpnt, Ch_LDONxpnt, Ch_RDONxpnt, &
                Ch_NH4xpnt, Ch_NO3xpnt)    !BK20240701
     deallocate(Ch_LPONxabs, Ch_RPONxabs, Ch_LDONxabs, Ch_RDONxabs, &
                Ch_NH4xabs, Ch_NO3xabs)    !BK20240729
     deallocate(Ch_LPONxdis, Ch_RPONxdis, Ch_LDONxdis, Ch_RDONxdis, &
                Ch_NH4xdis, Ch_NO3xdis)    !BK20240729
     deallocate(Ch_NSC, Ch_NMBError)

     ! --------------------------------------------------------------------
     ! Channel Phosphorus
     ! --------------------------------------------------------------------
     deallocate(Ch_ALGP, Ch_ZOOP, Ch_MBMP)
     deallocate(Ch_PO4, Ch_PIPA, Ch_PIPS)
     deallocate(Ch_LPOP, Ch_RPOP, Ch_LDOP, Ch_RDOP, Ch_POPDEPOSIT)
     deallocate(Ch_PO4uptk, Ch_PO4ADS, Ch_PIPAADS)
     deallocate(Ch_ALGPG, Ch_ALGPR, Ch_ALGPE, Ch_ALGPM)
     deallocate(Ch_ALGPZ, Ch_ALGPS)
     deallocate(Ch_ZOOPG, Ch_ZOOPR, Ch_ZOOPE, Ch_ZOOPM, Ch_ZOOPZ, Ch_ZOOPS)
     deallocate(Ch_ALGPdsch, Ch_ZOOPdsch, Ch_MBMPdsch)
     deallocate(Ch_PO4dsch, Ch_PIPAdsch, Ch_PIPSdsch)
     deallocate(Ch_LPOPdsch, Ch_RPOPdsch, Ch_LDOPdsch, Ch_RDOPdsch)
     deallocate(Ch_LPOPpnt, Ch_RPOPpnt, Ch_LDOPpnt, Ch_RDOPpnt, &
                Ch_PO4pnt)   !BK20250702
     deallocate(Ch_LPOPabs, Ch_RPOPabs, Ch_LDOPabs, Ch_RDOPabs, &
                Ch_PO4abs)   !BK20250702
     deallocate(Ch_LPOPdis, Ch_RPOPdis, Ch_LDOPdis, Ch_RDOPdis, &
                Ch_PO4dis)   !BK20250702
     deallocate(Ch_LPOPxpnt, Ch_RPOPxpnt, Ch_LDOPxpnt, Ch_RDOPxpnt, &
                Ch_PO4xpnt)    !BK20240701
     deallocate(Ch_LPOPxabs, Ch_RPOPxabs, Ch_LDOPxabs, Ch_RDOPxabs, &
                Ch_PO4xabs)    !BK20240729
     deallocate(Ch_LPOPxdis, Ch_RPOPxdis, Ch_LDOPxdis, Ch_RDOPxdis, &
                Ch_PO4xdis)    !BK20240729
     deallocate(Ch_PSC, Ch_PMBError)

     deallocate(SSA)

  end subroutine Deallocate_Temp

end module module_AlloDeallocate
