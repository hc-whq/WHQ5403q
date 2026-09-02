!  Program Name:
!  Author(s)/Contact(s): Y.Kwon
!  Abstract:
!  History Log:
!
!  Usage:
!  Parameters: <Specify typical arguments passed>
!  Input Files:
!        <list file names and briefly describe the data they include>
!  Output Files:
!        <list file names and briefly describe the information they include>
!
!  Condition codes:
!        <list exit condition or error codes returned >
!        If appropriate, descriptive troubleshooting instructions or
!        likely causes for failures could be mentioned here with the
!        appropriate error code
!
!  User controllable options: <if applicable>

module module_SedCNPvariables

  implicit none

  !domain
  type domainSedCNP
     integer              :: ix, jx       !for LSM output
     integer              :: ixrt, jxrt   !for terrain RT output
     integer              :: nch          !number of channel
     integer              :: nbasin       !number of gw basin
     integer              :: ntime_sedcnp !Y.Kwon20230316
     integer              :: nsl          !number of soil layers (default nsl = 4)
     !integer, allocatable :: sbid_cell(:,:) !subbasin ID (sbid) for grid cell (i,j)  !BK20231011
     integer, allocatable :: gwid_ch(:)   !gw basin ID (gwid) for a given channel ID     !BK20231011 !BK20240726
     integer, allocatable :: chid_i(:)   ! Grid i-index for each channel ID  !BK20250701
     integer, allocatable :: chid_j(:)   ! Grid j-index for each channel ID  !BK20250701
     !integer, allocatable :: x_ch(:)      !x-coord. of the representative cell for a given channel  !BK20231029  !BK20240509
     !integer, allocatable :: y_ch(:)      !y-coord. of the representative cell for a given channel  !BK20231029  !BK20240509
     real(8)              :: areaxy       !cell area (m2)
     real(8)              :: DT           !time-step (sec)
     integer              :: ISWATER = -999 !ISWATER parameter, read from MPTABLE.TBL !BK20251208
     
  end type
  type(domainSedCNP) domain

  !date
  type dateSedCNPforHydroOutput
     character(len=19)  :: olddate, newdate
  end type
  type(dateSedCNPforHydroOutput) dateSedCNP 

  !Hydrology outputs
  type hydro_output_forSedCNP
     real(8), allocatable :: rainfall(:,:)     ! [mm/s]
     real(8), allocatable :: runsrf(:,:)       ! [mm/s]
     real(8), allocatable :: runsrfRT(:,:)     ! [mm]
     real(8), allocatable :: olwdepth(:,:)     ! [mm]
     real(8), allocatable :: irridepth(:,:)    ! [mm]
     real(8), allocatable :: infltrate(:,:)    ! [mm/s]
     real(8), allocatable :: evapwater(:,:)    ! [mm/s]
     real(8), allocatable :: slope(:,:)        ! land slope [-] (from terrain RT output)
     real(8), allocatable :: FVEG(:,:)         ! Green vegetation fraction [-]
     real(8), allocatable :: CanHt(:,:)        ! vegetatio height [m]
     real(8), allocatable :: vgtyp(:,:)        ! dominant vegetation category
     real(8), allocatable :: sltyp(:,:)        ! dominant soil category !BK20240205
     !integer(8), allocatable :: vgtyp(:,:)        ! dominant vegetation category   !BK20240915
     !integer(8), allocatable :: sltyp(:,:)        ! dominant soil category !BK20240205 !BK20240915
     real(8), allocatable :: chVa(:)           ! channel velocity [m/s]

     real(8), allocatable :: t_soil(:,:,:)     ! soil layer temperature [K]
     real(8), allocatable :: theta(:,:,:)      ! volumetric soil moisture [m3 m-3]
     real(8), allocatable :: smcmax(:,:,:)     ! soil water content at saturation [m3 m-3]
     real(8), allocatable :: smcref(:,:,:)     ! soil water content at field capacity [m3 m-3]
     real(8), allocatable :: percolate(:,:,:)  ! percolation between soil layers [mm s-1]
     real(8), allocatable :: ETRANLayer(:,:,:) ! vegetation water uptake from soil layers [mm s-1]
     real(8), allocatable :: WA(:,:)           ! water storage in aquifer [mm]
     real(8), allocatable :: GWDischarge(:,:)  ! groundwater discharge [mm s-1]
     real(8), allocatable :: EvapCanopy(:,:)   ! canopy water evaporation [mm/s]
     real(8), allocatable :: EvapSoil(:,:)     ! evaporation from soil [mm/s]
     real(8), allocatable :: Etran(:,:)        ! transpiration [mm/s]
     real(8), allocatable :: T2m(:,:)          ! 2m Air Temperature [K]    !Y.Kwon 20240510


     !--- better be  moved to domainSedCNP??? --- !BK20240509
     integer, allocatable :: subbasinID(:,:)   ! subbasin ID for grid cells
     !integer, allocatable :: Basin(:)          ! subbasin ID (1...n)     !BK20240509
     real(8), allocatable :: AREAbasin(:)      ! subbasin area [km2]
     real(8), allocatable :: Bucket_max(:)     ! subbasin bucket height [mm]     !BK20240509
     real(8), allocatable :: Bucket_ini(:)     ! initial height of water in the bucket [mm]     !BK20240509
     integer, allocatable :: linkID(:)         ! channel link ID for reach-based links        !Y.Kwon 20230604
     integer, allocatable :: dwnstrm_linkID(:) ! downstream (channel) ID !Y.Kwon 20230604
     integer, allocatable :: linkID_grid(:,:)  !channel link ID for grid cells, from Fulldom_hires.nc  !BK20240616
     integer, allocatable :: gwbasin(:,:)    !gw basin ID for LSM grids (ix,jx)                !BK20240613
     integer, allocatable :: gwbasin_ch(:,:) !gw basin ID  for channel grids (ixrt,jxrt)       !BK20240613
     !--- better be moved to domainSedCNP??? --- !BK20240509



     real(8), allocatable :: head(:)           ! river stage [m]
     !real(8), allocatable :: qGwBucket(:)      ! flux from gw bucket [m3/s]
     real(8), allocatable :: q_lateral(:)      ! lateral influx to channel (= q_gwch + q_surf) [m3/s]  !BK20240726
     real(8), allocatable :: z_gwbas(:)         ! depth in gw bucket [m]  !BK20260305
     real(8), allocatable :: q_gwch(:)         ! flux from gw bucket to channel [m3/s]  !BK20240620 !BK20240726
     !real(8), allocatable :: q_base(:)         ! flux from gw bucket to channel [m3/s]  !BK20240620
     real(8), allocatable :: q_surf(:)         ! flux from surface runoff to channel [m3/s]
     real(8), allocatable :: q_dsch(:)         ! channel discharge [m3/s]  !Y.Kwon 20230604
     !real(8), allocatable :: q_aqgw(:,:)       ! flux from aquifer to gw bucket [mm] !Y.Kwon 20230501
     real(8), allocatable :: q_sogw(:,:)       ! flux from soil to gw bucket [mm] !BK20260313
     real(8), allocatable :: q_intf(:,:)       ! flux from inter flow (subsurface flow --> exfiltration) [mm s-1] !Y.Kwon20230604
     real(8), allocatable :: q_usch0(:)        ! channel inflow from upstream (previous time step) [m3/s]
     real(8), allocatable :: Wstg_st(:)        ! stream water storage at current time step [m3]
     !real(8), allocatable :: Wstg_st0(:)       ! stream water storage at previous time step [m3]

     !real(8), allocatable :: Aq_WStorage(:,:)                  ! (m3)
     real(8), allocatable :: Gw_WStorage(:)    ! (m3)   !BK
     
     real(8), allocatable :: Qpnt(:,:)        ! water influx to channel from point sources (m3 s-1)   !BK20231016
     !real(8), allocatable :: Qabs_max(:,:)     ! water abstraction capacity (m3 s-1)   !BK20240703 !BK20241016
     real(8), allocatable :: Qabs(:,:)         ! water outflux from channel by abstraction (m3 s-1) !BK20240703
     real(8), allocatable :: Qdis(:,:)         ! water influx to channel by discharge (m3 s-1) !BK20240703
     real(8), allocatable :: swdown(:,:)       ! short-wave downward solar radiation (W m-2)   !BK20240129
     real(8), allocatable :: SFCTMPavg(:)      ! basin-averaged surface air temperature (K)  !Y.Kwon 20240510
     real(8), allocatable :: SWDOWNavg(:)      ! basin-averaged downward solar radiation (W m-2)  !Y.Kwon 20240830
     integer              :: time(1)           ! minutes since 1970-01-01 00:00:00 UTC   !Y.Kwon20241125
     real(8)              :: dx(1)             ! grid spacing (in m) of the water outputs from LSMs  !Y.Kwon20250621
  end type hydro_output_forSedCNP
  type(hydro_output_forSedCNP) SedCNP_hydro

  !channel inputs
  type channel_input_forSedCNP
     real(8), allocatable :: Wst(:)         !bottom width of channel [m]
     real(8), allocatable :: SideSlopch(:)  !Channel side slope
     !real(8), allocatable :: TopWdth(:)     !Top width of channel [m]        !BK20251117
     !real(8), allocatable :: TopWdthCC(:)   !Compound Channel Top Width [m]  !BK20251117
     real(8), allocatable :: Lst(:)         !stream length [m]
     real(8), allocatable :: Sf(:)          !channel slope
  end type channel_input_forSedCNP
  type(channel_input_forSedCNP) channel_input

  !for sediment
  type overlandSed_struct
     !real(8), allocatable :: Sow0(:,:,:)   !overland water sediment storage [kg] (grain_size,lon,lat)  !BK20240618 !BK20250925
     !real(8), allocatable :: Sol0(:,:,:)   !overland sediment storage [kg] (grain_size,lon,lat)  !BK20240618 !BK20250925
     real(8), allocatable :: Ssurf0(:,:,:) !overland sediment transport [kg/s] (grain_size,lon,lat)  !BK20240618 
     real(8), allocatable :: Sow(:,:,:)    !overland water sediment storage [kg] (grain_size,lon,lat)
     real(8), allocatable :: Sol(:,:,:)    !overland sediment storage [kg] (grain_size,lon,lat)
     real(8), allocatable :: Solcp(:,:,:)  !sediment transport capacity [kg/s] (grain_size,lon,lat)
     real(8), allocatable :: Ssurf(:,:,:)  !overland sediment transport [kg/s] (grain_size,lon,lat)
     real(8), allocatable :: Srdet(:,:,:)  !sediment detachment by rainfall and irrigation water [kg/s] (grain_size,lon,lat)
     real(8), allocatable :: Soler(:,:,:)  !sediment erosion [kg/s] (grain_size,lon,lat)
     real(8), allocatable :: Soldp(:,:,:)  !sediment deposition [kg/s] (grain_size,lon,lat
     real(8), allocatable :: Colsd(:,:,:)  !concentration of sediment particles in the overland flow [kg/m3] (grain_size,lon,lat)
     real(8), allocatable :: SOILPSF(:,:)  !soil particle size fractions [-]  !BK20240205 
  end type overlandSed_struct
  type(overlandSed_struct) overSed

  type channelSed_struct
     real(8), allocatable :: Scw0(:,:)    !channel water sediment storage (grain_size,nch)  !BK20240618
     real(8), allocatable :: Sch0(:,:)    !channel sediment storage (grain_size,nch)  !BK20240618
     real(8), allocatable :: Susch0(:,:)  !upstream discharge at previous time step [kg] (grain_size,nch)  !BK20240618
     real(8), allocatable :: Scw(:,:)     !channel water sediment storage (grain_size,nch)
     real(8), allocatable :: Sch(:,:)     !channel sediment storage (grain_size,nch)
     !real(8), allocatable :: Susch(:,:)   !upstream discharge [kg] (grain_size,nch)  !BK20250925
     real(8), allocatable :: Spnt(:,:,:) !sediment influx from point sources [kg s-1] (grain_size,nch,ntime)  !BK20231016
     real(8), allocatable :: Sxpnt(:,:)   !sediment influx from point sources [kg] (grain_size,nch)  !BK20240701
     real(8), allocatable :: Sabs(:,:,:)  !sediment outflux by abstraction [kg s-1] (nps,nch,ntime) !BK20240703
     real(8), allocatable :: Sxabs(:,:)   !sediment outflux by abstraction [kg] (grain_size,nch)  !BK20240729
     real(8), allocatable :: Sdis(:,:,:)  !sediment influx by discharge [kg s-1] (nps,nch,ntime) !BK20240703
     real(8), allocatable :: Sxdis(:,:)   !sediment influx by discharge [kg] (grain_size,nch)  !BK20240729
     !real(8), allocatable :: Solch(:,:)   !overland to channel [kg] (grain_size,nch)
     real(8), allocatable :: Solch0(:,:)   !overland to channel [kg] (grain_size,nch)  !BK20250925
     !real(8), allocatable :: Sirri(:,:)    !irrigation [kg/s] (grain_size,nch)  !BK20260510 removed Sirri
     real(8), allocatable :: Sinput(:,:)  !total sediment input [kg] (grain_size,nch)
     real(8), allocatable :: Sdsch(:,:)   !downstream discharge [kg/s] (grain_size,nch)
     real(8), allocatable :: Scher(:,:)   !channel erosion [kg/s] (grain_size,nch)
     real(8), allocatable :: Schdp(:,:)   !channel deposition [kg/s] (grain_size,nch)
     real(8), allocatable :: Schcp(:,:)   !channel capacity [kg/s] (grain_size,nch)
     real(8), allocatable :: Cchsd(:,:)   !channel sediment concentration [kg/m3] (grain_size,nch)
  end type channelSed_struct
  type(channelSed_struct) channelSed

  !for CNP
  !type CNP_struct
     integer                :: ntm    ! number of temperature multiplier factors (4)
     integer                :: nalg   ! number of algae species (3)
     integer                :: nps    ! number of sediment particle size classes
     integer                :: nst    ! number of soil texture classes
     integer                :: nlc    ! number of land cover types (default 27)

     ! --- input paramaters  !#BK
     real(8), allocatable   :: KLIT(:), KRES(:), KEXC(:), KMAN(:), KLDOM(:), KRPOM(:), KRDOM(:)
     real(8), allocatable   :: KMBM(:), KDEN(:), KNIT(:), KVOL(:)
     real(8), allocatable   :: KAGMAX(:), KARMAX(:), KAEMAX(:), KAMMAX(:)
     real(8), allocatable   :: FROOT(:,:)
     real(8), allocatable   :: TMBM(:), MMBM(:)
     real(8), allocatable   :: TZOO(:), MZOO(:)
     real(8), allocatable   :: TALG(:,:), MALG(:,:)
     !real(8), allocatable   :: THETA_S(:), THETA_F(:)  ! (-)        !Y.Kwon   !BK20231102
     real(8), allocatable   :: DWATER(:)                                ! (m)  !BK
     !real(8), allocatable   :: SOLRAD(:,:)                                ! (W m-2)  !BK20240129
     real(8), allocatable   :: DSOIL(:)                                  ! (m)  !BK
     real(8), allocatable   :: TWATER(:)                ! (deg-C)  !Y.Kwon
     real(8), allocatable   :: Qsogw(:,:)   ! = percolate(m,n,4), mm sec-1  !WHQ 20260313
     !integer, allocatable   :: st(:,:)    ! soil texture code !BK20240205

     real(8)                :: FBIO, FMET, FDEC, FREF   ! (FBIO + FMET + FDEC = 1,  0 < FREF < 1)  !BK
     real(8)                :: FMBMC, FMBMN, FMBMP  !#BK
     real(8)                :: FALGC, FALGN, FALGP  !#BK
     real(8)                :: FZOOC, FZOON, FZOOP  !#BK
     !real(8)                :: CNRMBM, CNRALG, CNRZOO, CPRMBM, CPRALG, CPRZOO  !BK20240114
     real(8)                :: KTB, KTC, KTHETA  !#BK
     real(8)                :: RK4SOL            !Y.Kwon
     integer                :: JDAY    !BK20231018
     real(8)                :: KLPOMW, KLDOMW, KRPOMW, KRDOMW, KMBMW, KDENW, KNITW, KVOLW  !BK
     real(8)                :: KADSPO4, KADSPIPA  !#BK
     real(8)                :: SIGMAPO4, SIGMAPIPA  !#BK
     real(8)                :: ALPHA, ALPHADOC, ALPHACHL, ALPHATSS  !#BK
     real(8), dimension(12) :: ALBEDO  !#BK
     real(8)                :: SOLRADMAX, MONOD_N, MONOD_P  !#BK
     real(8)                :: EZI, KZIMAX, ZLOW, ZHALF, KZRMAX, KZMMAX  !#BK
     real(8)                :: OMEGAPOM, OMEGAALG, OMEGAZOO  !#BK
     real(8)                :: TAQUIFER  !#BK
   !   real(8)                :: Q90POMsurf   ! Q90POMdsch  !#BK
   !   real(8)                :: Q90DOMsurf, Q90NH4surf, Q90NO3surf, Q90PO4surf  !#BK
   !   real(8)                :: Q90DOMintf, Q90NH4intf, Q90NO3intf, Q90PO4intf  !#BK
   !   real(8)                :: Q90DOMperc, Q90NH4perc, Q90NO3perc, Q90PO4perc  !#BK
   !   real(8)                :: Q90DOMaqgw, Q90NH4aqgw, Q90NO3aqgw, Q90PO4aqgw  !#BK
   !   real(8)                :: Q90DOMgwch, Q90NH4gwch, Q90NO3gwch, Q90PO4gwch  !#BK
     !BK20260302
     real(8)                :: etaPOMsurf, etaDOMsurf 
     real(8)                :: etaNH4surf, etaNO3surf, etaPO4surf  
     real(8)                :: etaDOMintf, etaNH4intf, etaNO3intf, etaPO4intf  
     real(8)                :: etaDOMperc, etaNH4perc, etaNO3perc, etaPO4perc  
     real(8)                :: etaDOMsogw, etaNH4sogw, etaNO3sogw, etaPO4sogw  
     real(8)                :: etaDOMgwch, etaNH4gwch, etaNO3gwch, etaPO4gwch 
     real(8)                :: etaPOMdsch, etaDOMdsch 
     real(8)                :: etaNH4dsch, etaNO3dsch, etaPO4dsch 
     real(8)                :: Dmix = 0.01  ! effective mixing soil depth (m) !BK20260219

     ! --- SOIL-C
     real(8), allocatable :: So_LPOCLIT(:,:), So_LPOCRES(:,:,:), So_LPOCEXC(:,:), So_LPOCMAN(:,:)
     real(8), allocatable :: So_RPOC(:,:,:), So_LDOC(:,:,:), So_RDOC(:,:,:), So_MBMC(:,:,:), So_CO2(:,:,:)
     real(8), allocatable :: So_DIC(:,:,:) !BK20251215
     real(8), allocatable :: So_LPOCLIT0(:,:), So_LPOCRES0(:,:,:), So_LPOCEXC0(:,:), So_LPOCMAN0(:,:)
     real(8), allocatable :: So_RPOC0(:,:,:), So_LDOC0(:,:,:), So_RDOC0(:,:,:), So_MBMC0(:,:,:)
     real(8), allocatable :: So_DIC0(:,:,:) !BK20251215
     real(8), allocatable :: So_LPOCLITinso(:,:), So_LPOCRESinso(:,:), So_LPOCEXCinso(:,:), So_LPOCMANinso(:,:)
     real(8), allocatable :: So_LPOCLITin(:,:), So_LPOCRESin(:,:,:), So_LPOCEXCin(:,:), So_LPOCMANin(:,:)  !BK20240630
     real(8), allocatable :: So_LDOCperc0(:,:,:), So_RDOCperc0(:,:,:)
     real(8), allocatable :: So_DICperc0(:,:,:) !BK20251215
     real(8), allocatable :: So_LDOCperc(:,:,:), So_RDOCperc(:,:,:)
     real(8), allocatable :: So_DICperc(:,:,:) !BK20251215
     real(8), allocatable :: So_LPOCsurf(:,:), So_RPOCsurf(:,:), So_LDOCsurf(:,:), So_RDOCsurf(:,:), So_MBMCsurf(:,:)
     real(8), allocatable :: So_DICsurf(:,:) !BK20251215
     real(8), allocatable :: So_LPOCsurf0(:,:), So_RPOCsurf0(:,:), So_LDOCsurf0(:,:), So_RDOCsurf0(:,:), So_MBMCsurf0(:,:)
     real(8), allocatable :: So_DICsurf0(:,:) !BK20251215
     real(8), allocatable :: So_LDOCintf(:,:,:), So_RDOCintf(:,:,:)
     real(8), allocatable :: So_DICintf(:,:,:) !BK20251215
     real(8), allocatable :: So_LDOCintf0(:,:,:), So_RDOCintf0(:,:,:)
     real(8), allocatable :: So_DICintf0(:,:,:) !BK20251215
     real(8), allocatable :: So_LDOCsogw(:,:), So_RDOCsogw(:,:) !BK20260316
     real(8), allocatable :: So_DICsogw(:,:) !BK20251215 !BK20260316
     real(8), allocatable :: So_LDOCsogw0(:,:), So_RDOCsogw0(:,:) !BK20260316
     real(8), allocatable :: So_DICsogw0(:,:) !BK20251215 !BK20260316
     real(8), allocatable :: So_LDOCgwso(:,:), So_RDOCgwso(:,:)  !BK20260316
     real(8), allocatable :: So_LDOCgwso0(:,:), So_RDOCgwso0(:,:)  !BK20260316
     real(8), allocatable :: So_CSC(:,:,:), So_CMBError(:,:,:) 

     ! --- SOIL-N
     real(8), allocatable :: So_LPONLIT(:,:), So_LPONRES(:,:,:), So_LPONEXC(:,:), So_LPONMAN(:,:)
     real(8), allocatable :: So_RPON(:,:,:), So_LDON(:,:,:), So_RDON(:,:,:), So_MBMN(:,:,:)
     real(8), allocatable :: So_NH4(:,:,:), So_NO3(:,:,:)
     real(8), allocatable :: So_LPONLIT0(:,:), So_LPONRES0(:,:,:), So_LPONEXC0(:,:), So_LPONMAN0(:,:)
     real(8), allocatable :: So_RPON0(:,:,:), So_LDON0(:,:,:), So_RDON0(:,:,:), So_MBMN0(:,:,:)
     real(8), allocatable :: So_NH40(:,:,:), So_NO30(:,:,:)
     real(8), allocatable :: So_LPONLITinso(:,:), So_LPONRESinso(:,:), So_LPONEXCinso(:,:), So_LPONMANinso(:,:)
     real(8), allocatable :: So_NH4fert(:,:), So_NO3fert(:,:)
     real(8), allocatable :: So_LPONLITin(:,:), So_LPONRESin(:,:,:), So_LPONEXCin(:,:), So_LPONMANin(:,:)  !BK20240630
     real(8), allocatable :: So_NH4fxin(:,:), So_NO3fxin(:,:)  !BK20240630
     real(8), allocatable :: So_NITRI(:,:,:), So_DENIT(:,:,:), So_NH3VOL(:,:,:), So_NH4uptk(:,:,:), So_NO3uptk(:,:,:)
     real(8), allocatable :: So_LDONperc0(:,:,:), So_RDONperc0(:,:,:), So_NH4perc0(:,:,:),So_NO3perc0(:,:,:)
     real(8), allocatable :: So_LDONperc(:,:,:), So_RDONperc(:,:,:), So_NH4perc(:,:,:),So_NO3perc(:,:,:)
     real(8), allocatable :: So_LPONsurf(:,:), So_RPONsurf(:,:), So_LDONsurf(:,:), So_RDONsurf(:,:), So_MBMNsurf(:,:)
     real(8), allocatable :: So_NH4surf(:,:), So_NO3surf(:,:)
     real(8), allocatable :: So_LPONsurf0(:,:), So_RPONsurf0(:,:), So_LDONsurf0(:,:), So_RDONsurf0(:,:), So_MBMNsurf0(:,:)
     real(8), allocatable :: So_NH4surf0(:,:), So_NO3surf0(:,:)
     real(8), allocatable :: So_LDONintf(:,:,:), So_RDONintf(:,:,:), So_NH4intf(:,:,:), So_NO3intf(:,:,:)
     real(8), allocatable :: So_LDONintf0(:,:,:), So_RDONintf0(:,:,:), So_NH4intf0(:,:,:), So_NO3intf0(:,:,:)
     real(8), allocatable :: So_LDONsogw(:,:), So_RDONsogw(:,:), So_NH4sogw(:,:), So_NO3sogw(:,:) !BK20260316
     real(8), allocatable :: So_LDONsogw0(:,:), So_RDONsogw0(:,:), So_NH4sogw0(:,:), So_NO3sogw0(:,:) !BK20260316
     real(8), allocatable :: So_LDONgwso(:,:), So_RDONgwso(:,:), So_NH4gwso(:,:), So_NO3gwso(:,:) !BK20260316
     real(8), allocatable :: So_LDONgwso0(:,:), So_RDONgwso0(:,:), So_NH4gwso0(:,:), So_NO3gwso0(:,:) !BK20260316
     real(8), allocatable :: So_NSC(:,:,:), So_NMBError(:,:,:)

     ! --- SOIL-P
     real(8), allocatable :: So_LPOPLIT(:,:), So_LPOPRES(:,:,:), So_LPOPEXC(:,:), So_LPOPMAN(:,:)
     real(8), allocatable :: So_RPOP(:,:,:), So_LDOP(:,:,:), So_RDOP(:,:,:), So_MBMP(:,:,:)
     real(8), allocatable :: So_PO4(:,:,:), So_PIPA(:,:,:), So_PIPS(:,:,:)
     real(8), allocatable :: So_LPOPLIT0(:,:), So_LPOPRES0(:,:,:), So_LPOPEXC0(:,:), So_LPOPMAN0(:,:)
     real(8), allocatable :: So_RPOP0(:,:,:), So_LDOP0(:,:,:), So_RDOP0(:,:,:), So_MBMP0(:,:,:)
     !real(8), allocatable :: So_PH40(:,:)  !Y.Kwon  !BK20241013 commented out
     real(8), allocatable :: So_PO40(:,:,:), So_PIPA0(:,:,:), So_PIPS0(:,:,:)
     real(8), allocatable :: So_LPOPLITinso(:,:), So_LPOPRESinso(:,:), So_LPOPEXCinso(:,:), So_LPOPMANinso(:,:)
     real(8), allocatable :: So_PO4fert(:,:)
     real(8), allocatable :: So_LPOPLITin(:,:), So_LPOPRESin(:,:,:), So_LPOPEXCin(:,:), So_LPOPMANin(:,:)  !BK20240630
     real(8), allocatable :: So_PO4fxin(:,:)  !BK20240630
     real(8), allocatable :: So_PO4uptk(:,:,:), So_PO4ADS(:,:,:), So_PIPAADS(:,:,:)
     real(8), allocatable :: So_LDOPperc0(:,:,:), So_RDOPperc0(:,:,:), So_PO4perc0(:,:,:)
     real(8), allocatable :: So_LDOPperc(:,:,:), So_RDOPperc(:,:,:), So_PO4perc(:,:,:)
     real(8), allocatable :: So_LPOPsurf(:,:), So_RPOPsurf(:,:), So_LDOPsurf(:,:), So_RDOPsurf(:,:), So_MBMPsurf(:,:)
     real(8), allocatable :: So_PO4surf(:,:), So_PIPAsurf(:,:), So_PIPSsurf(:,:)
     real(8), allocatable :: So_LPOPsurf0(:,:), So_RPOPsurf0(:,:), So_LDOPsurf0(:,:), So_RDOPsurf0(:,:), So_MBMPsurf0(:,:)
     real(8), allocatable :: So_PO4surf0(:,:), So_PIPAsurf0(:,:), So_PIPSsurf0(:,:)
     real(8), allocatable :: So_LDOPintf(:,:,:), So_RDOPintf(:,:,:), So_PO4intf(:,:,:)
     real(8), allocatable :: So_LDOPintf0(:,:,:), So_RDOPintf0(:,:,:), So_PO4intf0(:,:,:)
     real(8), allocatable :: So_LDOPsogw(:,:), So_RDOPsogw(:,:), So_PO4sogw(:,:)  !BK20260316
     real(8), allocatable :: So_LDOPsogw0(:,:), So_RDOPsogw0(:,:), So_PO4sogw0(:,:) !BK20260316
     real(8), allocatable :: So_LDOPgwso(:,:), So_RDOPgwso(:,:), So_PO4gwso(:,:)  !BK20260316
     real(8), allocatable :: So_LDOPgwso0(:,:), So_RDOPgwso0(:,:), So_PO4gwso0(:,:) !BK20260316
     real(8), allocatable :: So_PSC(:,:,:), So_PMBError(:,:,:)

   !   ! --- AQUIFER-C  !#BK
   !   real(8), allocatable :: Aq_LDOC(:,:), Aq_RDOC(:,:)
   !   real(8), allocatable :: Aq_LDOC0(:,:), Aq_RDOC0(:,:)
   !   real(8), allocatable :: Aq_LDOCaqso(:,:), Aq_RDOCaqso(:,:)
   !   real(8), allocatable :: Aq_LDOCaqso0(:,:), Aq_RDOCaqso0(:,:)
   !   real(8), allocatable :: Aq_LDOCaqgw(:,:), Aq_RDOCaqgw(:,:)
   !   real(8), allocatable :: Aq_LDOCaqgw0(:,:), Aq_RDOCaqgw0(:,:)
   !   real(8), allocatable :: Aq_CSC(:,:), Aq_CMBError(:,:)

   !   ! --- AQUIFER-N  !#BK
   !   real(8), allocatable :: Aq_LDON(:,:), Aq_RDON(:,:), Aq_NH4(:,:), Aq_NO3(:,:)
   !   real(8), allocatable :: Aq_LDON0(:,:), Aq_RDON0(:,:), Aq_NH40(:,:), Aq_NO30(:,:)
   !   real(8), allocatable :: Aq_LDONaqso(:,:), Aq_RDONaqso(:,:), Aq_NH4aqso(:,:), Aq_NO3aqso(:,:)
   !   real(8), allocatable :: Aq_LDONaqso0(:,:), Aq_RDONaqso0(:,:), Aq_NH4aqso0(:,:), Aq_NO3aqso0(:,:)
   !   real(8), allocatable :: Aq_LDONaqgw(:,:), Aq_RDONaqgw(:,:), Aq_NH4aqgw(:,:), Aq_NO3aqgw(:,:)
   !   real(8), allocatable :: Aq_LDONaqgw0(:,:), Aq_RDONaqgw0(:,:), Aq_NH4aqgw0(:,:), Aq_NO3aqgw0(:,:)
   !   real(8), allocatable :: Aq_NSC(:,:), Aq_NMBError(:,:)

   !   ! --- AQUIFER-P  !#BK
   !   real(8), allocatable :: Aq_LDOP(:,:), Aq_RDOP(:,:), Aq_PO4(:,:), Aq_PIPA(:,:), Aq_PIPS(:,:)
   !   real(8), allocatable :: Aq_LDOP0(:,:), Aq_RDOP0(:,:), Aq_PO40(:,:), Aq_PIPA0(:,:), Aq_PIPS0(:,:)
   !   real(8), allocatable :: Aq_PO4ADS(:,:), Aq_PIPAADS(:,:)  !BK20241216
   !   real(8), allocatable :: Aq_LDOPaqso(:,:), Aq_RDOPaqso(:,:), Aq_PO4aqso(:,:)
   !   real(8), allocatable :: Aq_LDOPaqso0(:,:), Aq_RDOPaqso0(:,:), Aq_PO4aqso0(:,:)
   !   real(8), allocatable :: Aq_LDOPaqgw(:,:), Aq_RDOPaqgw(:,:), Aq_PO4aqgw(:,:)
   !   real(8), allocatable :: Aq_LDOPaqgw0(:,:), Aq_RDOPaqgw0(:,:), Aq_PO4aqgw0(:,:)
   !   real(8), allocatable :: Aq_PSC(:,:), Aq_PMBError(:,:)

     ! --- GROUNDWATER-C
     real(8), allocatable :: Gw_LDOC(:), Gw_RDOC(:)
     real(8), allocatable :: Gw_LDOC0(:), Gw_RDOC0(:)
     real(8), allocatable :: Gw_LDOCsogw0(:), Gw_RDOCsogw0(:)  !BK20240513 !BK20240617 !BK20250310 !BK20260316
     real(8), allocatable :: Gw_LDOCgwso(:), Gw_RDOCgwso(:)    !BK20260320
     real(8), allocatable :: Gw_LDOCgwch(:), Gw_RDOCgwch(:)
     real(8), allocatable :: Gw_LDOCgwch0(:), Gw_RDOCgwch0(:)
     real(8), allocatable :: Gw_CSC(:), Gw_CMBError(:)

     ! --- GROUNDWATER-N  !#BK
     real(8), allocatable :: Gw_LDON(:), Gw_RDON(:), Gw_NH4(:), Gw_NO3(:)
     real(8), allocatable :: Gw_LDON0(:), Gw_RDON0(:), Gw_NH40(:), Gw_NO30(:)
     real(8), allocatable :: Gw_LDONsogw0(:), Gw_RDONsogw0(:), Gw_NH4sogw0(:), Gw_NO3sogw0(:)  !BK20240513 !BK20240617 !BK20250310 !BK20260316
     real(8), allocatable :: Gw_LDONgwso(:), Gw_RDONgwso(:), Gw_NH4gwso(:), Gw_NO3gwso(:)    !BK20260320
     real(8), allocatable :: Gw_LDONgwch(:), Gw_RDONgwch(:), Gw_NH4gwch(:), Gw_NO3gwch(:)
     real(8), allocatable :: Gw_LDONgwch0(:), Gw_RDONgwch0(:), Gw_NH4gwch0(:), Gw_NO3gwch0(:)
     real(8), allocatable :: Gw_NSC(:), Gw_NMBError(:)

     ! --- GROUNDWATER-P  !#BK  !BK20250809
     real(8), allocatable :: Gw_LDOP(:), Gw_RDOP(:), Gw_PO4(:)  
     real(8), allocatable :: Gw_LDOP0(:), Gw_RDOP0(:), Gw_PO40(:)
     real(8), allocatable :: Gw_LDOPsogw0(:), Gw_RDOPsogw0(:), Gw_PO4sogw0(:)  !BK20240513 !BK20240617 !BK20250310 !BK20260316
     real(8), allocatable :: Gw_LDOPgwso(:), Gw_RDOPgwso(:), Gw_PO4gwso(:)    !BK20260320
     real(8), allocatable :: Gw_LDOPgwch(:), Gw_RDOPgwch(:), Gw_PO4gwch(:)
     real(8), allocatable :: Gw_LDOPgwch0(:), Gw_RDOPgwch0(:), Gw_PO4gwch0(:)
     real(8), allocatable :: Gw_PSC(:), Gw_PMBError(:)

     ! --- CHANNEL-COMMON  
     !integer, allocatable :: abs2dis_ich(:,:)  !BK20240710 !BK20241018 !BK20250705
     !integer, allocatable :: abs2dis_lag(:,:)  !BK20240710 !BK20241018 !BK20250705
     
     ! --- CHANNEL-C  !#BK
     real(8), allocatable :: Ch_ALGC(:,:), Ch_ZOOC(:), Ch_MBMC(:), Ch_DIC(:)
     real(8), allocatable :: Ch_CO2EVAS(:)  !BK20251126
     real(8), allocatable :: Ch_LPOC(:), Ch_RPOC(:), Ch_LDOC(:), Ch_RDOC(:), Ch_POCDEPOSIT(:)
     real(8), allocatable :: Ch_ALGC0(:,:), Ch_ZOOC0(:), Ch_MBMC0(:), Ch_DIC0(:)  !BK20251120
     real(8), allocatable :: Ch_LPOC0(:), Ch_RPOC0(:), Ch_LDOC0(:), Ch_RDOC0(:), Ch_POCDEPOSIT0(:)
     real(8), allocatable :: Ch_ALGCusch0(:,:), Ch_ZOOCusch0(:), Ch_MBMCusch0(:), Ch_DICusch0(:)  !BK20251120
     real(8), allocatable :: Ch_LPOCusch0(:), Ch_RPOCusch0(:), Ch_LDOCusch0(:), Ch_RDOCusch0(:)
     real(8), allocatable :: Ch_LPOCsurf0(:), Ch_RPOCsurf0(:), Ch_LDOCsurf0(:), Ch_RDOCsurf0(:), Ch_MBMCsurf0(:)  !BK20240513 !BK20240617 !BK20250310
     real(8), allocatable :: Ch_LDOCintf0(:), Ch_RDOCintf0(:)  !BK20240513 !BK20240617 !BK20250310
     real(8), allocatable :: Ch_LDOCgwch0(:), Ch_RDOCgwch0(:)  !BK20250310
     real(8), allocatable :: Ch_ALGCG(:,:), Ch_ALGCR(:,:), Ch_ALGCE(:,:), Ch_ALGCM(:,:), Ch_ALGCZ(:,:)
     real(8), allocatable :: Ch_ALGCS(:,:)
     real(8), allocatable :: Ch_ZOOCG(:), Ch_ZOOCR(:), Ch_ZOOCE(:), Ch_ZOOCM(:), Ch_ZOOCZ(:), Ch_ZOOCS(:)
     real(8), allocatable :: Ch_ALGCdsch(:,:), Ch_ZOOCdsch(:), Ch_MBMCdsch(:), Ch_DICdsch(:)  !BK20251120
     real(8), allocatable :: Ch_LPOCdsch(:), Ch_RPOCdsch(:), Ch_LDOCdsch(:), Ch_RDOCdsch(:)
     real(8), allocatable :: Ch_LPOCpnt(:,:), Ch_RPOCpnt(:,:), Ch_LDOCpnt(:,:), Ch_RDOCpnt(:,:)   ![kg s-1]   !BK20231016
     real(8), allocatable :: Ch_LPOCxpnt(:), Ch_RPOCxpnt(:), Ch_LDOCxpnt(:), Ch_RDOCxpnt(:)           ![kg]       !BK20240701   
     real(8), allocatable :: Ch_LPOCabs(:,:), Ch_RPOCabs(:,:), Ch_LDOCabs(:,:), Ch_RDOCabs(:,:)       ![kg s-1]   !BK20240703
     real(8), allocatable :: Ch_LPOCxabs(:), Ch_RPOCxabs(:), Ch_LDOCxabs(:), Ch_RDOCxabs(:)           ![kg]       !BK20240729   
     real(8), allocatable :: Ch_LPOCdis(:,:), Ch_RPOCdis(:,:), Ch_LDOCdis(:,:), Ch_RDOCdis(:,:)       ![kg s-1]   !BK20240703   
     real(8), allocatable :: Ch_LPOCxdis(:), Ch_RPOCxdis(:), Ch_LDOCxdis(:), Ch_RDOCxdis(:)           ![kg]       !BK20240729   
     !real(8), allocatable :: Ch_LPOCsewr(:), Ch_RPOCsewr(:), Ch_LDOCsewr(:), Ch_RDOCsewr(:)
     real(8), allocatable :: Ch_CSC(:), Ch_CMBError(:)
     real(8), allocatable :: KAG(:,:), KAR(:,:), KAE(:,:), KAM(:,:)            !BK20231120  !BK20240115
     real(8), allocatable :: KZG(:), KZR(:), KZM(:) !, KZE(:)   !Y.Kwon  !BK20231120 !BK20251127

     ! --- CHANNEL-N
     real(8), allocatable :: Ch_ALGN(:,:), Ch_ZOON(:), Ch_MBMN(:), Ch_NH4(:), Ch_NO3(:)
     real(8), allocatable :: Ch_LPON(:), Ch_RPON(:), Ch_LDON(:), Ch_RDON(:), Ch_PONDEPOSIT(:)
     real(8), allocatable :: Ch_ALGN0(:,:), Ch_ZOON0(:), Ch_MBMN0(:), Ch_NH40(:), Ch_NO30(:)
     real(8), allocatable :: Ch_LPON0(:), Ch_RPON0(:), Ch_LDON0(:), Ch_RDON0(:), Ch_PONDEPOSIT0(:)
     real(8), allocatable :: Ch_NITRI(:), Ch_DENIT(:), Ch_NH3VOL(:), Ch_NH4uptk(:), Ch_NO3uptk(:)
     real(8), allocatable :: Ch_ALGNusch0(:,:), Ch_ZOONusch0(:), Ch_MBMNusch0(:), Ch_NH4usch0(:), Ch_NO3usch0(:)
     real(8), allocatable :: Ch_LPONusch0(:), Ch_RPONusch0(:), Ch_LDONusch0(:), Ch_RDONusch0(:)
     real(8), allocatable :: Ch_LPONsurf0(:), Ch_RPONsurf0(:), Ch_LDONsurf0(:), Ch_RDONsurf0(:)  !BK20240513 !BK20240617 !BK20250310
     real(8), allocatable :: Ch_MBMNsurf0(:), Ch_NH4surf0(:), Ch_NO3surf0(:)  !BK20240513 !BK20240617 !BK20250310
     real(8), allocatable :: Ch_LDONintf0(:), Ch_RDONintf0(:), Ch_NH4intf0(:), Ch_NO3intf0(:)  !BK20240513 !BK20240617 !BK20250310
     real(8), allocatable :: Ch_LDONgwch0(:), Ch_RDONgwch0(:), Ch_NH4gwch0(:), Ch_NO3gwch0(:)  !BK20250310
     real(8), allocatable :: Ch_ALGNG(:,:), Ch_ALGNR(:,:), Ch_ALGNE(:,:), Ch_ALGNM(:,:), Ch_ALGNZ(:,:)
     real(8), allocatable :: Ch_ALGNS(:,:)
     real(8), allocatable :: Ch_ZOONG(:), Ch_ZOONR(:), Ch_ZOONE(:), Ch_ZOONM(:), Ch_ZOONZ(:), Ch_ZOONS(:)
     real(8), allocatable :: Ch_ALGNdsch(:,:), Ch_ZOONdsch(:), Ch_MBMNdsch(:), Ch_NH4dsch(:), Ch_NO3dsch(:)
     real(8), allocatable :: Ch_LPONdsch(:), Ch_RPONdsch(:), Ch_LDONdsch(:), Ch_RDONdsch(:)
     real(8), allocatable :: Ch_LPONpnt(:,:), Ch_RPONpnt(:,:), Ch_LDONpnt(:,:), Ch_RDONpnt(:,:), &
                             Ch_NH4pnt(:,:), Ch_NO3pnt(:,:)  ![kg s-1]  !BK20231016 
     real(8), allocatable :: Ch_LPONxpnt(:), Ch_RPONxpnt(:), Ch_LDONxpnt(:), Ch_RDONxpnt(:), &
                             Ch_NH4xpnt(:), Ch_NO3xpnt(:)      ![kg]      !BK20240701 
     real(8), allocatable :: Ch_LPONabs(:,:), Ch_RPONabs(:,:), Ch_LDONabs(:,:), Ch_RDONabs(:,:), &
                             Ch_NH4abs(:,:), Ch_NO3abs(:,:)    ![kg s-1]  !BK20240703 
     real(8), allocatable :: Ch_LPONxabs(:), Ch_RPONxabs(:), Ch_LDONxabs(:), Ch_RDONxabs(:), &
                             Ch_NH4xabs(:), Ch_NO3xabs(:)      ![kg]      !BK20240729 
     real(8), allocatable :: Ch_LPONdis(:,:), Ch_RPONdis(:,:), Ch_LDONdis(:,:), Ch_RDONdis(:,:), &
                             Ch_NH4dis(:,:), Ch_NO3dis(:,:)    ![kg s-1]  !BK20240703 
     real(8), allocatable :: Ch_LPONxdis(:), Ch_RPONxdis(:), Ch_LDONxdis(:), Ch_RDONxdis(:), &
                             Ch_NH4xdis(:), Ch_NO3xdis(:)      ![kg]      !BK20240729 
     !real(8), allocatable :: Ch_LPONsewr(:), Ch_RPONsewr(:), Ch_LDONsewr(:), Ch_RDONsewr(:), Ch_NH4sewr(:), Ch_NO3sewr(:)
     real(8), allocatable :: Ch_NSC(:), Ch_NMBError(:)

     ! --- CHANNEL-P
     real(8), allocatable :: Ch_ALGP(:,:), Ch_ZOOP(:), Ch_MBMP(:)
     real(8), allocatable :: Ch_PO4(:), Ch_PIPA(:), Ch_PIPS(:)
     real(8), allocatable :: Ch_LPOP(:), Ch_RPOP(:), Ch_LDOP(:), Ch_RDOP(:), Ch_POPDEPOSIT(:)
     real(8), allocatable :: Ch_ALGP0(:,:), Ch_ZOOP0(:), Ch_MBMP0(:)
     real(8), allocatable :: Ch_PO40(:), Ch_PIPA0(:), Ch_PIPS0(:)
     real(8), allocatable :: Ch_PO4uptk(:), Ch_PO4ADS(:), Ch_PIPAADS(:)
     real(8), allocatable :: Ch_LPOP0(:), Ch_RPOP0(:), Ch_LDOP0(:), Ch_RDOP0(:), Ch_POPDEPOSIT0(:)
     real(8), allocatable :: Ch_ALGPusch0(:,:), Ch_ZOOPusch0(:), Ch_MBMPusch0(:)
     real(8), allocatable :: Ch_PO4usch0(:), Ch_PIPAusch0(:), Ch_PIPSusch0(:)
     real(8), allocatable :: Ch_LPOPusch0(:), Ch_RPOPusch0(:), Ch_LDOPusch0(:), Ch_RDOPusch0(:)
     real(8), allocatable :: Ch_LPOPsurf0(:), Ch_RPOPsurf0(:), Ch_LDOPsurf0(:), Ch_RDOPsurf0(:) !BK20240617 !BK20250310
     real(8), allocatable :: Ch_MBMPsurf0(:), Ch_PO4surf0(:), Ch_PIPAsurf0(:), Ch_PIPSsurf0(:) !BK20240617 !BK20250310
     real(8), allocatable :: Ch_LDOPintf0(:), Ch_RDOPintf0(:), Ch_PO4intf0(:) !BK20240617 !BK20250310
     real(8), allocatable :: Ch_LDOPgwch0(:), Ch_RDOPgwch0(:), Ch_PO4gwch0(:)   !BK20250310
     real(8), allocatable :: Ch_ALGPG(:,:), Ch_ALGPR(:,:), Ch_ALGPE(:,:), Ch_ALGPM(:,:)
     real(8), allocatable :: Ch_ALGPZ(:,:), Ch_ALGPS(:,:)
     real(8), allocatable :: Ch_ZOOPG(:), Ch_ZOOPR(:), Ch_ZOOPE(:), Ch_ZOOPM(:), Ch_ZOOPZ(:), Ch_ZOOPS(:)
     real(8), allocatable :: Ch_ALGPdsch(:,:), Ch_ZOOPdsch(:), Ch_MBMPdsch(:)
     real(8), allocatable :: Ch_PO4dsch(:), Ch_PIPAdsch(:), Ch_PIPSdsch(:)
     real(8), allocatable :: Ch_LPOPdsch(:), Ch_RPOPdsch(:), Ch_LDOPdsch(:), Ch_RDOPdsch(:)
     real(8), allocatable :: Ch_LPOPpnt(:,:), Ch_RPOPpnt(:,:), Ch_LDOPpnt(:,:), Ch_RDOPpnt(:,:), &
                             Ch_PO4pnt(:,:)  ![kg s-1]  !BK20231016 
     real(8), allocatable :: Ch_LPOPxpnt(:), Ch_RPOPxpnt(:), Ch_LDOPxpnt(:), Ch_RDOPxpnt(:), &
                             Ch_PO4xpnt(:)    ![kg]      !BK20240701 
     real(8), allocatable :: Ch_LPOPabs(:,:), Ch_RPOPabs(:,:), Ch_LDOPabs(:,:), Ch_RDOPabs(:,:), &
                             Ch_PO4abs(:,:)   ![kg s-1]  !BK20240703 
     real(8), allocatable :: Ch_LPOPxabs(:), Ch_RPOPxabs(:), Ch_LDOPxabs(:), Ch_RDOPxabs(:), &
                             Ch_PO4xabs(:)    ![kg]      !BK20240729 
     real(8), allocatable :: Ch_LPOPdis(:,:), Ch_RPOPdis(:,:), Ch_LDOPdis(:,:), Ch_RDOPdis(:,:), &
                             Ch_PO4dis(:,:)   ![kg s-1]  !BK20240703 
     real(8), allocatable :: Ch_LPOPxdis(:), Ch_RPOPxdis(:), Ch_LDOPxdis(:), Ch_RDOPxdis(:), &
                             Ch_PO4xdis(:)    ![kg]      !BK20240729 
     !real(8), allocatable :: Ch_LPOPsewr(:), Ch_RPOPsewr(:), Ch_LDOPsewr(:), Ch_RDOPsewr(:), Ch_PO4sewr(:)
     real(8), allocatable :: Ch_PSC(:), Ch_PMBError(:)
     real(8), allocatable :: SSA(:)     ! specific surface area per unit mass of sediment, m2 kg-1   !Y.Kwon

  !end type CNP_struct
  !type(CNP_struct) cnp



end module module_SedCNPvariables
