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

module module_SedCNPmodel_driver

   use module_hydro_stop, only: hydro_stop
   contains

   subroutine SedCNP_driver_ini(NTIME_out)

      use module_SedCNPvariables
!=====||__WHQ5403__||=====!
!
      ! SedCNPmodel already visible via ReadCNPini/module_SedCNP_in below; a redundant
      ! direct "use config_base" here triggers a gfortran diamond-import type mismatch.
      !use config_base, only: SedCNPmodel
!
!=====||__WHQ5403__||=====!
      use CNPparams
      !use CNPvariables
      use ReadCNPini
      use module_SedCNP_in  !BK20240205

      implicit none

      integer, intent(out) :: NTIME_out
      character(len=256)   :: geo_finegrid_flnm
      !logical              :: restart_flag  !BK20251209 commented out - not being used
      real                 :: dts    !timestep
      integer              :: khour
      integer              :: N_TIME
      integer              :: istat
      character(len=256)   :: filename_nch, filename_nbasin   !Y.Kwon(20250621)
      character(len=256)   :: filename_ixjxnsl                !Y.Kwon(20250621)
      character(len=256)   :: filename_geo                    !Y.Kwon(20250621)
      

      !get simulation period and timestep
      dts = real(SedCNPmodel%SedCNP_timestep)
      khour = SedCNPmodel%SedCNP_khour

      if ((khour < 0) .and. (SedCNPmodel%SedCNP_kday < 0)) then
         write(*, '("FATAL ERROR: In module_SedCNPmodel_driver SedCNP_driver_ini() - "// &
               "Namelist error: Either KHOUR or KDAY must be defined.")')
         call hydro_stop("FATAL ERROR: In SedCNP_driver_ini() - KHOUR or KDAY must be defined.")
      else if (( khour < 0 ) .and. (SedCNPmodel%SedCNP_kday > 0)) then
         khour = SedCNPmodel%SedCNP_kday * 24
      else if ((khour > 0) .and. (SedCNPmodel%SedCNP_kday > 0)) then
         write(*, '("WARNING: In SedCNP_driver_ini() - KHOUR and KDAY both defined. Using KHOUR.")')
      endif
      N_TIME = khour*3600./nint(dts)
      NTIME_out = N_TIME
      
      domain%ntime_sedcnp = NTIME_out  !Y.Kwon 20230316 (need to be fixed)
      domain%DT = dts  !BK20240711

      !--------------------------Y.Kwon(20250621)
      !ix, jx
      filename_ixjxnsl = "./DOMAIN/soil_properties.nc"
      call get_ixjxnsl(domain%ix,domain%jx,domain%nsl,trim(filename_ixjxnsl))
      domain%ixrt = domain%ix
      domain%jxrt = domain%jx

      !ISWATER
      call get_global_iswater(trim(filename_ixjxnsl), domain%ISWATER) !BK20251208
      !If failed (still default -999), default to 16 and warn
      if (domain%ISWATER == -999) then
         print*, "WARNING: ISWATER not found in input file, using default 16"
         domain%ISWATER = 16
      endif

      !nch
      filename_nch = "./DOMAIN/Route_Link.nc"
      call get_feature_id(domain%nch,trim(filename_nch))

      !nbasin
      filename_nbasin = "./DOMAIN/GWBUCKPARM.nc"
      call get_feature_id(domain%nbasin,trim(filename_nbasin))
      !--------------------------Y.Kwon(20250621)
   
      ntm  = 4    ! number of temperature multiplier factors (4)  !Y.Kwon
      !nalg = 3    ! number of algae species !Y.Kwon
      nps  = 4    ! number of sediment particle size classes
      nst  = 19   ! number of soil texture classes               *** HARDCODED ***  !BK20240128
      nlc  = 27   ! number of land cover types (default 27)      *** HARDCODED ***
      nalg = 3    ! number of algae species (default 3)          *** HARDCODED ***
      !geo_finegrid_flnm = "./DOMAIN/Fulldom_hires.nc"
      !restart_flag = .false.  !BK20251209

      !get soil property file name and directory
      !SedCNPspatial_filename = trim(SedCNPmodel%SedCNP_SPATIAL_FILENAME_

      !read land surface inputs required for SedCNPmodel
      !call SedCNPmodel_surface_input(SedCNPspatial_filename)

      !domain
      !allocate(domain%linkID_grid(domain%ixrt,domain%jxrt))  !BK20240619

      !channel input
      allocate(channel_input%Wst(domain%nch))
      !allocate(channel_input%TopWdth(domain%nch))     !BK20251117
      !allocate(channel_input%TopWdthCC(domain%nch))   !BK20251117
      allocate(channel_input%SideSlopch(domain%nch))
      allocate(channel_input%Lst(domain%nch))
      allocate(channel_input%Sf(domain%nch))

      allocate(SedCNP_hydro%subbasinID(domain%ixrt,domain%jxrt)) !Y.Kwon 20230422
      allocate(SedCNP_hydro%smcmax(domain%ix,domain%jx,4))
      allocate(SedCNP_hydro%smcref(domain%ix,domain%jx,4))
      allocate(SedCNP_hydro%linkID_grid(domain%ixrt,domain%jxrt))  !BK20240619
      allocate(SedCNP_hydro%gwbasin(domain%ix,domain%jx))         !BK20231011  !BK20240619
      allocate(SedCNP_hydro%gwbasin_ch(domain%ixrt,domain%jxrt))  !BK20231011  !BK20240619

      !mapping arrays
      allocate(domain%gwid_ch(domain%nch), stat=istat)    !BK20250701
      allocate(domain%chid_i(domain%nch), stat=istat)     !BK20250701
      allocate(domain%chid_j(domain%nch), stat=istat)     !BK20250701 

      !overland
      allocate(overSed%Sol(4,domain%ix,domain%jx))
      !allocate(overSed%Sol0(4,domain%ix,domain%jx))  !BK20240619  !BK20250925
      allocate(overSed%Sow(4,domain%ix,domain%jx))
      !allocate(overSed%Sow0(4,domain%ix,domain%jx))  !BK20240619  !BK20250925
      allocate(overSed%Solcp(4,domain%ix,domain%jx))
      allocate(overSed%Ssurf(4,domain%ix,domain%jx))
      allocate(overSed%Ssurf0(4,domain%ix,domain%jx))  !BK20240619
      allocate(overSed%Srdet(4,domain%ix,domain%jx)) 
      allocate(overSed%Soler(4,domain%ix,domain%jx))
      allocate(overSed%Soldp(4,domain%ix,domain%jx))
      allocate(overSed%Colsd(4,domain%ix,domain%jx))
      allocate(overSed%SOILPSF(nst,nps))                !BK20240205
      
      !channel
      allocate(channelSed%Sch(4,domain%nch))
      allocate(channelSed%Sch0(4,domain%nch))  !BK20240619
      allocate(channelSed%Scw(4,domain%nch))
      allocate(channelSed%Scw0(4,domain%nch))  !BK20240619
      !allocate(channelSed%Susch(4,domain%nch))  !BK20250925
      allocate(channelSed%Susch0(4,domain%nch))
      allocate(channelSed%Spnt(4,domain%nch,domain%ntime_sedcnp))
      allocate(channelSed%Sxpnt(4,domain%nch))   !BK20240701
      allocate(channelSed%Sabs(4,domain%nch,domain%ntime_sedcnp)) !BK20240710
      allocate(channelSed%Sxabs(4,domain%nch))   !BK20240729
      allocate(channelSed%Sdis(4,domain%nch,domain%ntime_sedcnp)) !BK20240710
      allocate(channelSed%Sxdis(4,domain%nch))   !BK20240729
      !allocate(channelSed%Solch(4,domain%nch))  
      allocate(channelSed%Solch0(4,domain%nch))  !BK20250925
      !allocate(channelSed%Sirri(4,domain%nch))  !BK20260510 removed Sirri
      allocate(channelSed%Sinput(4,domain%nch))
      allocate(channelSed%Sdsch(4,domain%nch)) 
      allocate(channelSed%Scher(4,domain%nch))
      allocate(channelSed%Schdp(4,domain%nch))
      allocate(channelSed%Schcp(4,domain%nch))
      allocate(channelSed%Cchsd(4,domain%nch))

      !CNP
      allocate(KLIT(nlc), KRES(nlc), KEXC(nlc), &
               KMAN(nlc), KLDOM(nlc), KRPOM(nlc), KRDOM(nlc), stat=istat)
      allocate(KMBM(nlc), KDEN(nlc), KNIT(nlc), KVOL(nlc), stat=istat)
      allocate(KAGMAX(nalg), KARMAX(nalg), KAEMAX(nalg), KAMMAX(nalg), &
            stat=istat)
      allocate(FROOT(nlc,domain%nsl),stat=istat)
      allocate(TMBM(ntm), MMBM(ntm),stat=istat)
      allocate(TZOO(ntm), MZOO(ntm),stat=istat)
      allocate(TALG(nalg,ntm), MALG(nalg,ntm),stat=istat)

      ! --- SOIL-C
      allocate(So_LPOCLIT0(domain%ix,domain%jx), &
               So_LPOCRES0(domain%ix,domain%nsl,domain%jx), &
               So_LPOCEXC0(domain%ix,domain%jx), &
               So_LPOCMAN0(domain%ix,domain%jx),stat=istat)
      allocate(So_RPOC0(domain%ix,domain%nsl,domain%jx), &
               So_LDOC0(domain%ix,domain%nsl,domain%jx), &
               So_RDOC0(domain%ix,domain%nsl,domain%jx), &
               So_MBMC0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LPOCsurf0(domain%ix,domain%jx), &
               So_RPOCsurf0(domain%ix,domain%jx), &
               So_LDOCsurf0(domain%ix,domain%jx), &
               So_RDOCsurf0(domain%ix,domain%jx), &
               So_MBMCsurf0(domain%ix,domain%jx),stat=istat)
      allocate(So_LDOCintf0(domain%ix,domain%nsl,domain%jx), &
               So_RDOCintf0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LDOCperc0(domain%ix,domain%nsl,domain%jx), &
               So_RDOCperc0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LDOCgwso0(domain%ix,domain%jx), &
               So_RDOCgwso0(domain%ix,domain%jx),stat=istat)  !BK20260322
      allocate(So_LDOCsogw0(domain%ix,domain%jx), &
               So_RDOCsogw0(domain%ix,domain%jx),stat=istat)  !BK20260318

      ! --- SOIL-N
      allocate(So_LPONLIT0(domain%ix,domain%jx), &
               So_LPONRES0(domain%ix,domain%nsl,domain%jx), &
               So_LPONEXC0(domain%ix,domain%jx), &
               So_LPONMAN0(domain%ix,domain%jx),stat=istat)
      allocate(So_RPON0(domain%ix,domain%nsl,domain%jx), &
               So_LDON0(domain%ix,domain%nsl,domain%jx), &
               So_RDON0(domain%ix,domain%nsl,domain%jx), &
               So_MBMN0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_NH40(domain%ix,domain%nsl,domain%jx), &
               So_NO30(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LDONperc0(domain%ix,domain%nsl,domain%jx), &
               So_RDONperc0(domain%ix,domain%nsl,domain%jx), &
               So_NH4perc0(domain%ix,domain%nsl,domain%jx), &
               So_NO3perc0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LPONsurf0(domain%ix,domain%jx), &
               So_RPONsurf0(domain%ix,domain%jx), &
               So_LDONsurf0(domain%ix,domain%jx), &
               So_RDONsurf0(domain%ix,domain%jx), &
               So_MBMNsurf0(domain%ix,domain%jx),stat=istat)
      allocate(So_NH4surf0(domain%ix,domain%jx), &
               So_NO3surf0(domain%ix,domain%jx),stat=istat)
      allocate(So_LDONintf0(domain%ix,domain%nsl,domain%jx), &
               So_RDONintf0(domain%ix,domain%nsl,domain%jx),&
               So_NH4intf0(domain%ix,domain%nsl,domain%jx), &
               So_NO3intf0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LDONgwso0(domain%ix,domain%jx), &
               So_RDONgwso0(domain%ix,domain%jx), &
               So_NH4gwso0(domain%ix,domain%jx), &
               So_NO3gwso0(domain%ix,domain%jx),stat=istat)    !BK20260322
      allocate(So_LDONsogw0(domain%ix,domain%jx), &
               So_RDONsogw0(domain%ix,domain%jx), &    !BK20260318
               So_NH4sogw0(domain%ix,domain%jx), &
               So_NO3sogw0(domain%ix,domain%jx),stat=istat)    !BK20260318

      ! --- SOIL-P
      allocate(So_LPOPLIT0(domain%ix,domain%jx), &
               So_LPOPRES0(domain%ix,domain%nsl,domain%jx), &
               So_LPOPEXC0(domain%ix,domain%jx), &
               So_LPOPMAN0(domain%ix,domain%jx),stat=istat)
      allocate(So_RPOP0(domain%ix,domain%nsl,domain%jx), &
               So_LDOP0(domain%ix,domain%nsl,domain%jx), &
               So_RDOP0(domain%ix,domain%nsl,domain%jx), &
               So_MBMP0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_PO40(domain%ix,domain%nsl,domain%jx), &
               So_PIPA0(domain%ix,domain%nsl,domain%jx), &
               So_PIPS0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LDOPperc0(domain%ix,domain%nsl,domain%jx), &
               So_RDOPperc0(domain%ix,domain%nsl,domain%jx), &
               So_PO4perc0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LPOPsurf0(domain%ix,domain%jx), &
               So_RPOPsurf0(domain%ix,domain%jx), &
               So_LDOPsurf0(domain%ix,domain%jx), &
               So_RDOPsurf0(domain%ix,domain%jx), &
               So_MBMPsurf0(domain%ix,domain%jx),stat=istat)
      allocate(So_PO4surf0(domain%ix,domain%jx), &
               So_PIPAsurf0(domain%ix,domain%jx), &
               So_PIPSsurf0(domain%ix,domain%jx),stat=istat)
      allocate(So_LDOPintf0(domain%ix,domain%nsl,domain%jx), &
               So_RDOPintf0(domain%ix,domain%nsl,domain%jx), &
               So_PO4intf0(domain%ix,domain%nsl,domain%jx),stat=istat)
      allocate(So_LDOPgwso0(domain%ix,domain%jx), &
               So_RDOPgwso0(domain%ix,domain%jx), &
               So_PO4gwso0(domain%ix,domain%jx),stat=istat)  !BK20260322
      allocate(So_LDOPsogw0(domain%ix,domain%jx), &
               So_RDOPsogw0(domain%ix,domain%jx), &
               So_PO4sogw0(domain%ix,domain%jx),stat=istat)  !BK20260318

      ! --- GROUNDWATER-C
      allocate(Gw_LDOC0(domain%nbasin), &
               Gw_RDOC0(domain%nbasin),stat=istat)
      allocate(Gw_LDOCsogw0(domain%nbasin), &
               Gw_RDOCsogw0(domain%nbasin),stat=istat)  !BK20240610  !BK20250310 !BK20260318
      allocate(Gw_LDOCgwch0(domain%nbasin), &
               Gw_RDOCgwch0(domain%nbasin),stat=istat)  !BK20240610

      ! --- GROUNDWATER-N  !#BK
      allocate(Gw_LDON0(domain%nbasin), &
               Gw_RDON0(domain%nbasin), &
               Gw_NH40(domain%nbasin), &
               Gw_NO30(domain%nbasin),stat=istat)
      allocate(Gw_LDONsogw0(domain%nbasin), &
               Gw_RDONsogw0(domain%nbasin), &
               Gw_NH4sogw0(domain%nbasin), & 
               Gw_NO3sogw0(domain%nbasin),stat=istat)  !BK20240610  !BK20250310 !BK20260318
      allocate(Gw_LDONgwch0(domain%nbasin), &
               Gw_RDONgwch0(domain%nbasin), &
               Gw_NH4gwch0(domain%nbasin), & 
               Gw_NO3gwch0(domain%nbasin),stat=istat)  !BK20240610

      ! --- GROUNDWATER-P  !#BK
      allocate(Gw_LDOP0(domain%nbasin), &
               Gw_RDOP0(domain%nbasin), &
               Gw_PO40(domain%nbasin), stat=istat)
      allocate(Gw_LDOPsogw0(domain%nbasin), &
               Gw_RDOPsogw0(domain%nbasin), &  
               Gw_PO4sogw0(domain%nbasin),stat=istat)  !BK20240610  !BK20250310 !BK20260318
      allocate(Gw_LDOPgwch0(domain%nbasin), &
               Gw_RDOPgwch0(domain%nbasin), & 
               Gw_PO4gwch0(domain%nbasin),stat=istat)  !BK20240610

      ! --- CHANNEL-C  !#BK
      allocate(Ch_ALGC0(domain%nch,nalg), & 
               Ch_ZOOC0(domain%nch), &
               Ch_MBMC0(domain%nch), &
               Ch_DIC0(domain%nch), &
               Ch_CO2EVAS(domain%nch),stat=istat)  !BK20251120 !BK20251126
      allocate(Ch_LPOC0(domain%nch), &
               Ch_RPOC0(domain%nch), &
               Ch_LDOC0(domain%nch), &
               Ch_RDOC0(domain%nch), &
               Ch_POCDEPOSIT0(domain%nch),stat=istat)
      allocate(Ch_ALGCusch0(domain%nch,nalg), &
               Ch_ZOOCusch0(domain%nch), &
               Ch_MBMCusch0(domain%nch), &
               Ch_DICusch0(domain%nch),stat=istat) !BK20251120 
      allocate(Ch_LPOCusch0(domain%nch), &
               Ch_RPOCusch0(domain%nch), &
               Ch_LDOCusch0(domain%nch), &
               Ch_RDOCusch0(domain%nch),stat=istat)
      allocate(Ch_LPOCsurf0(domain%nch), &
               Ch_RPOCsurf0(domain%nch), &
               Ch_LDOCsurf0(domain%nch), &
               Ch_RDOCsurf0(domain%nch), &
               Ch_MBMCsurf0(domain%nch),stat=istat)  !BK20240513  !BK20240610  !BK20250310
      allocate(Ch_LDOCintf0(domain%nch), &
               Ch_RDOCintf0(domain%nch),stat=istat)  !BK20240513  !BK20240610  !BK20250310
      allocate(Ch_LDOCgwch0(domain%nch), &
               Ch_RDOCgwch0(domain%nch),stat=istat)  !BK20250310              
      allocate(Ch_ALGCdsch(domain%nch,nalg), &
               Ch_ZOOCdsch(domain%nch), &
               Ch_MBMCdsch(domain%nch), &
               Ch_DICdsch(domain%nch),stat=istat) !BK20250917 !BK20251120
      allocate(Ch_LPOCdsch(domain%nch), &
               Ch_RPOCdsch(domain%nch), &
               Ch_LDOCdsch(domain%nch), &
               Ch_RDOCdsch(domain%nch),stat=istat)  !BK20250917

      ! --- CHANNEL-N
      allocate(Ch_ALGN0(domain%nch,nalg), &
               Ch_ZOON0(domain%nch), &
               Ch_MBMN0(domain%nch), &
               Ch_NH40(domain%nch), &
               Ch_NO30(domain%nch),stat=istat)
      allocate(Ch_LPON0(domain%nch), &
               Ch_RPON0(domain%nch), &
               Ch_LDON0(domain%nch), &
               Ch_RDON0(domain%nch), &
               Ch_PONDEPOSIT0(domain%nch),stat=istat)
      allocate(Ch_ALGNusch0(domain%nch,nalg), &
               Ch_ZOONusch0(domain%nch), &
               Ch_MBMNusch0(domain%nch), &
               Ch_NH4usch0(domain%nch), &
               Ch_NO3usch0(domain%nch),stat=istat)
      allocate(Ch_LPONusch0(domain%nch), &
               Ch_RPONusch0(domain%nch), &
               Ch_LDONusch0(domain%nch), &
               Ch_RDONusch0(domain%nch),stat=istat)
      allocate(Ch_LPONsurf0(domain%nch), &
               Ch_RPONsurf0(domain%nch), &
               Ch_LDONsurf0(domain%nch), &
               Ch_RDONsurf0(domain%nch),stat=istat)  !BK20240610 !BK20250310
      allocate(Ch_MBMNsurf0(domain%nch), &
               Ch_NH4surf0(domain%nch), &
               Ch_NO3surf0(domain%nch),stat=istat)  !BK20240610 !BK20250310
      allocate(Ch_LDONintf0(domain%nch), &
               Ch_RDONintf0(domain%nch), &
               Ch_NH4intf0(domain%nch), &
               Ch_NO3intf0(domain%nch),stat=istat)  !BK20240610 !BK20250310
      allocate(Ch_LDONgwch0(domain%nch), &
               Ch_RDONgwch0(domain%nch), &
               Ch_NH4gwch0(domain%nch), &
               Ch_NO3gwch0(domain%nch),stat=istat)  !BK20250310
      allocate(Ch_ALGNdsch(domain%nch,nalg), &
               Ch_ZOONdsch(domain%nch), &
               Ch_MBMNdsch(domain%nch), &
               Ch_NH4dsch(domain%nch), &
               Ch_NO3dsch(domain%nch),stat=istat)  !BK20250917
      allocate(Ch_LPONdsch(domain%nch), &
               Ch_RPONdsch(domain%nch), &
               Ch_LDONdsch(domain%nch), &
               Ch_RDONdsch(domain%nch),stat=istat)  !BK20250917

      ! --- CHANNEL-P
      allocate(Ch_ALGP0(domain%nch,nalg), &
               Ch_ZOOP0(domain%nch), &
               Ch_MBMP0(domain%nch),stat=istat)
      allocate(Ch_PO40(domain%nch), &
               Ch_PIPA0(domain%nch), &
               Ch_PIPS0(domain%nch),stat=istat)
      allocate(Ch_LPOP0(domain%nch), &
               Ch_RPOP0(domain%nch), &
               Ch_LDOP0(domain%nch), &
               Ch_RDOP0(domain%nch), &
               Ch_POPDEPOSIT0(domain%nch),stat=istat)
      allocate(Ch_ALGPusch0(domain%nch,nalg), &
               Ch_ZOOPusch0(domain%nch), &
               Ch_MBMPusch0(domain%nch),stat=istat)
      allocate(Ch_PO4usch0(domain%nch), &
               Ch_PIPAusch0(domain%nch), &
               Ch_PIPSusch0(domain%nch),stat=istat)
      allocate(Ch_LPOPusch0(domain%nch), &
               Ch_RPOPusch0(domain%nch), &
               Ch_LDOPusch0(domain%nch), &
               Ch_RDOPusch0(domain%nch),stat=istat)
      allocate(Ch_LPOPsurf0(domain%nch), &
               Ch_RPOPsurf0(domain%nch), &
               Ch_LDOPsurf0(domain%nch), &
               Ch_RDOPsurf0(domain%nch),stat=istat)  !BK20240610 !BK20250310
      allocate(Ch_MBMPsurf0(domain%nch), &
               Ch_PO4surf0(domain%nch), &
               Ch_PIPAsurf0(domain%nch), &
               Ch_PIPSsurf0(domain%nch),stat=istat)  !BK20240610 !BK20250310
      allocate(Ch_LDOPintf0(domain%nch), &
               Ch_RDOPintf0(domain%nch), &
               Ch_PO4intf0(domain%nch),stat=istat)  !BK20240610 !BK20250310
      allocate(Ch_LDOPgwch0(domain%nch), &
               Ch_RDOPgwch0(domain%nch), &
               Ch_PO4gwch0(domain%nch),stat=istat)   !BK20250310              
      allocate(Ch_ALGPdsch(domain%nch,nalg), &
               Ch_ZOOPdsch(domain%nch), &
               Ch_MBMPdsch(domain%nch),stat=istat) !BK20250917
      allocate(Ch_PO4dsch(domain%nch), &
               Ch_PIPAdsch(domain%nch), &
               Ch_PIPSdsch(domain%nch),stat=istat) !BK20250917
      allocate(Ch_LPOPdsch(domain%nch), &
               Ch_RPOPdsch(domain%nch), &
               Ch_LDOPdsch(domain%nch), &
               Ch_RDOPdsch(domain%nch),stat=istat) !BK20250917

      !initialize
      overSed%Sol       = 0.0
      overSed%Sow       = 0.0
      channelSed%Sch    = 0.0
      channelSed%Scw    = 0.0 
      channelSed%Spnt   = 0.0
      channelSed%Susch0 = 0.0
      channelSed%Solch0 = 0.0   !BK20250925

      !-----------------------------------Y.Kwon(20250621)
      !overSed%Sol0      = 0.0  !BK20251203
      !overSed%Sow0      = 0.0  !BK20251203
      overSed%Solcp     = 0.0
      overSed%Ssurf     = 0.0
      overSed%Ssurf0    = 0.0
      overSed%Srdet     = 0.0
      overSed%Soler     = 0.0
      overSed%Soldp     = 0.0
      overSed%Colsd     = 0.0
      overSed%SOILPSF   = 0.0
      channelSed%Sch0   = 0.0
      channelSed%Scw0   = 0.0
      channelSed%Spnt   = 0.0
      channelSed%Sxpnt  = 0.0
      channelSed%Sabs   = 0.0
      channelSed%Sdis   = 0.0
      channelSed%Sxabs  = 0.0
      channelSed%Sxdis  = 0.0
      !channelSed%Sirri  = 0.0  !BK20260510 removed Sirri
      channelSed%Sinput = 0.0
      channelSed%Sdsch  = 0.0
      channelSed%Scher  = 0.0
      channelSed%Schdp  = 0.0
      channelSed%Schcp  = 0.0
      channelSed%Cchsd  = 0.0

      KLIT = 0.0
      KRES = 0.0
      KEXC = 0.0
      KMAN = 0.0
      KLDOM = 0.0
      KRPOM = 0.0
      KRDOM = 0.0
      KMBM = 0.0
      KDEN = 0.0
      KNIT = 0.0
      KVOL = 0.0
      KAGMAX = 0.0
      KARMAX = 0.0
      KAEMAX = 0.0
      KAMMAX = 0.0
      FROOT = 0.0
      TMBM = 0.0
      MMBM = 0.0
      TZOO = 0.0
      MZOO = 0.0
      TALG = 0.0
      MALG = 0.0
      !-----------------------------------Y.Kwon(20250621)

      !endif  !BK20251209 commented out - restart_flag not being used

      channel_input%Wst        = 0.0
      channel_input%SideSlopch = 0.0
      !channel_input%TopWdth    = 0.0  !BK20251117
      !channel_input%TopWdthCC  = 0.0  !BK20251117
      channel_input%Lst        = 0.0
      channel_input%Sf         = 0.0
      SedCNP_hydro%subbasinID  = -999
      SedCNP_hydro%smcmax      = 0.0
      SedCNP_hydro%smcref      = 0.0

      ! --- Read in soil particle size fractions for soil texture classes
      call ReadSOILPSF()  !BK20240205

      ! --- Initialize all CNP variables to zero. !BK20250906
      ! C
      So_LPOCLIT0      = 0.0; So_LPOCRES0      = 0.0
      So_LPOCEXC0      = 0.0; So_LPOCMAN0      = 0.0
      So_RPOC0         = 0.0; So_LDOC0         = 0.0
      So_RDOC0         = 0.0; So_MBMC0         = 0.0
      So_LDOCperc0     = 0.0; So_RDOCperc0     = 0.0 !BK20251209
      So_LPOCsurf0     = 0.0; So_RPOCsurf0     = 0.0
      So_LDOCsurf0     = 0.0; So_RDOCsurf0     = 0.0
      So_MBMCsurf0     = 0.0; So_LDOCintf0     = 0.0
      So_RDOCintf0     = 0.0
      So_LDOCgwso0     = 0.0; So_RDOCgwso0     = 0.0  !BK20260316
      So_LDOCsogw0     = 0.0; So_RDOCsogw0     = 0.0  !BK20260316

      Gw_LDOC0         = 0.0; Gw_RDOC0         = 0.0
      Gw_LDOCsogw0     = 0.0; Gw_RDOCsogw0     = 0.0 !BK20251209  !BK20260316
      Gw_LDOCgwch0     = 0.0; Gw_RDOCgwch0     = 0.0

      Ch_ALGC0         = 0.0; Ch_ZOOC0         = 0.0 
      Ch_MBMC0         = 0.0; Ch_DIC0          = 0.0  !BK20251120
      Ch_CO2EVAS       = 0.0                           !BK20251126
      Ch_LPOC0         = 0.0; Ch_RPOC0         = 0.0
      Ch_LDOC0         = 0.0; Ch_RDOC0         = 0.0
      Ch_POCDEPOSIT0   = 0.0; Ch_ALGCusch0     = 0.0
      Ch_ZOOCusch0     = 0.0; Ch_MBMCusch0     = 0.0  
      Ch_DICusch0      = 0.0; Ch_LPOCusch0     = 0.0
      Ch_RPOCusch0     = 0.0; Ch_LDOCusch0     = 0.0  !BK20251120
      Ch_RDOCusch0     = 0.0
      Ch_LPOCsurf0     = 0.0; Ch_RPOCsurf0     = 0.0
      Ch_LDOCsurf0     = 0.0; Ch_RDOCsurf0     = 0.0 !BK20251209
      Ch_MBMCsurf0     = 0.0
      Ch_LDOCintf0     = 0.0; Ch_RDOCintf0     = 0.0
      Ch_LDOCgwch0     = 0.0; Ch_RDOCgwch0     = 0.0 !BK20251209

      ! N
      So_LPONLIT0      = 0.0; So_LPONRES0      = 0.0
      So_LPONEXC0      = 0.0; So_LPONMAN0      = 0.0
      So_RPON0         = 0.0; So_LDON0         = 0.0
      So_RDON0         = 0.0; So_MBMN0         = 0.0
      So_NH40          = 0.0; So_NO30          = 0.0
      So_LDONperc0     = 0.0; So_RDONperc0     = 0.0
      So_NH4perc0      = 0.0; So_NO3perc0      = 0.0 !BK20251209
      So_LPONsurf0     = 0.0; So_RPONsurf0     = 0.0
      So_LDONsurf0     = 0.0; So_RDONsurf0     = 0.0
      So_MBMNsurf0     = 0.0; So_NH4surf0      = 0.0
      So_NO3surf0      = 0.0
      So_LDONintf0     = 0.0; So_RDONintf0     = 0.0
      So_NH4intf0      = 0.0; So_NO3intf0      = 0.0
      So_LDONgwso0     = 0.0; So_RDONgwso0     = 0.0  !BK20260316
      So_NH4gwso0      = 0.0; So_NO3gwso0      = 0.0  !BK20260316
      So_LDONsogw0     = 0.0; So_RDONsogw0     = 0.0  !BK20260316
      So_NH4sogw0      = 0.0; So_NO3sogw0      = 0.0  !BK20260316

      Gw_LDON0         = 0.0; Gw_RDON0         = 0.0
      Gw_NH40          = 0.0; Gw_NO30          = 0.0
      Gw_LDONsogw0     = 0.0; Gw_RDONsogw0     = 0.0  !BK20260316
      Gw_NH4sogw0      = 0.0; Gw_NO3sogw0      = 0.0 !BK20251209  !BK20260316
      Gw_LDONgwch0     = 0.0; Gw_RDONgwch0     = 0.0 
      Gw_NH4gwch0      = 0.0; Gw_NO3gwch0      = 0.0
      
      Ch_ALGN0         = 0.0; Ch_ZOON0         = 0.0
      Ch_MBMN0         = 0.0; Ch_LPON0         = 0.0
      Ch_RPON0         = 0.0; Ch_LDON0         = 0.0
      Ch_RDON0         = 0.0; Ch_NH40          = 0.0
      Ch_NO30          = 0.0; Ch_PONDEPOSIT0   = 0.0
      Ch_ALGNusch0     = 0.0; Ch_ZOONusch0     = 0.0
      Ch_MBMNusch0     = 0.0
      Ch_LPONusch0     = 0.0; Ch_RPONusch0     = 0.0
      Ch_LDONusch0     = 0.0; Ch_RDONusch0     = 0.0
      Ch_NH4usch0      = 0.0; Ch_NO3usch0      = 0.0
      Ch_LPONsurf0     = 0.0; Ch_RPONsurf0     = 0.0
      Ch_LDONsurf0     = 0.0; Ch_RDONsurf0     = 0.0 !BK20251209
      Ch_MBMNsurf0     = 0.0; Ch_NH4surf0      = 0.0
      Ch_NO3surf0      = 0.0                           !BK20251209
      Ch_LDONintf0     = 0.0; Ch_RDONintf0     = 0.0
      Ch_NH4intf0      = 0.0; Ch_NO3intf0      = 0.0 !BK20251209
      Ch_LDONgwch0     = 0.0; Ch_RDONgwch0     = 0.0
      Ch_NH4gwch0      = 0.0; Ch_NO3gwch0      = 0.0 !BK20251209

      ! P
      So_LPOPLIT0      = 0.0; So_LPOPRES0      = 0.0
      So_LPOPEXC0      = 0.0; So_LPOPMAN0      = 0.0
      So_RPOP0         = 0.0; So_LDOP0         = 0.0
      So_RDOP0         = 0.0; So_MBMP0         = 0.0
      So_PO40          = 0.0; So_PIPA0         = 0.0
      So_PIPS0         = 0.0
      So_LDOPperc0     = 0.0; So_RDOPperc0     = 0.0
      So_PO4perc0      = 0.0                           !BK20251209
      So_LPOPsurf0     = 0.0; So_RPOPsurf0     = 0.0
      So_LDOPsurf0     = 0.0; So_RDOPsurf0     = 0.0
      So_MBMPsurf0     = 0.0; So_PIPAsurf0     = 0.0
      So_PIPSsurf0     = 0.0; So_PO4surf0      = 0.0
      So_LDOPintf0     = 0.0; So_RDOPintf0     = 0.0
      So_PO4intf0      = 0.0
      So_LDOPgwso0     = 0.0; So_RDOPgwso0     = 0.0  !BK20260316
      So_PO4gwso0      = 0.0                          !BK20260316
      So_LDOPsogw0     = 0.0; So_RDOPsogw0     = 0.0  !BK20260316
      So_PO4sogw0      = 0.0                          !BK20260316  

      Gw_LDOP0         = 0.0; Gw_RDOP0         = 0.0
      Gw_PO40          = 0.0
      Gw_LDOPsogw0     = 0.0; Gw_RDOPsogw0     = 0.0  !BK20260316
      Gw_PO4sogw0      = 0.0                           !BK20251209  !BK20260316
      Gw_LDOPgwch0     = 0.0; Gw_RDOPgwch0     = 0.0
      Gw_PO4gwch0      = 0.0

      Ch_ALGP0         = 0.0; Ch_ZOOP0         = 0.0
      Ch_MBMP0         = 0.0; Ch_LPOP0         = 0.0
      Ch_RPOP0         = 0.0; Ch_LDOP0         = 0.0
      Ch_RDOP0         = 0.0; Ch_PIPA0         = 0.0
      Ch_PIPS0         = 0.0; Ch_PO40          = 0.0
      Ch_POPDEPOSIT0   = 0.0
      Ch_ALGPusch0     = 0.0; Ch_ZOOPusch0     = 0.0
      Ch_MBMPusch0     = 0.0
      Ch_LPOPusch0     = 0.0; Ch_RPOPusch0     = 0.0
      Ch_LDOPusch0     = 0.0; Ch_RDOPusch0     = 0.0
      Ch_PO4usch0      = 0.0; Ch_PIPAusch0     = 0.0
      Ch_PIPSusch0     = 0.0
      Ch_LPOPsurf0     = 0.0; Ch_RPOPsurf0     = 0.0
      Ch_LDOPsurf0     = 0.0; Ch_RDOPsurf0     = 0.0 !BK20251209
      Ch_MBMPsurf0     = 0.0; Ch_PO4surf0      = 0.0
      Ch_PIPAsurf0     = 0.0; Ch_PIPSsurf0     = 0.0 !BK20251209
      Ch_LDOPintf0     = 0.0; Ch_RDOPintf0     = 0.0
      Ch_PO4intf0      = 0.0                           !BK20251209
      Ch_LDOPgwch0     = 0.0; Ch_RDOPgwch0     = 0.0
      Ch_PO4gwch0      = 0.0                           !BK20251209

   end subroutine SedCNP_driver_ini


   subroutine SedCNPmodel_surface_input()

      use module_SedCNPvariables
      use module_SedCNP_in
      use module_SedCNP_out
!=====||__WHQ5403__||=====!
!
      ! SedCNPmodel already visible via module_SedCNP_in/module_SedCNP_out above.
      !use config_base, only: SedCNPmodel
!
!=====||__WHQ5403__||=====!

      implicit none

      character(len=256)                   :: geo_finegrid_flnm
      character(len=256)                   :: var_name
      character(len=256)                   :: filename_CH, filename_soil
      character(len=256)                   :: filename_subbasinID
      real, dimension(domain%ix,domain%jx) :: lat, lon
      integer                              :: status

      !read surface inputs
      !var_name = "LATITUDE"
      !call readRT2d_real(var_name,lat,domain%ixrt,domain%jxrt,&
      !                  trim(geo_finegrid_flnm))
      !var_name = "LONGITUDE"
      !call readRT2d_real(var_name,lon,domain%ixrt,domain%jxrt,&
      !                  trim(geo_finegrid_flnm)) 

      filename_CH = './DOMAIN/Route_Link.nc'
      !bottom width of channel [m]
      status = get1d_ch_real("BtmWdth",channel_input%Wst,domain%nch,trim(filename_CH))

      ! compound channel NOT implemented for WHQ yet   !BK20251117
      !if (SedCNPmodel%SedCNP_compound_channel) then
      !   !Top width of channel [m]
      !   status = get1d_ch_real("TopWdth",channel_input%TopWdth,domain%nch,trim(filename_CH))
      !   !Compound Channel Top Width [m]
      !   status = get1d_ch_real("TopWdthCC",channel_input%TopWdthCC,domain%nch,trim(filename_CH))
      !endif

      ! Channel side slope
      status = get1d_ch_real("ChSlp", channel_input%SideSlopch, &
                        domain%nch, trim(filename_CH))
      !stream length [m]
      status = get1d_ch_real("Length", channel_input%Lst, &
                        domain%nch, trim(filename_CH))
      !channel slope
      status = get1d_ch_real("So", channel_input%Sf, &
                        domain%nch, trim(filename_CH))

      filename_subbasinID = './DOMAIN/GWBASINS.nc'
      !subbasin ID (from GWBASINS.nc)
      status = get2d_int("BASIN", SedCNP_hydro%subbasinID, &
                        domain%ixrt, domain%jxrt, trim(filename_subbasinID))

      geo_finegrid_flnm = './DOMAIN/Fulldom_hires.nc'
      !link ID of grid cells - excluding "connect-only" links (from Fulldom_hires.nc)
      status = get2d_int("LINKID", SedCNP_hydro%linkID_grid, &
                        domain%ixrt, domain%jxrt, trim(geo_finegrid_flnm))   !BK20240616

      !=====||__WHQ5403__||=====!
      ! Fulldom_hires.nc and GWBASINS.nc are stored north-up (y decreasing), whereas the LSM/RTOUT
      ! outputs read by SedCNP are south-up (the model grid, as read by the hydro model). Without
      ! this flip, the gw basins and stream pixels were overlaid on the mirrored LSM cells.
      if (is_north_up(trim(filename_subbasinID))) then
         SedCNP_hydro%subbasinID = SedCNP_hydro%subbasinID(:, domain%jxrt:1:-1)
         write(6,*) 'INFO: SedCNP: ', trim(filename_subbasinID), ' is north-up; flipped to the model grid'
      endif
      if (is_north_up(trim(geo_finegrid_flnm))) then
         SedCNP_hydro%linkID_grid = SedCNP_hydro%linkID_grid(:, domain%jxrt:1:-1)
         write(6,*) 'INFO: SedCNP: ', trim(geo_finegrid_flnm), ' is north-up; flipped to the model grid'
      endif
      !=====||__WHQ5403__||=====!

      !---WHQ
      filename_soil= './DOMAIN/soil_properties.nc'
      !soil water content at saturation [m3 m-3]
      status = get4d_soil_real("smcmax", SedCNP_hydro%smcmax, &
                        domain%ix, domain%jx, trim(filename_soil))
      !soil water content at field capacity [m3 m-3]
      status = get4d_soil_real("smcref", SedCNP_hydro%smcref, &
                        domain%ix, domain%jx, trim(filename_soil))
      !---WHQ

   end subroutine SedCNPmodel_surface_input


   subroutine SedCNPmodel_hydro_input(itime)
      
      use module_SedCNPvariables
      !use module_HYDRO_io,          only: get2d_lsm_real
      use module_SedCNP_in,         only: get3d_lsm_real, get4d_lsm_real, &
                                          get4d_soil_real, get1d_ch_real, &
                                          get1d_ch_int, get2d_int, &
                                          get3d_lsm_int, get_lsm_time, &  !BK20231005 !Y.Kwon20241125 
                                          get_lsm_dx                      !Y.Kwon(20250621)
      use module_SedCNP_in,         only: get1d_bas_real   !WHQ5403 per-gw-basin variables by basin ID
      use module_SedCNP_out
!=====||__WHQ5403__||=====!
!
      ! SedCNPmodel already visible via module_SedCNP_out above.
      !use config_base,              only: SedCNPmodel
!
!=====||__WHQ5403__||=====!
      use Module_Date_utilities_rt, only: geth_newdate

      implicit none

      integer, intent(in) :: itime  ! timestep loop
      character(len=256)  :: filename, hydro_outdir, hydro_indir  !BK20251207
      character(len=256)  :: filename_RT, filename_CHRT, filename_subbasinAREA !Y.Kwon 20230422
      character(len=256)  :: filename_link !Y.Kwon 20230604
      !character(len=256)  :: filename_forcing  !BK20231029  !BK20240129
      character(len=256)  :: filename_GW   !BK20240620
      integer             :: ich  !, chid, gwid   !BK20240620 !BK20240726 !BK20241102
      integer             :: i,j    !kyh_debug  
      integer             :: status
      real                :: dt    !timestep
      !real(8)             :: buf1(domain%nch), buf2(domain%nch), buf3(domain%nch)  !BK20240620
      !real(8)             :: buf4(domain%nch), buf5(domain%nch)    !BK20240620

      ! Temporary buffers to hold data read from NetCDF files before re-indexing !BK20250907
      !real(8), allocatable :: q_lateral_in(:), q_dsch_in(:), head_in(:), chVa_in(:)
      !real(8), allocatable :: Qpnt_in(:), Qabs_in(:), Qdis_in(:)
      !integer              :: chid

      character(len=4)    :: yyyy
      character(len=2)    :: mm,dd,hr,mn

      dt = real(SedCNPmodel%SedCNP_timestep)
      hydro_outdir = trim(SedCNPmodel%LDASOUT_dir)
      hydro_indir = trim(".")  !BK20251207 hardcoded

         if(itime == 1) then
            write(unit=yyyy, fmt='(i4.4)') SedCNPmodel%SedCNP_start_year
            write(unit=mm,   fmt='(i2.2)') SedCNPmodel%SedCNP_start_month
            write(unit=dd,   fmt='(i2.2)') SedCNPmodel%SedCNP_start_day
            write(unit=hr,   fmt='(i2.2)') SedCNPmodel%SedCNP_start_hour
            write(unit=mn,   fmt='(i2.2)') SedCNPmodel%SedCNP_start_min

            dateSedCNP%olddate = trim(yyyy)//'-'//trim(mm)//'-'//trim(dd)//&
                              '-'//trim(hr)//':'//trim(mn)//':00'

            !get newdate for next timestep !Advance dateSedCNP%olddate to the NEXT timestep  
            call geth_newdate(dateSedCNP%newdate, dateSedCNP%olddate, nint(dt))  !BK20250717 
            dateSedCNP%olddate = dateSedCNP%newdate                              !BK20250717 
         endif

      !create filename
      filename = trim(hydro_outdir)//'/'//&
                  trim(dateSedCNP%olddate(1:4))//&
                  trim(dateSedCNP%olddate(6:7))//&
                  trim(dateSedCNP%olddate(9:10))//&
                  trim(dateSedCNP%olddate(12:13))//&
                  trim(dateSedCNP%olddate(15:16))//&
                  '.LDASOUT_DOMAIN1'

      filename_RT = trim(hydro_outdir)//'/'//&
                  trim(dateSedCNP%olddate(1:4))//&
                  trim(dateSedCNP%olddate(6:7))//&
                  trim(dateSedCNP%olddate(9:10))//&
                  trim(dateSedCNP%olddate(12:13))//&
                  trim(dateSedCNP%olddate(15:16))//&
                  '.RTOUT_DOMAIN1'

      filename_CHRT = trim(hydro_outdir)//'/'//&
                  trim(dateSedCNP%olddate(1:4))//&
                  trim(dateSedCNP%olddate(6:7))//&
                  trim(dateSedCNP%olddate(9:10))//&
                  trim(dateSedCNP%olddate(12:13))//&
                  trim(dateSedCNP%olddate(15:16))//&
                  '.CHRTOUT_DOMAIN1'

      filename_GW = trim(hydro_outdir)//'/'//&
                  trim(dateSedCNP%olddate(1:4))//&
                  trim(dateSedCNP%olddate(6:7))//&
                  trim(dateSedCNP%olddate(9:10))//&
                  trim(dateSedCNP%olddate(12:13))//&
                  trim(dateSedCNP%olddate(15:16))//&
                  '.GWOUT_DOMAIN1'

      !filename_subbasinID = trim(hydro_outdir)//'/subbasinID.nc'   !Y.Kwon 20230422   !BK20231119
      !filename_subbasinAREA = trim(hydro_outdir)//'/DOMAIN/GWBUCKPARM.nc'  !Y.Kwon 20230604
      filename_subbasinAREA = trim(hydro_indir)//'/DOMAIN/GWBUCKPARM.nc'  !Y.Kwon 20230604 !BK20251207

      !filename_link = trim(hydro_outdir)//'/DOMAIN/Route_Link.nc'  !Y.Kwon 20230604
      filename_link = trim(hydro_indir)//'/DOMAIN/Route_Link.nc'  !Y.Kwon 20230604 !BK20251207
      !*** NOTE: Route_Link.nc include "connect-only" links *** !BK20240616


      !filename_forcing = trim(hydro_outdir)//'FORCING/'//trim(dateSedCNP%olddate(1:4))//&   !BK20240129
      !           trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
      !           trim(dateSedCNP%olddate(12:13))//&
      !           '.LDASIN_DOMAIN1'
      
      !time [minutes since 1970-01-01 00:00:00 UTC]
      status = get_lsm_time("time", SedCNP_hydro%time, trim(filename))  !Y.Kwon 20241125
      !dx [m]
      !status = get_lsm_dx("x",SedCNP_hydro%dx,trim(filename))  !Y.Kwon(20250621)
      SedCNP_hydro%dx(1) = get_lsm_dx("x", trim(filename))  !BK20251210
      !rainfall [mm/s]
      status = get3d_lsm_real("QRAIN", SedCNP_hydro%rainfall, &
                        domain%ix, domain%jx, trim(filename))
      !surface runoff [mm/s]
      status = get3d_lsm_real("RUNSRF", SedCNP_hydro%runsrf, &
                        domain%ix, domain%jx, trim(filename))
      !overland water depth [mm]
      status = get3d_lsm_real("OLWDepth", SedCNP_hydro%olwdepth, &
                        domain%ix, domain%jx, trim(filename))
      !irrigation water depth [mm]
      status = get3d_lsm_real("IRRIDepth", SedCNP_hydro%irridepth,&
                        domain%ix, domain%jx, trim(filename))
      !infiltration [mm/s]
      status = get3d_lsm_real("INFLTrate", SedCNP_hydro%infltrate, &
                        domain%ix, domain%jx, trim(filename))
      !evaporation from overland water [mm/s]
      status = get3d_lsm_real("EvapWater", SedCNP_hydro%evapwater, &
                        domain%ix, domain%jx, trim(filename))
      !green vegetation fraction [-]
      status = get3d_lsm_real("FVEG", SedCNP_hydro%FVEG, &
                        domain%ix, domain%jx, trim(filename))
      !vegetation height [m]
      status = get3d_lsm_real("CanHeight", SedCNP_hydro%CanHt, &
                        domain%ix, domain%jx, trim(filename))
      !land slope (from RT output)
      status = get3d_lsm_real("SO8LD_Vmax", SedCNP_hydro%slope, &
                        domain%ixrt, domain%jxrt, trim(filename_RT))
      !subbasin ID (from RT output)   !Y.Kwon 20230422
      !status = get3d_lsm_real("subbasinID", SedCNP_hydro%subbasinID, &
      !                  domain%ixrt, domain%jxrt, trim(filename_subbasinID))
      !flux from soil to gw gucket [mm] 
      status = get3d_lsm_real("q_sogw", SedCNP_hydro%q_sogw, &
                        domain%ixrt, domain%jxrt, trim(filename_RT))  !BK20260313
      !inter flow [mm] !Y.Kwon 20230501
      status = get3d_lsm_real("q_intf", SedCNP_hydro%q_intf, &
                        domain%ixrt, domain%jxrt, trim(filename_RT))
      !dominant vegetation category
      status = get3d_lsm_real("IVGTYP", SedCNP_hydro%vgtyp, &
                        domain%ix, domain%jx, trim(filename))
      !dominant soil texture category
      status = get3d_lsm_real("ISLTYP", SedCNP_hydro%sltyp, &
                        domain%ix, domain%jx, trim(filename))  !BK20240205 
      !! ----only for CNP-----------------------------------


      !solar radiation [W m-2]
      !status = get3d_lsm_real("SWDOWN", SOLRAD, &
      !                  domain%ix, domain%jx, trim(filename_forcing))  
      !SOLRAD = 300.0 !BK20231029   hardcoded for debugging - to be fixed 
      !status = get3d_lsm_real("SWDOWN", SedCNP_hydro%swdown, &
      !                  domain%ix,domain%jx,trim(filename))  !BK20240129 
      !status = get3d_lsm_real("SWFORC", SedCNP_hydro%swdown, &
      !                  domain%ix,domain%jx,trim(filename))  !BK20240129
      
      !soil layer temperature [K]
      status = get4d_lsm_real("SOIL_T", SedCNP_hydro%t_soil, &
                        domain%ix, domain%jx, trim(filename))
      !volumetric soil moisture [m3 m-3]
      status = get4d_lsm_real("SOIL_M", SedCNP_hydro%theta, &
                        domain%ix, domain%jx, trim(filename))
      !percolation between soil layers [mm s-1]
      status = get4d_lsm_real("Percolate", SedCNP_hydro%percolate, &
                        domain%ix, domain%jx, trim(filename))
      !vegetation water uptake from soil layers [mm s-1]
      status = get4d_lsm_real("ETRANLayer", SedCNP_hydro%ETRANLayer, &
                        domain%ix, domain%jx, trim(filename))
      !water storage in aquifer [mm]
      status = get3d_lsm_real("WA", SedCNP_hydro%WA, &
                        domain%ix, domain%jx, trim(filename))
      !groundwater discharge [mm/s]
      status = get3d_lsm_real("GWDischarge", SedCNP_hydro%GWDischarge, &
                        domain%ix, domain%jx, trim(filename))
      !channel inflow by surface (terrain) runoff [mm]
      status = get3d_lsm_real("QSTRMVOLRT", SedCNP_hydro%runsrfRT, &
                        domain%ixrt, domain%jxrt, trim(filename_RT))   !Y.Kwon 20230427
      !-----------------------------------------------------

      !canopy water evaporation [mm/s]
      status = get3d_lsm_real("ECAN", SedCNP_hydro%EvapCanopy, &
                        domain%ix, domain%jx, trim(filename))
      !evaporation from soil [mm/s]
      status = get3d_lsm_real("EvapSoil", SedCNP_hydro%EvapSoil, &
                        domain%ix, domain%jx, trim(filename))
      !transpiration [mm/s]
      status = get3d_lsm_real("ETRAN", SedCNP_hydro%Etran, &
                        domain%ix, domain%jx, trim(filename))
      !flux from gw bucket outflow to channel - q_gwch [m3/s]
      !status = get1d_ch_real("outflow", SedCNP_hydro%q_gwch, &
      !                  domain%nbasin, trim(filename_GW)) !BK20240620
      !WHQ5403 by gw basin ID (GWOUT has no Basin = -9999 fill row of GWBUCKPARM, so nbasin differs)
      status = get1d_bas_real("outflow", "feature_id", SedCNP_hydro%q_gwch, &
                        domain%nbasin, trim(filename_GW))
      if (status /= 0) call hydro_stop("SedCNPmodel_hydro_input: failed to read outflow in "//trim(filename_GW))
      !2m Air Temp [K]   (Y.Kwon 20240510)

      !depth in gw bucket [mm in GWOUT] converted to [m] for SedCNP !BK20260508
      !status = get1d_ch_real("depth", SedCNP_hydro%z_gwbas, &
      !                  domain%nbasin, trim(filename_GW))
      status = get1d_bas_real("depth", "feature_id", SedCNP_hydro%z_gwbas, &   !WHQ5403 by gw basin ID
                        domain%nbasin, trim(filename_GW))
      if (status /= 0) call hydro_stop("SedCNPmodel_hydro_input: failed to read depth in "//trim(filename_GW))
      if (status == 0) then   !---BK20260508
         SedCNP_hydro%z_gwbas = SedCNP_hydro%z_gwbas / 1000.0  ! mm -> m
      endif                   !---BK20260508

      status = get3d_lsm_real("T2M", SedCNP_hydro%T2m, &
                        domain%ix, domain%jx, trim(filename))
      !shortwave downward solar radiation [W/s2]  (Y.Kwon 20240830)
      status = get3d_lsm_real("SWDOWN", SedCNP_hydro%swdown, & 
                        domain%ix, domain%jx, trim(filename))

      !-------- CHRTOUT variables
      !flux from surface runoff to channel - q_surf (m3/s)  !BK20240620 !BK20240726 commented out
      !status = get1d_ch_real("qSfcLatRunoff",SedCNP_hydro%q_surf,domain%nch,trim(filename_CHRT))
      
      !flux from gw bucket to channel - q_lateral [m3/s]  
      status = get1d_ch_real("q_lateral", SedCNP_hydro%q_lateral, &
                        domain%nch, trim(filename_CHRT)) !BK20240620 qBucket replaced with q_lateral 
      !q_lateral: m3/s --> (m3/s) / (subbasin area)

      !channel discharge [m3/s]
      status = get1d_ch_real("streamflow",SedCNP_hydro%q_dsch,domain%nch,trim(filename_CHRT))
      !status = get1d_ch_real("streamflow",q_dsch_in,domain%nch,trim(filename_CHRT)) !BK20250907
      
      !River stage [m]
      status = get1d_ch_real("Head",SedCNP_hydro%head,domain%nch,trim(filename_CHRT))
      !status = get1d_ch_real("Head",head_in,domain%nch,trim(filename_CHRT)) !BK20250907
      
      !channel velocity [m/s]
      status = get1d_ch_real("velocity",SedCNP_hydro%chVa,domain%nch,trim(filename_CHRT))
      !status = get1d_ch_real("velocity",chVa_in,domain%nch,trim(filename_CHRT)) !BK20250907

      !WHQ5403 links inside lakes (lake_option > 0) have fill values (-9999) for streamflow and
      ! velocity in CHRTOUT; treat them as no flow so that the fill values are neither used in the
      ! channel calculations nor added to the upstream inflow of the downstream link
      where (SedCNP_hydro%q_dsch < 0.0) SedCNP_hydro%q_dsch = 0.0
      where (SedCNP_hydro%chVa < 0.0) SedCNP_hydro%chVa = 0.0
   
      ! !Qpnt [m3 s-1]
      ! status = get1d_ch_real("Qpnt",buf1,domain%nch,trim(filename_CHRT))  !BK20241028
      ! SedCNP_hydro%Qpnt(:,itime) = buf1(:)  !BK20241028

      ! !Qabs [m3 s-1]
      ! status = get1d_ch_real("Qabs",buf2,domain%nch,trim(filename_CHRT))  !BK20241028
      ! SedCNP_hydro%Qabs(:,itime) = buf2(:)  !BK20241028

      ! !Qdis [m3 s-1]
      ! status = get1d_ch_real("Qdis",buf3,domain%nch,trim(filename_CHRT))  !BK20241028
      ! SedCNP_hydro%Qdis(:,itime) = buf3(:)  !BK20241028

      !Qpnt [m3 s-1]
      status = get1d_ch_real("Qpnt", SedCNP_hydro%Qpnt(:,itime), &
                        domain%nch, trim(filename_CHRT))  !BK20241028 !BK20241102
      !Qabs [m3 s-1]
      status = get1d_ch_real("Qabs", SedCNP_hydro%Qabs(:,itime),&
                        domain%nch, trim(filename_CHRT))  !BK20241028 !BK20241102
      !Qdis [m3 s-1]
      status = get1d_ch_real("Qdis", SedCNP_hydro%Qdis(:,itime), &
                        domain%nch, trim(filename_CHRT))  !BK20241028 !BK20241102
      !basin area [km2] !Y.Kwon 20230604
      !status = get1d_ch_real("Area_sqkm", SedCNP_hydro%AREAbasin, &
      !                  domain%nbasin, trim(filename_subbasinAREA))  !kyh_debug (need to be modified)
      status = get1d_bas_real("Area_sqkm", "Basin", SedCNP_hydro%AREAbasin, &   !WHQ5403 by gw basin ID
                        domain%nbasin, trim(filename_subbasinAREA))
      !subbasin bucket height [mm]
      !status = get1d_ch_real("Zmax", SedCNP_hydro%Bucket_max, &
      !                  domain%nbasin, trim(filename_subbasinAREA))  !BK20240510
      status = get1d_bas_real("Zmax", "Basin", SedCNP_hydro%Bucket_max, &   !WHQ5403 by gw basin ID
                        domain%nbasin, trim(filename_subbasinAREA))
      !initial height of water in the bucket [mm]
      !status = get1d_ch_real("Zinit", SedCNP_hydro%Bucket_ini, &
      !                  domain%nbasin, trim(filename_subbasinAREA))  !BK20240510
      status = get1d_bas_real("Zinit", "Basin", SedCNP_hydro%Bucket_ini, &   !WHQ5403 by gw basin ID
                        domain%nbasin, trim(filename_subbasinAREA))
      !channel link ID
      status = get1d_ch_int("link", SedCNP_hydro%linkID, &
                        domain%nch, trim(filename_link)) !Y.Kwon 20230604
      !downstream (channel) ID
      status = get1d_ch_int("to", SedCNP_hydro%dwnstrm_linkID, &
                        domain%nch, trim(filename_link))  !Y.Kwon 20230604

      !!get newdate for next timestep
      !call geth_newdate(dateSedCNP%newdate, dateSedCNP%olddate, nint(dt))
      !dateSedCNP%olddate = dateSedCNP%newdate

   end subroutine SedCNPmodel_hydro_input


   subroutine SedCNPmodel_output(itime)

      use module_SedCNPvariables
!=====||__WHQ5403__||=====!
!
      ! SedCNPmodel already visible via module_SedCNP_out below.
      !use config_base, only: SedCNPmodel
!
!=====||__WHQ5403__||=====!
      use Module_Date_utilities_rt, only: geth_newdate
      use module_SedCNP_out
      use module_lakeSedCNP, only: write_Lake_Sed_output, write_Lake_C_output, &   !WHQ5403 lakes/reservoirs
                                   write_Lake_N_output, write_Lake_P_output

      implicit none

      integer, intent(in) :: itime  ! timestep loop
      real                :: dt     !timestep
      character(len=256)  :: sedCNP_outdir
      character(len=256)  :: filename_SEDOUT, filename_CHSEDOUT
      character(len=256)  :: filename_LAKE   !WHQ5403 lakes/reservoirs
      character(len=256)  :: filename_SOCOUT, filename_AQCOUT, filename_GWCOUT 
      character(len=256)  :: filename_SONOUT, filename_AQNOUT
      character(len=256)  :: filename_WSoCOUT, filename_CHCOUT
      character(len=256)  :: filename_SOPOUT, filename_AQPOUT
      character(len=256)  :: filename_CHNOUT, filename_GWNOUT
      character(len=256)  :: filename_CHPOUT, filename_GWPOUT
      !integer                :: lcover     !land code (Y.Kwon)  !BK20231005
      ! real                :: lcover     !land code (Y.Kwon)

      integer :: mkdir_status

      dt = real(SedCNPmodel%SedCNP_timestep)
      sedCNP_outdir = trim(SedCNPmodel%WHQOUT_dir)

      ! --- Check if the output directory is specified
      if (len_trim(sedCNP_outdir) == 0) then
         write(*,*) 'FATAL ERROR: WHQOUT_dir is not specified in the hydro.namelist.'
         call hydro_stop('FATAL ERROR: WHQOUT_dir not specified.')
      endif

      call system('mkdir -p '//trim(sedCNP_outdir), status=mkdir_status)
      if (mkdir_status /= 0) then
         write(*,*) 'FATAL ERROR: Could not create output directory: ', trim(sedCNP_outdir)
         call hydro_stop('FATAL ERROR: Could not create output directory.')
      endif

      !create SEDOUT filename
      filename_SEDOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.SEDOUT_DOMAIN1'
      call write_Sed_output(filename_SEDOUT)

      !create CHSEDOUT filename
      filename_CHSEDOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                           '.CHSEDOUT_DOMAIN1'
      call write_Ch_Sed_output(filename_CHSEDOUT)

      !create SOCOUT filename
      filename_SOCOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.SOCOUT_DOMAIN1'
      call write_So_C_output(filename_SOCOUT)  !BK20250928

      !create AQCOUT filename
      ! filename_AQCOUT = trim(sedCNP_outdir)//'/'//&
      !                   trim(dateSedCNP%olddate(1:4))//&
      !                   trim(dateSedCNP%olddate(6:7))//&
      !                   trim(dateSedCNP%olddate(9:10))//&
      !                   trim(dateSedCNP%olddate(12:13))//&
      !                   trim(dateSedCNP%olddate(15:16))//&
      !                      '.AQCOUT_DOMAIN1'
      !call write_Aq_C_output(filename_AQCOUT)

      !create GWCOUT filename
      filename_GWCOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.GWCOUT_DOMAIN1'
      call write_Gw_C_output(filename_GWCOUT)

      !create CHCOUT filename
      filename_CHCOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.CHCOUT_DOMAIN1'
      call write_Ch_C_output(filename_CHCOUT)   !kyh_debug

   !-------------------------------------------------------------------------------Y.Kwon 20240714
      !create SONOUT filename
      filename_SONOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.SONOUT_DOMAIN1'
      call write_So_N_output(filename_SONOUT)  !BK20250928

      !create AQNOUT filename
      ! filename_AQNOUT = trim(sedCNP_outdir)//'/'//&
      !                   trim(dateSedCNP%olddate(1:4))//&
      !                   trim(dateSedCNP%olddate(6:7))//&
      !                   trim(dateSedCNP%olddate(9:10))//&
      !                   trim(dateSedCNP%olddate(12:13))//&
      !                   trim(dateSedCNP%olddate(15:16))//&
      !                   '.AQNOUT_DOMAIN1'
      !call write_Aq_N_output(filename_AQNOUT)

      !create GWNOUT filename
      filename_GWNOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.GWNOUT_DOMAIN1'
      call write_Gw_N_output(filename_GWNOUT)
   !-------------------------------------------------------------------------------Y.Kwon 20240714

      !create CHNOUT filename
      filename_CHNOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.CHNOUT_DOMAIN1'
      call write_Ch_N_output(filename_CHNOUT)

   !-------------------------------------------------------------------------------Y.Kwon 20240714
      !create SOPOUT filename
      filename_SOPOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.SOPOUT_DOMAIN1'
      call write_So_P_output(filename_SOPOUT)  !BK20250928

      !create AQPOUT filename
      ! filename_AQPOUT = trim(sedCNP_outdir)//'/'//&
      !                   trim(dateSedCNP%olddate(1:4))//&
      !                   trim(dateSedCNP%olddate(6:7))//&
      !                   trim(dateSedCNP%olddate(9:10))//&
      !                   trim(dateSedCNP%olddate(12:13))//&
      !                   trim(dateSedCNP%olddate(15:16))//&
      !                   '.AQPOUT_DOMAIN1'
      !call write_Aq_P_output(filename_AQPOUT)

      !create GWPOUT filename
      filename_GWPOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.GWPOUT_DOMAIN1'
      call write_Gw_P_output(filename_GWPOUT)
   !-------------------------------------------------------------------------------Y.Kwon 20240714

      !create CHPOUT filename
      filename_CHPOUT = trim(sedCNP_outdir)//'/'//&
                        trim(dateSedCNP%olddate(1:4))//&
                        trim(dateSedCNP%olddate(6:7))//&
                        trim(dateSedCNP%olddate(9:10))//&
                        trim(dateSedCNP%olddate(12:13))//&
                        trim(dateSedCNP%olddate(15:16))//&
                        '.CHPOUT_DOMAIN1'
      call write_Ch_P_output(filename_CHPOUT)

      !=====||__WHQ5403__||=====!
      !create LAKESEDOUT, LAKECOUT, LAKENOUT, LAKEPOUT filenames (lakes/reservoirs: lake water balance
      !and the CHSEDOUT/CHCOUT/CHNOUT/CHPOUT variables at the lake outlet links)
      if (lake%active) then
         filename_LAKE = trim(sedCNP_outdir)//'/'//&
                         trim(dateSedCNP%olddate(1:4))//&
                         trim(dateSedCNP%olddate(6:7))//&
                         trim(dateSedCNP%olddate(9:10))//&
                         trim(dateSedCNP%olddate(12:13))//&
                         trim(dateSedCNP%olddate(15:16))
         call write_Lake_Sed_output(trim(filename_LAKE)//'.LAKESEDOUT_DOMAIN1')
         call write_Lake_C_output(trim(filename_LAKE)//'.LAKECOUT_DOMAIN1')
         call write_Lake_N_output(trim(filename_LAKE)//'.LAKENOUT_DOMAIN1')
         call write_Lake_P_output(trim(filename_LAKE)//'.LAKEPOUT_DOMAIN1')
      endif
      !=====||__WHQ5403__||=====!

      !get newdate for next timestep
      call geth_newdate(dateSedCNP%newdate, dateSedCNP%olddate, nint(dt))
      dateSedCNP%olddate = dateSedCNP%newdate

      !call write_C_output(itime)  !BK                             !BK20240703
      !call write_N_output(itime)  !#BK                            !BK20240703
      !call write_P_output(itime)  !BK   !kyh_debug  !BK20230311   !BK20240703

   end subroutine SedCNPmodel_output


   subroutine SedCNP_driver_exe(itime)

!=====||__WHQ5403__||=====!
!
      ! SedCNPmodel is already visible via ReadCNPini/WriteCNPini/CNPmain/module_SedCNP_in
      ! below. A direct "use config_base" here, combined with any sibling module that
      ! itself uses config_base, triggers a gfortran diamond-import derived-type mismatch.
      !use config_base,           only: SedCNPmodel
!
!=====||__WHQ5403__||=====!
      use module_SedCNPvariables
      use module_overlandSed,    only: overlandSedTransport
      use module_channelSed,     only: channelSedTransport
      use ReadCNPini
      use WriteCNPini
      use CNPmain
      use module_SedCNP_in
      !use CNPvariables
      use CNPparams      !Y.Kwon 20230312
      use module_AlloDeallocate
      use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT, int64
      use hashtable, only: hash_t   !WHQ5403 link ID -> channel index mapping
      use module_lakeSedCNP, only: Lake_Init, Lake_UpdateHydro   !WHQ5403 lakes/reservoirs
      use module_lateralTransport, only: Lateral_Transport_Init, Lateral_Loads   !WHQ5403 lateral loads to channels
#ifdef MPP_LAND
      use module_mpp_land, only: mpp_land_sync, my_id, io_id
#endif


      implicit none

      integer, intent(in) :: itime  ! timestep loop
      integer             :: i, j, k
      !integer             :: sbid   !subbasin id    !BK20231010
      integer             :: gwid   !gw basin id    !BK20240613
      !integer             :: ibasin    !index for gw basins  !BK20301030
      integer             :: chid   !channel link id
      integer             :: ich    !index for channels  !BK20231029
      real                :: prg    !BK20230402
      logical             :: first_call  !BK20231114


      !show the current simulation timestep on the screen   !BK20230402 !BK20250616
      prg = real(itime) / real(domain%ntime_sedcnp) * 100.0
#ifdef MPP_LAND
      if (my_id .eq. io_id) then
#endif
         write(ERROR_UNIT, '(A,F6.2,A,I0,A,I0,A)', advance='no') &
               " WRF-HydroQual running (water quality) ... ", &
               prg, "% (", itime, " / ", domain%ntime_sedcnp, ")"
         if (itime == domain%ntime_sedcnp) then  !---BK20250807
            ! On the final timestep, print a newline, then the success message.
            write(ERROR_UNIT,*)
            write(ERROR_UNIT,*) "The model finished successfully"
         else
            write(ERROR_UNIT, '(A)', advance='no') char(13)
         endif
         call flush(ERROR_UNIT)    !---BK20250807
#ifdef MPP_LAND
      endif
#endif
      
      if(SedCNPmodel%SedCNP_option == 1) then
         !sediment only; need water outputs
         !will be developed later
         !stop "FATAL ERROR: SedCNP_option = 1 does not work in the current version"
         if (my_id .eq. io_id) then  !BK20250803
            write(ERROR_UNIT,*) "WARNING: SedCNP_option = 1 is not yet implemented. Skipping."
         endif
         return
      elseif(SedCNPmodel%SedCNP_option == 2) then
         !CNP only; need water and sediment outputs
         !will be developed later
         !stop "FATAL ERROR: SedCNP_option = 2 does not work in the current version"
         if (my_id .eq. io_id) then  !BK20250803
            write(ERROR_UNIT,*) "WARNING: SedCNP_option = 2 is not yet implemented. Skipping."
         endif
         return            
      elseif(SedCNPmodel%SedCNP_option == 3) then

         if (itime == 1) then
            call Allocate_Temp()   !20231018   !BK20231017 better move to SedCNP_driver_ini???

            !domain%areaxy = SedCNPmodel%SedCNP_dx * SedCNPmodel%SedCNP_dx  !BK20230319 !Y.Kwon
            
            !gwbasin ID for rt grid cells   BK20240616
            SedCNP_hydro%gwbasin = SedCNP_hydro%subbasinID

            !gwbasin_ch: gwbasin IDs for channel grids   !BK20240616
            SedCNP_hydro%gwbasin_ch = -9999
            do j = 1, domain%jxrt
               do i = 1, domain%ixrt
                  if(SedCNP_hydro%linkID_grid(i,j) .ge. 1.) then
                     SedCNP_hydro%gwbasin_ch(i,j) = SedCNP_hydro%gwbasin(i,j)  !! mapping
                  endif
               end do
            end do
         endif

         !sediment + CNP; need water outputs
         !get water outputs
         call SedCNPmodel_hydro_input(itime)

         if (itime == 1) then
            !--- BK20250701
            ! create mapping from channel index (ich) to a single, (firstly-identified) 'representative' grid
            ! coordinates (i,j) and groundwater basin ID (gwid).
            ! Initialize with -9999 to identify unmapped links later
            !=====||__WHQ5403__||=====!
            ! Link IDs in Route_Link.nc need not be contiguous (1..nch), so they cannot be used as
            ! array indices. As in the hydro routing (module_RT.F90), grid LINKID values are converted
            ! to the channel index (position in Route_Link.nc) via a hash table of link IDs.
            ! Must follow SedCNPmodel_hydro_input, which reads SedCNP_hydro%linkID.
            !=====||__WHQ5403__||=====!
            domain%chid_i = -9999
            domain%chid_j = -9999
            domain%gwid_ch = -9999

            !=====||__WHQ5403__||=====! original code (link ID used as array index; moved from before SedCNPmodel_hydro_input)
            !do j = 1, domain%jxrt
            !   do i = 1, domain%ixrt
            !      if(SedCNP_hydro%linkID_grid(i,j) >= 1) then
            !         chid = SedCNP_hydro%linkID_grid(i,j)
            !         !--- BK20250701
            !         if ((chid > 0) .and. &
            !            (chid <= domain%nch) .and. &
            !            (domain%chid_i(chid) == -9999)) then
            !            domain%gwid_ch(chid) = SedCNP_hydro%gwbasin(i,j)
            !            domain%chid_i(chid) = i
            !            domain%chid_j(chid) = j
            !         endif
            !         !--- BK20250701
            !      endif
            !   end do
            !end do
            !=====||__WHQ5403__||=====!

            block
               type(hash_t) :: hash_table
               integer(kind=int64) :: val
               logical :: found

               call hash_table%set_all_idx(int(SedCNP_hydro%linkID, int64), domain%nch)
               do j = 1, domain%jxrt
                  do i = 1, domain%ixrt
                     if(SedCNP_hydro%linkID_grid(i,j) >= 1) then
                        call hash_table%get(int(SedCNP_hydro%linkID_grid(i,j), int64), val, found)
                        if (found) then
                           ich = int(val)
                           if (domain%chid_i(ich) == -9999) then
                              domain%gwid_ch(ich) = SedCNP_hydro%gwbasin(i,j)
                              domain%chid_i(ich) = i
                              domain%chid_j(ich) = j
                           endif
                        endif
                     endif
                  end do
               end do
               call hash_table%clear()
            end block
         endif

         !=====||__WHQ5403__||=====!
         ! lakes/reservoirs: identify lake links (first time-step) and set the lake water
         ! balance and the outlet-link hydraulics (every time-step, after the hydro inputs)
         if (itime == 1) call Lake_Init()
         if (itime == 1) call Lateral_Map_Init()   ! lateral inflows from gw basins to channels (after Lake_Init)
         if (itime == 1) call Lateral_Transport_Init()
         call Lake_UpdateHydro(itime)
         !=====||__WHQ5403__||=====!
         
         ! Calculate basin-wide averages once per timestep. This is called after
         ! hydro inputs are updated and before any CNP processing for the timestep.
         call Calculate_Basin_Averages()

         !initialisation ------ better move to SedCNP_driver_ini !BK20250917 restored
         if (itime == 1) then    !------------- BK20241013 
            !--- Sed
         
            !Soil
            So_LPOCLIT = 0.0
            So_LPOCRES = 0.0
            So_LPOCEXC = 0.0
            So_LPOCMAN = 0.0
            So_RPOC = 0.0
            So_LDOC = 0.0
            So_RDOC = 0.0
            So_MBMC = 0.0
            So_CO2 = 0.0
            So_LPOCLITinso = 0.0
            So_LPOCRESinso = 0.0
            So_LPOCEXCinso = 0.0
            So_LPOCMANinso = 0.0
            So_LPOCLITin = 0.0
            So_LPOCRESin = 0.0
            So_LPOCEXCin = 0.0
            So_LPOCMANin = 0.0
            So_LDOCperc = 0.0
            So_RDOCperc = 0.0
            So_LPOCsurf = 0.0
            So_RPOCsurf = 0.0
            So_LDOCsurf = 0.0
            So_RDOCsurf = 0.0
            So_MBMCsurf = 0.0
            So_LDOCintf = 0.0
            So_RDOCintf = 0.0
            So_LDOCsogw = 0.0  !BK20260316
            So_RDOCsogw = 0.0  !BK20260316
            So_LDOCgwso = 0.0  !BK20260316
            So_RDOCgwso = 0.0  !BK20260316
            So_CSC = 0.0
            So_CMBError = 0.0
                  
            !Groundwater
            Gw_LDOC = 0.0
            Gw_RDOC = 0.0
            Gw_LDOCgwso = 0.0  !BK20260318
            Gw_RDOCgwso = 0.0  !BK20260318              
            Gw_LDOCgwch = 0.0
            Gw_RDOCgwch = 0.0
            Gw_CSC = 0.0
            Gw_CMBError = 0.0
         
            !Channel
            Ch_ALGC = 0.0
            Ch_ZOOC = 0.0
            Ch_MBMC = 0.0
            Ch_DIC = 0.0
            Ch_CO2EVAS = 0.0  !BK20251126
            Ch_LPOC = 0.0
            Ch_RPOC = 0.0
            Ch_LDOC = 0.0
            Ch_RDOC = 0.0
            Ch_POCDEPOSIT = 0.0
            Ch_ALGCG = 0.0
            Ch_ALGCR = 0.0
            Ch_ALGCE = 0.0
            Ch_ALGCM = 0.0
            Ch_ALGCZ = 0.0
            Ch_ALGCS = 0.0
            Ch_ZOOCG = 0.0
            Ch_ZOOCR = 0.0
            Ch_ZOOCE = 0.0
            Ch_ZOOCM = 0.0
            Ch_ZOOCZ = 0.0
            Ch_ZOOCS = 0.0
            Ch_ALGCdsch = 0.0
            Ch_ZOOCdsch = 0.0
            Ch_MBMCdsch = 0.0
            Ch_DICdsch = 0.0  !BK20251121
            Ch_LPOCdsch = 0.0
            Ch_RPOCdsch = 0.0
            Ch_LDOCdsch = 0.0
            Ch_RDOCdsch = 0.0
            Ch_LPOCpnt = 0.0
            Ch_RPOCpnt = 0.0
            Ch_LDOCpnt = 0.0
            Ch_RDOCpnt = 0.0
            Ch_LPOCxpnt = 0.0
            Ch_RPOCxpnt = 0.0
            Ch_LDOCxpnt = 0.0
            Ch_RDOCxpnt = 0.0
            Ch_LPOCabs = 0.0
            Ch_RPOCabs = 0.0
            Ch_LDOCabs = 0.0
            Ch_RDOCabs = 0.0
            Ch_LPOCxabs = 0.0
            Ch_RPOCxabs = 0.0
            Ch_LDOCxabs = 0.0
            Ch_RDOCxabs = 0.0
            Ch_LPOCdis = 0.0
            Ch_RPOCdis = 0.0
            Ch_LDOCdis = 0.0
            Ch_RDOCdis = 0.0
            Ch_LPOCxdis = 0.0
            Ch_RPOCxdis = 0.0
            Ch_LDOCxdis = 0.0
            Ch_RDOCxdis = 0.0
            Ch_CSC = 0.0
            Ch_CMBError = 0.0
         
            !--- N
            !Soil
            So_LPONLIT = 0.0
            So_LPONRES = 0.0
            So_LPONEXC = 0.0
            So_LPONMAN = 0.0
            So_RPON = 0.0
            So_LDON = 0.0
            So_RDON = 0.0
            So_MBMN = 0.0
            So_NH4 = 0.0
            So_NO3 = 0.0
            So_LPONLITinso = 0.0
            So_LPONRESinso = 0.0
            So_LPONEXCinso = 0.0
            So_LPONMANinso = 0.0
            So_LPONLITin = 0.0
            So_LPONRESin = 0.0
            So_LPONEXCin = 0.0
            So_LPONMANin = 0.0
            So_NH4fert = 0.0
            So_NO3fert = 0.0
            So_NH4fxin = 0.0
            So_NO3fxin = 0.0
            So_NITRI = 0.0
            So_NH3VOL = 0.0
            So_DENIT = 0.0
            So_NH4uptk = 0.0
            So_NO3uptk = 0.0
            So_LDONperc = 0.0
            So_RDONperc = 0.0
            So_NH4perc = 0.0
            So_NO3perc = 0.0
            So_LPONsurf = 0.0
            So_RPONsurf = 0.0
            So_LDONsurf = 0.0
            So_RDONsurf = 0.0
            So_MBMNsurf = 0.0
            So_NH4surf = 0.0
            So_NO3surf = 0.0
            So_LDONintf = 0.0
            So_RDONintf = 0.0
            So_NH4intf = 0.0
            So_NO3intf = 0.0
            So_LDONsogw = 0.0  !BK20260316
            So_RDONsogw = 0.0  !BK20260316
            So_NH4sogw = 0.0  !BK20260316
            So_NO3sogw = 0.0  !BK20260316
            So_LDONgwso = 0.0  !BK20260316
            So_RDONgwso = 0.0  !BK20260316
            So_NH4gwso = 0.0  !BK20260316
            So_NO3gwso = 0.0  !BK20260316
            So_NSC = 0.0
            So_NMBError = 0.0

            !Groundwater
            Gw_LDON = 0.0
            Gw_RDON = 0.0
            Gw_NH4 = 0.0
            Gw_NO3 = 0.0
            Gw_LDONgwso = 0.0  !BK20260318
            Gw_RDONgwso = 0.0  !BK20260318
            Gw_NH4gwso = 0.0  !BK20260318
            Gw_NO3gwso = 0.0  !BK20260318
            Gw_LDONgwch = 0.0
            Gw_RDONgwch = 0.0
            Gw_NH4gwch = 0.0
            Gw_NO3gwch = 0.0
            Gw_NSC = 0.0
            Gw_NMBError = 0.0

            !Channel
            Ch_ALGN = 0.0
            Ch_ZOON = 0.0
            Ch_MBMN = 0.0
            Ch_LPON = 0.0
            Ch_RPON = 0.0
            Ch_LDON = 0.0
            Ch_RDON = 0.0
            Ch_NH4 = 0.0
            Ch_NO3 = 0.0
            Ch_PONDEPOSIT = 0.0
            Ch_ALGNG = 0.0
            Ch_ALGNR = 0.0
            Ch_ALGNE = 0.0
            Ch_ALGNM = 0.0
            Ch_ALGNZ = 0.0
            Ch_ALGNS = 0.0
            Ch_ZOONG = 0.0
            Ch_ZOONR = 0.0
            Ch_ZOONE = 0.0
            Ch_ZOONM = 0.0
            Ch_ZOONZ = 0.0
            Ch_ZOONS = 0.0
            Ch_NITRI = 0.0
            Ch_DENIT = 0.0
            Ch_NH3VOL = 0.0
            Ch_NH4uptk = 0.0
            Ch_NO3uptk = 0.0
            Ch_ALGNdsch = 0.0
            Ch_ZOONdsch = 0.0
            Ch_MBMNdsch = 0.0
            Ch_LPONdsch = 0.0
            Ch_RPONdsch = 0.0
            Ch_LDONdsch = 0.0
            Ch_RDONdsch = 0.0
            Ch_NH4dsch = 0.0
            Ch_NO3dsch = 0.0
            Ch_LPONpnt = 0.0
            Ch_RPONpnt = 0.0
            Ch_LDONpnt = 0.0
            Ch_RDONpnt = 0.0
            Ch_NH4pnt = 0.0
            Ch_NO3pnt = 0.0
            Ch_LPONxpnt = 0.0
            Ch_RPONxpnt = 0.0
            Ch_LDONxpnt = 0.0
            Ch_RDONxpnt = 0.0
            Ch_NH4xpnt = 0.0
            Ch_NO3xpnt = 0.0
            Ch_LPONabs = 0.0
            Ch_RPONabs = 0.0
            Ch_LDONabs = 0.0
            Ch_RDONabs = 0.0
            Ch_NH4abs = 0.0
            Ch_NO3abs = 0.0
            Ch_LPONxabs = 0.0
            Ch_RPONxabs = 0.0
            Ch_LDONxabs = 0.0
            Ch_RDONxabs = 0.0
            Ch_NH4xabs = 0.0
            Ch_NO3xabs = 0.0
            Ch_LPONdis = 0.0
            Ch_RPONdis = 0.0
            Ch_LDONdis = 0.0
            Ch_RDONdis = 0.0
            Ch_NH4dis = 0.0
            Ch_NO3dis = 0.0
            Ch_LPONxdis = 0.0
            Ch_RPONxdis = 0.0
            Ch_LDONxdis = 0.0
            Ch_RDONxdis = 0.0
            Ch_NH4xdis = 0.0
            Ch_NO3xdis = 0.0
            Ch_NSC = 0.0
            Ch_NMBError = 0.0

            !--- P
            !Soil
            So_LPOPLIT = 0.0
            So_LPOPRES = 0.0
            So_LPOPEXC = 0.0
            So_LPOPMAN = 0.0
            So_RPOP = 0.0
            So_LDOP = 0.0
            So_RDOP = 0.0
            So_MBMP = 0.0
            So_PIPA = 0.0
            So_PIPS = 0.0
            So_PO4 = 0.0
            So_LPOPLITinso = 0.0
            So_LPOPRESinso = 0.0
            So_LPOPEXCinso = 0.0
            So_LPOPMANinso = 0.0
            So_LPOPLITin = 0.0
            So_LPOPRESin = 0.0
            So_LPOPEXCin = 0.0
            So_LPOPMANin = 0.0
            So_PO4fert = 0.0
            So_PO4fxin = 0.0
            So_PO4uptk = 0.0
            So_PO4ADS = 0.0
            So_PIPAADS = 0.0
            So_LDOPperc = 0.0
            So_RDOPperc = 0.0
            So_PO4perc = 0.0
            So_LPOPsurf = 0.0
            So_RPOPsurf = 0.0
            So_LDOPsurf = 0.0
            So_RDOPsurf = 0.0
            So_MBMPsurf = 0.0
            So_PIPAsurf = 0.0
            So_PIPSsurf = 0.0
            So_PO4surf = 0.0
            So_LDOPintf = 0.0
            So_RDOPintf = 0.0
            So_PO4intf = 0.0
            So_LDOPsogw = 0.0  !BK20260316
            So_RDOPsogw = 0.0  !BK20260316
            So_PO4sogw = 0.0   !BK20260316
            So_LDOPgwso = 0.0  !BK20260316
            So_RDOPgwso = 0.0  !BK20260316
            So_PO4gwso = 0.0   !BK20260316
            So_PSC = 0.0
            So_PMBError = 0.0

            !Groundwater
            Gw_LDOP = 0.0
            Gw_RDOP = 0.0
            Gw_PO4 = 0.0
            Gw_LDOPgwso = 0.0  !BK20260318
            Gw_RDOPgwso = 0.0  !BK20260318
            Gw_PO4gwso = 0.0  !BK20260318
            Gw_LDOPgwch = 0.0
            Gw_RDOPgwch = 0.0
            Gw_PO4gwch = 0.0
            Gw_PSC = 0.0
            Gw_PMBError = 0.0

            !Channel
            Ch_ALGP = 0.0
            Ch_ZOOP = 0.0
            Ch_MBMP = 0.0
            Ch_LPOP = 0.0
            Ch_RPOP = 0.0
            Ch_LDOP = 0.0
            Ch_RDOP = 0.0
            Ch_PIPA = 0.0
            Ch_PIPS = 0.0
            Ch_PO4 = 0.0
            Ch_POPDEPOSIT = 0.0
            Ch_ALGPG = 0.0
            Ch_ALGPR = 0.0
            Ch_ALGPE = 0.0
            Ch_ALGPM = 0.0
            Ch_ALGPZ = 0.0
            Ch_ALGPS = 0.0
            Ch_ZOOPG = 0.0
            Ch_ZOOPR = 0.0
            Ch_ZOOPE = 0.0
            Ch_ZOOPM = 0.0
            Ch_ZOOPZ = 0.0
            Ch_ZOOPS = 0.0
            Ch_PO4ADS = 0.0
            Ch_PIPAADS = 0.0
            Ch_PO4uptk = 0.0
            Ch_ALGPdsch = 0.0
            Ch_ZOOPdsch = 0.0
            Ch_MBMPdsch = 0.0
            Ch_LPOPdsch = 0.0
            Ch_RPOPdsch = 0.0
            Ch_LDOPdsch = 0.0
            Ch_RDOPdsch = 0.0
            Ch_PIPAdsch = 0.0
            Ch_PIPSdsch = 0.0
            Ch_PO4dsch = 0.0
            Ch_LPOPpnt = 0.0
            Ch_RPOPpnt = 0.0
            Ch_LDOPpnt = 0.0
            Ch_RDOPpnt = 0.0
            Ch_PO4pnt = 0.0
            Ch_LPOPxpnt = 0.0
            Ch_RPOPxpnt = 0.0
            Ch_LDOPxpnt = 0.0
            Ch_RDOPxpnt = 0.0
            Ch_PO4xpnt = 0.0
            Ch_LPOPabs = 0.0
            Ch_RPOPabs = 0.0
            Ch_LDOPabs = 0.0
            Ch_RDOPabs = 0.0
            Ch_PO4abs = 0.0
            Ch_LPOPxabs = 0.0
            Ch_RPOPxabs = 0.0
            Ch_LDOPxabs = 0.0
            Ch_RDOPxabs = 0.0
            Ch_PO4xabs = 0.0
            Ch_LPOPdis = 0.0
            Ch_RPOPdis = 0.0
            Ch_LDOPdis = 0.0
            Ch_RDOPdis = 0.0
            Ch_PO4dis = 0.0
            Ch_LPOPxdis = 0.0
            Ch_RPOPxdis = 0.0
            Ch_LDOPxdis = 0.0
            Ch_RDOPxdis = 0.0
            Ch_PO4xdis = 0.0
            Ch_PSC = 0.0
            Ch_PMBError = 0.0

         endif                   !------------- BK20241013

         !read in CNP storages and inputs at the beginning of simulation   !BK20230319  
         if (itime == 1) then
            domain%areaxy = SedCNP_hydro%dx(1) * SedCNP_hydro%dx(1)         !Y.Kwon(20250621)
            !areaxy = SedCNPmodel%SedCNP_dx * SedCNPmodel%SedCNP_dx     
            call ReadCNPini_So ()  !BK20230319
            !call ReadCNPini_Aq () !BK20230319
            call ReadCNPini_Gw () !BK20230319
            call ReadCNPini_Ch () !BK20230319

            ! --- initialize current state with initial conditions  
            ! --- to prevent mass balance error at first step   !BK20251216
            ! C
            So_LPOCLIT = So_LPOCLIT0
            So_LPOCRES = So_LPOCRES0
            So_LPOCEXC = So_LPOCEXC0
            So_LPOCMAN = So_LPOCMAN0
            So_RPOC = So_RPOC0
            So_LDOC = So_LDOC0
            So_RDOC = So_RDOC0
            So_MBMC = So_MBMC0
            Gw_LDOC = Gw_LDOC0
            Gw_RDOC = Gw_RDOC0
            Ch_ALGC = Ch_ALGC0
            Ch_ZOOC = Ch_ZOOC0
            Ch_MBMC = Ch_MBMC0
            Ch_DIC  = Ch_DIC0
            Ch_LPOC = Ch_LPOC0
            Ch_RPOC = Ch_RPOC0
            Ch_LDOC = Ch_LDOC0
            Ch_RDOC = Ch_RDOC0
            Ch_POCDEPOSIT = Ch_POCDEPOSIT0
            ! N
            So_LPONLIT = So_LPONLIT0
            So_LPONRES = So_LPONRES0
            So_LPONEXC = So_LPONEXC0
            So_LPONMAN = So_LPONMAN0
            So_RPON = So_RPON0
            So_LDON = So_LDON0
            So_RDON = So_RDON0
            So_MBMN = So_MBMN0
            So_NH4  = So_NH40
            So_NO3  = So_NO30
            Gw_LDON = Gw_LDON0
            Gw_RDON = Gw_RDON0
            Gw_NH4  = Gw_NH40
            Gw_NO3  = Gw_NO30
            Ch_ALGN = Ch_ALGN0
            Ch_ZOON = Ch_ZOON0
            Ch_MBMN = Ch_MBMN0
            Ch_LPON = Ch_LPON0
            Ch_RPON = Ch_RPON0
            Ch_LDON = Ch_LDON0
            Ch_RDON = Ch_RDON0
            Ch_NH4  = Ch_NH40
            Ch_NO3  = Ch_NO30
            Ch_PONDEPOSIT = Ch_PONDEPOSIT0
            ! P
            So_LPOPLIT = So_LPOPLIT0
            So_LPOPRES = So_LPOPRES0
            So_LPOPEXC = So_LPOPEXC0
            So_LPOPMAN = So_LPOPMAN0
            So_RPOP = So_RPOP0
            So_LDOP = So_LDOP0
            So_RDOP = So_RDOP0
            So_MBMP = So_MBMP0
            So_PO4  = So_PO40
            So_PIPA = So_PIPA0
            So_PIPS = So_PIPS0
            Gw_LDOP = Gw_LDOP0
            Gw_RDOP = Gw_RDOP0
            Gw_PO4  = Gw_PO40
            Ch_ALGP = Ch_ALGP0
            Ch_ZOOP = Ch_ZOOP0
            Ch_MBMP = Ch_MBMP0
            Ch_LPOP = Ch_LPOP0
            Ch_RPOP = Ch_RPOP0
            Ch_LDOP = Ch_LDOP0
            Ch_RDOP = Ch_RDOP0
            Ch_PO4  = Ch_PO40
            Ch_PIPA = Ch_PIPA0
            Ch_PIPS = Ch_PIPS0
            Ch_POPDEPOSIT = Ch_POPDEPOSIT0
            ! --- End of initialization !BK20251216


            !KZG = 0.01 !Y.Kwon  !commented out BK20231027

            !call ReadCNPInputs_LC (itime)   !BK20230313
            call ReadCNPinputs () !BK20230313 BK20231024 
            !call ReadCNPInputs_PS (itime)   !kyh_debug(20230911)
            !call ReadCNPInputs_PS ()   !kyh_debug(20230911)   !BK20231016
            call ReadCH_PNTSRC ()   !kyh_debug(20230911)   !BK20231016 !BK20240703
            call ReadCH_ABSDIS ()   !BK20240703
            call ReadCNPparams ()    !Y.Kwon 20230312 
         endif

         !overland
         do i = 1,domain%ix
            do j = 1,domain%jx
               !WHQ5403 skip cells without LSM outputs (fill value -9999 in LDASOUT, e.g. outside the domain)
               if (SedCNP_hydro%vgtyp(i,j) < 0.0) cycle
               !overland sediment
               call overlandSedTransport(i,j)
               !Soil CNP
               call RunCNP_So(i,j,itime)
            enddo
         enddo

         !groundwater  
         !gwid = monotonic subbasin ID as in GWBUCKPARM.nc !BK20240509  
         do gwid = 1,domain%nbasin   !BK20240613 
            call RunCNP_Gw(gwid)  !BK20240613
         enddo

         !WHQ5403 lateral (surface runoff, interflow) loads of this time-step to the channel links
         call Lateral_Loads(itime)

         !channel
         if (itime == 1) then  !BK20231114
            !call ReadSBID_CH()   !BK20231011  !BK20240525
            first_call = .true.

            !--- write ich-chid-gwid table  !BK20250110
            call system('mkdir -p ./debug')
            open (5110, file='./debug/ich_chid_gwid.csv', status='replace')
            write (5110,'(A12,",",A12,",",A12)') 'ich','chid','gwid'
            do ich = 1, domain%nch
               chid = SedCNP_hydro%linkID(ich)
               !gwid = domain%gwid_ch(chid)
               gwid = domain%gwid_ch(ich)   !WHQ5403 indexed by channel index, not link ID
               write (5110,'(I12,",",I12,",",I12)') ich,chid,gwid
            enddo
            close (5110)
            !--- write ich-chid-gwid table  !BK20250110
         else
            first_call = .false.
         endif

         !do ich = 0, domain%nch - 1  !1,domain%nch    !BK20231010  !BK20231029 !BK20240523
         do ich = 1, domain%nch               !BK20240610
            !chid = SedCNP_hydro%linkID(ich)  !BK20240512  !BK20241101
            !sbid = domain%sbid_ch(chid)   !BK20240512
            !!--- skip 'connect-only' channels  !BK20240111
            !if (sbid == 0) then
            !   cycle  ! NO SedCNP simulation for obs'connect-only' channels
            !endif
            !--- BK20240111  

            !--- BK20250701      
            chid = SedCNP_hydro%linkID(ich)
            ! Add a bounds check to prevent out-of-bounds access on chid
            !if (chid < 1 .or. chid > domain%nch) then    !--- BK20250702
            !   cycle
            !endif                                        !--- BK20250702

            ! Check if this channel has a valid grid mapping. If not, it's likely
            ! a "connect-only" link with no spatial representation. Skip it
            ! to prevent out-of-bounds errors downstream.
            !if (domain%chid_i(chid) < 1) then  !--- BK20250924 commented out 
            !   cycle
            !endif                              !--- BK20250924 commented out 
            !--- BK20250701

            !WHQ5403 internal lake links are not simulated (the lake is simulated at its outlet link)
            if (lake%typel(ich) == 2) cycle

            !channel sediment
            call channelSedTransport(ich,itime)    !BK20231016  !BK20231029  !BK20240512 !BK20241101

            !channel CNP
            call RunCNP_Ch(ich,itime,first_call)              !BK20231029  !BK20231114  !BK20240512 !BK20241101
            first_call = .false.

         enddo

         !--- update storages and variables for the next time-step
         !Soil
         do i = 1,domain%ix
            do j = 1,domain%jx
               call UpdateCNP_So(i,j)  !BK20231027
               !call UpdateCNP_Aq(i,j)  !BK20231027
            enddo
         enddo
         !Groundwater
         !do ibasin = 1,domain%nbasin    !BK20240509
         !do sbid = 1,domain%nbasin  !1,domain%nbasin !BK20240108 (sbid = 0 for unintentionally-generated obsolete channels) 
         do gwid = 1,domain%nbasin  !BK20240613 
            !sbid = ibasin  !BK20240509  monotonic subbasin ID as in GWBUCKPARM.nc
            !call UpdateCNP_Gw(sbid)  !BK20231027
            call UpdateCNP_Gw(gwid)  !BK20240613
         enddo
         !Channel
         !do ich = 0, domain%nch - 1  !1,domain%nch    !BK20240523
         do ich = 1, domain%nch  !BK20240613
            !chid = SedCNP_hydro%linkID(ich)  !BK20240512 !BK20241101
            if (lake%typel(ich) == 2) cycle   !WHQ5403 internal lake links are not simulated
            call UpdateCNP_Ch(ich)
         enddo  

         !write output
         call SedCNPmodel_output(itime) 

         !write CNP storages at the end of simulation   !BK20230316
         if (itime == domain%ntime_sedcnp) then   !Y.Kwon 20230316
            !domain%areaxy = SedCNPmodel%SedCNP_dx * SedCNPmodel%SedCNP_dx  !BK20230319
            domain%areaxy = SedCNP_hydro%dx(1) * SedCNP_hydro%dx(1)        !Y.Kwon(20250621)
            call WriteCNPini_So()
            !call WriteCNPini_Aq()
            call WriteCNPini_Gw()
            call WriteCNPini_Ch()
            call Deallocate_Temp()      !20231018
         endif

      elseif(SedCNPmodel%SedCNP_option == 4) then
         !water + sediment + CNP
         !will be developed later
         !stop "FATAL ERROR: SedCNP_option = 4 does not work in the current version"
         if (my_id .eq. io_id) then  !BK20250803
            write(ERROR_UNIT,*) "WARNING: SedCNP_option = 4 is not yet implemented. Skipping."
         endif
         return            
      endif

   end subroutine SedCNP_driver_exe

!=====||__WHQ5403__||=====!
!
   !> @brief Builds the distribution of the lateral inflows (surface runoff, interflow, groundwater)
   !! from the gw basins to the channels as in the hydro model (simp_gw_buck, module_GW_baseflow):
   !! the flux of gw basin b is divided uniformly among its stream pixels, so link k receives the
   !! fraction w(k,b) = npix(k,b)/npix(b). Stream pixels of internal lake links are counted for the
   !! lake outlet link. Also lists the LSM cells of each gw basin and writes diagnostics.
   !! Called once, after the channel-to-grid mapping and Lake_Init.
   subroutine Lateral_Map_Init()

      use module_SedCNPvariables
      use module_SedCNP_in, only: get2d_int, is_north_up
      use SedCNP_config,    only: SedCNPmodel
      use hashtable,        only: hash_t
      use ISO_FORTRAN_ENV,  only: int64

      implicit none

      integer, allocatable :: npix_b(:)          ! stream pixels of each gw basin
      integer, allocatable :: nbas_k(:)          ! gw basins of each link
      integer, allocatable :: nlnk_b(:)          ! links of each gw basin
      integer, allocatable :: pk(:), pb(:)       ! (link, basin) of each stream pixel (pk = 0: link not found)
      integer, allocatable :: fill(:), msk(:,:)
      integer              :: i, j, k, b, n, m, npix, nnz, status, ncell
      integer              :: nb_nostream, nb_nostream_in, ncell_nostream, nlink_multi, nbas_multi
      integer              :: nlink_nobas, npix_nolink
      logical              :: found_b, has_msk
      real(8)              :: wsum_min, wsum_max, wsum
      type(hash_t)         :: hash_table
      integer(kind=int64)  :: val
      logical              :: found

      ! --- stream pixels: (link, gw basin)
      npix = count(SedCNP_hydro%linkID_grid >= 1)
      allocate(pk(npix), pb(npix), npix_b(domain%nbasin))
      npix_b = 0
      npix_nolink = 0
      call hash_table%set_all_idx(int(SedCNP_hydro%linkID, int64), domain%nch)
      n = 0
      do j = 1, domain%jxrt
         do i = 1, domain%ixrt
            if (SedCNP_hydro%linkID_grid(i,j) >= 1) then
               n = n + 1
               pb(n) = SedCNP_hydro%gwbasin(i,j)   ! AGGFACTRT = 1 (routing grid = LSM grid)
               pk(n) = 0
               call hash_table%get(int(SedCNP_hydro%linkID_grid(i,j), int64), val, found)
               if (found) pk(n) = lake_ich(int(val))
               if (pb(n) >= 1 .and. pb(n) <= domain%nbasin) then
                  npix_b(pb(n)) = npix_b(pb(n)) + 1
                  if (pk(n) == 0) npix_nolink = npix_nolink + 1
               else
                  pb(n) = 0
               endif
            endif
         enddo
      enddo
      call hash_table%clear()

      ! --- CSR: basins and weights of each link
      allocate(nbas_k(domain%nch), nlnk_b(domain%nbasin), lat_ptr(domain%nch+1), fill(domain%nch))
      allocate(lat_gw(npix), lat_w(npix))
      nbas_k = 0
      do n = 1, npix
         if (pk(n) > 0 .and. pb(n) > 0) nbas_k(pk(n)) = nbas_k(pk(n)) + 1   ! upper bound (duplicates)
      enddo
      lat_ptr(1) = 1
      do k = 1, domain%nch
         lat_ptr(k+1) = lat_ptr(k) + nbas_k(k)
      enddo
      fill = 0
      lat_w = 0.0d0
      do n = 1, npix
         k = pk(n); b = pb(n)
         if (k == 0 .or. b == 0) cycle
         found_b = .false.
         do m = lat_ptr(k), lat_ptr(k) + fill(k) - 1
            if (lat_gw(m) == b) then
               lat_w(m) = lat_w(m) + 1.0d0
               found_b = .true.
               exit
            endif
         enddo
         if (.not. found_b) then
            m = lat_ptr(k) + fill(k)
            lat_gw(m) = b
            lat_w(m) = 1.0d0
            fill(k) = fill(k) + 1
         endif
      enddo
      ! compact (remove the duplicate slots) and normalise: w = npix(k,b)/npix(b)
      nnz = 0
      do k = 1, domain%nch
         m = lat_ptr(k)
         lat_ptr(k) = nnz + 1
         do n = m, m + fill(k) - 1
            nnz = nnz + 1
            lat_gw(nnz) = lat_gw(n)
            lat_w(nnz)  = lat_w(n) / real(npix_b(lat_gw(n)), 8)
         enddo
      enddo
      lat_ptr(domain%nch+1) = nnz + 1

      ! --- SedCNP_lateral_option = 0 (previous method): each link takes its own gw basin gwid_ch
      !     (weight 1, whole basin), and a lake outlet link also the gw basins of the internal links
      if (SedCNPmodel%SedCNP_lateral_option == 0) then
         nnz = 0
         do k = 1, domain%nch
            lat_ptr(k) = nnz + 1
            if (lake%active) then
               if (lake%typel(k) == 2) cycle
            endif
            if (domain%gwid_ch(k) >= 1 .and. domain%gwid_ch(k) <= domain%nbasin) then
               nnz = nnz + 1
               lat_gw(nnz) = domain%gwid_ch(k)
               lat_w(nnz)  = 1.0d0
            endif
            if (lake%active) then
               if (lake%typel(k) == 1) then
                  do m = 1, domain%nch
                     if (lake%typel(m) /= 2 .or. lake%lake_of_ich(m) /= lake%lake_of_ich(k)) cycle
                     b = domain%gwid_ch(m)
                     if (b < 1 .or. b > domain%nbasin .or. b == domain%gwid_ch(k)) cycle
                     if (any(lat_gw(lat_ptr(k):nnz) == b)) cycle
                     nnz = nnz + 1
                     lat_gw(nnz) = b
                     lat_w(nnz)  = 1.0d0
                  enddo
               endif
            endif
         enddo
         lat_ptr(domain%nch+1) = nnz + 1
      endif

      ! --- LSM cells of each gw basin
      allocate(bas_ptr(domain%nbasin+1))
      ncell = count(SedCNP_hydro%gwbasin >= 1 .and. SedCNP_hydro%gwbasin <= domain%nbasin)
      allocate(bas_ci(ncell), bas_cj(ncell))
      nlnk_b = 0   ! temporarily: cells per basin
      do j = 1, domain%jx
         do i = 1, domain%ix
            b = SedCNP_hydro%gwbasin(i,j)
            if (b >= 1 .and. b <= domain%nbasin) nlnk_b(b) = nlnk_b(b) + 1
         enddo
      enddo
      bas_ptr(1) = 1
      do b = 1, domain%nbasin
         bas_ptr(b+1) = bas_ptr(b) + nlnk_b(b)
      enddo
      nlnk_b = 0
      do j = 1, domain%jx
         do i = 1, domain%ix
            b = SedCNP_hydro%gwbasin(i,j)
            if (b >= 1 .and. b <= domain%nbasin) then
               m = bas_ptr(b) + nlnk_b(b)
               bas_ci(m) = i
               bas_cj(m) = j
               nlnk_b(b) = nlnk_b(b) + 1
            endif
         enddo
      enddo

      ! --- diagnostics
      allocate(msk(domain%ixrt,domain%jxrt))
      msk = 1
      status = get2d_int("basn_msk", msk, domain%ixrt, domain%jxrt, './DOMAIN/Fulldom_hires.nc')
      has_msk = (status == 0)
      if (has_msk .and. is_north_up('./DOMAIN/Fulldom_hires.nc')) msk = msk(:, domain%jxrt:1:-1)   ! to the model grid (south-up)

      nlnk_b = 0
      do k = 1, domain%nch
         do m = lat_ptr(k), lat_ptr(k+1) - 1
            nlnk_b(lat_gw(m)) = nlnk_b(lat_gw(m)) + 1
         enddo
      enddo
      nlink_multi = 0; nlink_nobas = 0
      do k = 1, domain%nch
         if (lake%active) then
            if (lake%typel(k) == 2) cycle
         endif
         n = lat_ptr(k+1) - lat_ptr(k)
         if (n > 1)  nlink_multi = nlink_multi + 1
         if (n == 0) nlink_nobas = nlink_nobas + 1
      enddo
      nbas_multi = count(nlnk_b > 1)
      nb_nostream = 0; nb_nostream_in = 0; ncell_nostream = 0
      wsum_min = huge(1.0d0); wsum_max = 0.0d0
      do b = 1, domain%nbasin
         if (bas_ptr(b+1) == bas_ptr(b)) cycle   ! basin without cells (e.g. fill row of GWBUCKPARM)
         if (npix_b(b) == 0) then
            nb_nostream = nb_nostream + 1
            ncell_nostream = ncell_nostream + (bas_ptr(b+1) - bas_ptr(b))
            if (has_msk) then
               do m = bas_ptr(b), bas_ptr(b+1) - 1
                  if (msk(bas_ci(m), bas_cj(m)) > 0) then
                     nb_nostream_in = nb_nostream_in + 1
                     write(6,*) 'WARNING: SedCNP gw basin ', b, ' is inside the basin mask but has no stream pixel;', &
                                ' its lateral inflows reach no channel (as in the hydro model)'
                     exit
                  endif
               enddo
            endif
         else
            wsum = 0.0d0
            do k = 1, domain%nch
               do m = lat_ptr(k), lat_ptr(k+1) - 1
                  if (lat_gw(m) == b) wsum = wsum + lat_w(m)
               enddo
            enddo
            wsum_min = min(wsum_min, wsum); wsum_max = max(wsum_max, wsum)
         endif
      enddo

      if (SedCNPmodel%SedCNP_lateral_option == 1) then
         write(6,'(A)')      ' INFO: SedCNP lateral inflows (SedCNP_lateral_option = 1): groundwater from gw basins with the'
         write(6,'(A)')      '       hydro-model weights; surface runoff and interflow along the hydro-model cell paths'
      else
         write(6,'(A)')      ' INFO: SedCNP lateral inflows (SedCNP_lateral_option = 0): sum over the gw basin of each link'
      endif
      write(6,'(A)')         ' INFO: SedCNP gw basins and channel links (stream pixels):'
      write(6,'(A,I8,A,I8)') '   gw basins with cells              : ', count(bas_ptr(2:) > bas_ptr(:domain%nbasin)), &
                             '   channel links: ', domain%nch
      write(6,'(A,I8)')      '   gw basins draining to >1 link     : ', nbas_multi
      write(6,'(A,I8)')      '   links receiving from >1 gw basin  : ', nlink_multi
      write(6,'(A,I8)')      '   links without gw basin            : ', nlink_nobas
      write(6,'(A,I8,A,I8,A)') '   gw basins without stream pixels   : ', nb_nostream, '  (', ncell_nostream, ' cells)'
      if (has_msk) then
         write(6,'(A,I8,A,I8)') '     outside the basin mask (basn_msk): ', nb_nostream - nb_nostream_in, &
                                '   inside: ', nb_nostream_in
      else
         write(6,'(A)')      '     basn_msk not found in Fulldom_hires.nc: mask check skipped'
      endif
      if (npix_nolink > 0) write(6,'(A,I8)') '   stream pixels of links not in Route_Link: ', npix_nolink
      if (wsum_max > 0.0d0 .and. SedCNPmodel%SedCNP_lateral_option == 1) &
         write(6,'(A,F10.6,A,F10.6)') '   sum of groundwater weights per basin (min/max): ', wsum_min, ' / ', wsum_max

      deallocate(pk, pb, npix_b, nbas_k, nlnk_b, fill, msk)

   end subroutine Lateral_Map_Init
!
!=====||__WHQ5403__||=====!

end module module_SedCNPmodel_driver
