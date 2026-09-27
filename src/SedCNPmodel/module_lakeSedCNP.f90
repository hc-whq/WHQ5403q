!=====||__WHQ5403__||=====!
!
!> @brief Sediment and CNP in lakes/reservoirs (level-pool lakes of the hydro model).
!! @details
!! The hydro model (reach-based routing, lake_option = 1) classifies the links of a lake
!! (Route_Link NHDWaterbodyComID = LAKEPARM lake_id) as
!!   TYPEL 1: lake outlet link = the level-pool reservoir (streamflow = lake outflow),
!!   TYPEL 2: internal lake links (not routed; fill values in CHRTOUT),
!!   TYPEL 3: links flowing into the lake (their 'to' is redirected to the outlet link),
!! and the lake inflow is the TYPEL 3 streamflow plus the lateral inflow of all lake links.
!!
!! Here each lake is simulated as a completely mixed reactor by reusing the channel
!! sediment/CNP equations at its outlet link, with
!!   - geometry: vertically-walled square of the lake area (as the level-pool scheme),
!!     i.e. bottom width = length = sqrt(LkArea), vertical side walls,
!!   - water depth = water_sfc_elev (LAKEOUT) - BottomE (LAKEPARM), storage = LkArea * depth,
!!   - discharge = lake outflow (LAKEOUT), velocity = outflow / cross-section (no resuspension),
!!   - upstream inflow from the TYPEL 3 links ('to' redirected to the outlet, as in the hydro model),
!!   - lateral inflows of the internal links (surface runoff, interflow, groundwater, point
!!     sources, water distribution) added to the outlet link.
!! Internal lake links are excluded from the channel calculations.
!! Outputs: LAKESEDOUT, LAKECOUT, LAKENOUT, LAKEPOUT (lake water balance and the variables
!! of CHSEDOUT, CHCOUT, CHNOUT, CHPOUT at the lake outlet links).
!
module module_lakeSedCNP

   use module_SedCNPvariables
   use SedCNP_config,     only: SedCNPmodel
   use module_SedCNP_in,  only: get1d_ch_real, get1d_ch_int, get_feature_id
   use module_hydro_stop, only: hydro_stop

   implicit none

   character(len=*), parameter :: filename_lakeparm  = './DOMAIN/LAKEPARM.nc'
   character(len=*), parameter :: filename_routelink = './DOMAIN/Route_Link.nc'

   contains

   !> @brief Identifies lakes and lake links and sets the outlet-link geometry.
   !! Called once at the first time-step, after SedCNPmodel_hydro_input (link IDs) and the
   !! channel-to-grid mapping (domain%gwid_ch).
   subroutine Lake_Init()

      implicit none

      integer, allocatable :: wb(:)        ! NHDWaterbodyComID of each link
      integer, allocatable :: ich_dwn(:)   ! channel index of the downstream link (0: none)
      integer, allocatable :: nout(:)      ! number of outlet links of each lake
      real(8), allocatable :: buf(:)
      logical              :: file_exists
      integer              :: ich, k, il, status, gw_out, gw_int
      integer              :: nint_lk, ninf_lk

      ! --- default: no lakes
      lake%active = .false.
      lake%nlake  = 0
      allocate(lake%typel(domain%nch), lake%lake_of_ich(domain%nch))
      lake%typel       = 0
      lake%lake_of_ich = 0
      allocate(lake%gwbasin_lat(size(SedCNP_hydro%gwbasin,1), size(SedCNP_hydro%gwbasin,2)))
      lake%gwbasin_lat = SedCNP_hydro%gwbasin
      allocate(lake%gw_lake(domain%nbasin))
      lake%gw_lake = 0

      inquire(file=filename_lakeparm, exist=file_exists)
      if (.not. file_exists) then
         write(6,*) 'INFO: ', filename_lakeparm, ' not found. Lakes are not simulated in SedCNP.'
         return
      endif

      ! --- lake parameters
      call get_feature_id(lake%nlake, filename_lakeparm)
      if (lake%nlake <= 0) return
      allocate(lake%lake_id(lake%nlake), lake%outlet_ich(lake%nlake), &
               lake%area(lake%nlake), lake%bottomE(lake%nlake), &
               lake%wse(lake%nlake), lake%depth(lake%nlake), lake%storage(lake%nlake), &
               lake%inflow(lake%nlake), lake%outflow(lake%nlake), buf(lake%nlake), nout(lake%nlake))
      status = get1d_ch_int("lake_id", lake%lake_id, lake%nlake, filename_lakeparm)
      if (status /= 0) call hydro_stop("Lake_Init: failed to read lake_id in "//filename_lakeparm)
      status = get1d_ch_real("LkArea", buf, lake%nlake, filename_lakeparm)
      if (status /= 0) call hydro_stop("Lake_Init: failed to read LkArea in "//filename_lakeparm)
      lake%area = buf * 1.0d6   ! km2 -> m2
      status = get1d_ch_real("BottomE", lake%bottomE, lake%nlake, filename_lakeparm)
      if (status /= 0) call hydro_stop("Lake_Init: BottomE (lake bottom elevation [m]) is required in " &
                                       //filename_lakeparm//" to simulate lakes in SedCNP")

      ! --- lake links (same classification as nhdLakeMap in module_HYDRO_io.F90)
      allocate(wb(domain%nch), ich_dwn(domain%nch))
      status = get1d_ch_int("NHDWaterbodyComID", wb, domain%nch, filename_routelink)
      if (status /= 0) call hydro_stop("Lake_Init: failed to read NHDWaterbodyComID in "//filename_routelink)

      ich_dwn = 0
      do ich = 1, domain%nch
         do k = 1, domain%nch
            if (SedCNP_hydro%linkID(k) == SedCNP_hydro%dwnstrm_linkID(ich)) then
               ich_dwn(ich) = k
               exit
            endif
         enddo
      enddo

      do ich = 1, domain%nch
         if (wb(ich) > 0) then
            il = lake_index(wb(ich))
            if (il == 0) cycle   ! waterbody not in LAKEPARM: routed as a channel
            lake%lake_of_ich(ich) = il
            lake%typel(ich) = 2
            if (ich_dwn(ich) == 0) then
               lake%typel(ich) = 1
            elseif (wb(ich_dwn(ich)) /= wb(ich)) then
               lake%typel(ich) = 1
            endif
         endif
      enddo
      do ich = 1, domain%nch
         if (lake%typel(ich) == 0 .and. ich_dwn(ich) > 0) then
            if (lake%typel(ich_dwn(ich)) > 0) then
               lake%typel(ich) = 3
               lake%lake_of_ich(ich) = lake%lake_of_ich(ich_dwn(ich))
            endif
         endif
      enddo

      nout = 0
      lake%outlet_ich = 0
      do ich = 1, domain%nch
         if (lake%typel(ich) == 1) then
            il = lake%lake_of_ich(ich)
            nout(il) = nout(il) + 1
            lake%outlet_ich(il) = ich
         endif
      enddo
      do il = 1, lake%nlake
         if (nout(il) /= 1) then
            write(6,*) 'ERROR: Lake_Init: lake ', lake%lake_id(il), ' has ', nout(il), ' outlet links'
            call hydro_stop("Lake_Init: each lake must have exactly one outlet link")
         endif
      enddo

      ! --- outlet link: vertically-walled square of the lake area
      do il = 1, lake%nlake
         ich = lake%outlet_ich(il)
         channel_input%Wst(ich)        = sqrt(lake%area(il))   ! bottom width [m]
         channel_input%Lst(ich)        = sqrt(lake%area(il))   ! length [m]
         channel_input%SideSlopch(ich) = 1.0d6                 ! vertical side walls
      enddo

      ! --- lateral inflows of internal links -> lake outlet
      do ich = 1, domain%nch
         if (lake%typel(ich) == 2) then
            il = lake%lake_of_ich(ich)
            gw_out = domain%gwid_ch(lake%outlet_ich(il))
            gw_int = domain%gwid_ch(ich)
            if (gw_int > 0 .and. gw_int /= gw_out) then
               if (gw_out > 0) then
                  where (SedCNP_hydro%gwbasin == gw_int) lake%gwbasin_lat = gw_out
                  lake%gw_lake(gw_int) = il
               else
                  write(6,*) 'WARNING: Lake_Init: outlet link of lake ', lake%lake_id(il), &
                             ' has no gw basin; lateral inflows to internal link ', &
                             SedCNP_hydro%linkID(ich), ' are not simulated'
               endif
            endif
         endif
      enddo

      lake%wse = 0.0;    lake%depth = 0.0;   lake%storage = 0.0
      lake%inflow = 0.0; lake%outflow = 0.0

      lake%active = .true.

      write(6,*) 'INFO: SedCNP lakes: ', lake%nlake
      do il = 1, lake%nlake
         nint_lk = count(lake%typel == 2 .and. lake%lake_of_ich == il)
         ninf_lk = count(lake%typel == 3 .and. lake%lake_of_ich == il)
         write(6,'(A,I10,A,I10,A,F12.1,A,F10.3,A,I4,A,I4)') '   lake ', lake%lake_id(il), &
               '  outlet link ', SedCNP_hydro%linkID(lake%outlet_ich(il)), &
               '  area [m2] ', lake%area(il), '  BottomE [m] ', lake%bottomE(il), &
               '  internal links ', nint_lk, '  inflow links ', ninf_lk
      enddo

      deallocate(wb, ich_dwn, buf, nout)

   end subroutine Lake_Init


   !> @brief Sets the lake water balance and the outlet-link hydraulics for the time-step.
   !! Called every time-step after SedCNPmodel_hydro_input (which re-reads CHRTOUT and 'to').
   subroutine Lake_UpdateHydro(itime)

      implicit none

      integer, intent(in)  :: itime
      character(len=256)   :: filename_LAKEOUT
      integer, allocatable :: fid(:)
      real(8), allocatable :: wse_in(:), qin_in(:), qout_in(:)
      integer              :: n, k, il, ich, o, status
      logical              :: file_exists

      if (.not. lake%active) return

      filename_LAKEOUT = trim(SedCNPmodel%LDASOUT_dir)//'/'//&
                         trim(dateSedCNP%olddate(1:4))//&
                         trim(dateSedCNP%olddate(6:7))//&
                         trim(dateSedCNP%olddate(9:10))//&
                         trim(dateSedCNP%olddate(12:13))//&
                         trim(dateSedCNP%olddate(15:16))//&
                         '.LAKEOUT_DOMAIN1'
      inquire(file=trim(filename_LAKEOUT), exist=file_exists)
      if (.not. file_exists) then
         write(6,*) 'ERROR: Lake_UpdateHydro: ', trim(filename_LAKEOUT), ' not found'
         call hydro_stop("Lake_UpdateHydro: LAKEOUT files are required to simulate lakes (outlake = 1 in hydro.namelist)")
      endif

      n = 0
      call get_feature_id(n, trim(filename_LAKEOUT))
      allocate(fid(n), wse_in(n), qin_in(n), qout_in(n))
      status = get1d_ch_int("feature_id", fid, n, trim(filename_LAKEOUT))
      status = status + abs(get1d_ch_real("water_sfc_elev", wse_in, n, trim(filename_LAKEOUT)))
      status = status + abs(get1d_ch_real("inflow", qin_in, n, trim(filename_LAKEOUT)))
      status = status + abs(get1d_ch_real("outflow", qout_in, n, trim(filename_LAKEOUT)))
      if (status /= 0) call hydro_stop("Lake_UpdateHydro: failed to read "//trim(filename_LAKEOUT))

      do k = 1, n
         il = lake_index(fid(k))
         if (il == 0) cycle
         lake%wse(il)     = wse_in(k)
         lake%inflow(il)  = max(0.0d0, qin_in(k))
         lake%outflow(il) = max(0.0d0, qout_in(k))
      enddo
      deallocate(fid, wse_in, qin_in, qout_in)

      do il = 1, lake%nlake
         lake%depth(il)   = max(0.0d0, lake%wse(il) - lake%bottomE(il))
         lake%storage(il) = lake%area(il) * lake%depth(il)
         ! outlet link = lake (head = mean depth; CHRTOUT has Head = -999, velocity = 0 for the outlet)
         o = lake%outlet_ich(il)
         SedCNP_hydro%head(o)   = lake%depth(il)
         SedCNP_hydro%q_dsch(o) = lake%outflow(il)
         if (lake%depth(il) > 1.0d-6) then
            SedCNP_hydro%chVa(o) = lake%outflow(il) / (sqrt(lake%area(il)) * lake%depth(il))
         else
            SedCNP_hydro%chVa(o) = 0.0
         endif
      enddo

      do ich = 1, domain%nch
         if (lake%typel(ich) == 3) then
            ! inflow links drain to the lake outlet link (as TO_NODE in the hydro model)
            SedCNP_hydro%dwnstrm_linkID(ich) = SedCNP_hydro%linkID(lake%outlet_ich(lake%lake_of_ich(ich)))
         elseif (lake%typel(ich) == 2) then
            ! internal lake links are not simulated
            SedCNP_hydro%dwnstrm_linkID(ich) = -9999
            SedCNP_hydro%q_dsch(ich) = 0.0
            SedCNP_hydro%chVa(ich)   = 0.0
            ! water abstraction from an internal link -> from the lake
            o = lake%outlet_ich(lake%lake_of_ich(ich))
            SedCNP_hydro%Qabs(o,itime)   = SedCNP_hydro%Qabs(o,itime) + SedCNP_hydro%Qabs(ich,itime)
            SedCNP_hydro%Qabs(ich,itime) = 0.0
         endif
      enddo

   end subroutine Lake_UpdateHydro


   !> @brief Writes the lake sediment output file (LAKESEDOUT): lake water balance and the
   !! sediment variables of the channel output (write_Ch_Sed_output) at the lake outlet links.
   subroutine write_Lake_Sed_output(output_flnm)

      implicit none

      character(len=*), intent(in) :: output_flnm
      integer, parameter           :: nv = 13
      character(len=64)            :: vn(nv), vu(nv)
      character(len=128)           :: vl(nv)
      integer                      :: o(lake%nlake)
      real(8)                      :: v(nps,lake%nlake,nv)

      if (.not. lake%active) return
      o = lake%outlet_ich

      vn = [character(len=64) :: &
         "Sch", "Scw", "Susch0", "Sxpnt", "Solch0", "Sinput", "Sdsch", "Scher", "Schdp", "Schcp", &
         "Cchsd", "Sxabs", "Sxdis"]
      vl = [character(len=128) :: &
         "channel sediment storage", "channel water sediment storage", "upstream discharge", &
         "sed. in-flux from point source", "overland to channel", "total sediment input", &
         "downstream discharge", "channel erosion", "channel deposition", "channel capacity", &
         "channel sediment concentration", "sed. out-flux by abstraction", &
         "sed. in-flux by discharge"]
      vu = [character(len=64) :: &
         "kg", "kg", "kg s-1", "kg", "kg s-1", "kg", "kg s-1", "kg s-1", "kg s-1", "kg", "kg m-3", &
         "kg", "kg"]

      v(:,:,1) = channelSed%Sch(:,o)
      v(:,:,2) = channelSed%Scw(:,o)
      v(:,:,3) = channelSed%Susch0(:,o)
      v(:,:,4) = channelSed%Sxpnt(:,o)
      v(:,:,5) = channelSed%Solch0(:,o)
      v(:,:,6) = channelSed%Sinput(:,o)
      v(:,:,7) = channelSed%Sdsch(:,o)
      v(:,:,8) = channelSed%Scher(:,o)
      v(:,:,9) = channelSed%Schdp(:,o)
      v(:,:,10) = channelSed%Schcp(:,o)
      v(:,:,11) = channelSed%Cchsd(:,o)
      v(:,:,12) = channelSed%Sxabs(:,o)
      v(:,:,13) = channelSed%Sxdis(:,o)

      call write_lake_nc(output_flnm, 'LAKESEDOUT', nv, vn, vl, vu, v2=v)

   end subroutine write_Lake_Sed_output

   !> @brief Writes the lake carbon output file (LAKECOUT): lake water balance and the
   !! carbon variables of the channel output (write_Ch_C_output) at the lake outlet links.
   subroutine write_Lake_C_output(output_flnm)

      implicit none

      character(len=*), intent(in) :: output_flnm
      integer, parameter           :: nv = 54
      character(len=64)            :: vn(nv), vu(nv)
      character(len=128)           :: vl(nv)
      integer                      :: o(lake%nlake)
      real(8)                      :: v(lake%nlake,nv)

      if (.not. lake%active) return
      o = lake%outlet_ich

      vn = [character(len=64) :: &
         "Ch_ALGCusch0_1", "Ch_ALGCusch0_2", "Ch_ALGCusch0_3", "Ch_ZOOCusch0", "Ch_LPOCusch0", &
         "Ch_RPOCusch0", "Ch_LDOCusch0", "Ch_RDOCusch0", "Ch_MBMCusch0", "Ch_DICusch0", &
         "Ch_LPOCxpnt", "Ch_RPOCxpnt", "Ch_LDOCxpnt", "Ch_RDOCxpnt", "Ch_LPOCxabs", "Ch_RPOCxabs", &
         "Ch_LDOCxabs", "Ch_RDOCxabs", "Ch_LPOCxdis", "Ch_RPOCxdis", "Ch_LDOCxdis", "Ch_RDOCxdis", &
         "Ch_ALGCdsch_1", "Ch_ALGCdsch_2", "Ch_ALGCdsch_3", "Ch_ZOOCdsch", "Ch_LPOCdsch", &
         "Ch_RPOCdsch", "Ch_LDOCdsch", "Ch_RDOCdsch", "Ch_MBMCdsch", "Ch_DICdsch", "Ch_ALGC_1", &
         "Ch_ALGC_2", "Ch_ALGC_3", "Ch_ZOOC", "Ch_LPOC", "Ch_RPOC", "Ch_LDOC", "Ch_RDOC", &
         "Ch_MBMC", "Ch_DIC", "Ch_CSC", "Ch_CMBError", "Ch_LPOCsurf0", "Ch_RPOCsurf0", &
         "Ch_LDOCsurf0", "Ch_RDOCsurf0", "Ch_MBMCsurf0", "Ch_LDOCintf0", "Ch_RDOCintf0", &
         "Ch_LDOCgwch0", "Ch_RDOCgwch0", "Twater"]
      vl = [character(len=128) :: &
         "Ch_ALGCusch0_1", "Ch_ALGCusch0_2", "Ch_ALGCusch0_3", "Ch_ZOOCusch0", "Ch_LPOCusch0", &
         "Ch_RPOCusch0", "Ch_LDOCusch0", "Ch_RDOCusch0", "Ch_MBMCusch0", "Ch_DICusch0", &
         "Ch_LPOCxpnt", "Ch_RPOCxpnt", "Ch_LDOCxpnt", "Ch_RDOCxpnt", "Ch_LPOCxabs", "Ch_RPOCxabs", &
         "Ch_LDOCxabs", "Ch_RDOCxabs", "Ch_LPOCxdis", "Ch_RPOCxdis", "Ch_LDOCxdis", "Ch_RDOCxdis", &
         "Ch_ALGCdsch_1", "Ch_ALGCdsch_2", "Ch_ALGCdsch_3", "Ch_ZOOCdsch", "Ch_LPOCdsch", &
         "Ch_RPOCdsch", "Ch_LDOCdsch", "Ch_RDOCdsch", "Ch_MBMCdsch", "Ch_DICdsch", "Ch_ALGC_1", &
         "Ch_ALGC_2", "Ch_ALGC_3", "Ch_ZOOC", "Ch_LPOC", "Ch_RPOC", "Ch_LDOC", "Ch_RDOC", &
         "Ch_MBMC", "Ch_DIC", "Ch_CSC", "Ch_CMBError", "Ch_LPOCsurf0", "Ch_RPOCsurf0", &
         "Ch_LDOCsurf0", "Ch_RDOCsurf0", "Ch_MBMCsurf0", "Ch_LDOCintf0", "Ch_RDOCintf0", &
         "Ch_LDOCgwch0", "Ch_RDOCgwch0", "Twater"]
      vu = [character(len=64) :: &
         "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", &
         "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", &
         "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", &
         "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", &
         "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC", "kgC", "kgC", "kgC", "kgC", "kgC", &
         "kgC", "kgC", "kgC", "kgC", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", &
         "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "kgC dt-1", "deg-C"]

      v(:,1) = Ch_ALGCusch0(o,1)
      v(:,2) = Ch_ALGCusch0(o,2)
      v(:,3) = Ch_ALGCusch0(o,3)
      v(:,4) = Ch_ZOOCusch0(o)
      v(:,5) = Ch_LPOCusch0(o)
      v(:,6) = Ch_RPOCusch0(o)
      v(:,7) = Ch_LDOCusch0(o)
      v(:,8) = Ch_RDOCusch0(o)
      v(:,9) = Ch_MBMCusch0(o)
      v(:,10) = Ch_DICusch0(o)
      v(:,11) = Ch_LPOCxpnt(o)
      v(:,12) = Ch_RPOCxpnt(o)
      v(:,13) = Ch_LDOCxpnt(o)
      v(:,14) = Ch_RDOCxpnt(o)
      v(:,15) = Ch_LPOCxabs(o)
      v(:,16) = Ch_RPOCxabs(o)
      v(:,17) = Ch_LDOCxabs(o)
      v(:,18) = Ch_RDOCxabs(o)
      v(:,19) = Ch_LPOCxdis(o)
      v(:,20) = Ch_RPOCxdis(o)
      v(:,21) = Ch_LDOCxdis(o)
      v(:,22) = Ch_RDOCxdis(o)
      v(:,23) = Ch_ALGCdsch(o,1)
      v(:,24) = Ch_ALGCdsch(o,2)
      v(:,25) = Ch_ALGCdsch(o,3)
      v(:,26) = Ch_ZOOCdsch(o)
      v(:,27) = Ch_LPOCdsch(o)
      v(:,28) = Ch_RPOCdsch(o)
      v(:,29) = Ch_LDOCdsch(o)
      v(:,30) = Ch_RDOCdsch(o)
      v(:,31) = Ch_MBMCdsch(o)
      v(:,32) = Ch_DICdsch(o)
      v(:,33) = Ch_ALGC(o,1)
      v(:,34) = Ch_ALGC(o,2)
      v(:,35) = Ch_ALGC(o,3)
      v(:,36) = Ch_ZOOC(o)
      v(:,37) = Ch_LPOC(o)
      v(:,38) = Ch_RPOC(o)
      v(:,39) = Ch_LDOC(o)
      v(:,40) = Ch_RDOC(o)
      v(:,41) = Ch_MBMC(o)
      v(:,42) = Ch_DIC(o)
      v(:,43) = Ch_CSC(o)
      v(:,44) = Ch_CMBError(o)
      v(:,45) = Ch_LPOCsurf0(o)
      v(:,46) = Ch_RPOCsurf0(o)
      v(:,47) = Ch_LDOCsurf0(o)
      v(:,48) = Ch_RDOCsurf0(o)
      v(:,49) = Ch_MBMCsurf0(o)
      v(:,50) = Ch_LDOCintf0(o)
      v(:,51) = Ch_RDOCintf0(o)
      v(:,52) = Ch_LDOCgwch0(o)
      v(:,53) = Ch_RDOCgwch0(o)
      v(:,54) = TWATER(o)

      call write_lake_nc(output_flnm, 'LAKECOUT', nv, vn, vl, vu, v1=v)

   end subroutine write_Lake_C_output

   !> @brief Writes the lake nitrogen output file (LAKENOUT): lake water balance and the
   !! nitrogen variables of the channel output (write_Ch_N_output) at the lake outlet links.
   subroutine write_Lake_N_output(output_flnm)

      implicit none

      character(len=*), intent(in) :: output_flnm
      integer, parameter           :: nv = 70
      character(len=64)            :: vn(nv), vu(nv)
      character(len=128)           :: vl(nv)
      integer                      :: o(lake%nlake)
      real(8)                      :: v(lake%nlake,nv)

      if (.not. lake%active) return
      o = lake%outlet_ich

      vn = [character(len=64) :: &
         "Ch_ALGNusch0_1", "Ch_ALGNusch0_2", "Ch_ALGNusch0_3", "Ch_ZOONusch0", "Ch_LPONusch0", &
         "Ch_RPONusch0", "Ch_LDONusch0", "Ch_RDONusch0", "Ch_MBMNusch0", "Ch_NH4usch0", &
         "Ch_NO3usch0", "Ch_LPONxpnt", "Ch_RPONxpnt", "Ch_LDONxpnt", "Ch_RDONxpnt", "Ch_NH4xpnt", &
         "Ch_NO3xpnt", "Ch_LPONxabs", "Ch_RPONxabs", "Ch_LDONxabs", "Ch_RDONxabs", "Ch_NH4xabs", &
         "Ch_NO3xabs", "Ch_LPONxdis", "Ch_RPONxdis", "Ch_LDONxdis", "Ch_RDONxdis", "Ch_NH4xdis", &
         "Ch_NO3xdis", "Ch_ALGNdsch_1", "Ch_ALGNdsch_2", "Ch_ALGNdsch_3", "Ch_ZOONdsch", &
         "Ch_LPONdsch", "Ch_RPONdsch", "Ch_LDONdsch", "Ch_RDONdsch", "Ch_MBMNdsch", "Ch_NH4dsch", &
         "Ch_NO3dsch", "Ch_NH3VOL", "Ch_DENIT", "Ch_ALGN_1", "Ch_ALGN_2", "Ch_ALGN_3", "Ch_ZOON", &
         "Ch_LPON", "Ch_RPON", "Ch_LDON", "Ch_RDON", "Ch_MBMN", "Ch_NH4", "Ch_NO3", "Ch_NSC", &
         "Ch_NMBError", "Ch_LPONsurf0", "Ch_RPONsurf0", "Ch_LDONsurf0", "Ch_RDONsurf0", &
         "Ch_MBMNsurf0", "Ch_NH4surf0", "Ch_NO3surf0", "Ch_LDONintf0", "Ch_RDONintf0", &
         "Ch_NH4intf0", "Ch_NO3intf0", "Ch_LDONgwch0", "Ch_RDONgwch0", "Ch_NH4gwch0", &
         "Ch_NO3gwch0"]
      vl = [character(len=128) :: &
         "Ch_ALGNusch0_1", "Ch_ALGNusch0_2", "Ch_ALGNusch0_3", "Ch_ZOONusch0", "Ch_LPONusch0", &
         "Ch_RPONusch0", "Ch_LDONusch0", "Ch_RDONusch0", "Ch_MBMNusch0", "Ch_NH4usch0", &
         "Ch_NO3usch0", "Ch_LPONxpnt", "Ch_RPONxpnt", "Ch_LDONxpnt", "Ch_RDONxpnt", "Ch_NH4xpnt", &
         "Ch_NO3xpnt", "Ch_LPONxabs", "Ch_RPONxabs", "Ch_LDONxabs", "Ch_RDONxabs", "Ch_NH4xabs", &
         "Ch_NO3xabs", "Ch_LPONxdis", "Ch_RPONxdis", "Ch_LDONxdis", "Ch_RDONxdis", "Ch_NH4xdis", &
         "Ch_NO3xdis", "Ch_ALGNdsch_1", "Ch_ALGNdsch_2", "Ch_ALGNdsch_3", "Ch_ZOONdsch", &
         "Ch_LPONdsch", "Ch_RPONdsch", "Ch_LDONdsch", "Ch_RDONdsch", "Ch_MBMNdsch", "Ch_NH4dsch", &
         "Ch_NO3dsch", "Ch_NH3VOL", "Ch_DENIT", "Ch_ALGN_1", "Ch_ALGN_2", "Ch_ALGN_3", "Ch_ZOON", &
         "Ch_LPON", "Ch_RPON", "Ch_LDON", "Ch_RDON", "Ch_MBMN", "Ch_NH4", "Ch_NO3", "Ch_NSC", &
         "Ch_NMBError", "Ch_LPONsurf0", "Ch_RPONsurf0", "Ch_LDONsurf0", "Ch_RDONsurf0", &
         "Ch_MBMNsurf0", "Ch_NH4surf0", "Ch_NO3surf0", "Ch_LDONintf0", "Ch_RDONintf0", &
         "Ch_NH4intf0", "Ch_NO3intf0", "Ch_LDONgwch0", "Ch_RDONgwch0", "Ch_NH4gwch0", &
         "Ch_NO3gwch0"]
      vu = [character(len=64) :: &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN", "kgN", "kgN", "kgN", "kgN", "kgN", "kgN", "kgN", "kgN", "kgN", "kgN", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", "kgN dt-1", &
         "kgN dt-1", "kgN dt-1"]

      v(:,1) = Ch_ALGNusch0(o,1)
      v(:,2) = Ch_ALGNusch0(o,2)
      v(:,3) = Ch_ALGNusch0(o,3)
      v(:,4) = Ch_ZOONusch0(o)
      v(:,5) = Ch_LPONusch0(o)
      v(:,6) = Ch_RPONusch0(o)
      v(:,7) = Ch_LDONusch0(o)
      v(:,8) = Ch_RDONusch0(o)
      v(:,9) = Ch_MBMNusch0(o)
      v(:,10) = Ch_NH4usch0(o)
      v(:,11) = Ch_NO3usch0(o)
      v(:,12) = Ch_LPONxpnt(o)
      v(:,13) = Ch_RPONxpnt(o)
      v(:,14) = Ch_LDONxpnt(o)
      v(:,15) = Ch_RDONxpnt(o)
      v(:,16) = Ch_NH4xpnt(o)
      v(:,17) = Ch_NO3xpnt(o)
      v(:,18) = Ch_LPONxabs(o)
      v(:,19) = Ch_RPONxabs(o)
      v(:,20) = Ch_LDONxabs(o)
      v(:,21) = Ch_RDONxabs(o)
      v(:,22) = Ch_NH4xabs(o)
      v(:,23) = Ch_NO3xabs(o)
      v(:,24) = Ch_LPONxdis(o)
      v(:,25) = Ch_RPONxdis(o)
      v(:,26) = Ch_LDONxdis(o)
      v(:,27) = Ch_RDONxdis(o)
      v(:,28) = Ch_NH4xdis(o)
      v(:,29) = Ch_NO3xdis(o)
      v(:,30) = Ch_ALGNdsch(o,1)
      v(:,31) = Ch_ALGNdsch(o,2)
      v(:,32) = Ch_ALGNdsch(o,3)
      v(:,33) = Ch_ZOONdsch(o)
      v(:,34) = Ch_LPONdsch(o)
      v(:,35) = Ch_RPONdsch(o)
      v(:,36) = Ch_LDONdsch(o)
      v(:,37) = Ch_RDONdsch(o)
      v(:,38) = Ch_MBMNdsch(o)
      v(:,39) = Ch_NH4dsch(o)
      v(:,40) = Ch_NO3dsch(o)
      v(:,41) = Ch_NH3VOL(o)
      v(:,42) = Ch_DENIT(o)
      v(:,43) = Ch_ALGN(o,1)
      v(:,44) = Ch_ALGN(o,2)
      v(:,45) = Ch_ALGN(o,3)
      v(:,46) = Ch_ZOON(o)
      v(:,47) = Ch_LPON(o)
      v(:,48) = Ch_RPON(o)
      v(:,49) = Ch_LDON(o)
      v(:,50) = Ch_RDON(o)
      v(:,51) = Ch_MBMN(o)
      v(:,52) = Ch_NH4(o)
      v(:,53) = Ch_NO3(o)
      v(:,54) = Ch_NSC(o)
      v(:,55) = Ch_NMBError(o)
      v(:,56) = Ch_LPONsurf0(o)
      v(:,57) = Ch_RPONsurf0(o)
      v(:,58) = Ch_LDONsurf0(o)
      v(:,59) = Ch_RDONsurf0(o)
      v(:,60) = Ch_MBMNsurf0(o)
      v(:,61) = Ch_NH4surf0(o)
      v(:,62) = Ch_NO3surf0(o)
      v(:,63) = Ch_LDONintf0(o)
      v(:,64) = Ch_RDONintf0(o)
      v(:,65) = Ch_NH4intf0(o)
      v(:,66) = Ch_NO3intf0(o)
      v(:,67) = Ch_LDONgwch0(o)
      v(:,68) = Ch_RDONgwch0(o)
      v(:,69) = Ch_NH4gwch0(o)
      v(:,70) = Ch_NO3gwch0(o)

      call write_lake_nc(output_flnm, 'LAKENOUT', nv, vn, vl, vu, v1=v)

   end subroutine write_Lake_N_output

   !> @brief Writes the lake phosphorus output file (LAKEPOUT): lake water balance and the
   !! phosphorus variables of the channel output (write_Ch_P_output) at the lake outlet links.
   subroutine write_Lake_P_output(output_flnm)

      implicit none

      character(len=*), intent(in) :: output_flnm
      integer, parameter           :: nv = 67
      character(len=64)            :: vn(nv), vu(nv)
      character(len=128)           :: vl(nv)
      integer                      :: o(lake%nlake)
      real(8)                      :: v(lake%nlake,nv)

      if (.not. lake%active) return
      o = lake%outlet_ich

      vn = [character(len=64) :: &
         "Ch_ALGPusch0_1", "Ch_ALGPusch0_2", "Ch_ALGPusch0_3", "Ch_ZOOPusch0", "Ch_LPOPusch0", &
         "Ch_RPOPusch0", "Ch_LDOPusch0", "Ch_RDOPusch0", "Ch_MBMPusch0", "Ch_PO4usch0", &
         "Ch_PIPAusch0", "Ch_PIPSusch0", "Ch_LPOPxpnt", "Ch_RPOPxpnt", "Ch_LDOPxpnt", &
         "Ch_RDOPxpnt", "Ch_PO4xpnt", "Ch_LPOPxabs", "Ch_RPOPxabs", "Ch_LDOPxabs", "Ch_RDOPxabs", &
         "Ch_PO4xabs", "Ch_LPOPxdis", "Ch_RPOPxdis", "Ch_LDOPxdis", "Ch_RDOPxdis", "Ch_PO4xdis", &
         "Ch_ALGPdsch_1", "Ch_ALGPdsch_2", "Ch_ALGPdsch_3", "Ch_ZOOPdsch", "Ch_LPOPdsch", &
         "Ch_RPOPdsch", "Ch_LDOPdsch", "Ch_RDOPdsch", "Ch_MBMPdsch", "Ch_PO4dsch", "Ch_PIPAdsch", &
         "Ch_PIPSdsch", "Ch_ALGP_1", "Ch_ALGP_2", "Ch_ALGP_3", "Ch_ZOOP", "Ch_LPOP", "Ch_RPOP", &
         "Ch_LDOP", "Ch_RDOP", "Ch_MBMP", "Ch_PO4", "Ch_PIPA", "Ch_PIPS", "Ch_PSC", "Ch_PMBError", &
         "Ch_LPOPsurf0", "Ch_RPOPsurf0", "Ch_LDOPsurf0", "Ch_RDOPsurf0", "Ch_MBMPsurf0", &
         "Ch_PO4surf0", "Ch_PIPAsurf0", "Ch_PIPSsurf0", "Ch_LDOPintf0", "Ch_RDOPintf0", &
         "Ch_PO4intf0", "Ch_LDOPgwch0", "Ch_RDOPgwch0", "Ch_PO4gwch0"]
      vl = [character(len=128) :: &
         "Ch_ALGPusch0_1", "Ch_ALGPusch0_2", "Ch_ALGPusch0_3", "Ch_ZOOPusch0", "Ch_LPOPusch0", &
         "Ch_RPOPusch0", "Ch_LDOPusch0", "Ch_RDOPusch0", "Ch_MBMPusch0", "Ch_PO4usch0", &
         "Ch_PIPAusch0", "Ch_PIPSusch0", "Ch_LPOPxpnt", "Ch_RPOPxpnt", "Ch_LDOPxpnt", &
         "Ch_RDOPxpnt", "Ch_PO4xpnt", "Ch_LPOPxabs", "Ch_RPOPxabs", "Ch_LDOPxabs", "Ch_RDOPxabs", &
         "Ch_PO4xabs", "Ch_LPOPxdis", "Ch_RPOPxdis", "Ch_LDOPxdis", "Ch_RDOPxdis", "Ch_PO4xdis", &
         "Ch_ALGPdsch_1", "Ch_ALGPdsch_2", "Ch_ALGPdsch_3", "Ch_ZOOPdsch", "Ch_LPOPdsch", &
         "Ch_RPOPdsch", "Ch_LDOPdsch", "Ch_RDOPdsch", "Ch_MBMPdsch", "Ch_PO4dsch", "Ch_PIPAdsch", &
         "Ch_PIPSdsch", "Ch_ALGP_1", "Ch_ALGP_2", "Ch_ALGP_3", "Ch_ZOOP", "Ch_LPOP", "Ch_RPOP", &
         "Ch_LDOP", "Ch_RDOP", "Ch_MBMP", "Ch_PO4", "Ch_PIPA", "Ch_PIPS", "Ch_PSC", "Ch_PMBError", &
         "Ch_LPOPsurf0", "Ch_RPOPsurf0", "Ch_LDOPsurf0", "Ch_RDOPsurf0", "Ch_MBMPsurf0", &
         "Ch_PO4surf0", "Ch_PIPAsurf0", "Ch_PIPSsurf0", "Ch_LDOPintf0", "Ch_RDOPintf0", &
         "Ch_PO4intf0", "Ch_LDOPgwch0", "Ch_RDOPgwch0", "Ch_PO4gwch0"]
      vu = [character(len=64) :: &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP", "kgP", "kgP", "kgP", "kgP", "kgP", &
         "kgP", "kgP", "kgP", "kgP", "kgP", "kgP", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", &
         "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1", "kgP dt-1"]

      v(:,1) = Ch_ALGPusch0(o,1)
      v(:,2) = Ch_ALGPusch0(o,2)
      v(:,3) = Ch_ALGPusch0(o,3)
      v(:,4) = Ch_ZOOPusch0(o)
      v(:,5) = Ch_LPOPusch0(o)
      v(:,6) = Ch_RPOPusch0(o)
      v(:,7) = Ch_LDOPusch0(o)
      v(:,8) = Ch_RDOPusch0(o)
      v(:,9) = Ch_MBMPusch0(o)
      v(:,10) = Ch_PO4usch0(o)
      v(:,11) = Ch_PIPAusch0(o)
      v(:,12) = Ch_PIPSusch0(o)
      v(:,13) = Ch_LPOPxpnt(o)
      v(:,14) = Ch_RPOPxpnt(o)
      v(:,15) = Ch_LDOPxpnt(o)
      v(:,16) = Ch_RDOPxpnt(o)
      v(:,17) = Ch_PO4xpnt(o)
      v(:,18) = Ch_LPOPxabs(o)
      v(:,19) = Ch_RPOPxabs(o)
      v(:,20) = Ch_LDOPxabs(o)
      v(:,21) = Ch_RDOPxabs(o)
      v(:,22) = Ch_PO4xabs(o)
      v(:,23) = Ch_LPOPxdis(o)
      v(:,24) = Ch_RPOPxdis(o)
      v(:,25) = Ch_LDOPxdis(o)
      v(:,26) = Ch_RDOPxdis(o)
      v(:,27) = Ch_PO4xdis(o)
      v(:,28) = Ch_ALGPdsch(o,1)
      v(:,29) = Ch_ALGPdsch(o,2)
      v(:,30) = Ch_ALGPdsch(o,3)
      v(:,31) = Ch_ZOOPdsch(o)
      v(:,32) = Ch_LPOPdsch(o)
      v(:,33) = Ch_RPOPdsch(o)
      v(:,34) = Ch_LDOPdsch(o)
      v(:,35) = Ch_RDOPdsch(o)
      v(:,36) = Ch_MBMPdsch(o)
      v(:,37) = Ch_PO4dsch(o)
      v(:,38) = Ch_PIPAdsch(o)
      v(:,39) = Ch_PIPSdsch(o)
      v(:,40) = Ch_ALGP(o,1)
      v(:,41) = Ch_ALGP(o,2)
      v(:,42) = Ch_ALGP(o,3)
      v(:,43) = Ch_ZOOP(o)
      v(:,44) = Ch_LPOP(o)
      v(:,45) = Ch_RPOP(o)
      v(:,46) = Ch_LDOP(o)
      v(:,47) = Ch_RDOP(o)
      v(:,48) = Ch_MBMP(o)
      v(:,49) = Ch_PO4(o)
      v(:,50) = Ch_PIPA(o)
      v(:,51) = Ch_PIPS(o)
      v(:,52) = Ch_PSC(o)
      v(:,53) = Ch_PMBError(o)
      v(:,54) = Ch_LPOPsurf0(o)
      v(:,55) = Ch_RPOPsurf0(o)
      v(:,56) = Ch_LDOPsurf0(o)
      v(:,57) = Ch_RDOPsurf0(o)
      v(:,58) = Ch_MBMPsurf0(o)
      v(:,59) = Ch_PO4surf0(o)
      v(:,60) = Ch_PIPAsurf0(o)
      v(:,61) = Ch_PIPSsurf0(o)
      v(:,62) = Ch_LDOPintf0(o)
      v(:,63) = Ch_RDOPintf0(o)
      v(:,64) = Ch_PO4intf0(o)
      v(:,65) = Ch_LDOPgwch0(o)
      v(:,66) = Ch_RDOPgwch0(o)
      v(:,67) = Ch_PO4gwch0(o)

      call write_lake_nc(output_flnm, 'LAKEPOUT', nv, vn, vl, vu, v1=v)

   end subroutine write_Lake_P_output

   !> @brief Writes a LAKE*OUT NetCDF file: lake ID, outlet link ID, lake water balance
   !! (water surface elevation, depth, inflow, storage, outflow, area, bottom elevation) and
   !! the material variables, either per lake (v1) or per grain size and lake (v2).
   subroutine write_lake_nc(output_flnm, ftype, nv, vn, vl, vu, v1, v2)

      use netcdf
      use module_NWM_io_dict

      implicit none

      character(len=*),  intent(in)           :: output_flnm, ftype
      integer,           intent(in)           :: nv
      character(len=64), intent(in)           :: vn(nv), vu(nv)
      character(len=128),intent(in)           :: vl(nv)
      real(8),           intent(in), optional :: v1(:,:)     ! (nlake,nv)
      real(8),           intent(in), optional :: v2(:,:,:)   ! (nps,nlake,nv)
      integer, parameter  :: nw = 7
      character(len=64)   :: wn(nw), wu(nw)
      character(len=128)  :: wl(nw)
      real(8)             :: w(lake%nlake,nw)
      integer             :: diagFlag, ftn, iret, varId, dimId(3), iv, il
      integer             :: outlet_link(lake%nlake)

      diagFlag = 0

      wn = [character(len=64) :: "water_sfc_elev", "depth", "inflow", "storage", "outflow", "area", "BottomE"]
      wl = [character(len=128) :: "lake water surface elevation", "lake mean water depth", &
            "lake water inflow", "lake water storage", "lake water outflow", &
            "lake surface area", "lake bottom elevation"]
      wu = [character(len=64) :: "m", "m", "m3 s-1", "m3", "m3 s-1", "m2", "m"]
      w(:,1) = lake%wse
      w(:,2) = lake%depth
      w(:,3) = lake%inflow
      w(:,4) = lake%storage
      w(:,5) = lake%outflow
      w(:,6) = lake%area
      w(:,7) = lake%bottomE
      do il = 1, lake%nlake
         outlet_link(il) = SedCNP_hydro%linkID(lake%outlet_ich(il))
      enddo

      iret = nf90_create(trim(output_flnm), cmode=NF90_NETCDF4, ncid=ftn)
      call nwmCheck(diagFlag,iret,'ERROR: Unable to create '//ftype//' NetCDF file.')
      iret = nf90_def_dim(ftn,'time',NF90_UNLIMITED,dimId(1))
      call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
      iret = nf90_def_dim(ftn,'nlake',lake%nlake,dimId(2))
      call nwmCheck(diagFlag,iret,'ERROR: Unable to define nlake dimension')
      if (present(v2)) then
         iret = nf90_def_dim(ftn,'grain_size',nps,dimId(3))
         call nwmCheck(diagFlag,iret,'ERROR: Unable to define grain_size dimension')
      endif

      iret = nf90_def_var(ftn,"time",nf90_float,(/dimId(1)/),varId)
      call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: time")
      iret = nf90_put_att(ftn,varId,'long_name','valid output time')
      iret = nf90_put_att(ftn,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
      iret = nf90_def_var(ftn,"lake_id",nf90_int,(/dimId(2)/),varId)
      call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: lake_id")
      iret = nf90_put_att(ftn,varId,'long_name','lake ID (LAKEPARM lake_id)')
      iret = nf90_def_var(ftn,"outlet_link",nf90_int,(/dimId(2)/),varId)
      call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: outlet_link")
      iret = nf90_put_att(ftn,varId,'long_name','link ID of the lake outlet (Route_Link link)')
      do iv = 1, nw
         iret = nf90_def_var(ftn,trim(wn(iv)),nf90_double,(/dimId(2),dimId(1)/),varId)
         call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(wn(iv)))
         iret = nf90_put_att(ftn,varId,'long_name',trim(wl(iv)))
         iret = nf90_put_att(ftn,varId,'units',trim(wu(iv)))
      enddo
      do iv = 1, nv
         if (present(v2)) then
            iret = nf90_def_var(ftn,trim(vn(iv)),nf90_double,(/dimId(3),dimId(2),dimId(1)/),varId)
         else
            iret = nf90_def_var(ftn,trim(vn(iv)),nf90_double,(/dimId(2),dimId(1)/),varId)
         endif
         call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(vn(iv)))
         iret = nf90_put_att(ftn,varId,'long_name',trim(vl(iv)))
         iret = nf90_put_att(ftn,varId,'units',trim(vu(iv)))
      enddo
      iret = nf90_enddef(ftn)
      call nwmCheck(diagFlag,iret,'ERROR: Unable to take '//ftype//' file out of definition mode')

      iret = nf90_inq_varid(ftn,'time',varId)
      iret = nf90_put_var(ftn,varId,SedCNP_hydro%time,(/1/),(/1/))
      call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')
      iret = nf90_inq_varid(ftn,'lake_id',varId)
      iret = nf90_put_var(ftn,varId,lake%lake_id)
      call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: lake_id')
      iret = nf90_inq_varid(ftn,'outlet_link',varId)
      iret = nf90_put_var(ftn,varId,outlet_link)
      call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: outlet_link')
      do iv = 1, nw
         iret = nf90_inq_varid(ftn,trim(wn(iv)),varId)
         iret = nf90_put_var(ftn,varId,w(:,iv),(/1,1/),(/lake%nlake,1/))
         call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: '//trim(wn(iv)))
      enddo
      do iv = 1, nv
         iret = nf90_inq_varid(ftn,trim(vn(iv)),varId)
         if (present(v2)) then
            iret = nf90_put_var(ftn,varId,v2(:,:,iv),(/1,1,1/),(/nps,lake%nlake,1/))
         else
            iret = nf90_put_var(ftn,varId,v1(:,iv),(/1,1/),(/lake%nlake,1/))
         endif
         call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: '//trim(vn(iv)))
      enddo
      iret = nf90_close(ftn)
      call nwmCheck(diagFlag,iret,'ERROR: Unable to close '//ftype//' file')

   end subroutine write_lake_nc


   !> @brief Lake index (1..nlake) of a lake ID, 0 if not found.
   integer function lake_index(id)
      implicit none
      integer, intent(in) :: id
      integer             :: il

      lake_index = 0
      do il = 1, lake%nlake
         if (lake%lake_id(il) == id) then
            lake_index = il
            return
         endif
      enddo
   end function lake_index

end module module_lakeSedCNP
!
!=====||__WHQ5403__||=====!
