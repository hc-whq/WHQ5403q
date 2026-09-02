module CNPmain

   use module_SedCNPvariables
   use Carbon
   use Nitrogen
   use Phosphorus
   use module_CNPupdates
!=====||__WHQ5403q__||=====!
!
   ! SedCNPmodel is already visible transitively via Carbon/Nitrogen/Phosphorus/
   ! module_CNPupdates above; a redundant direct "use config_base, only: SedCNPmodel"
   ! here triggers a gfortran diamond-import derived-type mismatch (namelist_rt_).
!
!=====||__WHQ5403q__||=====!
   !use NOAHMP_TABLES, only: ISWATER_TABLE  !BK20251208 Reverted

   implicit none

   logical, save :: basin_averages_calculated = .false.  !BK20250831
   integer, save :: last_cnp_timestep = -1               !BK20250831

   contains


   subroutine RunCNP_So(i,j,itime)

      implicit none
   
      ! --- Input Arguments
      integer :: i, j, itime

      ! --- Local Hydrological & Physical Variables
      real(8), dimension(4) :: TSOIL       ! soil layer temperature [C]
      real(8), dimension(4) :: THETA       ! volumetric soil moisture [m3/m3]
      real(8)               :: THETA_S     ! saturation soil moisture [m3/m3]
      real(8)               :: THETA_F     ! field capacity soil moisture [m3/m3]
      real(8)               :: Qsurf       ! surface runoff [mm/s]
      real(8)               :: Qintf       ! interflow [mm/s]
      real(8), dimension(4) :: Qperc       ! percolation between layers [mm/s]
      real(8), dimension(4) :: Quptk       ! vegetation water uptake [mm/s]
      real(8)               :: Qgwso       ! irrigation water from GW abstraction [mm/s]
      integer               :: lcover      ! land cover code
      real(8)               :: dx          ! grid cell size [m]
      real(8), dimension(4) :: So_WStorage ! soil layer water storage [m3]
      integer               :: isl

      ! Set model time-step and grid size
      domain%DT = real(SedCNPmodel%SedCNP_timestep)
      dx = SedCNP_hydro%dx(1)       
      domain%areaxy = dx * dx

      ! TODO: These should be read from an input file, not hardcoded.
      DSOIL(1) = 0.10
      DSOIL(2) = 0.30
      DSOIL(3) = 0.60
      DSOIL(4) = 1.00

      SSA(1) = 2264.2   ! clay (diameter 1 um)
      SSA(2) = 226.42   ! silt (diameter 10 um)
      SSA(3) = 22.642   ! fine sand (diameter 100 um)
      SSA(4) = 2.2642   ! coarse sand (diameter 1000 um)

      !Dmix = 0.01  ! effective mixing depth interacting with surface runoff [m] 

      ! --- Get hydrological inputs for the current grid cell
      TSOIL(:) = SedCNP_hydro%t_soil(i,:,j) - 273.15 ! K -> C
      THETA(:) = SedCNP_hydro%theta(i,:,j)
      THETA_S  = SedCNP_hydro%smcmax(i,j,1)
      THETA_F  = SedCNP_hydro%smcref(i,j,1)
      Qsurf    = SedCNP_hydro%runsrf(i,j)                ! [mm/s]
      Qintf    = SedCNP_hydro%q_intf(i,j) / domain%DT    ! [mm/s]
      Qperc(:) = SedCNP_hydro%percolate(i,:,j)           ! [mm/s]
      Quptk(:) = SedCNP_hydro%ETRANLayer(i,:,j)          ! [mm/s]
      Qsogw    = SedCNP_hydro%q_sogw(i,j) / domain%DT    ! [mm/s]  !BK20260313   
      Qgwso    = SedCNP_hydro%irridepth(i,j) / domain%DT ! [mm/s]  !BK20260313

      if (SedCNP_hydro%vgtyp(i,j) /= SedCNP_hydro%vgtyp(i,j)) then
         lcover = domain%ISWATER
      else
         lcover = int(SedCNP_hydro%vgtyp(i,j))
      endif

      do isl = 1, domain%nsl
         So_WStorage(isl) = SedCNP_hydro%theta(i,isl,j) * (dx*dx*DSOIL(isl))
      enddo

      ! --- Check for valid input values
      if(lcover.gt.0) then
         if((Qsurf<0.0).or.(Qsurf/=Qsurf)) Qsurf = 0.0
         if((Qgwso<0.0).or.(Qgwso/=Qgwso)) Qgwso = 0.0    !BK20260316
         if((Qintf<0.0).or.(Qintf/=Qintf)) Qintf = 0.0
      
         do isl = 1, domain%nsl
            if(THETA(isl)<0.0) THETA(isl) = 0.0
            if(THETA_S<0.0) THETA_S = 0.0
            if(THETA_F<0.0) THETA_F = 0.0
            if((Qperc(isl)<0.0).or.(Qperc(isl)/=Qperc(isl))) Qperc(isl) = 0.0
            if((Quptk(isl)<0.0).or.(Quptk(isl)/=Quptk(isl))) Quptk(isl) = 0.0
         enddo
      endif

      ! --- Run CNP cycle for land grid cells
      if(lcover.gt.0) then 

         ! Carbon cycle for soil 
         call SoilCarbonCycle(i,j,lcover,itime,TSOIL,THETA,THETA_S,THETA_F, &
            So_WStorage,Quptk,Qsurf,Qperc,Qintf)  
         call CheckCMB_So(i,j,lcover)

         ! Nitrogen cycle for soil 
         call SoilNitrogenCycle(i,j,lcover,itime,TSOIL,THETA,THETA_S,THETA_F, &
            So_WStorage,Quptk,Qsurf,Qperc,Qintf)  
         call CheckNMB_So(i,j,lcover)

         ! Phosphorus cycle for soil 
         call SoilPhosphorusCycle(i,j,lcover,itime,TSOIL,THETA,THETA_S,THETA_F, &
            So_WStorage,Quptk,Qsurf,Qperc,Qintf)
         call CheckPMB_So(i,j,lcover)

      endif

   end subroutine RunCNP_So


   subroutine RunCNP_Gw(gwid)  !BK20260316

      implicit none
      
      ! --- Input Arguments
      integer, intent(in) :: gwid

      ! --- Local Variables
      integer :: i, j
      real(8) :: Gw_WS       ! groundwater storage [m3]
      real(8) :: Qgwch       ! groundwater discharge to channel [mm s-1]
      real(8) :: basin_area  ! subbasin area [m2]

      ! Calculate total basin area  !BK20260320
      basin_area = 0.0
      DO i = 1, domain%ix
         DO j = 1, domain%jx
            IF (SedCNP_hydro%gwbasin(i,j) == gwid) THEN
               basin_area = basin_area + domain%areaxy
            END IF
         END DO
      END DO

      ! groundwater storage
      SedCNP_hydro%Gw_WStorage(gwid) = SedCNP_hydro%z_gwbas(gwid) * basin_area !(m3) !BK20260305
      Gw_WS = SedCNP_hydro%Gw_WStorage(gwid)

      ! Calculate groundwater discharge flux to the channel
      if (SedCNP_hydro%q_gwch(gwid).lt.0.0) then
         Qgwch = 0.0
      else
         Qgwch = (SedCNP_hydro%q_gwch(gwid) / domain%areaxy) * 1000.0  ! [mm s-1]
      endif
      
      ! --- Carbon cycle for groundwater
      call GroundwaterCarbonCycle(gwid, Gw_WS, Qgwch)  !BK20260316
      call CheckCMB_Gw(gwid)

      ! --- Nitrogen cycle for groundwater
      call GroundwaterNitrogenCycle(gwid, Gw_WS, Qgwch)  !BK20260316
      call CheckNMB_Gw(gwid)

      ! --- Phosphorus cycle for groundwater
      call GroundwaterPhosphorusCycle(gwid, Gw_WS, Qgwch)  !BK20260316
      call CheckPMB_Gw(gwid)

   end subroutine RunCNP_Gw


   subroutine RunCNP_Ch(ich,itime,first_call)
           
      implicit none
      
      ! --- Input Arguments
      integer, intent(in) :: ich, itime
      logical, intent(in) :: first_call

      ! --- Local Variables
      real(8) :: Qdsch, WStorage, Dch_thres, BtmWdth, SideSlopch
      real(8) :: WtopWdth, Axsect, Lst
      real(8) :: dt
      integer :: chid, gwid

      dt = real(SedCNPmodel%SedCNP_timestep)
      chid  = SedCNP_hydro%linkID(ich)
      gwid = domain%gwid_ch(chid)

      !BK20260213
      SSA(1) = 2264.2   ! clay (diameter 1 um)
      SSA(2) = 226.42   ! silt (diameter 10 um)
      SSA(3) = 22.642   ! fine sand (diameter 100 um)
      SSA(4) = 2.2642   ! coarse sand (diameter 1000 um)

      ! Get channel hydraulic parameters
      DWATER(ich)  = SedCNP_hydro%head(ich)
      SideSlopch = channel_input%SideSlopch(ich)
      BtmWdth    = channel_input%Wst(ich)
      Lst        = channel_input%Lst(ich)

      ! Prevent invalid geometry values
      if (SideSlopch <= 0.0) SideSlopch = 1.0E-6
      if (DWATER(ich) < 0.0) DWATER(ich) = 0.0
      if (BtmWdth < 0.0) BtmWdth = 0.0
      if (Lst < 0.0) Lst = 0.0

      ! Compute cross-section area and water storage
      WtopWdth = 2.0 * (DWATER(ich) / SideSlopch) + BtmWdth
      Axsect = (BtmWdth + WtopWdth) * DWATER(ich) / 2.0
      WStorage = Axsect * Lst
      if (WStorage < 0.0) then
         WStorage = 0.0
         write(*,'(A, I5, A)') 'WARNING: Channel Water Storage below ZERO at channel (', chid, ')'
      endif
      SedCNP_hydro%Wstg_st(ich) = WStorage
      
      Qdsch = SedCNP_hydro%q_dsch(ich) ! [m3/s]

      ! Estimate channel water temperature from areally-averaged air temp
      if (gwid > 0) then
         ! Regression for Korean rivers (Hydrocore in-house, 2023)
         TWATER(ich) = -3. + (36.1 + 3.) / (1. + exp(0.1 * (13.3 - &
            (SedCNP_hydro%SFCTMPavg(gwid) - 273.15))))
      else
         ! Use average of basin 1 as a fallback
         TWATER(ich) = -3. + (36.1 + 3.) / (1. + exp(0.1 * (13.3 - &
            (SedCNP_hydro%SFCTMPavg(1) - 273.15))))
      endif

      ! Run CNP cycles if there is water in the channel
      !if (WStorage > 0.0) then  !BK20260318 commented out obsolete check 

         ! Carbon cycle for the channel
         call ChannelCarbonCycle(ich,itime,Qdsch,WtopWdth,WStorage)
         call CheckCMB_Ch(ich)

         ! Nitrogen cycle for the channel
         call ChannelNitrogenCycle(ich,itime,Qdsch,WtopWdth,WStorage)
         call CheckNMB_Ch(ich)

         ! Phosphorus cycle for the channel
         call ChannelPhosphorusCycle(ich,itime,Qdsch,WtopWdth,WStorage)
         call CheckPMB_Ch(ich)

      !endif

   end subroutine RunCNP_Ch


   subroutine Calculate_Basin_Averages()
      implicit none
      integer :: gwid, ix, jx
      real(8), dimension(domain%nbasin) :: count_t2m, count_swdown
  
      ! Initialize arrays
      SedCNP_hydro%SFCTMPavg = 0.0
      SedCNP_hydro%SWDOWNavg = 0.0
      count_t2m = 0.0
      count_swdown = 0.0
  
      ! Loop over grid cells to accumulate values for each basin
      do ix = 1, domain%ix
         do jx = 1, domain%jx
            gwid = SedCNP_hydro%gwbasin(ix,jx)
            if (gwid > 0 .and. gwid <= domain%nbasin) then
               if(SedCNP_hydro%T2m(ix,jx) > 0.0) then
                  SedCNP_hydro%SFCTMPavg(gwid) = &
                     SedCNP_hydro%SFCTMPavg(gwid) + SedCNP_hydro%T2m(ix,jx)
                  count_t2m(gwid) = count_t2m(gwid) + 1.0
               endif
               if(SedCNP_hydro%swdown(ix,jx) > 0.0) then
                  SedCNP_hydro%SWDOWNavg(gwid) = &
                     SedCNP_hydro%SWDOWNavg(gwid) + SedCNP_hydro%swdown(ix,jx)
                  count_swdown(gwid) = count_swdown(gwid) + 1.0
               endif
            endif
         enddo
      enddo

      ! Loop over basins to compute the final average
      do gwid = 1, domain%nbasin 
          if(count_t2m(gwid) > 0.0) then
              SedCNP_hydro%SFCTMPavg(gwid) = &
                 SedCNP_hydro%SFCTMPavg(gwid) / count_t2m(gwid)
          else
              SedCNP_hydro%SFCTMPavg(gwid) = -999.0
          endif
  
          if(count_swdown(gwid) > 0.0) then
              SedCNP_hydro%SWDOWNavg(gwid) = &
                 SedCNP_hydro%SWDOWNavg(gwid) / count_swdown(gwid)
          else
              SedCNP_hydro%SWDOWNavg(gwid) = 0.0
          endif
      enddo
   end subroutine Calculate_Basin_Averages

end module CNPmain