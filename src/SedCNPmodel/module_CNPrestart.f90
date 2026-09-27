!=====||__WHQ5403__||=====!
!
!> @brief SedCNP restart (initial condition) files in NetCDF.
!! @details
!! Write: every WHQ_RESTART_DT hours of model time and at the end of the simulation, nine files
!!   <WHQOUT_dir>/RESTART_<kind>.YYYYMMDDHH_DOMAIN1, kind = SOC, SON, SOP (soil C/N/P),
!!   GWC, GWN, GWP (groundwater), CHC, CHN, CHP (channels); YYYYMMDDHH is the valid time of the
!!   state (end of the time-step). Units: soil kg ha-1, groundwater and channels kg.
!! Read: at the first time-step, from RESTART_FILENAME_<kind> in whq.namelist (the next run starts
!!   with SedCNP_START_* = the valid time of the restart files). A kind without a file name falls
!!   back to the previous text initial condition files (Cini*_file) or a cold start.
!! The variables are those of the previous text files (module_WriteCNPini / module_ReadCNPini).
!
module module_CNPrestart

   use netcdf
   use module_SedCNPvariables
   use SedCNP_config,     only: SedCNPmodel
   use module_hydro_stop, only: hydro_stop

   implicit none

   private
   public :: Write_CNP_Restart, Read_CNP_Restart, CNP_Restart_Due

   character(len=*), parameter :: kinds(9) = ['SOC','SON','SOP','GWC','GWN','GWP','CHC','CHN','CHP']

   contains

   !> @brief .true. if restart files are written at the end of time-step itime.
   logical function CNP_Restart_Due(itime)
      implicit none
      integer, intent(in) :: itime
      integer(8)          :: tsec, rsec

      CNP_Restart_Due = (itime == domain%ntime_sedcnp)
      if (SedCNPmodel%whq_restart_dt > 0) then
         tsec = int(itime, 8) * int(SedCNPmodel%SedCNP_timestep, 8)
         !rsec = int(SedCNPmodel%whq_restart_dt, 8) * 60_8
         rsec = int(SedCNPmodel%whq_restart_dt, 8) * 3600_8   ! WHQ_RESTART_DT [hours]
         if (mod(tsec, rsec) == 0) CNP_Restart_Due = .true.
      endif
   end function CNP_Restart_Due


   !> @brief Writes the nine restart files for the current state (valid time dateSedCNP%olddate).
   subroutine Write_CNP_Restart()

      implicit none

      character(len=256) :: fname(9), stamp
      integer            :: ncid, k
      real(8)            :: ha

      ha = domain%areaxy / 10000.0d0   ! cell area [ha]
      call system('mkdir -p '//trim(SedCNPmodel%WHQOUT_dir))
      stamp = dateSedCNP%olddate(1:4)//dateSedCNP%olddate(6:7)//dateSedCNP%olddate(9:10)// &
              dateSedCNP%olddate(12:13)
      do k = 1, 9
         !fname(k) = trim(SedCNPmodel%WHQOUT_dir)//'/RESTART_'//kinds(k)//'.'//trim(stamp)//'_DOMAIN1.nc'
         fname(k) = trim(SedCNPmodel%WHQOUT_dir)//'/RESTART_'//kinds(k)//'.'//trim(stamp)//'_DOMAIN1'   !WHQ5403 no .nc
      enddo

      ! --- soil C
      call create_soil(fname(1), ncid)
      call put2(ncid, 'So_LPOCLIT', 'soil litter labile particulate organic C (layer 1)', 'kg-C ha-1', So_LPOCLIT/ha)
      call put2(ncid, 'So_LPOCEXC', 'soil excreta labile particulate organic C (layer 1)', 'kg-C ha-1', So_LPOCEXC/ha)
      call put2(ncid, 'So_LPOCMAN', 'soil manure labile particulate organic C (layer 1)', 'kg-C ha-1', So_LPOCMAN/ha)
      call put3(ncid, 'So_LPOCRES', 'soil residue labile particulate organic C', 'kg-C ha-1', So_LPOCRES/ha)
      call put3(ncid, 'So_RPOC', 'soil refractory particulate organic C', 'kg-C ha-1', So_RPOC/ha)
      call put3(ncid, 'So_LDOC', 'soil labile dissolved organic C', 'kg-C ha-1', So_LDOC/ha)
      call put3(ncid, 'So_RDOC', 'soil refractory dissolved organic C', 'kg-C ha-1', So_RDOC/ha)
      call put3(ncid, 'So_MBMC', 'soil microbial biomass C', 'kg-C ha-1', So_MBMC/ha)
      call nc_close(ncid, fname(1))
      ! --- soil N
      call create_soil(fname(2), ncid)
      call put2(ncid, 'So_LPONLIT', 'soil litter labile particulate organic N (layer 1)', 'kg-N ha-1', So_LPONLIT/ha)
      call put2(ncid, 'So_LPONEXC', 'soil excreta labile particulate organic N (layer 1)', 'kg-N ha-1', So_LPONEXC/ha)
      call put2(ncid, 'So_LPONMAN', 'soil manure labile particulate organic N (layer 1)', 'kg-N ha-1', So_LPONMAN/ha)
      call put3(ncid, 'So_LPONRES', 'soil residue labile particulate organic N', 'kg-N ha-1', So_LPONRES/ha)
      call put3(ncid, 'So_RPON', 'soil refractory particulate organic N', 'kg-N ha-1', So_RPON/ha)
      call put3(ncid, 'So_LDON', 'soil labile dissolved organic N', 'kg-N ha-1', So_LDON/ha)
      call put3(ncid, 'So_RDON', 'soil refractory dissolved organic N', 'kg-N ha-1', So_RDON/ha)
      call put3(ncid, 'So_MBMN', 'soil microbial biomass N', 'kg-N ha-1', So_MBMN/ha)
      call put3(ncid, 'So_NH4', 'soil ammonium', 'kg-N ha-1', So_NH4/ha)
      call put3(ncid, 'So_NO3', 'soil nitrate', 'kg-N ha-1', So_NO3/ha)
      call nc_close(ncid, fname(2))
      ! --- soil P
      call create_soil(fname(3), ncid)
      call put2(ncid, 'So_LPOPLIT', 'soil litter labile particulate organic P (layer 1)', 'kg-P ha-1', So_LPOPLIT/ha)
      call put2(ncid, 'So_LPOPEXC', 'soil excreta labile particulate organic P (layer 1)', 'kg-P ha-1', So_LPOPEXC/ha)
      call put2(ncid, 'So_LPOPMAN', 'soil manure labile particulate organic P (layer 1)', 'kg-P ha-1', So_LPOPMAN/ha)
      call put3(ncid, 'So_LPOPRES', 'soil residue labile particulate organic P', 'kg-P ha-1', So_LPOPRES/ha)
      call put3(ncid, 'So_RPOP', 'soil refractory particulate organic P', 'kg-P ha-1', So_RPOP/ha)
      call put3(ncid, 'So_LDOP', 'soil labile dissolved organic P', 'kg-P ha-1', So_LDOP/ha)
      call put3(ncid, 'So_RDOP', 'soil refractory dissolved organic P', 'kg-P ha-1', So_RDOP/ha)
      call put3(ncid, 'So_MBMP', 'soil microbial biomass P', 'kg-P ha-1', So_MBMP/ha)
      call put3(ncid, 'So_PO4', 'soil phosphate', 'kg-P ha-1', So_PO4/ha)
      call put3(ncid, 'So_PIPA', 'soil active particulate inorganic P', 'kg-P ha-1', So_PIPA/ha)
      call put3(ncid, 'So_PIPS', 'soil stable particulate inorganic P', 'kg-P ha-1', So_PIPS/ha)
      call nc_close(ncid, fname(3))

      ! --- groundwater
      call create_gw(fname(4), ncid)
      call put1g(ncid, 'Gw_LDOC', 'groundwater labile dissolved organic C', 'kg-C', Gw_LDOC)
      call put1g(ncid, 'Gw_RDOC', 'groundwater refractory dissolved organic C', 'kg-C', Gw_RDOC)
      call nc_close(ncid, fname(4))
      call create_gw(fname(5), ncid)
      call put1g(ncid, 'Gw_LDON', 'groundwater labile dissolved organic N', 'kg-N', Gw_LDON)
      call put1g(ncid, 'Gw_RDON', 'groundwater refractory dissolved organic N', 'kg-N', Gw_RDON)
      call put1g(ncid, 'Gw_NH4', 'groundwater ammonium', 'kg-N', Gw_NH4)
      call put1g(ncid, 'Gw_NO3', 'groundwater nitrate', 'kg-N', Gw_NO3)
      call nc_close(ncid, fname(5))
      call create_gw(fname(6), ncid)
      call put1g(ncid, 'Gw_LDOP', 'groundwater labile dissolved organic P', 'kg-P', Gw_LDOP)
      call put1g(ncid, 'Gw_RDOP', 'groundwater refractory dissolved organic P', 'kg-P', Gw_RDOP)
      call put1g(ncid, 'Gw_PO4', 'groundwater phosphate', 'kg-P', Gw_PO4)
      call nc_close(ncid, fname(6))

      ! --- channels
      call create_ch(fname(7), ncid)
      call put2c(ncid, 'Ch_ALGC', 'channel algae C (3 species)', 'kg-C', Ch_ALGC)
      call put1c(ncid, 'Ch_ZOOC', 'channel zooplankton C', 'kg-C', Ch_ZOOC)
      call put1c(ncid, 'Ch_MBMC', 'channel microbial biomass C', 'kg-C', Ch_MBMC)
      call put1c(ncid, 'Ch_DIC', 'channel dissolved inorganic C', 'kg-C', Ch_DIC)
      call put1c(ncid, 'Ch_LPOC', 'channel labile particulate organic C', 'kg-C', Ch_LPOC)
      call put1c(ncid, 'Ch_RPOC', 'channel refractory particulate organic C', 'kg-C', Ch_RPOC)
      call put1c(ncid, 'Ch_LDOC', 'channel labile dissolved organic C', 'kg-C', Ch_LDOC)
      call put1c(ncid, 'Ch_RDOC', 'channel refractory dissolved organic C', 'kg-C', Ch_RDOC)
      call put1c(ncid, 'Ch_POCDEPOSIT', 'channel deposited particulate organic C', 'kg-C', Ch_POCDEPOSIT)
      call nc_close(ncid, fname(7))
      call create_ch(fname(8), ncid)
      call put2c(ncid, 'Ch_ALGN', 'channel algae N (3 species)', 'kg-N', Ch_ALGN)
      call put1c(ncid, 'Ch_ZOON', 'channel zooplankton N', 'kg-N', Ch_ZOON)
      call put1c(ncid, 'Ch_MBMN', 'channel microbial biomass N', 'kg-N', Ch_MBMN)
      call put1c(ncid, 'Ch_LPON', 'channel labile particulate organic N', 'kg-N', Ch_LPON)
      call put1c(ncid, 'Ch_RPON', 'channel refractory particulate organic N', 'kg-N', Ch_RPON)
      call put1c(ncid, 'Ch_LDON', 'channel labile dissolved organic N', 'kg-N', Ch_LDON)
      call put1c(ncid, 'Ch_RDON', 'channel refractory dissolved organic N', 'kg-N', Ch_RDON)
      call put1c(ncid, 'Ch_NH4', 'channel ammonium', 'kg-N', Ch_NH4)
      call put1c(ncid, 'Ch_NO3', 'channel nitrate', 'kg-N', Ch_NO3)
      call put1c(ncid, 'Ch_PONDEPOSIT', 'channel deposited particulate organic N', 'kg-N', Ch_PONDEPOSIT)
      call nc_close(ncid, fname(8))
      call create_ch(fname(9), ncid)
      call put2c(ncid, 'Ch_ALGP', 'channel algae P (3 species)', 'kg-P', Ch_ALGP)
      call put1c(ncid, 'Ch_ZOOP', 'channel zooplankton P', 'kg-P', Ch_ZOOP)
      call put1c(ncid, 'Ch_MBMP', 'channel microbial biomass P', 'kg-P', Ch_MBMP)
      call put1c(ncid, 'Ch_LPOP', 'channel labile particulate organic P', 'kg-P', Ch_LPOP)
      call put1c(ncid, 'Ch_RPOP', 'channel refractory particulate organic P', 'kg-P', Ch_RPOP)
      call put1c(ncid, 'Ch_LDOP', 'channel labile dissolved organic P', 'kg-P', Ch_LDOP)
      call put1c(ncid, 'Ch_RDOP', 'channel refractory dissolved organic P', 'kg-P', Ch_RDOP)
      call put1c(ncid, 'Ch_PO4', 'channel phosphate', 'kg-P', Ch_PO4)
      call put1c(ncid, 'Ch_PIPA', 'channel active particulate inorganic P', 'kg-P', Ch_PIPA)
      call put1c(ncid, 'Ch_PIPS', 'channel stable particulate inorganic P', 'kg-P', Ch_PIPS)
      call put1c(ncid, 'Ch_POPDEPOSIT', 'channel deposited particulate organic P', 'kg-P', Ch_POPDEPOSIT)
      call nc_close(ncid, fname(9))

      !WHQ5403 not printed (user request)
      !call progress_clear()   !WHQ5403
      !write(6,*) 'INFO: SedCNP restart files written: ', trim(SedCNPmodel%WHQOUT_dir)//'/RESTART_*.'// &
      !           trim(stamp)//'_DOMAIN1'

   end subroutine Write_CNP_Restart


   !> @brief Reads the restart files given in whq.namelist into the initial state (*0 arrays).
   !! Called at the first time-step after ReadCNPini_So/Gw/Ch (text files / cold start).
   subroutine Read_CNP_Restart()

      implicit none

      character(len=256) :: f
      integer            :: ncid
      real(8)            :: ha

      ha = domain%areaxy / 10000.0d0

      f = SedCNPmodel%restart_filename_soc
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'soil')
         call get2(ncid, f, 'So_LPOCLIT', So_LPOCLIT0, ha); call get2(ncid, f, 'So_LPOCEXC', So_LPOCEXC0, ha)
         call get2(ncid, f, 'So_LPOCMAN', So_LPOCMAN0, ha); call get3(ncid, f, 'So_LPOCRES', So_LPOCRES0, ha)
         call get3(ncid, f, 'So_RPOC', So_RPOC0, ha);       call get3(ncid, f, 'So_LDOC', So_LDOC0, ha)
         call get3(ncid, f, 'So_RDOC', So_RDOC0, ha);       call get3(ncid, f, 'So_MBMC', So_MBMC0, ha)
         call nc_close(ncid, f)
      endif
      f = SedCNPmodel%restart_filename_son
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'soil')
         call get2(ncid, f, 'So_LPONLIT', So_LPONLIT0, ha); call get2(ncid, f, 'So_LPONEXC', So_LPONEXC0, ha)
         call get2(ncid, f, 'So_LPONMAN', So_LPONMAN0, ha); call get3(ncid, f, 'So_LPONRES', So_LPONRES0, ha)
         call get3(ncid, f, 'So_RPON', So_RPON0, ha);       call get3(ncid, f, 'So_LDON', So_LDON0, ha)
         call get3(ncid, f, 'So_RDON', So_RDON0, ha);       call get3(ncid, f, 'So_MBMN', So_MBMN0, ha)
         call get3(ncid, f, 'So_NH4', So_NH40, ha);         call get3(ncid, f, 'So_NO3', So_NO30, ha)
         call nc_close(ncid, f)
      endif
      f = SedCNPmodel%restart_filename_sop
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'soil')
         call get2(ncid, f, 'So_LPOPLIT', So_LPOPLIT0, ha); call get2(ncid, f, 'So_LPOPEXC', So_LPOPEXC0, ha)
         call get2(ncid, f, 'So_LPOPMAN', So_LPOPMAN0, ha); call get3(ncid, f, 'So_LPOPRES', So_LPOPRES0, ha)
         call get3(ncid, f, 'So_RPOP', So_RPOP0, ha);       call get3(ncid, f, 'So_LDOP', So_LDOP0, ha)
         call get3(ncid, f, 'So_RDOP', So_RDOP0, ha);       call get3(ncid, f, 'So_MBMP', So_MBMP0, ha)
         call get3(ncid, f, 'So_PO4', So_PO40, ha);         call get3(ncid, f, 'So_PIPA', So_PIPA0, ha)
         call get3(ncid, f, 'So_PIPS', So_PIPS0, ha)
         call nc_close(ncid, f)
      endif

      f = SedCNPmodel%restart_filename_gwc
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'gw')
         call get1(ncid, f, 'Gw_LDOC', Gw_LDOC0); call get1(ncid, f, 'Gw_RDOC', Gw_RDOC0)
         call nc_close(ncid, f)
      endif
      f = SedCNPmodel%restart_filename_gwn
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'gw')
         call get1(ncid, f, 'Gw_LDON', Gw_LDON0); call get1(ncid, f, 'Gw_RDON', Gw_RDON0)
         call get1(ncid, f, 'Gw_NH4', Gw_NH40);   call get1(ncid, f, 'Gw_NO3', Gw_NO30)
         call nc_close(ncid, f)
      endif
      f = SedCNPmodel%restart_filename_gwp
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'gw')
         call get1(ncid, f, 'Gw_LDOP', Gw_LDOP0); call get1(ncid, f, 'Gw_RDOP', Gw_RDOP0)
         call get1(ncid, f, 'Gw_PO4', Gw_PO40)
         call nc_close(ncid, f)
      endif

      ! channels (minimum algae/zooplankton as in ReadCNPini_Ch)
      f = SedCNPmodel%restart_filename_chc
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'ch')
         call get2c(ncid, f, 'Ch_ALGC', Ch_ALGC0); call get1(ncid, f, 'Ch_ZOOC', Ch_ZOOC0)
         call get1(ncid, f, 'Ch_MBMC', Ch_MBMC0);  call get1(ncid, f, 'Ch_DIC', Ch_DIC0)
         call get1(ncid, f, 'Ch_LPOC', Ch_LPOC0);  call get1(ncid, f, 'Ch_RPOC', Ch_RPOC0)
         call get1(ncid, f, 'Ch_LDOC', Ch_LDOC0);  call get1(ncid, f, 'Ch_RDOC', Ch_RDOC0)
         call get1(ncid, f, 'Ch_POCDEPOSIT', Ch_POCDEPOSIT0)
         call nc_close(ncid, f)
         Ch_ALGC0 = max(Ch_ALGC0, 0.00106d0); Ch_ZOOC0 = max(Ch_ZOOC0, 0.00106d0)
      endif
      f = SedCNPmodel%restart_filename_chn
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'ch')
         call get2c(ncid, f, 'Ch_ALGN', Ch_ALGN0); call get1(ncid, f, 'Ch_ZOON', Ch_ZOON0)
         call get1(ncid, f, 'Ch_MBMN', Ch_MBMN0);  call get1(ncid, f, 'Ch_LPON', Ch_LPON0)
         call get1(ncid, f, 'Ch_RPON', Ch_RPON0);  call get1(ncid, f, 'Ch_LDON', Ch_LDON0)
         call get1(ncid, f, 'Ch_RDON', Ch_RDON0);  call get1(ncid, f, 'Ch_NH4', Ch_NH40)
         call get1(ncid, f, 'Ch_NO3', Ch_NO30);    call get1(ncid, f, 'Ch_PONDEPOSIT', Ch_PONDEPOSIT0)
         call nc_close(ncid, f)
         Ch_ALGN0 = max(Ch_ALGN0, 0.00016d0); Ch_ZOON0 = max(Ch_ZOON0, 0.00016d0)
      endif
      f = SedCNPmodel%restart_filename_chp
      if (len_trim(f) > 0) then
         call open_check(f, ncid, 'ch')
         call get2c(ncid, f, 'Ch_ALGP', Ch_ALGP0); call get1(ncid, f, 'Ch_ZOOP', Ch_ZOOP0)
         call get1(ncid, f, 'Ch_MBMP', Ch_MBMP0);  call get1(ncid, f, 'Ch_LPOP', Ch_LPOP0)
         call get1(ncid, f, 'Ch_RPOP', Ch_RPOP0);  call get1(ncid, f, 'Ch_LDOP', Ch_LDOP0)
         call get1(ncid, f, 'Ch_RDOP', Ch_RDOP0);  call get1(ncid, f, 'Ch_PO4', Ch_PO40)
         call get1(ncid, f, 'Ch_PIPA', Ch_PIPA0);  call get1(ncid, f, 'Ch_PIPS', Ch_PIPS0)
         call get1(ncid, f, 'Ch_POPDEPOSIT', Ch_POPDEPOSIT0)
         call nc_close(ncid, f)
         Ch_ALGP0 = max(Ch_ALGP0, 0.00001d0); Ch_ZOOP0 = max(Ch_ZOOP0, 0.00001d0)
      endif

   end subroutine Read_CNP_Restart


   ! ------------------------------------------------------------------ file creation
   subroutine create_common(fname, ncid)
      character(len=*), intent(in) :: fname
      integer, intent(out)         :: ncid
      call chk(nf90_create(trim(fname), NF90_NETCDF4, ncid), 'create', fname)
      call chk(nf90_put_att(ncid, NF90_GLOBAL, 'title', 'WRF-HydroQual SedCNP restart'), 'att', fname)
      call chk(nf90_put_att(ncid, NF90_GLOBAL, 'valid_time', trim(dateSedCNP%olddate)), 'att', fname)
   end subroutine create_common

   subroutine create_soil(fname, ncid)
      character(len=*), intent(in) :: fname
      integer, intent(out)         :: ncid
      integer                      :: d
      call create_common(fname, ncid)
      call chk(nf90_def_dim(ncid, 'x', domain%ix, d), 'dim x', fname)
      call chk(nf90_def_dim(ncid, 'soil_layers', domain%nsl, d), 'dim soil_layers', fname)
      call chk(nf90_def_dim(ncid, 'y', domain%jx, d), 'dim y', fname)
      call chk(nf90_enddef(ncid), 'enddef', fname)
   end subroutine create_soil

   subroutine create_gw(fname, ncid)
      character(len=*), intent(in) :: fname
      integer, intent(out)         :: ncid
      integer                      :: d, v, b
      call create_common(fname, ncid)
      call chk(nf90_def_dim(ncid, 'nbasin', domain%nbasin, d), 'dim nbasin', fname)
      call chk(nf90_def_var(ncid, 'gwid', NF90_INT, (/d/), v), 'var gwid', fname)
      call chk(nf90_put_att(ncid, v, 'long_name', 'gw basin ID'), 'att', fname)
      call chk(nf90_enddef(ncid), 'enddef', fname)
      call chk(nf90_put_var(ncid, v, (/(b, b = 1, domain%nbasin)/)), 'put gwid', fname)
   end subroutine create_gw

   subroutine create_ch(fname, ncid)
      character(len=*), intent(in) :: fname
      integer, intent(out)         :: ncid
      integer                      :: d, a, v
      call create_common(fname, ncid)
      call chk(nf90_def_dim(ncid, 'nch', domain%nch, d), 'dim nch', fname)
      call chk(nf90_def_dim(ncid, 'nalg', nalg, a), 'dim nalg', fname)
      call chk(nf90_def_var(ncid, 'link', NF90_INT, (/d/), v), 'var link', fname)
      call chk(nf90_put_att(ncid, v, 'long_name', 'link ID (Route_Link.nc)'), 'att', fname)
      call chk(nf90_enddef(ncid), 'enddef', fname)
      call chk(nf90_put_var(ncid, v, SedCNP_hydro%linkID), 'put link', fname)
   end subroutine create_ch

   ! ------------------------------------------------------------------ write variables
   subroutine def_var(ncid, fname, name, lname, units, dnames, v)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: fname, name, lname, units, dnames(:)
      integer, intent(out)         :: v
      integer                      :: dids(size(dnames)), n
      do n = 1, size(dnames)
         call chk(nf90_inq_dimid(ncid, trim(dnames(n)), dids(n)), 'dim '//trim(dnames(n)), fname)
      enddo
      call chk(nf90_redef(ncid), 'redef', fname)
      call chk(nf90_def_var(ncid, name, NF90_DOUBLE, dids, v), 'def '//name, fname)
      call chk(nf90_put_att(ncid, v, 'long_name', lname), 'att', fname)
      call chk(nf90_put_att(ncid, v, 'units', units), 'att', fname)
      call chk(nf90_enddef(ncid), 'enddef', fname)
   end subroutine def_var

   subroutine put2(ncid, name, lname, units, a)       ! (ix,jx)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: name, lname, units
      real(8), intent(in)          :: a(:,:)
      integer                      :: v
      call def_var(ncid, 'restart', name, lname, units, [character(len=16) :: 'x', 'y'], v)
      call chk(nf90_put_var(ncid, v, a), 'put '//name, 'restart')
   end subroutine put2

   subroutine put3(ncid, name, lname, units, a)       ! (ix,nsl,jx)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: name, lname, units
      real(8), intent(in)          :: a(:,:,:)
      integer                      :: v
      call def_var(ncid, 'restart', name, lname, units, [character(len=16) :: 'x', 'soil_layers', 'y'], v)
      call chk(nf90_put_var(ncid, v, a), 'put '//name, 'restart')
   end subroutine put3

   subroutine put1g(ncid, name, lname, units, a)      ! (nbasin)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: name, lname, units
      real(8), intent(in)          :: a(:)
      integer                      :: v
      call def_var(ncid, 'restart', name, lname, units, [character(len=16) :: 'nbasin'], v)
      call chk(nf90_put_var(ncid, v, a), 'put '//name, 'restart')
   end subroutine put1g

   subroutine put1c(ncid, name, lname, units, a)      ! (nch)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: name, lname, units
      real(8), intent(in)          :: a(:)
      integer                      :: v
      call def_var(ncid, 'restart', name, lname, units, [character(len=16) :: 'nch'], v)
      call chk(nf90_put_var(ncid, v, a), 'put '//name, 'restart')
   end subroutine put1c

   subroutine put2c(ncid, name, lname, units, a)      ! (nch,nalg)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: name, lname, units
      real(8), intent(in)          :: a(:,:)
      integer                      :: v
      call def_var(ncid, 'restart', name, lname, units, [character(len=16) :: 'nch', 'nalg'], v)
      call chk(nf90_put_var(ncid, v, a), 'put '//name, 'restart')
   end subroutine put2c

   ! ------------------------------------------------------------------ read
   subroutine open_check(fname, ncid, kind)
      character(len=*), intent(in) :: fname, kind
      integer, intent(out)         :: ncid
      integer, allocatable         :: link(:)
      integer                      :: v

      if (nf90_open(trim(fname), NF90_NOWRITE, ncid) /= 0) &
         call hydro_stop('Read_CNP_Restart: cannot open restart file '//trim(fname))
      select case (kind)
      case ('soil')
         call check_dim(ncid, fname, 'x', domain%ix)
         call check_dim(ncid, fname, 'y', domain%jx)
         call check_dim(ncid, fname, 'soil_layers', domain%nsl)
      case ('gw')
         call check_dim(ncid, fname, 'nbasin', domain%nbasin)
      case ('ch')
         call check_dim(ncid, fname, 'nch', domain%nch)
         allocate(link(domain%nch))
         call chk(nf90_inq_varid(ncid, 'link', v), 'var link', fname)
         call chk(nf90_get_var(ncid, v, link), 'get link', fname)
         if (any(link /= SedCNP_hydro%linkID)) &
            call hydro_stop('Read_CNP_Restart: link IDs in '//trim(fname)//' differ from Route_Link.nc')
         deallocate(link)
      end select
      call progress_clear()   !WHQ5403
      write(6,*) 'INFO: SedCNP initial conditions read from ', trim(fname)
   end subroutine open_check

   subroutine check_dim(ncid, fname, name, n)
      integer, intent(in)          :: ncid, n
      character(len=*), intent(in) :: fname, name
      integer                      :: d, len
      call chk(nf90_inq_dimid(ncid, name, d), 'dim '//name, fname)
      call chk(nf90_inquire_dimension(ncid, d, len=len), 'dim '//name, fname)
      if (len /= n) call hydro_stop('Read_CNP_Restart: dimension '//name//' of '//trim(fname)// &
                                    ' does not match the model domain')
   end subroutine check_dim

   subroutine get2(ncid, fname, name, a, fac)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: fname, name
      real(8), intent(out)         :: a(:,:)
      real(8), intent(in)          :: fac
      integer                      :: v
      call chk(nf90_inq_varid(ncid, name, v), 'var '//name, fname)
      call chk(nf90_get_var(ncid, v, a), 'get '//name, fname)
      a = a * fac
   end subroutine get2

   subroutine get3(ncid, fname, name, a, fac)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: fname, name
      real(8), intent(out)         :: a(:,:,:)
      real(8), intent(in)          :: fac
      integer                      :: v
      call chk(nf90_inq_varid(ncid, name, v), 'var '//name, fname)
      call chk(nf90_get_var(ncid, v, a), 'get '//name, fname)
      a = a * fac
   end subroutine get3

   subroutine get1(ncid, fname, name, a)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: fname, name
      real(8), intent(out)         :: a(:)
      integer                      :: v
      call chk(nf90_inq_varid(ncid, name, v), 'var '//name, fname)
      call chk(nf90_get_var(ncid, v, a), 'get '//name, fname)
   end subroutine get1

   subroutine get2c(ncid, fname, name, a)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: fname, name
      real(8), intent(out)         :: a(:,:)
      integer                      :: v
      call chk(nf90_inq_varid(ncid, name, v), 'var '//name, fname)
      call chk(nf90_get_var(ncid, v, a), 'get '//name, fname)
   end subroutine get2c

   ! ------------------------------------------------------------------ utilities
   subroutine nc_close(ncid, fname)
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: fname
      call chk(nf90_close(ncid), 'close', fname)
   end subroutine nc_close

   subroutine chk(status, what, fname)
      integer, intent(in)          :: status
      character(len=*), intent(in) :: what, fname
      if (status /= NF90_NOERR) then
         call progress_clear()   !WHQ5403
         write(6,*) 'ERROR: SedCNP restart: ', trim(what), ' (', trim(fname), '): ', trim(nf90_strerror(status))
         call hydro_stop('SedCNP restart file error')
      endif
   end subroutine chk

end module module_CNPrestart
!
!=====||__WHQ5403__||=====!
