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

module module_overlandSed

   !use module_SedCNPvariables, only: overSed, SedCNP_hydro
   use module_SedCNPvariables  !BK20240619
!=====||__WHQ5403q__||=====!
!
   use SedCNP_config,          only: SedCNPmodel
!
!=====||__WHQ5403q__||=====!
   use module_hydro_stop,      only: HYDRO_stop !BK20251130
   !use NOAHMP_TABLES, only: ISWATER_TABLE   !BK20251208 Reverted


   !integer             :: headerPrint = 1          !hlcho header print 여부 확인 !BK20251016 commented out
   !integer             :: headerPrint_vege = 1     !hlcho header print 여부 확인 !BK20251016 commented out
   
   ! Parameters to be read from file !BK20250831
   integer, parameter :: num_lcover_types = 27  !*** hardcoded
   !real(8), dimension(:), allocatable :: crVegeCover_param  !BK20260605
   !real(8), dimension(:), allocatable :: crCropHeight_param !BK20260605
   !real(8), dimension(:), allocatable :: Sf_param  !BK20260525
   real(8), dimension(:), allocatable :: areaP_Ratio_param
   real(8), dimension(:), allocatable :: manningsn_param
   logical :: params_read = .false.

   !save headerPrint        !BK20251016 commented out
   !save headerPrint_vege   !BK20251016 commented out
   !save crVegeCover_param, crCropHeight_param, Sf_param, areaP_Ratio_param, manningsn_param, params_read  !BK20250831 !BK20251016 commented out
   real(8), dimension(4) :: k_hs_dstc_param !BK20251130
   real(8), dimension(4) :: k_setV_param    !BK20251130
   logical :: paddy_params_read = .false.   !BK20251130

   real(8)               :: Gs = 2.65  !specific gravity of sediment !BK20260218
   real(8)               :: g = 9.8    !gravitational acceleration [m/s2]
   real(8)               :: v = 1.0E-6 !kinematic viscosity of water at 20 deg-C [m2/s] 

   save params_read, k_hs_dstc_param, k_setV_param, paddy_params_read  !BK20251016 !BK20251130

   contains

   subroutine overlandSedTransport(i,j)

      implicit none

      integer, intent(in) :: i, j 
      real(8)                :: dt               !time-step [sec]
      real(8)                :: dx               !grid cell size [m]
      real(8)                :: OLWDepth         !overland water depth [mm]   !BK20240421
      real(8)                :: IRRIDepth        !irrigation water depth [mm]
      real(8)                :: Qrain            !rainfall rate [mm/s]    !BK20240421
      real(8)                :: Qsurf            !surface runoff [mm/s]   !BK20240421
      real(8)                :: Qinfl            !infiltration [mm/s]  --> need to be updated (kyh_debug)  !BK20240421
      real(8)                :: Qevap            !evaporation [mm/s] --> need to be updated (kyh_debug)   !BK20240421
      real(8)                :: crVegeCover
      real(8)                :: crCropHeight     !vegetation height
      real(8)                :: Sf               !overland slope
      real(8)                :: areaP_Ratio      !impervious areas ratio [0~1]	 
      integer             :: ips
      real(8), dimension(4)  :: dp               !particle diameter [μm]
      real(8), dimension(4)  :: vc               !critical transport velocity for each grid size (m/s)
      real(8), dimension(4)  :: k_hs_dstc        !coefficient of disturbance for clay, silt, fine sand and coarse sand in the paddy field (-)
      real(8), dimension(4)  :: k_setV           !coefficient of settling velocity for clay, silt, fine sand and coarse sand in the paddy field (-)
      !real(8), dimension(4)  :: topSoilFractions 
      !real(8)                :: vgtyp           !vegetation type  !BK20240205
      integer                :: lcover           !land cover code   !BK20231006
      logical                :: is_paddy_cell    ! Flag for paddy cells !BK20260517
      real(8)                :: manningsn        ! Manning's n
      
      !particle diameter
      !hard-coded; should be relocated or moved to config file
      dp(1) = 1.      !clay [μm]
      dp(2) = 10.     !silt [μm]
      dp(3) = 100.    !fine sand [μm]
      dp(4) = 1000.   !coarse sand  [μm]

      !hard-coded; should be relocated or moved to config file hlcho 
      !critical velocity
      vc(1) = 0.09 / 100.0     !clay 0.01 to 0.1 cm/s -> (m/s)
      vc(2) = 0.40 / 100.0     !silt 0.1 to 1 cm/s
      vc(3) = 1.50 / 100.0     !fine sand 1 to 10 cm/s
      vc(4) = 3.00 / 100.0     !coarse sand 10 to 100 cm/s
      
      if (.not. paddy_params_read) then  !BK20251130  

            !call read_config_paddy()   !BK20251203 commented out

            !Default hardcoded values  !BK20251203
            k_hs_dstc_param(1) = 0.100 / 150.0 !clay
            k_hs_dstc_param(2) = 0.200 / 50.0  !silt
            k_hs_dstc_param(3) = 0.300 / 15.0  !fine sand
            k_hs_dstc_param(4) = 0.400 / 5.0   !coarse sand

            !Default hardcoded values  !BK20251203
            k_setV_param(1) = 1.0 !clay
            k_setV_param(2) = 1.0 !silt
            k_setV_param(3) = 1.0 !fine sand
            k_setV_param(4) = 1.0 !coarse sand

            paddy_params_read = .true.  !BK20260605

      endif

      k_hs_dstc = k_hs_dstc_param
      k_setV = k_setV_param  

      !hard-coded; should be relocated or moved to config file
      !토양 구성비율
      !topSoilFractions(1) = 0.30 !clay  !BK20240205
      !topSoilFractions(2) = 0.35 !silt
      !topSoilFractions(3) = 0.20 !fine sand
      !topSoilFractions(4) = 0.15 !coarse sand

      !get model time-step 
      dt = real(SedCNPmodel%SedCNP_timestep)

      !get model grid size
      !dx = SedCNPmodel%SedCNP_dx
      dx = SedCNP_hydro%dx(1)   !Y.Kwon(20250621)

      !get surface inputs for a grid
      !Sf = 0.01 !overland slope --> temporary coded, need to be fixed (kyh_debug)

      !get inputs for a grid
      Qrain        = SedCNP_hydro%rainfall(i,j)        ! rainfall [mm/s]   !BK20240421
      OLWDepth     = SedCNP_hydro%olwdepth(i,j)        ! overland water depth in paddy fields [mm]   !BK20240421
      IRRIDepth    = SedCNP_hydro%irridepth(i,j)       ! irrigation water depth [mm] 
      Qsurf        = SedCNP_hydro%runsrf(i,j)          ! surface runoff [mm/s]  !BK20240420 BK20240421 BK20240821
      Qinfl        = SedCNP_hydro%infltrate(i,j)       ! infiltration [mm/s]   !BK20240421
      Qevap        = SedCNP_hydro%evapwater(i,j)       ! evaporation from overland water surface [mm/s]  !BK20240421
      crVegeCover  = SedCNP_hydro%FVEG(i,j)            ! Green vegetation fraction 
      crCropHeight = SedCNP_hydro%CanHt(i,j)           ! vegetation height
      Sf = SedCNP_hydro%slope(i,j)  !BK20260525

      !--- BK20240917
      !if (isnan(SedCNP_hydro%vgtyp(i,j))) then
      if (SedCNP_hydro%vgtyp(i,j) /= SedCNP_hydro%vgtyp(i,j)) then  !BK20250720
         lcover = domain%ISWATER   !BK20251208  Read from parameter file
      else
         lcover = int(SedCNP_hydro%vgtyp(i,j))               ! land cover code   !BK20240205
      endif
      !--- BK20240917

      !--- BK20250831
      if (.not. params_read) then
      call read_SEDparams()
      endif

      if (lcover > 0 .and. lcover <= num_lcover_types) then
         !crVegeCover = crVegeCover_param(lcover)     !BK20260605 commented out (already populated above)
         !crCropHeight = crCropHeight_param(lcover)   !BK20260605 commented out (already populated above)
         !Sf = SedCNP_hydro%slope(i,j)                !BK20260605 commented out (already populated above)
         areaP_Ratio = areaP_Ratio_param(lcover)
         manningsn = manningsn_param(lcover)
      else
            ! Default values for out-of-bounds lcover
            !crVegeCover = 0.2      !BK20260605 commented out (already populated above)
            !crCropHeight = 0.15    !BK20260605 commented out (already populated above)
            !Sf = 0.002             !BK20260605 commented out (already populated above)
            areaP_Ratio = 0.0
            manningsn = 0.030
      endif
      !--- BK20250831

      if (lcover == 7) then !---BK20260517
         is_paddy_cell = .true.
      else
         is_paddy_cell = .false.
      endif !---BK20260517
   
      Qinfl = max(Qinfl, 0.)  !BK20240421

      do ips = 1,4
         !1: clay
         !2: silt
         !3: fine sand
         !4: coarse sand

         !call splashDetachment(i,j,ips,dt,dx,rainRate,OLWDepth,Runsrf,IRRIDepth,&
         !                      topSoilFractions,crVegeCover,crCropHeight,areaP_Ratio,vgtyp)
         call splashDetachment(i,j,ips,dt,dx,Qrain,OLWDepth,Qsurf,IRRIDepth,&   !BK20240421 BK20240421
                                 crVegeCover,crCropHeight,areaP_Ratio,lcover)  !BK20240205

         !if(OLWDepth > 0) then
         !   !call paddySedTransport(i,j,ips,dt,dx,rainRate,OLWDepth,IRRIDepth,&   
         !   !                       Runsrf,infl,evap,k_hs_dstc,k_setV,dp)  
         if (is_paddy_cell) then  !BK20260517
            ! This is a paddy cell. Use the transport logic for ponded conditions,
            ! which is valid for both wet and dry states. 
            call paddySedTransport(i,j,ips,dt,dx,Qrain,OLWDepth,IRRIDepth,&   !BK20240421
                                    Qsurf,Qinfl,Qevap,k_hs_dstc,k_setV,dp)  !BK20240421
         else
            !Julien and Simons, 1985; Rustomji and Prosser, 2000 hlcho OLWDepth 추가
            !call overlandSedTransportJS(i,j,ips,dt,dx,Runsrf,Sf,areaP_Ratio,&
            !                            dp,vc,OLWDepth,manningsn,vgtyp)
            call overlandSedTransportJS(i,j,ips,dt,dx,Qsurf,Sf,areaP_Ratio,&   !BK20240421
                                          dp,vc,OLWDepth,manningsn,lcover)  !BK20240205
            !hlcho_debug
            !call paddyPrint(i,j,ips,dt,dx,rainRate,OLWDepth,IRRIDepth,&
            !                Runsrf,infl,evap,k_hs_dstc,k_setV)  !BK20240420 commented out
         endif
         !call sedimentAggregation(i,j,ips,dx,dt,OLWDepth, Runsrf) !hlcho Runsrf 추가
         call sedimentAggregation(i,j,ips,dx,dt,OLWDepth,Qsurf)   !BK20240421
      enddo
   
   end subroutine overlandSedTransport

   subroutine read_SEDparams()   !BK20250831
      implicit none
      character(len=256) :: file_path
      integer :: iu, i, lcover_idx, ierr
      logical :: fexist

      !BK20251130 Check for user-specified file
      if (len_trim(SedCNPmodel%SEDparams_file) > 0) then
         inquire(file=trim(SedCNPmodel%SEDparams_file), exist=fexist)
         if (fexist) then
            file_path = trim(SedCNPmodel%SEDparams_file)
         else
            fexist = .false.
         endif
      else
         fexist = .false.
      endif

      if (.not. fexist) then
         file_path = './WHQIN/SEDparams.dat'
         inquire(file=trim(file_path), exist=fexist)
         if (.not. fexist) then
               call wrf_error_fatal('SEDparams.dat not found in ./WHQIN/')
         endif
      endif

      !allocate(crVegeCover_param(num_lcover_types), stat=ierr) !BK20260605
      !allocate(crCropHeight_param(num_lcover_types), stat=ierr) !BK20260605
      !allocate(Sf_param(num_lcover_types), stat=ierr)  !BK20260525
      allocate(areaP_Ratio_param(num_lcover_types), stat=ierr)
      allocate(manningsn_param(num_lcover_types), stat=ierr)

      open(newunit=iu, file=trim(file_path), status='old', action='read')

      ! Skip header lines (assuming 2 header lines)
      read(iu,*,iostat=ierr)
      if (ierr /= 0) call wrf_error_fatal('Error reading header from SEDparams.dat')
      read(iu,*,iostat=ierr)
      if (ierr /= 0) call wrf_error_fatal('Error reading header from SEDparams.dat')

      do i = 1, num_lcover_types
         !read(iu, *, iostat=ierr) lcover_idx, crVegeCover_param(i), crCropHeight_param(i), &
         !                           Sf_param(i), areaP_Ratio_param(i), manningsn_param(i)
         !read(iu, *, iostat=ierr) lcover_idx, crVegeCover_param(i), crCropHeight_param(i), &
         !                         areaP_Ratio_param(i), manningsn_param(i)  !BK20260525
         read(iu, *, iostat=ierr) lcover_idx, areaP_Ratio_param(i), manningsn_param(i)  !BK20260605
         if (ierr < 0) then
               call wrf_error_fatal('End of file reached prematurely in SEDparams.dat')
         else if (ierr > 0) then
               call wrf_error_fatal('Format error reading SEDparams.dat')
         endif
      end do

      close(iu)
      params_read = .true.

   end subroutine read_SEDparams  !BK20250831

   ! subroutine read_config_paddy()   !BK20251130
   !    implicit none
   !    character(len=256) :: file_path
   !    integer :: iu, i, ierr
   !    logical :: fexist

   !    ! Default hardcoded values
   !    k_hs_dstc_param(1) = 0.100 / 10.0  !clay
   !    k_hs_dstc_param(2) = 0.200 / 10.0  !silt
   !    k_hs_dstc_param(3) = 0.300 / 10.0  !fine sand
   !    k_hs_dstc_param(4) = 0.400 / 10.0  !coarse sand

   !    k_setV_param(1) = 1.0 !clay
   !    k_setV_param(2) = 1.0 !silt
   !    k_setV_param(3) = 1.0 !fine sand
   !    k_setV_param(4) = 1.0 !coarse sand

   !    ! Check for user-specified file
   !    if (len_trim(SedCNPmodel%config_paddy_file) > 0) then
   !       inquire(file=trim(SedCNPmodel%config_paddy_file), exist=fexist)
   !       if (fexist) then
   !          file_path = trim(SedCNPmodel%config_paddy_file)
   !       else
   !          fexist = .false.
   !       endif
   !    else
   !       fexist = .false.
   !    endif

   !    if (.not. fexist) then
   !       file_path = './WHQIN/config_paddy.dat'
   !       inquire(file=trim(file_path), exist=fexist)
   !    endif

   !    if (fexist) then
   !       open(newunit=iu, file=trim(file_path), status='old', action='read', iostat=ierr)
   !       if (ierr == 0) then
   !             ! Assuming format: 
   !             ! Header
   !             ! k_hs_dstc(1) k_hs_dstc(2) k_hs_dstc(3) k_hs_dstc(4)
   !             ! k_setV(1) k_setV(2) k_setV(3) k_setV(4)
   !             read(iu, *, iostat=ierr) ! Skip header
   !             read(iu, *, iostat=ierr) k_hs_dstc_param(1:4)
   !             read(iu, *, iostat=ierr) k_setV_param(1:4)
   !             close(iu)
   !       endif
   !    endif

   !    paddy_params_read = .true.

   ! end subroutine read_config_paddy


   
   !subroutine splashDetachment(i,j,ips,dt,dx,rainRate,OLWDepth,Runsrf,IRRIDepth,&
   !                            topSoilFractions,crVegeCover,crCropHeight,areaP_Ratio,vgtyp)
   subroutine splashDetachment(i,j,ips,dt,dx,Qrain,OLWDepth,Qsurf,IRRIDepth,&   !BK20240421
                                 crVegeCover,crCropHeight,areaP_Ratio,lcover)  !BK20240205

      implicit none

      integer,            intent(in) :: i, j
      integer,            intent(in) :: ips
      real(8),               intent(in) :: dt              !time-step [sec]
      real(8),               intent(in) :: dx              !grid cell size [m]
      !real(8),               intent(in) :: rainRate        !rainfall rate [mm/s]  
      real(8),               intent(in) :: Qrain           !rainfall rate [mm/s]   !BK20240421
      real(8),               intent(in) :: OLWDepth        !overland water depth [mm]
      !real(8),               intent(in) :: Runsrf          !surface runoff [mm/s]	 
      real(8),               intent(in) :: Qsurf           !surface runoff [mm/s]  !BK20240421
      real(8),               intent(in) :: IRRIDepth       !irrigation water depth [mm]
      !real(8), dimension(4), intent(in) :: topSoilFractions     !BK20240205
      real(8),               intent(in) :: crVegeCover
      real(8),               intent(out):: crCropHeight    !vegetation height
      real(8),               intent(in) :: areaP_Ratio     !pervious areas ration (0~1) 
      !real(8),               intent(in) :: vgtyp           !landuse !BK20240205
      integer,            intent(in) :: lcover          !land cover code   !BK20231006 !BK20240205
      integer                        :: st              !soil texture code  !BK20240205  
      !real(8)                           :: rain            !rainfall during the time-step [mm]  !BK20240421
      real(8)                           :: Abund_ha
      real(8)                           :: areaIrri_Ha
      real(8)                           :: areaP_Ha        !pervious area [ha]
      real(8)                           :: Tbund_now
      real(8)                           :: crRainIC        !canopy-intercepted rainfall [mm]
      real(8)                           :: det
      real(8)                           :: Hbund_now
      !real(8)                           :: irrigationIntensity  !BK20240517
      real(8)                           :: irriheight  
      real(8)                           :: kDet            !rainfall detachment coefficient
      real(8)                           :: KEr             !kinetic energy of rainfall [J/m^2]
      real(8)                           :: KEl             !kinetic energy of leaf drainage [J/m^2]
      real(8)                           :: KEi             !kinetic energy of irrigation [J/m^2]
      real(8)                           :: rainIntensity   !mean rainfall intensity [mm/hour]
      real(8)                           :: rainSedDet      !sediment detachment by rainfall [kg s-1]
      real(8)                           :: rainSedDet_bund !sediment detachment at the bund by rainfall [kg s-1]
      real(8)                           :: irriSedDet      !sediment detachment by irrigation water [kg s-1]
      !real(8), dimension(4)             :: SSirri          !need to be checked
      !real(8), dimension(4)             :: topSoilFractions_forest !임시 !BK20240205 
      !real(8)                           :: crVegeCover_tmp  ! 토지피복별 crVegeCover (임시) 입력자료로 구성
   
      crRainIC     = 2.    !temporary hard-coded
      irriheight   = 1.    !hard-coded
      kDet         = 100. !hard-coded
      !SSirri       = 0.0   !temporary hard-coded (kyh)
      det          = 0.25   !Soil detachability (g/J) temporary hard-coded (kyh)
      !areaP_Ha     = 10000. !temporary hard-coded (kyh)
      Tbund_now = 1. !100. !temporary hard-coded (kyh) hlcho  !BK20240213 real number
      
      areaP_Ha = (dx * dx * areaP_Ratio) / 10000.  ![ha]  (/ 10000.: m2 -> ha)
      Abund_ha = (dx - Tbund_now) * Tbund_now * 4. / 10000. ! (/ 10000.: m2 -> ha) !BK20240213 real number
      areaIrri_Ha = areaP_Ha  * 0.0001 

      !1-1) splash detachment by rain drops (kg)
      !rainfall (mm)
      !if(rainRate >= 0.0) then
      !   rain = rainRate * dt ! rainRate mm/s
      !else
      !   rain = 0.0
      !endif    !BK20240421

      !mean rainfall intensity (mm/hour)     
      !mm/s --> mm/hr
      !if(rainRate >= 0.0) then
      !   rainIntensity = rainRate * 3600.
      if(Qrain >= 0.0) then   !BK20240421
         rainIntensity = Qrain * 3600.   !BK20240421
      else
         rainIntensity = 0.0
      endif

   
      !estimate kinetic energy of the rain (J/m^2) - from EUROSEM model
      !rain: rainfall during the time-step (mm)
      !lcover = int(vgtyp)    !BK20231006  !BK20240205
      !if(vgtyp >= 15 .AND. vgtyp <= 17) then ! 식생의 경우, 바닥에 낙엽 등이 있어 강우에너지 감소 필요    !BK20231006
      if(lcover >= 15 .AND. lcover <= 17) then ! 식생의 경우, 바닥에 낙엽 등이 있어 강우에너지 감소 필요    !BK20231006
         crCropHeight = 0.2
         !topSoilFractions_forest(1) = 0.30 !clay  !BK20240205
         !topSoilFractions_forest(2) = 0.35 !silt
         !topSoilFractions_forest(3) = 0.20 !fine sand
         !topSoilFractions_forest(4) = 0.15 !coarse sand
      endif
   
      ! kE 'J/m2/mm' -> rain / dt로 => 'J/m2/s'로 변경
      !if (rain > 0.0) then
      !   if (rain > crRainIC) then
      !      KEr = max(8.95 + 8.44 * log(rainIntensity), 0.0) * (rain - crRainIC) / dt
      !      KEl = max(max(15.8 * (crCropHeight**0.5) - 5.87, 0.0) * crRainIC, 0.0) / dt
      !   else
      !      KEr = 0.0
      !      KEl = max(max(15.8 * (crCropHeight**0.5) - 5.87, 0.0) * rain, 0.0) /dt	
      !   endif   
      !else
      !   KEr = 0.
      !   KEl = 0.
      !endif

      ! kinetic energy [J/m^2] of raindrops and leaf drips (Brandt, 1990)
      if (Qrain > 0.0) then    
            KEr = max(8.95 + 8.44 * log(rainIntensity), 0.0) * Qrain   !BK20240421
            if (Qrain*dt > crRainIC) then           
               KEl = max(max(15.8 * (crCropHeight**0.5) - 5.87, 0.0) * (Qrain*dt - crRainIC), 0.0) / dt   !BK20240421
            else           
               KEl = 0.0
            endif   
      else
         KEr = 0.
         KEl = 0.
      endif
      
      !if (vgtyp >=  1 .AND. vgtyp <= 6) then ! 도시지역, 이후 토지피복별 매개변수로 변경   !BK20231006
      if (lcover >=  1 .AND. lcover <= 6) then ! 도시지역, 이후 토지피복별 매개변수로 변경   !BK20231006
            det = 0.1  !HARDCODED
      endif
   
   ! 아래 식에서 OLWDepth => 이전 타임스텝에서의 overland flow depth를 나타냄
   ! 현재 Noah-MP에서 해당값이 없어 if 분기문을 만들고 Runsrf 구분 추가. 
   ! 단 Runsrf 단위가 m/s여서 불합리하며 overland flow depth로 이후 코드 변명 필요	 
      ! splash detachment (Morgan et al., 1998)
      If (Hbund_now > 0.0) then
         if (OLWDepth > 0.0) then
            kDet = 100.
            rainSedDet = 10. * det * (KEr * (1. - crVegeCover) &
                        + KEl * crVegeCover) &
                        !* exp(-kDet * OLWDepth) * (areaP_Ha - Abund_ha) * dt
                        * exp(-kDet * OLWDepth/1000.) * (areaP_Ha - Abund_ha) !BK20240421
         else
            kDet = 0.1 !20000.0
            !if (Runsrf > 0.0) then ! Runsrf 초기값 -9999
            if (Qsurf > 0.0) then    !BK20240421
               rainSedDet = 10. * det * (KEr * (1. - crVegeCover) &
                        + KEl * crVegeCover) &
                        !* exp(-kDet * Runsrf) * (areaP_Ha - Abund_ha) * dt
                        * exp(-kDet * Qsurf/1000.) * (areaP_Ha - Abund_ha) !BK20240421
            else
               rainSedDet = 10. * det * (KEr * (1. - crVegeCover) &
                        + KEl * crVegeCover) * (areaP_Ha - Abund_ha)      
            endif
         endif
      else  
         if (OLWDepth > 0.0) then
            kDet = 100.
            rainSedDet = 10. * det * (KEr * (1. - crVegeCover) &
                        + KEl * crVegeCover) &
                        !* exp(-kDet * OLWDepth) * areaP_Ha * dt
                        * exp(-kDet * OLWDepth/1000.) * areaP_Ha  !BK20240421
         else
            kDet = 0.1 !20000.0
            !if(Runsrf > 0.0) then
            if(Qsurf > 0.0) then   !BK20240421
               rainSedDet = 10. * det * (KEr * (1. - crVegeCover) &
                        + KEl * crVegeCover) &
                        !* exp(-kDet * Runsrf) * areaP_Ha * dt
                        * exp(-kDet * Qsurf/1000.) * areaP_Ha  !BK20240421
            else
               rainSedDet = 10. * det * (KEr * (1. - crVegeCover) &
                        + KEl * crVegeCover) * areaP_Ha 
            endif
         endif
      endif


   
      !1-2) bund splash detachment by rain drops (kg)
      if (Hbund_now > 0.0) then
         !rainSedDet_bund = 10. * det * KEr * Abund_ha
         rainSedDet_bund = 10. * det * KEr * Abund_ha  !BK20240517
      else
         rainSedDet_bund = 0.
      endif

      !2) splash detachment by irrigation water (kg)
      !irrigation(mm)
      !irri = Qirri / (areaIrri_Ha * 10000) * 1000  !need to be checked

      !mean irrigation intensity (mm/hr) : m3 -> mm/hr
      !irrigationIntensity = irri / dt * 3600  !need to be checked
      !irrigationIntensity = IRRIDepth / dt * 3600.   !BK20240517  OBSOLETE - not being used!

      !if (Qirri > 0.0) then
      if (IRRIDepth > 0.0) then
         !KEi = max(max(15.8 * (irriheight**0.5) - 5.87, 0.0) * irri, 0.0)  !need to be checked
         !KEi = max(max(15.8 * (irriheight**0.5) - 5.87, 0.0) * IRRIDepth, 0.0) 
         KEi = max(max(15.8 * (irriheight**0.5) - 5.87, 0.0) * IRRIDepth, 0.0) / dt  !BK20240517
      else
         KEi = 0.0
      endif

      !irriSedDet = 10. * det * KEi * exp(-kDet * OLWDepth) * areaIrri_Ha
      !irriSedDet = 10. * det * KEi * exp(-kDet * OLWDepth/1000.) * areaIrri_Ha  !BK20240421
      !--- BK20240517
      if (OLWDepth > 0.0) then
         irriSedDet = 10. * det * KEi * exp(-kDet * OLWDepth/1000.) * areaIrri_Ha 
      else
         irriSedDet = 10. * det * KEi * areaIrri_Ha 
      endif
      !--- BK20240517


      !assign the total amount of splash detachment to each sediment class (kg)
      !if(vgtyp >=  15 .AND. vgtyp <=17) then ! 식생의 경우, top soil의 토성을 다르게 주기 위해, 이후 토양 또는 토지이용 정보에 따른 매개변수로 입력 받도록 함   !BK20231006
      !if(lcover >=  15 .AND. lcover <=17) then ! 식생의 경우, top soil의 토성을 다르게 주기 위해, 이후 토양 또는 토지이용 정보에 따른 매개변수로 입력 받도록 함   !BK20231006
      !   overSed%Srdet(ips,i,j) = tmpSedDet * topSoilFractions_forest(ips) !* areaP_Ha ! [kg/cellarea]
      !else
      !   overSed%Srdet(ips,i,j) = tmpSedDet * topSoilFractions(ips) !* areaP_Ha ! [kg/cellarea]	 
      !endif

         !--- BK20240914
         !if (isnan(SedCNP_hydro%sltyp(i,j))) then
         if (SedCNP_hydro%sltyp(i,j) /= SedCNP_hydro%sltyp(i,j)) then  !BK20250720
         st = 14  !hardcoded
         else
         st = int(SedCNP_hydro%sltyp(i,j))               !BK20240205
         endif
         !--- BK20240914

      !--- BK20240915
      if (st .gt. 0) then 
         overSed%Srdet(ips,i,j) = (rainSedDet + rainSedDet_bund + irriSedDet) &
                                  * overSed%SOILPSF(st,ips)  !BK20240205
      else
         overSed%Srdet(ips,i,j) = 0.
      endif
      !--- BK20240915
      
      !update overland sediment storage (kg)
      overSed%Sol(ips,i,j) = overSed%Sol(ips,i,j) + overSed%Srdet(ips,i,j) * dt !+ SSirri(ips)  !BK20260219   

   end subroutine splashDetachment


   !subroutine paddySedTransport(i,j,ips,dt,dx,rainRate,OLWDepth,IRRIDepth,&
   !                             Runsrf,infl,evap,k_hs_dstc,k_setV,dp)
   subroutine paddySedTransport(i,j,ips,dt,dx,Qrain,OLWDepth,IRRIDepth,&   !BK20240421
                                 Qsurf,Qinfl,Qevap,k_hs_dstc,k_setV,dp)   !BK20240421

      implicit none

      integer,              intent(in) :: i, j
      integer,              intent(in) :: ips
      real(8),                 intent(in) :: dt          !time-step [sec]
      real(8),                 intent(in) :: dx          !grid cell size [m]
      !real(8),                 intent(in) :: rainRate    !rainfall rate [mm/s]  
      real(8),                 intent(in) :: Qrain       !rainfall rate [mm/s]   !BK20240421
      !real(8),                 intent(in) :: OLWDepth    !overland water depth [m]
      real(8),                 intent(in) :: OLWDepth    !overland water depth [mm]   !BK20240421
      real(8),                 intent(in) :: IRRIDepth   !irrigation water depth [mm]
      !real(8),                 intent(in) :: Runsrf      !surface runoff [mm/s]
      real(8),                 intent(in) :: Qsurf       !surface runoff [mm/s]   !BK20240421
      !real(8),                 intent(in) :: infl        !infiltration [mm/s]
      real(8),                 intent(in) :: Qinfl       !infiltration [mm/s]   !BK20240421
      !real(8),                 intent(in) :: evap        !evaporation [mm/s] --> need to be updated (kyh_debug)
      real(8),                 intent(in) :: Qevap       !evaporation [mm/s] --> need to be updated (kyh_debug)  !BK20240421
      real(8),   dimension(4), intent(in) :: k_hs_dstc
      real(8),   dimension(4), intent(in) :: k_setV
      real(8),   dimension(4)             :: setV        !values need to be defined or calculated
      real(8),   dimension(4), intent(in) :: dp        !particle diameter [μm]	 
      !real(8)                             :: rain        ![mm]  !BK20240421
      !real(8)                             :: irri        ![mm]
      !real(8)                             :: surf        ![mm]  !BK20240421
      real(8)                             :: k_sink
      real(8)                             :: Fstg
      !real(8)                             :: Dstc_RI
      !real(8)                             :: Qirri
      real(8)                             :: areaP       !pervious area in paddy fields [m2]
      !real(8)                             :: OLWStorageP !overland water storage in paddy fields [m3]
      real(8)                             :: OLWStorage !overland water storage [m3]   !BK20240421
      !real(8)                             :: QsurfP, Qinfl, QevapW   !BK20240421
      real(8)                             :: Dstc
      real(8)                             :: d         !dimensionless particle diameter [dimensionless]
      real(8)                             :: a, m      !experimentally determied constant
      real(8)                             :: df        !median floc diameter [um]
      !real(8)                             :: g         !gravitational acceleration [m/s2]
      !real(8)                             :: Gs        !specific gravity of soil
      !real(8)                             :: v         !kinematic viscosity of water [m2/s]
      real(8)                             :: Sow_local        !BK20260605 beginning-of-step Sow [kg]
      real(8)                             :: Soler_local      !BK20260605 erosion rate [kg/s]
      real(8)                             :: k_dep            !BK20260605 settling rate constant [s^-1]
      real(8)                             :: k_run            !BK20260605 runoff rate constant [s^-1]
      real(8)                             :: k_total_paddy    !BK20260605 [s^-1]
      real(8)                             :: Sow_ss_paddy     !BK20260605 CSTR steady-state Sow [kg]
      real(8)                             :: frac_rem_paddy   !BK20260605
      real(8)                             :: int_Sow_paddy    !BK20260605 time-integrated Sow [kg·s]
      real(8)                             :: dep_mass_paddy   !BK20260605 deposition mass [kg]
      real(8)                             :: run_mass_paddy   !BK20260605 runoff sediment mass [kg]
      real(8)                             :: Sow_new_paddy    !BK20260605 end-of-step Sow [kg]


      a  = 8.4E-3  !BK20260127 restored original value
      !a  = 8.4 * (10.**(-20.)) ! 침강속도 검토필요, 원 식으로 산정하면, size별로 8.40E-05 (=>8.40E-22), 6.82E-19, 6.72E-17, 6.56E-15 나옴
      m  = 0.024
      !g  = 9.8
      !Gs = 2.65
      !v  = 0.0001


      !areaP = dx * dx * 0.9  ![m2]   !BK20240421 arbitrary ratio applied
      areaP = dx * dx   ![m2]   !BK20240421 


      !if(rainRate >= 0.0) then
      !   rain = rainRate * dt  ![mm]
      !else
      !   rain = 0.0
      !endif   !BK20240421

      !irri = Qirri / areaP * 1000 ![mm]  --> need to be checked

      !surf = QsurfP / areaP * 1000 ![mm] --> need to be checked
      !if(Runsrf >= 0.0) then
      !   surf = Runsrf * dt ![mm]            --> need to be checked
      !else
      !   surf = 0.0
      !endif  !BK20240420

      !QsurfP = (surf/1000.) * areaP      ![m3] --> need to be checked (especially unit)  !BK20240421
      !OLWStorageP = OLWDepth * areaP    ![m3]                                            !BK20240421
      OLWStorage = OLWDepth/1000. * areaP    ![m3]                                        !BK20240421
      !Qinfl  = (infl*dt/1000.) * areaP   ![m3] --> need to be checked (kyh_debug)        !BK20240421 
      !QevapW = (evap*dt/1000.) * areaP   ![m3] --> need to be checked (kyh_debug)        !BK20240421

      !k_sink = exp(-(rain + irri + surf)) 
      !k_sink = exp(-(rain + IRRIDepth + surf))
      k_sink = exp(-(Qrain*dt + IRRIDepth + Qsurf*dt))  !BK20240421
      !Fstg = min(1.0, (OLWDepth * 100.)**2.)
      Fstg = min(1.0, (OLWDepth / 10.)**2.)   !BK20240421

      !BK20260605: replaced sequential Euler (erosion→deposition→runoff) with a CSTR analytical
      !solution. The prior sequential approach had two logical flaws:
      !(1) Soldp was computed from post-erosion Sow (inflated), then Ssurf from further-reduced
      !    Sow — deposition and runoff competed for the same Sow sequentially instead of
      !    simultaneously, causing Ssurf to be systematically underestimated vs Soler-Soldp.
      !(2) The runoff concentration used TotalOLWater (OLWStorage + Qsurf*dt + Qinfl*dt + Qevap*dt)
      !    as the denominator, incorrectly diluting the concentration by non-sediment-carrying
      !    outflows (Qinfl, Qevap) and including Qsurf in its own denominator.
      !The CSTR treats all three processes simultaneously from the same beginning-of-step Sow:
      !  d(Sow)/dt = Soler - k_total*Sow,  k_total = k_dep + k_run
      !  k_dep = k_sink * k_setV * vs / (OLWDepth/1000)   [s^-1]  (settling)
      !  k_run = Fstg * Qsurf / OLWDepth                  [s^-1]  (runoff)

      !--- Step 1: compute erosion mass ---
      Dstc = min(overSed%Sol(ips,i,j), overSed%Srdet(ips,i,j) * dt * Fstg * & !BK20260219
                (1. - k_hs_dstc(ips) / (k_hs_dstc(ips) + OLWDepth/1000.)) &
                + overSed%Sol(ips,i,j) * 0.001 * (1. - k_sink))
      Dstc = max(Dstc, 0.0)
      Soler_local = Dstc / dt                                                   !BK20260605

      !--- Step 2: compute settling velocity ---
      !setV(ips) = 898333.3 * dp(ips) ^ 2
      if(ips.eq.1) then  !clay
         df = dp(ips)
         setV(ips)  = a * (df**m)  ![cm/s]
         setV(ips)  = setV(ips)  / 100.    ![cm/s] -> [m/s]
      else                 !sand, silt
         d = dp(ips)*(10.**-6.) * (((Gs - 1.) * g / (v**2.))**(1./3.))
         ! Cheng (1997) formula for settling velocity  !BK20260511 typo corrected
         setV(ips)  = v / (dp(ips)*(10.**-6.)) * (((25. + 1.2 * (d**2.))**0.5 - 5.)**1.5)  ![m/s]
      endif

      !--- Step 3: CSTR rate constants ---                                      !BK20260605
      Sow_local = overSed%Sow(ips,i,j)                                         !BK20260605
      k_dep = 0.0d0                                                             !BK20260605
      k_run = 0.0d0                                                             !BK20260605
      if (OLWDepth > 0.0d0) then                                               !BK20260605
         k_dep = k_sink * k_setV(ips) * setV(ips) / (OLWDepth/1000.d0)        !BK20260605
         k_run = Fstg * (Qsurf/1000.d0) / (OLWDepth/1000.d0)                  !BK20260605
      endif                                                                      !BK20260605
      k_total_paddy = k_dep + k_run                                             !BK20260605

      !--- Step 4: CSTR analytical solution ---                                 !BK20260605
      if (k_total_paddy > 0.0d0) then                                          !BK20260605
         Sow_ss_paddy   = Soler_local / k_total_paddy                          !BK20260605
         frac_rem_paddy = exp(-k_total_paddy * dt)                             !BK20260605
         int_Sow_paddy  = Sow_ss_paddy * dt &                                  !BK20260605
                          - (Sow_ss_paddy - Sow_local) &                       !BK20260605
                          * (1.0d0 - frac_rem_paddy) / k_total_paddy           !BK20260605
         int_Sow_paddy  = max(0.0d0, int_Sow_paddy)                            !BK20260605
         dep_mass_paddy = k_dep * int_Sow_paddy                                !BK20260605
         run_mass_paddy = k_run * int_Sow_paddy                                !BK20260605
         Sow_new_paddy  = Sow_ss_paddy + (Sow_local - Sow_ss_paddy) &         !BK20260605
                          * frac_rem_paddy                                      !BK20260605
      else                                                                      !BK20260605
         dep_mass_paddy = 0.0d0                                                !BK20260605
         run_mass_paddy = 0.0d0                                                !BK20260605
         Sow_new_paddy  = Sow_local + Soler_local * dt                        !BK20260605
      endif                                                                     !BK20260605

      !--- Step 5: update state variables ---                                   !BK20260605
      overSed%Soler(ips,i,j) = Soler_local                                     !BK20260605
      overSed%Sol(ips,i,j)   = max(0.0d0, overSed%Sol(ips,i,j) &              !BK20260605
                               - Dstc + dep_mass_paddy)                        !BK20260605
      overSed%Sow(ips,i,j)   = max(0.0d0, Sow_new_paddy)                      !BK20260605
      overSed%Soldp(ips,i,j) = dep_mass_paddy / dt                             !BK20260605
      overSed%Ssurf(ips,i,j) = run_mass_paddy / dt                             !BK20260605
      if (OLWStorage > 1.0d-12) then                                           !BK20260605
         overSed%Colsd(ips,i,j) = overSed%Sow(ips,i,j) / OLWStorage          !BK20260605
      else                                                                      !BK20260605
         overSed%Colsd(ips,i,j) = 0.0d0                                       !BK20260605
      endif                                                                     !BK20260605

   end subroutine paddySedTransport


   !subroutine overlandSedTransportJS(i,j,ips,dt,dx,Runsrf,Sf,areaP_Ratio, &
   !                                  dp,vc,OLWDepth,manningsn, vgtyp)
   subroutine overlandSedTransportJS(i,j,ips,dt,dx,Qsurf,Sf,areaP_Ratio, &   !BK20240421
                                       dp,vc,OLWDepth,manningsn,lcover)  !BK20240205

      implicit none

      integer,            intent(in) :: i, j
      integer,            intent(in) :: ips
      real(8),               intent(in) :: dt       !time-step [sec]
      real(8),               intent(in) :: dx       !grid cell size [m]
      !real(8),               intent(in) :: Runsrf   !surface runoff [mm/s]
      real(8),               intent(in) :: Qsurf    !surface runoff [mm/s]   !BK20240421
      real(8),               intent(in) :: Sf       !overland slope
      real(8),               intent(in) :: areaP_Ratio       !pervious areas ratio (0~1)	 
      real(8), dimension(4), intent(in) :: dp       !particle diameter [μm]
      real(8), dimension(4), intent(in) :: vc       !critical discharge for each grid size; user-defined?
      real(8),               intent(in) :: OLWDepth !BK20260605 added intent(in)
      real(8),               intent(in) :: manningsn 
      !real(8),               intent(in) :: vgtyp   !BK20240205
      integer,               intent(in) :: lcover   !BK20240205
      real(8)                           :: Jc       !sediment transport capacity areal flux [kg/m2/s]
      real(8)                           :: Je       !erosion flux [kg/m2/s]
      real(8)                           :: Jd       !deposition flux [kg/m2/s]
      real(8)                           :: va       !surface runoff velocity
      !real(8)                           :: h_flow   !flow depth [m]  !BK20260511 !BK20260515
      real(8)                           :: Sow_local    !BK20260523
      real(8)                           :: vs           !BK20260523
      real(8)                           :: k_settling   !BK20260523
      real(8)                           :: k_discharge  !BK20260523
      real(8)                           :: k_total      !BK20260523
      real(8)                           :: frac_remaining, frac_lost  !BK20260523
      real(8)                           :: deposition_mass            !BK20260523
      real(8)                           :: discharge_rate             !BK20260523
      real(8)                           :: OLWater_cstr               !BK20260523
      real(8)                           :: dd_cstr, a_cstr, m_cstr, df_cstr  !BK20260523
      real(8)                           :: OLWDepth_m  !Manning-computed flow depth [m] !BK20260605


      !call overlandSedTransportCapa(i,j,ips,dt,dx,Runsrf,Sf,areaP_Ratio,vc,&
      !                              Jc,va,manningsn,vgtyp)
      call overlandSedTransportCapa(i,j,ips,dt,dx,Qsurf,Sf,areaP_Ratio,vc,&   !BK20240421
                                    Jc,va,OLWDepth_m,manningsn,lcover)  !BK20240205  !BK20260511 !BK20260515 !BK20260605

      !call overlandSedRunoffErosion(i,j,ips,dt,dx,Runsrf,&
      !                             Jc,va, OLWDepth,areaP_Ratio,vgtyp) ! OLWDepth 추가
      call overlandSedRunoffErosion(i,j,ips,dt,dx,Qsurf,&   !BK20240421
                                    Jc,va, OLWDepth_m,areaP_Ratio,lcover) !BK20240205  !BK20260511 !BK20260515 !BK20260605

      !BK20260523: replaced sequential overlandSedDeposition + overlandSedDischarge with a
      !CSTR analytical solution that computes both simultaneously from the same beginning-of-step
      !Sow. The prior sequential approach was numerically unstable when vs*dt/OLWDepth > 1:
      !for fine/coarse sand all Sow deposited before discharge was computed → Ssurf=0 always.
      !The CSTR partitions Sow loss between settling (→Sol) and outflow (→Ssurf) by rate fractions.
                                    
      !call overlandSedDeposition(i,j,ips,dt,dx,Qsurf,dp,OLWDepth,areaP_Ratio, va, lcover)  !BK20260523 replaced
      !call overlandSedDischarge(i,j,ips,dt,dx,Qsurf,OLWDepth,areaP_Ratio)                  !BK20260523 replaced

      Sow_local = overSed%Sow(ips,i,j)                                       !BK20260523

      !--- settling velocity (replicated from overlandSedDeposition) ---
      if (ips == 1) then                                                       !BK20260523
         a_cstr  = 8.4E-3                                                     !BK20260523
         m_cstr  = 0.024                                                      !BK20260523
         df_cstr = dp(ips)                                                    !BK20260523
         vs = a_cstr * (df_cstr**m_cstr) / 100.d0                            !BK20260523
      else                                                                     !BK20260523
         dd_cstr = dp(ips)*(10.d0**-6.) * (((Gs - 1.) * g / (v**2.))**(1./3.))  !BK20260523
         vs = v / (dp(ips)*(10.d0**-6.)) &                                    !BK20260523
              * (((25. + 1.2 * (dd_cstr**2.))**0.5 - 5.)**1.5)               !BK20260523
      endif                                                                    !BK20260523

      !--- rate constants [s^-1] ---
      OLWater_cstr = OLWDepth_m * dx * dx * areaP_Ratio                         !BK20260523 !BK20260605
      if (OLWDepth_m > 1.0e-9) then                                             !BK20260523 !BK20260605
         k_settling = vs / OLWDepth_m                                           !BK20260523 !BK20260605
      else                                                                     !BK20260523
         k_settling = 0.0d0                                                   !BK20260523
      endif                                                                    !BK20260523
      if (OLWater_cstr > 1.0e-12) then                                        !BK20260523
         k_discharge = (Qsurf / 1000.d0 * dx * dx) / OLWater_cstr            !BK20260523
      else                                                                     !BK20260523
         k_discharge = 0.0d0                                                  !BK20260523
      endif                                                                    !BK20260523
      k_total = k_settling + k_discharge                                       !BK20260523

      !--- analytical decay; partition loss by rate fraction ---
      if (k_total > 0.0d0) then                                                !BK20260523
         frac_remaining  = exp(-k_total * dt)                                 !BK20260523
         frac_lost       = 1.0d0 - frac_remaining                             !BK20260523
         deposition_mass = Sow_local * (k_settling  / k_total) * frac_lost   !BK20260523
         discharge_rate  = Sow_local * (k_discharge / k_total) * frac_lost / dt  !BK20260523
      else                                                                     !BK20260523
         frac_remaining  = 1.0d0                                              !BK20260523
         deposition_mass = 0.0d0                                              !BK20260523
         discharge_rate  = 0.0d0                                              !BK20260523
      endif                                                                    !BK20260523

      !--- update state variables (mass conserved: deposition+discharge = Sow_local*frac_lost) ---
      overSed%Sow(ips,i,j)   = max(0.0d0, Sow_local * frac_remaining)        !BK20260523
      overSed%Sol(ips,i,j)   = overSed%Sol(ips,i,j) + deposition_mass        !BK20260523
      overSed%Soldp(ips,i,j) = deposition_mass / dt                           !BK20260523
      overSed%Ssurf(ips,i,j) = discharge_rate                                 !BK20260523
      if (OLWater_cstr > 1.0e-12) then                                        !BK20260523
         overSed%Colsd(ips,i,j) = overSed%Sow(ips,i,j) / OLWater_cstr       !BK20260523
      else                                                                     !BK20260523
         overSed%Colsd(ips,i,j) = 0.0d0                                       !BK20260523
      endif                                                                    !BK20260523

   end subroutine overlandSedTransportJS


   subroutine overlandSedTransportCapa(i,j,ips,dt,dx,Qsurf,Sf,areaP_Ratio,vc, &  !in   !BK20240421
                                       Jc,va,OLWDepth,manningsn,lcover)    !out  !BK20240205  !BK20260511 !BK20260515

      implicit none

      integer,               intent(in)  :: i, j
      integer,               intent(in)  :: ips
      real(8),                  intent(in)  :: dt        !time-step [sec]
      real(8),                  intent(in)  :: dx        !grid cell size [m]
      !real(8),                  intent(in)  :: Runsrf    !surface runoff [mm/s]
      real(8),                  intent(in)  :: Qsurf     !surface runoff [mm/s]   !BK20240421
      real(8),                  intent(in)  :: Sf        !overland slope
      real(8),                  intent(in)  :: areaP_Ratio       !pervious areas ratio (0~1)	 
      real(8),    dimension(4), intent(in)  :: vc        !critical trasport velocity for each grid size
      real(8),                  intent(out) :: Jc        !sediment transport capacity areal flux [kg/m2/s]
      real(8),                  intent(out) :: va        !surface runoff velocity
      real(8),                  intent(out) :: OLWDepth    !flow depth [m]  !BK20260511 !BK20260515
      real(8),                  intent(in)  :: manningsn !manning's N [-] hlcho
      !real(8),                  intent(in)  :: vgtyp
      integer,               intent(in)  :: lcover        !land cover code   !BK20231006 !BK20240205
      !real(8)                               :: Runsrf_m  !surface runoff [m/s]  !BK20240420
      !real(8)                               :: Qsurf     !surface runoff [m3/s]  --> unit check    !BK20240421
      real(8)                               :: kt, beta_s, gamma_s    !calibration parameters; user-defined?
      real(8)                               :: Be        !width of eroding surface in flow direction [m]
      !real(8)                               :: Area_x    !cross section in the direction of flow [m2]
      real(8)                               :: q         !unit flow rate of water [m2/s]
      real(8)                               :: qs        !total sediment transport capacity [kg/m/s]
      real(8)                               :: hydraulic_r !hydraulic radius [m] hlcho
      real(8)                               :: hc        !overland flow depth at critical discharge(m)   !BK20260216
      real(8),    dimension(4)              :: qc        !critical discharge for each grid size [m2/s;
   
      !temporarily hard-coded
      ! hlcho kt, (*2000) erosion이 너무 일어나지 않음
      ! *** BK20260216: HARDWIRED - to be provided by the user ***
      ! The original value of 1.542E8 appears to be dimensionally inconsistent for SI units (kg/m/s for qs, m^2/s for q),
      ! likely stemming from a misapplication of a constant from an equation using different units (e.g., ppm and cfs/ft).
      ! This causes extreme over-prediction of transport capacity. A much smaller value is needed for physical realism.
      ! *** This value should be treated as a calibration parameter ***
      kt      = 1.542*(10.**3.) !Julien (1998, 2002) recommends a modified form of the Kilinc and Richardson (1973) relationship !BK20260511
      beta_s  = 1.4
      gamma_s = 1.4

      !width of eroding surface in flow direction [m]
      Be = dx

      !unit flow rate of water [m2/s] = meaning: volume flux (m3/s) per unit distance (m) in transverse direction
      !q = (Qsurf / 1000. * dx * dx) / Be    ! Qsurf in mm/s  !BK20240421
      q = (Qsurf / 1000. * dx * dx) / dx  !BK20260216

      ! Calculate actual flow depth and velocity from Manning's equation for a wide rectangular channel  !---BK20260511 !BK20260515
      if (Sf > 1.0e-6 .and. manningsn > 1.0e-6 .and. q > 0.0) then 
         ! Calculate normal flow depth from Manning's equation for a wide channel (R_h ~ depth)
         OLWDepth = (q * manningsn / (Sf**0.5))**(3./5.) !BK20260517 added missing line
         if (OLWDepth > 1.0e-6) then
            va = q / OLWDepth
         else
            va = 0.0
            OLWDepth = 0.0
         endif
      else
         OLWDepth = 0.0
         va = 0.0
      endif  !---BK20260511 !BK20260515

      !k = (manningsn*vc(ips)/(sf**0.5))**(3.0/2.0) ! k = 동수반경 = B*h/(B+2h)
      !qc(ips) = vc(ips)*(k*Be/(Be-2*k)) ! hc (수심) = k*Be/(Be-2*k)
      !BK20260216
      !hydraulic_r = (manningsn * vc(ips) / sf**0.5)**1.5
      !hc = hydraulic_r * dx / (dx - 2. * hydraulic_r)
      !qc(ips) = vc(ips) * hc
      if (Sf < 1.0e-8) then   !---BK20260511
         ! If slope is negligible, critical discharge is effectively infinite,
         ! so transport capacity is zero.
         qc(ips) = 1.0e36
      else
         ! The critical depth 'hc' can be calculated directly from Manning's equation
         ! using the critical velocity 'vc'. The previous formulation, which confused
         ! hydraulic radius and depth, was conceptually flawed and numerically unstable.
         hc = (manningsn * vc(ips) / Sf**0.5)**1.5
         qc(ips) = vc(ips) * hc
      endif   !---BK20260511

      !total sediment transport capacity [kg/m/s] hlcho else: qs = 0.0 
      if(q.gt.qc(ips)) then 
         qs = kt * ((q - qc(ips))**beta_s) * (Sf**gamma_s)
      else
         qs = 0.0
      endif

      !sediment transport capacity areal flux [kg/m2/s]
      !
      ! if (lcover >= 1 .AND. lcover <= 6) then       !BK20231006
      !    !*** BK20260216: WRONG equation - Be = 'slope length' in the flow direction!!! ***
      !    !Jc = qs / (Be * areaP_Ratio) ! 도시지역의 경우 유출 단면을 줄여줄 필요가 있음.  !BK20240421 도시지역에만 적용???   
      ! else
      !    Jc = qs / Be
      ! endif
      !BK20260216
      Jc = qs / Be   ! sediment transport capacity areal flux (kg/m2/s) = flux per unit surface area 

      !cross section in the direction of flow [m2]
      !Area_x = Be * OLWDepth !--> Flowdepth
   
      !BK20260216 *** WRONG equation: Be is 'length' along the flow direction!!! ***
      !Area_x = Be *  (Qsurf*dt/1000.)  !BK20240421 CHECK OUT THE EQUATION 
      !need to be checked 
      !hydraulic_r = Be * Runsrf_m / (Be + Runsrf_m) !hlcho
      !hydraulic_r = Be * (Qsurf*dt/1000.) / (Be + Qsurf*dt/1000.) !hlcho   !BK20240421 CHECK OUT THE EQUATION
      ! if(Area_x > 0.0) then
      !    !   va = Qsurf / Area_x   !--> Qsurf unit check is required
      !    va = 1.0 / manningsn * (hydraulic_r**(2.0/3.0)) * (Sf**0.5)
      ! else
      !    va = 0.0
      ! endif

      !sediment transport capacity [kg/s]
      !overSed%Solcp(ips,i,j) = Jc * Area_x * dt
      overSed%Solcp(ips,i,j) = Jc * dx * dx  !BK20260216  [kg s-1]

   end subroutine overlandSedTransportCapa


   subroutine overlandSedRunoffErosion(i,j,ips,dt,dx,Qsurf,&    !BK20240421
                                       Jc,va, OLWDepth, areaP_Ratio,lcover)  !BK20240205  !BK20260511 !BK20260515

      implicit none

      integer, intent(in)  :: i, j
      integer, intent(in)  :: ips
      real(8),    intent(in)  :: dt        !time-step [sec]
      real(8),    intent(in)  :: dx        !grid cell size [m]
      real(8),    intent(in)  :: Qsurf     !surface runoff [mm/s]   !BK20240421
      real(8),    intent(in)  :: Jc        !sediment transport capacity areal flux [kg/m2/s]
      real(8),    intent(in)  :: va        !surface runoff velocity
      real(8),    intent(in)  :: OLWDepth    !flow depth [m]  !BK20260511 !BK20260515
      real(8),    intent(in)  :: areaP_Ratio       !pervious areas ratio (0~1)
      integer, intent(in)  :: lcover        !land cover code   !BK20231006  !BK20240205
      real(8)                 :: Vmix      !bulk volume of the effective mixing depth of the soil [m3]
      real(8)                 :: vr        !resuspension (erosion) velocity [m/sec]
      real(8)                 :: Csb       !concentration of sediment at the bottom boundary
                                             !(in the bed) [kg/m3]
      real(8)                 :: Csw       !concentration of sediment for water 
      real(8)                 :: Je        !erosion flux [kg/m2/s]
      !real(8)                 :: Be        !width of eroding surface in flow direction [m]
      real(8)                 :: erosion_mass !sediment erosion mass [kg]  !BK20260511
      !real(8)                 :: Area_x    !cross section in the direction of flow [m2]
      !real(8)                 :: Runsrf_m  !surface runoff [m/s]  !BK20240420
      !real(8)                 :: kero      !erosion cal. coefficient 
      real(8)                 :: WaterVolume ! [m3]  !BK20260511


      !--- concentration of sediment at the bottom boundary (in the bed) [kg/m3]
      Vmix = dx * dx * Dmix * areaP_Ratio  
      if (Vmix > 0.) then  !BK20251031
         Csb = overSed%Sol(ips,i,j) / Vmix !hlcho  
      else
         Csb = 0.
      endif

      !--- concentration of sediment in the overland water and/or runoff flow [kg/m3]
      ! if (areaP_Ratio > 0.0) then   !BK20260218
      !    if ((OLWDepth + Qsurf*dt)/1000. > 0) then   !BK20240420   !BK20240421 CHECK OUT THE EQUATION
      !       Csw = overSed%Sow(ips,i,j) / (dx * dx * areaP_Ratio * (OLWDepth + Qsurf*dt)/1000.)  !BK20240420 !BK20240421 CHECK OUT THE EQUATION
      !    else
      !       Csw = 0.0
      !    endif
      WaterVolume = OLWDepth * dx * dx * areaP_Ratio  !---BK20260511 !BK20260515
      if (WaterVolume > 1.0e-9) then
         Csw = overSed%Sow(ips,i,j) / WaterVolume   !---BK20260511 !BK20260515
      else
         Csw = 0.0
      endif

      !--- resuspension (erosion) velocity [kg/m2/s]
      ! if (lcover >= 1 .AND. lcover <= 6) then   !BK20231006
      !    vr = (Jc - va / 3.0 * Csw) / 1300.0 ! 도시지역의 경우, 유속을 투수지표면으로 한정해야 함. 투수지표면 유속을 전체 유속의 1/3 가정, 
      ! else
      !    vr = (Jc - va * Csw) / 1300.0 ! (hlcho 검토 요망) 단위가 이상함 soil bulk density 1.3 g/cm³ * 1000 = 1300 kg/m³	 
      ! endif

      !BK20260218
      if (Jc > va * Csw) then
         vr = (Jc - va * Csw) / (Gs * 1000.0 * (1.0 - SedCNP_hydro%smcmax(i,j,1)))  
         vr = max(vr, 0.0)
      else
         vr = 0.0
      endif

      !erosion flux Je [kg/m2/s]
      Je = vr * Csb  !BK20260511

      
      !    ! arbitrary values hard-coded  !BK20260218 commented out
      !    if (ips == 1) then
      !       kero = 0.0009 / 2.0
      !    else if (ips == 2) then
      !       kero = 0.0008 / 3.2
      !    else if (ips == 3) then
      !       kero = 0.0007 / 4.0
      !    else if (ips == 4) then
      !       kero = 0.00045 / 4.0
      !    endif
      !   !    Je = vr * Csb * kero  !arbitrary

         !width of eroding surface in flow direction [m]
         !Be = dx

         !cross section in the direction of flow [m2]
         !Area_x = dx * dx * areaP_Ratio ! erosion은 바닥에서 발생해야 함 (hlcho 검토 요망) 
      
         !sediment erosion [kg/dt]
         !overSed%Soler(ips,i,j) = Je * Area_x * dt  
         !overSed%Soler(ips,i,j) = Je * dx * dx * dt  !BK20260218
         !overSed%Soler(ips,i,j) = min(overSed%Sol(ips,i,j), overSed%Soler(ips,i,j) ) 
         erosion_mass = Je * dx * dx * areaP_Ratio * dt   !BK20260511
         erosion_mass = min(overSed%Sol(ips,i,j), erosion_mass)  !BK20260511

         !sediment storage update
         !overSed%Sow(ips,i,j) = overSed%Sow(ips,i,j) + overSed%Soler(ips,i,j)  
         !overSed%Sol(ips,i,j) = overSed%Sol(ips,i,j) - overSed%Soler(ips,i,j) 
         overSed%Sow(ips,i,j) = overSed%Sow(ips,i,j) + erosion_mass  !BK20260511
         overSed%Sol(ips,i,j) = overSed%Sol(ips,i,j) - erosion_mass  !BK20260511

         if (dt > 0.) then  !---BK20260511
            overSed%Soler(ips,i,j) = erosion_mass / dt
         else
            overSed%Soler(ips,i,j) = 0.
         endif  !---BK20260511
   end subroutine overlandSedRunoffErosion


   !subroutine overlandSedDeposition(i,j,ips,dt,dx,Runsrf,dp,areaP_Ratio, va, vgtyp)
   subroutine overlandSedDeposition(i,j,ips,dt,dx,Qsurf,dp,OLWDepth,areaP_Ratio,va,lcover)  !BK20240205 !BK20240421 !BK20260219 !BK20260511 !BK20260515

      implicit none

      integer,            intent(in) :: i, j
      integer,            intent(in) :: ips
      real(8),               intent(in) :: dt        !time-step [sec]
      real(8),               intent(in) :: dx        !grid cell size [m]
      !real(8),               intent(in) :: Runsrf    !surface runoff [mm/s]
      real(8),               intent(in) :: Qsurf     !surface runoff [mm/s]   !BK20240421
      real(8), dimension(4), intent(in) :: dp        !particle diameter [μm]
      real(8),               intent(in) :: OLWDepth    !flow depth [m] !BK20260511 !BK20260515
      real(8),               intent(in) :: areaP_Ratio       !pervious areas ratio (0~1)
      real(8),               intent(in) :: va        !surface runoff velocity
      !real(8),               intent(in) :: vgtyp     !landuse 
      integer,            intent(in) :: lcover        !land cover code       !BK20231006  !BK20240205
      real(8)                           :: OLWater   !overland water storage [m3]
      !real(8)                           :: Qsurf     !surface runoff [m3/s]   !BK20240421
      real(8)                           :: vs
      real(8)                           :: Area_x    !cross section in the direction of flow [m2]
      real(8)                           :: a, m      !experimentally determied constant
      real(8)                           :: Be        !width of eroding surface in flow direction [m]
      real(8)                           :: dd        !dimensionless particle diameter [dimensionless]
      real(8)                           :: df        !median floc diameter [um]
      !real(8)                           :: g         !gravitational acceleration [m/s2]
      !real(8)                           :: Gs        !specific gravity of soil
      !real(8)                           :: v         !kinematic viscosity of water [m2/s]
      real(8)                           :: Jd        !deposition flux [kg/m2/s]  !BK20260511
      !real(8)                           :: Runsrf_m  !surface runoff [m/s]  !BK20240420
      real(8)                           :: Nf        !trap efficiency coefficient
      real(8)                           :: Nf_exp    !trap efficiency exponential coefficient
      real(8)                           :: Esedtrap   !trap efficiency
      real(8)                           :: deposition_mass !sediment deposition mass [kg]  !BK20260511


      OLWater = OLWDepth * dx * dx * areaP_Ratio ![m3] !BK20260511 !BK20260515

      !compute settling velocity
      if(ips == 1) then  
         !clay  ! Burban et al. (1990)
         a  = 8.4E-3  !BK20260127 restored original value
         m  = 0.024
         df = dp(ips)
         vs = a * (df**m)  ![cm/s]
         vs = vs / 100.    ![cm/s] -> [m/s] 
      else                
         !sand, silt  ! Cheng (1997)
         dd = dp(ips)*(10.**-6.) * (((Gs - 1.) * g / (v**2.))**(1./3.))
         vs = v / (dp(ips)*(10.**-6.)) * (((25. + 1.2 * (dd**2.))**0.5 - 5.)**1.5)  ![m/s]  !BK20260511
      endif
      
      !deposition flux  !---BK20260511 commented out  !BK20260515 restored
      if(OLWater > 1.0e-12) then  !BK20251204 prevent division by zero
        overSed%Colsd(ips,i,j) = overSed%Sow(ips,i,j) / OLWater
      else
        overSed%Colsd(ips,i,j) = 0.0
      endif
      Jd = vs * overSed%Colsd(ips,i,j)  !---BK20260511 commented out  !BK20260515 restored


      ! !========== SEDIMENT TRAP EFFICIENCY ========== !BK20260219
      ! ! THIS CODE IS TEMPORARILY COMMENTED OUT FOR TESTING PURPOSE.
      ! ! SHOULD BE REVIED AND IMPLEMENTED LATER.

      ! ! sediment trap, 도시지역의 경우 (낮은 매닝계수) 높은 유속으로 Esedtrap 값 매우 클 수 있음
      ! ! 따라서 도시지역은 불투수지표면이 아닌 투수지역 매닝계수를 부여해야 함
      ! Nf_exp = 0.87 ! 값 증가시 deposition 증가, 매개변수화 필요
      ! !if (va > 0 .AND. Runsrf_m > 0) then
      ! !        Nf = vs * (((dx*areaP_Ratio) / va) ** Nf_exp) / (va * Runsrf_m)
      ! if (va > 0 .AND. Qsurf/1000. > 0.) then  !BK20240420   !BK20240421
      !          Nf = vs * ((dx*areaP_Ratio/va)**Nf_exp) / (va * Qsurf/1000.)  !BK20240420   !BK20240421 CHECK OUT THE EQUATION
      ! else
      !          Nf = 0
      ! end if
   
      ! Esedtrap = 44.1 * (Nf**0.29) / 100 ! Tollner et al.(1976)
      ! Esedtrap = min(Esedtrap, 0.95)
   
      ! ! 도시지역 화단효과 반영
      ! !lcover = int(vgtyp)      !BK20231006
      ! !if (lcover >= 1 .AND. lcover <= 6) then   !BK20231006
      ! !   Esedtrap = min(Esedtrap*1.1, 0.95)
      ! !endif
   
      ! !if(Runsrf_m > 0.0) then        
      ! !   Jd = overSed%Sow(ips,i,j) * (1. - exp(-1. * vs * dt / Runsrf_m))
      ! !else
      ! !   Jd = 0.
      ! !endif
      Esedtrap = 0.0   !temporary setting
      ! !========== SEDIMENT TRAP EFFICIENCY ==========

      !width of eroding surface in flow direction [m]
      Be = dx
      !cross section in the direction of flow [m2]
      Area_x = Be * OLWDepth
      ! Area_x = Be * (Runsrf_m * dt)  !--> need to be checked hlcho
      ! Area_x = Be * Runsrf_m  !--> need to be checked hlcho hlcho (검토 요망)
      ! Area_x = dx * dx * areaP_Ratio !hlcho (검토 요망) deposition이 발생하는 flux 면

      !sediment deposition
      !overSed%Soldp(ips,i,j) = Jd * Area_x * dt + Esedtrap * overSed%Sow(ips,i,j) !BK20260218
      !overSed%Soldp(ips,i,j) = Jd * dx * dx * areaP_Ratio * dt + Esedtrap * overSed%Sow(ips,i,j) !BK20260218
      !overSed%Soldp(ips,i,j) = min(overSed%Sow(ips,i,j), overSed%Soldp(ips,i,j))/dt !BK20260517
      deposition_mass = Jd * dx * dx * areaP_Ratio * dt + Esedtrap * overSed%Sow(ips,i,j)  !BK20260511
      deposition_mass = min(overSed%Sow(ips,i,j), deposition_mass)   !BK20260511
      
      !---BK20260511  !BK20260515 commented out
      ! The original explicit Euler scheme (deposition_mass = Jd * Area * dt) is unstable
      ! when vs*dt/OLWDepth > 1. This can cause all suspended sediment to be deposited,
      ! zeroing out Sow before runoff concentration (Colsd) is calculated.
      ! Using the analytical solution of the first-order decay for deposition is more stable.
      ! if (OLWDepth > 1.0e-9) then
      !    deposition_mass = overSed%Sow(ips,i,j) * (1.0 - exp(-vs * dt / OLWDepth))
      ! else
      !    deposition_mass = 0.0
      ! endif
      ! deposition_mass = deposition_mass + Esedtrap * overSed%Sow(ips,i,j)
      ! ! The min() check is implicitly handled by the analytical solution, but kept for safety with Esedtrap
      ! deposition_mass = min(overSed%Sow(ips,i,j), deposition_mass)   
      !---BK20260511  !BK20260515 commented out
      
      ! Calculate deposition rate for output
      if (dt > 0.0) then
         overSed%Soldp(ips,i,j) = deposition_mass / dt
      else
         overSed%Soldp(ips,i,j) = 0.0
      endif
      
      !sediment storage update
      overSed%Sol(ips,i,j) = overSed%Sol(ips,i,j) + deposition_mass  !BK20260511 
      overSed%Sow(ips,i,j) = max(0.0, overSed%Sow(ips,i,j) - deposition_mass)  !BK20260511 

   end subroutine overlandSedDeposition


   !subroutine overlandSedDischarge(i,j,ips,dt,dx,Runsrf)
   !subroutine overlandSedDischarge(i,j,ips,dt,dx,Qsurf,OLWDepth)   !BK20240421 !BK20260219 !BK20260511 !BK20260515
   subroutine overlandSedDischarge(i,j,ips,dt,dx,Qsurf,OLWDepth,areaP_Ratio)   !BK20240421 !BK20260219 !BK20260511 !BK20260515

      implicit none

      integer, intent(in) :: i, j
      integer, intent(in) :: ips
      real(8),    intent(in) :: dt        !time-step [sec]
      real(8),    intent(in) :: dx        !grid cell size [m]
      real(8),    intent(in) :: Qsurf     !surface runoff [mm/s]   !BK20240421
      real(8),    intent(in) :: OLWDepth    !flow depth [m] !BK20260511 !BK20260515
      real(8),    intent(in) :: areaP_Ratio  !BK20260517
      !real(8)                :: Qsurf     !surface runoff [m3/s]   !BK20240421
      !real(8)                :: Runsrf_m  !surface runoff [m/s]  !BK20240421
      real(8)                :: Cs        !sediment concentration in water
      real(8)                :: OLWater   !overland water storage [m3]

      !!mm/s --> m/s
      !if(Runsrf >= 0.0) then
      !   Runsrf_m = Runsrf / 1000.0 !hlcho
      !!   Runsrf_m = Runsrf		
      !else
      !   Runsrf_m = 0.0
      !endif
      !Qsurf = max(Qsurf, 0.)  !BK20240420   !BK20240421

      !m/s --> m3/s
      !Qsurf = Runsrf_m * (dx * dx)
      !Qsurf = (Runsrf / 1000.) * (dx * dx)  !BK20240420   !BK20240421
      OLWater = OLWDepth * dx * dx * areaP_Ratio ![m3] !BK20260511 !BK20260515 !BK20260517

         if(OLWater > 1.0e-12) then  !BK20251204 prevent division by zero
            overSed%Colsd(ips,i,j) = overSed%Sow(ips,i,j) / OLWater 
         else
            overSed%Colsd(ips,i,j) = 0
         endif
      overSed%Ssurf(ips,i,j) = (Qsurf/1000.*dx*dx) * overSed%Colsd(ips,i,j) !BK20240421 !BK20260219

      !sediment storage update
      overSed%Ssurf(ips,i,j) = min(overSed%Sow(ips,i,j), overSed%Ssurf(ips,i,j)*dt) / dt !BK20260219
      overSed%Sow(ips,i,j) = max(overSed%Sow(ips,i,j) - overSed%Ssurf(ips,i,j)*dt, 0.0) ! hlcho min 추가 !BK20260219

   end subroutine overlandSedDischarge


   !subroutine sedimentAggregation(i,j,ips,dx,dt,OLWDepth, Runsrf)
   subroutine sedimentAggregation(i,j,ips,dx,dt,OLWDepth,Qsurf)  !BK20240421

      implicit none

      integer, intent(in) :: i, j
      integer, intent(in) :: ips
      real(8),    intent(in) :: dx          
      real(8),    intent(in) :: dt          !time-step [sec]
      real(8),    intent(in) :: OLWDepth    !overland water depth [mm]
      !real(8),    intent(in) :: Runsrf      !surface runoff [mm/s] ! HLCHO 추가
      real(8),    intent(in) :: Qsurf      !surface runoff [mm/s]  !BK20240421
      real(8)                :: E
      real(8)                :: Fstg
      real(8)                :: sedAggr_OL

      E = 2.7182818284590451
      !Fstg = exp(-OLWDepth * 1000) 
      !Fstg = exp(-OLWDepth * 1000)  ! [mm]  !conceptual estimation
   
      !if (Runsrf < 0.0) then  ! 초기값 -9999, 유속을 이용한 수심으로 변경해야 함
      if (Qsurf < 0.0) then  ! 초기값 -9999, 유속을 이용한 수심으로 변경해야 함   !BK20240421
         Fstg = exp(-OLWDepth)  ! [mm]  !conceptual estimation, assuming exponential decay with depth in mm  !BK20260511
      else
         !Fstg = exp(-(OLWDepth * 1000. + Runsrf * (dx *1000.) / 10.))  ! [mm]  유속 10 mm/sec으로 일단 가정하여 수식 적용	 
         Fstg = exp(-(OLWDepth + (Qsurf*dt) / 10.))  ! [mm]  유속 10 mm/sec으로 일단 가정하여 수식 적용 !BK20240421 CHECK OUT THE EQUATION !BK20260511 unit conversion	 
      endif 
      
      sedAggr_OL = overSed%Sol(ips,i,j) * Fstg / E * (dt / 24. / 3600.)


      !update overland sediment storage (kg)
      overSed%Sol(ips,i,j) = overSed%Sol(ips,i,j) - sedAggr_OL  
   
   end subroutine sedimentAggregation


end module module_overlandSed
