module module_channelSed

   use module_SedCNPvariables
!=====||__WHQ5403q__||=====!
!
   use SedCNP_config,          only: SedCNPmodel
!
!=====||__WHQ5403q__||=====!

   real(8)               :: Gs = 2.65  !specific gravity of sediment 
   real(8)               :: g = 9.8    !gravitational acceleration [m/s2]
   real(8)               :: v = 0.0001        !kinematic viscosity of water [m2/s]

   contains

   !subroutine channelSedTransport(ich)      
   subroutine channelSedTransport(ich,itime)      !BK20231029 !BK20240519(ich>chid) !BK20241101(chid>ich)

      implicit none

      integer, intent(in)             :: ich       !channel index
      integer, intent(in)             :: itime                   !BK20231016
      integer                         :: chid, gwid 
      integer                         :: ips
      real(8)                            :: Qdsch      !channel discharge (should be read from Noah-MP outpus) [m3/s]
      real(8),   dimension(4)            :: dp         !particle diameter [μm]
      real(8)                            :: va         !advective (flow) velocity (in the down-gradient direction) [m/s2]
      real(8)                            :: Sf         !channel (friction) slope [-]
      real(8)                            :: Lst        !channel length [m]
      real(8)                            :: dt         !time-step [sec]
      real(8)                            :: h          !channel flow depth [m]
      real(8)                            :: h_thres    !channel depth below the compound channel [m]
      real(8)                            :: BtmWdth    !bottom width of channel [m]
      real(8)                            :: SideSlopch !Channel side slope
      !real(8)                            :: TopWdth    !Top width of channel [m]          !BK20251117
      !real(8)                            :: TopWdthCC  !Compound Channel Top Width [m]    !BK20251117
      real(8)                            :: WtopWdth   !channel water top width [m]
      real(8)                            :: Axsect     !channel water flow cross section area [m2]
      !real(8)                            :: Axsect_1   !cross section area below compound channel [m2] !BK20251117
      !real(8)                            :: Axsect_2   !cross section area of compund channel [2]      !BK20251117
      real(8)                            :: WStorage   !water storage [m3]  !BK20240710

      chid = SedCNP_hydro%linkID(ich)         !BK20231029


      !particle diameter
      !hard-coded; should be relocated or moved to config file
      dp(1) = 1.      !clay [μm]
      dp(2) = 10.     !silt [μm]
      dp(3) = 100.    !fine sand [μm]
      dp(4) = 1000.   !coarse sand  [μm]

      !get channel geometric inputs
      BtmWdth    = channel_input%Wst(ich)            !bottom width of channel [m] 
      SideSlopch = channel_input%SideSlopch(ich)    !channel side slope   !BK20231029 (not working with the new nc file) !BK20251117
      !SideSlopch = 1.0          !BK20231029  !BK20241028 hardcoded for debugging
      
      ! compound channel NOT implemented for WHQ yet   !BK20251117
      !if (SedCNPmodel%SedCNP_compound_channel) then
      !   TopWdth    = channel_input%TopWdth(ich)     ! Top width of channel [m]
      !   TopWdthCC  = channel_input%TopWdthCC(ich)   ! Compound Channel Top Width [m]
      !endif

      Lst = channel_input%Lst(ich)                   !stream length [m]
      !Sf  = channel_input%Sf(ich)                    !channel slope
      Sf  = min(channel_input%Sf(ich), 0.01)          !channel slope  !BK20240620 for debugging

      !get model time-step
      dt = real(SedCNPmodel%SedCNP_timestep)

      !get hydrology inputs for a channel
      va         = SedCNP_hydro%chVa(ich)          !channel velocity [m/s]
                                                   !kyh_debug; need to be checked if chVa can be
                                                   !used for va (advective flow velocity [m/s2])
      Qdsch      = SedCNP_hydro%q_dsch(ich)        !channel discharge [m3/s]
      !Qusch      = SedCNP_hydro%q_usch0(ich)      !channel inflow from upstream [m3/s]
      !channelSed%Susch(:,ich) = channelSed%Susch0(:,ich)  !sediment transport from upstream [kg/dt]  !BK20240618 commented out
      h          = SedCNP_hydro%head(ich)          !channel flow depth [m]

      if (va < 0.0) then
         va = 0.0
      endif
      if (BtmWdth < 0.0) then
         BtmWdth = 0.0
      endif
      if (Lst < 0.0) then
         Lst = 0.0
      endif
      if (Sf < 0.0) then
         Sf = 0.0
      endif
      if (Qdsch < 0.0) then
         Qdsch = 0.0
      endif
      !if (Qusch < 0.0) then
      !   Qusch = 0.0
      !endif
      if (h < 0.0) then
         h = 0.0
      endif

      !Compute cross-section area of channel flow [m2]
      !WtopWdth = (2 * h * SideSlopch) + BtmWdth              !channel water top width [m]
      WtopWdth = (2 * h / SideSlopch) + BtmWdth              !channel water top width [m]  !BK20251117 corrected equation  

      ! compound channel NOT implemented for WHQ yet   !BK20251117
      !   if (SedCNPmodel%SedCNP_compound_channel) then
      !      if (abs(SideSlopch) > 1.0E-9) then   !BK20250810
      !         h_thres = (TopWdth - BtmWdth) / (2.0 * SideSlopch)  !channel depth below compound channel [m]
      !      else
      !         h_thres = 0.0
      !      endif
      !      if(h.le.h_thres) then
      !         Axsect = ((BtmWdth + WtopWdth) * h) / 2.0        !channel water flow cross section area
      !                                                            !isosceles trapezoid is assumed
      !      else
      !         Axsect_1 = ((BtmWdth + TopWdth) * h) / 2.0       !cross section area below the compound channel
      !         Axsect_2 = TopWdthCC * (h - h_thres)           !cross section area of compund channel
      !                                                            !rectangle is assumed
      !         Axsect = Axsect_1 + Axsect_2
      !      endif
      !   else
      !     Axsect = (BtmWdth + WtopWdth) * h / 2.0  !channel flow x-section area (isosceles trapezoid is assumed)
      !   endif

      Axsect = (BtmWdth + WtopWdth) * h / 2.0  !channel flow x-section area (isosceles trapezoid is assumed) !BK20251117


         !***** variables below are read from WHQ hydro sim. results ****************************** !BK20240713 
         !
         WStorage = SedCNP_hydro%Wstg_st(ich)
         !
         !Qabs (assumption: Qabs can't be greater than 50% of WStorage)
         !SedCNP_hydro%Qabs(ich,itime) = min(WStorage*0.5, SedCNP_hydro%Qabs_max(ich,itime)*dt) / dt
         !SedCNP_hydro%Qabs(ich,itime) = min(WStorage*0.5, SedCNP_hydro%Qabs(ich,itime)*dt) / dt  !BK20241016
         !************************************************************************************************************


         !call sedInput(ich,ips)     
         !call sedInput(ich,ips,itime)     !BK20231016 BK20231029 
         
         !--- BK20240726
         ! gwid = domain%gwid_ch(chid)  !BK20250928
         ! if (gwid .ge. 0) then  !'genuine' channel to which gwbasin drains !BK20250928
            do ips = 1, nps
               !call sedInput(ich,ips,itime,WStorage)
               call sedInput(ich,ips,itime,Qdsch,WStorage)  !BK20241105
            enddo     
         ! endif  !BK20250928

         do ips = 1, nps
            !call channelSedTransportEH(ich,ips,dt,Qdsch,Qusch,h,dp,va,Sf,BtmWdth,Lst,Axsect)  !Engelund and Hansen (1967)
            call channelSedTransportEH(ich,ips,dt,Qdsch,h,dp,va,Sf,BtmWdth,Lst,Axsect)  !Engelund and Hansen (1967)  !BK20231029
         enddo
         !--- BK20240726

   end subroutine channelSedTransport

   
   subroutine sedInput(ich,ips,itime,Qdsch,WStorage)    !BK20231016 Bk20231029 !BK20240519(ich>chid) !BK20240710 !BK20241101(chid>ich) !BK20241105

      implicit none

      integer, intent(in) :: ich  
      integer, intent(in) :: ips
      integer, intent(in) :: itime           !BK20231016
      real(8), intent(in) :: Qdsch         !BK20241105     
      real(8), intent(in) :: WStorage         !BK20240710             
      integer             :: i, j
      real(8)             :: dt              !BK20231016
      !integer             :: chid_dis, tlag  !BK20240710
      integer             :: ich_dis, lag_ts  !BK20241101
      integer             :: chid, gwid            !BK20240726
      real(8)             :: Cchsd_local
      real(8)             :: abs_ratio_local  !BK20260523

      dt = real(SedCNPmodel%SedCNP_timestep)   !BK20231016
      chid = SedCNP_hydro%linkID(ich)         !BK20231029
      gwid = domain%gwid_ch(chid)

      !initialisation
      channelSed%Solch0(ips,ich) = 0.  !BK20250925

      ! --- influx: overland surface runoff  !BK20240726
      do i = 1, domain%ix
         do j = 1, domain%jx
               if (gwid > 0) then  !BK20250928
                  if (SedCNP_hydro%gwbasin(i,j) .eq. gwid) then    !BK20231012  !BK20240610 !BK20240726
                     if (overSed%Ssurf0(ips,i,j) >= 0.0) then  
                        !channelSed%Solch(ips,ich) = channelSed%Solch(ips,ich) + overSed%Ssurf(ips,i,j)
                        channelSed%Solch0(ips,ich) = channelSed%Solch0(ips,ich) + overSed%Ssurf0(ips,i,j)  !BK20240618  !BK20250925
                     endif
                  endif
               else  !BK20250928
                  cycle  !BK20250928
               endif  !BK20250928
         enddo
      enddo

      !BK20240618 commented out
      !if (channelSed%Susch(ips,ich) < 0.0) then
      !   channelSed%Susch(ips,ich) = 0.0
      !endif
      !if (channelSed%Solch(ips,ich) < 0.0) then
      !   channelSed%Solch(ips,ich) = 0.0
      !endif
      !!if (channelSed%Spnt(ips,ich) < 0.0) then
      !!   channelSed%Spnt(ips,ich) = 0.0
      !if (channelSed%Spnt(ips,ich,itime) < 0.0) then   !BK20231016
      !   channelSed%Spnt(ips,ich,itime) = 0.0
      !endif
      !if (channelSed%Scw(ips,ich) < 0.0) then
      !   channelSed%Scw(ips,ich) = 0.0
      !endif

      !---BK20260523 commented out (old block)
      !    ! --- influx: point sources
      !    channelSed%Sxpnt(ips,ich) = channelSed%Spnt(ips,ich,itime) * dt  !BK20240701

      !    ! --- outflux: water distribution - abstraction !BK20240710
      !    if(SedCNP_hydro%Qabs(ich,itime) .gt. 0.) then
      !       ! channelSed%Sabs(ips,ich,itime) = channelSed%Scw(ips,ich) &
      !       !                                 * (SedCNP_hydro%Qabs(ich,itime)*dt / WStorage)
      !       !channelSed%Sxabs(ips,ich) = channelSed%Sabs(ips,ich,itime)*dt  !BK20240729    
      !       !-- BK20250705 modification for numerical stability
      !       !if ((WStorage + Qdsch*dt) > 1.0E-9) then
      !       !   channelSed%Sabs(ips,ich,itime) = channelSed%Scw(ips,ich) &
      !       !         * (SedCNP_hydro%Qabs(ich,itime) / (WStorage + Qdsch*dt))  !BK20241105
      !       !   channelSed%Sxabs(ips,ich) = channelSed%Sabs(ips,ich,itime)*dt  !BK20240729

      !       ! The calculation for sediment abstraction is harmonized with the sediment discharge
      !       ! calculation (S_out = C * Q_out) for consistency. The concentration is based on
      !       ! the water and sediment storage at the beginning of the timestep.
      !       if (WStorage > 1.0E-9) then
      !          Cchsd_local = channelSed%Scw(ips,ich) / WStorage  !BK20260521
      !       else
      !          !channelSed%Sabs(ips,ich,itime) = 0.0
      !          !channelSed%Sxabs(ips,ich) = 0.0
      !          Cchsd_local = 0.0  !BK20260521
      !       endif
      !       channelSed%Sabs(ips,ich,itime) = Cchsd_local * SedCNP_hydro%Qabs(ich,itime) !BK20260521
      !       channelSed%Sxabs(ips,ich) = channelSed%Sabs(ips,ich,itime) * dt !BK20260521
      !    endif
   
      !    !update storage
      !    !channelSed%Scw(ips,ich) = channelSed%Scw(ips,ich) - channelSed%Sabs(ips,ich,itime) * dt

      !    ! --- influx: water distribution - discharge !BK20240710 
      !    !--- BK20250705 commented out (to be restored later if the model can use MPI with time-lag simulations)
      !    !ich_dis = abs2dis_ich(ich,itime)
      !    !lag_ts = abs2dis_lag(ich,itime)
      !    ! if(SedCNP_hydro%Qabs(ich,itime) .gt. 0.) then
      !    !    if(ich_dis .ge. 1 .and. lag_ts .ge. 1) then
      !    !       if(itime+lag_ts .le. domain%ntime_sedcnp) then  !BK20240723
      !    !          channelSed%Sdis(ips,ich_dis,itime+lag_ts) = channelSed%Sdis(ips,ich_dis,itime+lag_ts) &
      !    !                                                    + channelSed%Sabs(ips,ich,itime) 
      !    !       endif
      !    !       !channelSed%Sxdis(ips,ich_dis) = channelSed%Sdis(ips,ich_dis,itime)*dt  !BK20240729                        
      !    !    endif
      !    ! endif
      !    !--- BK20250705 commented out

      !    channelSed%Sxdis(ips,ich) = channelSed%Sdis(ips,ich,itime)*dt  !BK20240705

      !    !update storage
      !    !channelSed%Scw(ips,ich) = channelSed%Scw(ips,ich) + channelSed%Sdis(ips,ich,itime) * dt

      ! !channelSed%Sinput(ips,ich) = channelSed%Susch(ips,ich) &       !BK20240618
      ! !channelSed%Sinput(ips,ich) = channelSed%Susch0(ips,ich) &       !BK20240618
      ! !                               + channelSed%Solch(ips,ich) &
      ! channelSed%Sinput(ips,ich) = channelSed%Susch0(ips,ich)*dt &       !BK20240618  !BK20250914
      !                                  + channelSed%Solch0(ips,ich)*dt &              !BK20250914  !BK20250925                  
      !                                  + channelSed%Sxpnt(ips,ich) &      !BK20240701
      !                                  - channelSed%Sxabs(ips,ich) &      !BK20240729
      !                                  + channelSed%Sxdis(ips,ich)        !BK20240729
                                       
      !                                  !- Sirri (should be included 
      !                                  !       once irrigation water extraction 
      !                                  !       from stream is represented in the model) 

      ! !update channel water sediment storage
      ! !channelSed%Scw(ips,ich) = channelSed%Scw(ips,ich) + channelSed%Sinput(ips,ich)
      ! channelSed%Scw(ips,ich) = max(0.0d0, channelSed%Scw(ips,ich) + channelSed%Sinput(ips,ich)) !BK20260522
      !---BK20260523 commented out (old block)


      !---BK20260523 new block
         ! --- influx: point sources
         channelSed%Sxpnt(ips,ich) = channelSed%Spnt(ips,ich,itime) * dt  !BK20240701

         ! --- influx: water distribution - discharge !BK20240710
         !--- BK20250705 commented out (to be restored later if the model can use MPI with time-lag simulations)
         !ich_dis = abs2dis_ich(ich,itime)
         !lag_ts = abs2dis_lag(ich,itime)
         ! if(SedCNP_hydro%Qabs(ich,itime) .gt. 0.) then
         !    if(ich_dis .ge. 1 .and. lag_ts .ge. 1) then
         !       if(itime+lag_ts .le. domain%ntime_sedcnp) then
         !          channelSed%Sdis(ips,ich_dis,itime+lag_ts) = ...
         !       endif
         !    endif
         ! endif
         !--- BK20250705 commented out
         channelSed%Sxdis(ips,ich) = channelSed%Sdis(ips,ich,itime)*dt  !BK20240705 !BK20260523 moved before abstraction

         !BK20260523: apply all non-abstraction inflows to Scw before computing abstraction,
         !mirroring the sequencing in module_Carbon/Nitrogen/Phosphorus. This ensures Sxabs
         !is nonzero even when Scw_old = 0 but sediment arrives from upstream/lateral in this step.
         channelSed%Scw(ips,ich) = max(0.0d0, channelSed%Scw(ips,ich)  &  !BK20260523
                                        + channelSed%Susch0(ips,ich)*dt &  !BK20260523
                                        + channelSed%Solch0(ips,ich)*dt &  !BK20260523
                                        + channelSed%Sxpnt(ips,ich)     &  !BK20260523
                                        + channelSed%Sxdis(ips,ich))       !BK20260523

         ! --- outflux: water distribution - abstraction !BK20240710
         !BK20260523: abs_ratio approach consistent with module_Carbon/Nitrogen/Phosphorus;
         !applied to Scw after non-abstraction inflows have already been added above.
         channelSed%Sxabs(ips,ich) = 0.0d0  !BK20260523
         if(SedCNP_hydro%Qabs(ich,itime) .gt. 0.) then
            if ((WStorage + Qdsch * dt) > 1.0E-9) then                                   !BK20260523
               abs_ratio_local = SedCNP_hydro%Qabs(ich,itime) / (WStorage + Qdsch * dt)  !BK20260523
            else                                                                           !BK20260523
               abs_ratio_local = 0.0d0                                                    !BK20260523
            endif                                                                          !BK20260523
            channelSed%Sabs(ips,ich,itime) = channelSed%Scw(ips,ich) * abs_ratio_local   !BK20260523
            channelSed%Sxabs(ips,ich)      = channelSed%Sabs(ips,ich,itime) * dt         !BK20260523
            channelSed%Scw(ips,ich) = max(0.0d0, channelSed%Scw(ips,ich) &               !BK20260523
                                            - channelSed%Sxabs(ips,ich))                  !BK20260523
         endif

         !BK20260523: Sinput is for output bookkeeping only — Scw is already fully updated
         !inline above (inflows added, abstraction subtracted). Do NOT add Sinput to Scw again.
         channelSed%Sinput(ips,ich) = channelSed%Susch0(ips,ich)*dt &      !BK20240618  !BK20250914
                                        + channelSed%Solch0(ips,ich)*dt &   !BK20250914  !BK20250925
                                        + channelSed%Sxpnt(ips,ich)     &   !BK20240701
                                        - channelSed%Sxabs(ips,ich)     &   !BK20240729
                                        + channelSed%Sxdis(ips,ich)         !BK20240729

         !BK20260523: commented out — Scw already updated inline above; applying Sinput here
         !would double-count all inflows.
         !channelSed%Scw(ips,ich) = max(0.0d0, channelSed%Scw(ips,ich) + channelSed%Sinput(ips,ich)) !BK20260522
      !---BK20260523 new block

   end subroutine sedInput

   !subroutine channelSedTransportEH(ich,ips,dt,Qdsch,Qusch,h,dp,va,Sf,BtmWdth,Lst,Area_x)
   subroutine channelSedTransportEH(ich,ips,dt,Qdsch,h,dp,va,Sf,BtmWdth,Lst,Area_x)  !BK20240519(ich>chid) !BK20241101(chid>ich)

      implicit none

      integer,               intent(in)  :: ich  
      integer,               intent(in)  :: ips
      real(8),               intent(in)  :: dt    !time-step [sec]
      real(8),               intent(in)  :: h     !channel flow depth [m]
      real(8),               intent(in)  :: Qdsch !channel discharge [m3/s]
      real(8), dimension(4), intent(in)  :: dp    !particle diameter [μm]
      real(8),               intent(in)  :: va    !advective (flow) velocity (in the down-gradient direction) [m/s2]
      real(8),               intent(in)  :: Sf    !channel (friction) slope [-]
      real(8),               intent(in)  :: BtmWdth   !channel width [m]
      real(8),               intent(in)  :: Lst   !channel length [m]
      real(8),               intent(in)  :: Area_x   !cross section in the direction of flow [m2]
      real(8)                            :: Jc           !sediment transport capacity areal flux [kg/m2/s]
      real(8)                            :: V            !BK20260523 bed layer volume [m3]
      real(8)                            :: Csb          !BK20260523 bed concentration [kg/m3]
      real(8)                            :: Csw          !BK20260523 water-column concentration [kg/m3]
      real(8)                            :: rho_b        !BK20260523 dry bulk density of bed [kg/m3]
      real(8)                            :: Je           !BK20260523 erosion flux [kg/m2/s]
      real(8)                            :: vr           !BK20260523 erosion velocity [m/s]
      real(8)                            :: bed_area     !BK20260523 channel bed area [m2]
      real(8)                            :: Scher_local  !BK20260523 erosion rate [kg/s]
      real(8)                            :: vs           !BK20260523 settling velocity [m/s]
      real(8)                            :: a_cstr, m_cstr, df_cstr, d_cstr  !BK20260523
      real(8)                            :: Wstg_local   !BK20260523 water storage [m3]
      real(8)                            :: k_settling   !BK20260523 [s^-1]
      real(8)                            :: k_discharge  !BK20260523 [s^-1]
      real(8)                            :: k_total      !BK20260523 [s^-1]
      real(8)                            :: Scw_ss       !BK20260523 CSTR steady-state Scw [kg]
      real(8)                            :: frac_remaining   !BK20260523
      real(8)                            :: integrated_Scw   !BK20260523 time-integrated Scw [kg·s]
      real(8)                            :: deposition_mass  !BK20260523 [kg]
      real(8)                            :: discharge_mass   !BK20260523 [kg]
      real(8)                            :: Scw_local    !BK20260523 beginning-of-step Scw [kg]
      real(8)                            :: Scw_new      !BK20260523 end-of-step Scw [kg]

      call channelSedTransportCapa(ich,ips,va,Sf,BtmWdth,h,dp,Area_x,Jc)

      !BK20260523: replaced sequential channelSedErosion + channelSedDeposition +
      !channelSedDischarge with a CSTR analytical solution. The prior sequential approach
      !overestimated channel deposition (~4x for clay/silt) and completely blocked downstream
      !sand transport (vs*dt/h >> 1 → cap removes all Scw before discharge runs).
      !The CSTR treats all three processes simultaneously from the same beginning-of-step Scw:
      !  d(Scw)/dt = Scher − k_total·Scw  (first-order linear ODE with constant source)
      !  Scw(t+dt) = Scw_ss + (Scw₀ − Scw_ss)·exp(−k_total·dt)
      !  where Scw_ss = Scher/k_total,  k_total = k_settling + k_discharge
      
      !call channelSedErosion(ich,ips,dt,Jc,va,BtmWdth,Lst,h,Area_x)       !BK20260523 replaced
      !call channelSedDeposition(ich,ips,dt,dp,BtmWdth,Lst,h,Qdsch,Area_x)  !BK20260523 replaced
      !call channelSedDischarge(ich,ips,Qdsch,dt)                            !BK20260523 replaced

      !--- Step 1: compute Scher (channelSedErosion logic) ---
      V     = BtmWdth * Lst * 0.01d0                                                !BK20260523
      rho_b = 1300.d0                                                                !BK20260523
      if (V > 1.0E-9) then                                                           !BK20260523
         channelSed%Sch(ips,ich) = min(channelSed%Sch(ips,ich), V * rho_b)         !BK20260523
         Csb = channelSed%Sch(ips,ich) / V                                          !BK20260523
      else                                                                           !BK20260523
         Csb = 0.0d0                                                                !BK20260523
      endif                                                                          !BK20260523
      if (SedCNP_hydro%Wstg_st(ich) > 1.0E-9) then                                 !BK20260523
         Csw = channelSed%Scw(ips,ich) / SedCNP_hydro%Wstg_st(ich)                !BK20260523
      else                                                                           !BK20260523
         Csw = 0.0d0                                                                !BK20260523
      endif                                                                          !BK20260523
      if (Jc > va * Csw) then                                                       !BK20260523
         vr = max(0.0d0, (Jc - va * Csw) / rho_b)                                  !BK20260523
      else                                                                           !BK20260523
         vr = 0.0d0                                                                 !BK20260523
      endif                                                                          !BK20260523
      Je          = vr * Csb                                                         !BK20260523
      bed_area    = BtmWdth * Lst                                                    !BK20260523
      Scher_local = Je * bed_area                                                    !BK20260523
      Scher_local = min(max(0.0d0, channelSed%Sch(ips,ich)), &                      !BK20260523
                        max(0.0d0, Scher_local * dt)) / dt                          !BK20260523
      channelSed%Scher(ips,ich) = Scher_local                                       !BK20260523

      !--- Step 2: compute settling velocity vs (channelSedDeposition logic) ---
      a_cstr = 8.4E-3                                                                !BK20260523
      m_cstr = 0.024d0                                                               !BK20260523
      if (ips == 1) then                                                             !BK20260523
         df_cstr = dp(ips)                                                           !BK20260523
         vs = a_cstr * (df_cstr**m_cstr) / 100.d0                                   !BK20260523
      else                                                                           !BK20260523
         d_cstr = dp(ips)*1.0E-6 * (((Gs - 1.) * g / (v**2.))**(1./3.))           !BK20260523
         vs = v / (dp(ips)*1.0E-6) &                                                !BK20260523
              * (((25. + 1.2 * (d_cstr**2.))**0.5 - 5.)**(1.5))                    !BK20260523
      endif                                                                          !BK20260523

      !--- Step 3: rate constants [s^-1] ---
      Wstg_local = SedCNP_hydro%Wstg_st(ich)                                        !BK20260523
      if (Wstg_local > 1.0E-12) then                                                !BK20260523
         k_settling  = vs * bed_area / Wstg_local                                   !BK20260523
         k_discharge = Qdsch / Wstg_local                                           !BK20260523
      else                                                                           !BK20260523
         k_settling  = 0.0d0                                                        !BK20260523
         k_discharge = 0.0d0                                                        !BK20260523
      endif                                                                          !BK20260523
      k_total = k_settling + k_discharge                                             !BK20260523

      !--- Step 4: CSTR analytical solution ---
      Scw_local = channelSed%Scw(ips,ich)                                           !BK20260523
      if (k_total > 0.0d0) then                                                     !BK20260523
         Scw_ss         = Scher_local / k_total                                     !BK20260523
         frac_remaining = exp(-k_total * dt)                                        !BK20260523
         integrated_Scw = Scw_ss * dt &                                             !BK20260523
                          - (Scw_ss - Scw_local) * (1.0d0 - frac_remaining) / k_total  !BK20260523
         integrated_Scw = max(0.0d0, integrated_Scw)                               !BK20260523
         deposition_mass = k_settling  * integrated_Scw                             !BK20260523
         discharge_mass  = k_discharge * integrated_Scw                             !BK20260523
         Scw_new         = Scw_ss + (Scw_local - Scw_ss) * frac_remaining          !BK20260523
      else                                                                           !BK20260523
         deposition_mass = 0.0d0                                                    !BK20260523
         discharge_mass  = 0.0d0                                                    !BK20260523
         Scw_new         = Scw_local + Scher_local * dt                             !BK20260523
      endif                                                                          !BK20260523

      !--- Step 5: update state variables and output fields ---
      !mass conserved: Scw_new = Scw_local + Scher*dt − deposition_mass − discharge_mass
      !Sch_new ≥ 0 guaranteed: Scher*dt ≤ Sch (from cap above), deposition_mass ≥ 0
      channelSed%Scw(ips,ich)   = max(0.0d0, Scw_new)                              !BK20260523
      channelSed%Sch(ips,ich)   = max(0.0d0, channelSed%Sch(ips,ich) &             !BK20260523
                                   - Scher_local * dt + deposition_mass)            !BK20260523
      channelSed%Schdp(ips,ich) = max(0.0d0, deposition_mass) / dt                 !BK20260523
      channelSed%Sdsch(ips,ich) = max(0.0d0, discharge_mass)  / dt                 !BK20260523
      if (Wstg_local > 1.0E-9) then                                                 !BK20260523
         channelSed%Cchsd(ips,ich) = channelSed%Scw(ips,ich) / Wstg_local         !BK20260523
      else                                                                           !BK20260523
         channelSed%Cchsd(ips,ich) = 0.0d0                                         !BK20260523
      endif                                                                          !BK20260523

   end subroutine channelSedTransportEH


   subroutine channelSedTransportCapa(ich,ips,va,Sf,BtmWdth,h,dp,Area_x,Jc)  !BK20240519(ich>chid) !BK20241101(chid>ich)

      implicit none

      integer,               intent(in)  :: ich      !channel index
      integer,               intent(in)  :: ips
      real(8),               intent(in)  :: va       !advective (flow) velocity (in the down-gradient direction) [m/s2]
      real(8),               intent(in)  :: Sf       !channel (friction) slope [-]
      real(8),               intent(in)  :: BtmWdth  !channel width [m]
      real(8),               intent(in)  :: h        !channel flow depth [m]
      real(8), dimension(4), intent(in)  :: dp       !particle diameter [μm]
      real(8),               intent(in)  :: Area_x   !cross section in the direction of flow [m2]
      real(8),               intent(out) :: Jc       !sediment transport capacity areal flux [kg/m2/s]
      real(8)                            :: Cw       !concentration of entrained sediment particles by weight at the transport capacity [-]
      real(8)                            :: Ct       !concentration of entrained sediment particles at the transport capacity 
      !real(8)                            :: Gs       !particle specific gravity [-]
      !real(8)                            :: g        !gravitation acceleration [m/s2]
      real(8)                            :: Rh       !hydraulic radius of flow
      real(8), dimension(4)              :: vc       !critical velocity for erosion


      !Gs = 2.65
      !g = 9.8

      !hlcho 20221108 아래 클래스별 vc값 설정 (문헌을 통해 해당값 도출 필요)
      !hard-coded; should be relocated or moved to config file
      vc(1) = 0.1      !clay [m/s]
      vc(2) = 0.2      !silt [m/s]
      vc(3) = 0.5      !fine sand [m/s]
      vc(4) = 1.0      !coarse sand  [m/s]

      !concentration of entrained sediment particles by weight at the transport capacity
      !Rh = BtmWdth * h / (BtmWdth + 2. * h)
      !BK20250705 modification for numerical stability
      if ((BtmWdth + 2. * h) > 1.0E-9) then
         Rh = BtmWdth * h / (BtmWdth + 2. * h)
      else
         Rh = 0.0
      endif

      if (va > vc(ips)) then
         Cw = 0.05 * (Gs / (Gs - 1.)) * &
               ((va - vc(ips)) * Sf) / (((Gs - 1.) * g * dp(ips)*1.0E-6)**0.5) * &  !BK20250702 corrected to (-6.)
               !((Rh * Sf) / (((Gs - 1.) * dp(ips)*1.0E-6)**0.5))   !BK20250702 corrected to (-6.) 
               ((Rh * Sf) / ((Gs - 1.) * dp(ips)*1.0E-6))**0.5  !BK20260518 Corrected to be dimensionless Shields parameter

         ! The empirical formula for Cw (a weight fraction) can yield unphysically large
         ! values (> 1.0), especially for small clay particles. This causes the denominator
         ! in the Ct calculation to become negative, resulting in a negative transport capacity.
         ! Cw is capped at a physical limit (e.g., 0.8) to prevent this.
         Cw = min(Cw, 0.8d0) !BK20260518

         if (abs(Gs + (1. - Gs) * Cw) > 1.0E-9) then  !BK20250810
            !Ct = 1.0E-6 * Gs * Cw / (Gs + (1. - Gs) * Cw)   ! g m-3
            !Ct = 1.0E-3 * Gs * Cw / (Gs + (1. - Gs) * Cw)   ! kg m-3  !BK20260219
            ! The conversion from weight fraction (Cw) to mass concentration (Ct, kg/m3)
            ! requires multiplying by water density (~1000 kg/m3). The previous factor was incorrect.
            Ct = 1000.d0 * Gs * Cw / (Gs + (1. - Gs) * Cw)   ! kg m-3 !BK20260518  
         else
            Ct = 0.0
         endif
      else
         Cw = 0.0
         Ct = 0.0
      endif

      !sediment transport capacity areal flux
      Jc = va * Ct   ! kg m-2 s-1

      !cross section in the direction of flow [m2]
      !Area_x = BtmWdth * h

      !sediment transport capacity
      channelSed%Schcp(ips,ich) = Jc * Area_x  ! kg s-1

   end subroutine channelSedTransportCapa


   subroutine channelSedErosion(ich,ips,dt,Jc,va,BtmWdth,Lst,h,Area_x)  !BK20240519(ich>chid) !BK20241101(chid>ich)

      implicit none

      integer, intent(in)    :: ich     !channel index
      integer, intent(in)    :: ips
      real(8),    intent(in) :: dt      !time-step [sec] (get this from config file)
      real(8),    intent(in) :: Jc      !sediment transport capacity areal flux [kg/m2/s]
      real(8),    intent(in) :: va      !advective (flow) velocity (in the down-gradient direction) [m/s] !BK20260518
      real(8),    intent(in) :: BtmWdth !channel width [m]
      real(8),    intent(in) :: Lst     !channel length [m]
      real(8),    intent(in) :: h       !channel flow depth [m]
      real(8),    intent(in) :: Area_x  !cross section in the direction of flow [m2]
      real(8)                :: V       !sediment storage [m3]
      real(8)                :: Csb     !concentration of sediment at the bottom boundary (in the bed) [kg/m3] !BK20260518
      real(8)                :: Csw     !concentration of sediment in water column [kg/m3]
      real(8)                :: WaterVolume ! [m3]
      real(8)                :: rho_b   ! bulk density of bed sediment [kg/m3]
      real(8)                :: Je      !erosion flux [kg/m2/s] !BK20260518
      real(8)                :: vr      !resuspension (erosion) velocity [m/s] !BK20260518
      real(8)                :: bed_area ! [m2]


      !channel sediment storage volume [m3]
      V = BtmWdth * Lst * 0.01

      !concentration of sediment at the bottom boundary (in the bed)
      !BK20250705 modification for numerical stability
      if (V > 1.0E-9) then
         ! Enforce physical bed capacity: Csb cannot exceed bulk dry density
         rho_b = 1300.d0  ! Dry bulk density of bed material [kg/m3], hardcoded for now. !BK20260522
         channelSed%Sch(ips,ich) = min(channelSed%Sch(ips,ich), V * rho_b) !BK20260522
         Csb = channelSed%Sch(ips,ich) / V
      else
         Csb = 0.0
      endif

      ! Concentration of sediment in the water column
      ! WaterVolume = Area_x * Lst !---BK20260518
      ! if (WaterVolume > 1.0E-9) then
      !    Csw = channelSed%Scw(ips,ich) / WaterVolume
      ! else
      !    Csw = 0.0
      ! endif 
      if (SedCNP_hydro%Wstg_st(ich) > 1.0E-9) then !---BK20260522
         Csw = channelSed%Scw(ips,ich) / SedCNP_hydro%Wstg_st(ich)
      else
         Csw = 0.0
      endif !---BK20260522

      !resuspension (erosion) velocity [m/s] !---BK20260518
      !vr = Jc / Csb - va
      !vr = max(vr, 0.0)
      ! The original equation was unstable (division by zero if Csb=0), causing NaN values.
      ! It is replaced with a more physically robust formulation based on the difference
      ! between transport capacity and advective flux.
      if (Jc > va * Csw) then
         vr = (Jc - va * Csw) / rho_b
         vr = max(vr, 0.0)
      else
         vr = 0.0
      endif !---BK20260518

      !erosion flux [kg/m2/s] 
      Je = vr * Csb

      ! The erosion area is the channel bed area. The original code incorrectly used
      ! the flow cross-sectional area (Area_x).
      bed_area = BtmWdth * Lst  !BK20260518

      !sediment erosion [kg/s]
      channelSed%Scher(ips,ich) = Je * bed_area !BK20260518
      !channelSed%Scher(ips,ich) = min(channelSed%Sch(ips,ich), channelSed%Scher(ips,ich) * dt) / dt
      channelSed%Scher(ips,ich) = min(max(0.0d0, channelSed%Sch(ips,ich)), &
                                max(0.0d0, channelSed%Scher(ips,ich) * dt)) / dt  !BK20260522

      !BK20260523: state update deferred to channelSedDeposition (net flux approach)
      !sediment storage update
      !channelSed%Scw(ips,ich) = channelSed%Scw(ips,ich) + channelSed%Scher(ips,ich) * dt
      !channelSed%Sch(ips,ich) = channelSed%Sch(ips,ich) - channelSed%Scher(ips,ich) * dt

   end subroutine channelSedErosion


   !subroutine channelSedDeposition(ich,ips,dt,dp,BtmWdth,h,Qdsch,Qusch,Area_x)
   subroutine channelSedDeposition(ich,ips,dt,dp,BtmWdth,Lst,h,Qdsch,Area_x)  !BK20240519(ich>chid) !BK20241101(chid>ich) !BK20260518

      implicit none

      integer,              intent(in)    :: ich    !channel index
      integer,              intent(in)    :: ips
      real(8),                 intent(in) :: dt      !time-step [sec] (get this from config file)
      real(8),   dimension(4), intent(in) :: dp      !particle diameter [μm]
      real(8),                 intent(in) :: BtmWdth     !channel width [m]
      real(8),                 intent(in) :: Lst     !channel length [m]
      real(8),                 intent(in) :: h       !channel flow depth [m]
      real(8),                 intent(in) :: Qdsch   !channel discharge [m3/s]
      real(8),                 intent(in) :: Area_x  !cross section in the direction of flow [m2]
      real(8)                             :: vs      !quiescent settling velocity [m/s]
      !real(8)                             :: Gs      !particle specific gravity
      !real(8)                             :: g       !gravitational acceleration [m/s2]
      !real(8)                             :: v       !kinematic viscosity of water [m2/s]
      real(8)                             :: a       !experimentally determined constant
      real(8)                             :: m       !experimentally determined constant
      real(8)                             :: d       !dimensionless particle diameter [-]
      real(8)                             :: df      !median floc diameter [um]
      real(8)                             :: Jd      !deposition flux [kg/m2/s]
      real(8)                             :: bed_area  ! [m2]


      !Gs = 2.65
      !g = 9.8
      !v = 0.0001
      a = 8.4E-3
      m = 0.024

      !temporary coded, need to be fixed (kyh_debug)
      !Wstg_st = 10. * 86400.
      !SedCNP_hydro%Wstg_st(ich) = SedCNP_hydro%Wstg_st0(ich) + Qusch - Qdsch
      !if (SedCNP_hydro%Wstg_st(ich) < 0.0) then
      !   SedCNP_hydro%Wstg_st(ich) = 0.0
      !endif
      ! WRF-Hydro에서 확인
      ! 채널에 남아있는 유량 (Wstg_st = Wstg_st0 + Qusch - Qdsch)

      !settling velocity 
      if(ips.eq.1) then  
         !clay (Burban et al., 1990)
         df = dp(ips) 
         vs = a * (df**m)  ![cm/s]
         vs = vs / 100.    ![m/s]
      else 
         !sand, silt (Cheng, 1997)
         ! 𝑑 = dimensionless particle diameter [dimensionless]
         d = dp(ips)*1.0E-6 * (((Gs - 1.) * g / (v**2.))**(1./3.))
         !vs = v / dp(ips)*1.0E-6 * (((25. + 1.2 * (d**2.))**0.5 - 5.)**(1.5))  ![m/s]
         vs = v / (dp(ips)*1.0E-6) * (((25. + 1.2 * (d**2.))**0.5 - 5.)**(1.5))  ![m/s]  !BK20260127 
      endif

      !deposition flux
      if (SedCNP_hydro%Wstg_st(ich) > 0.0) then
         channelSed%Cchsd(ips,ich) = channelSed%Scw(ips,ich) / SedCNP_hydro%Wstg_st(ich)
      else
         channelSed%Cchsd(ips,ich) = 0.0
      endif
      Jd = vs * channelSed%Cchsd(ips,ich)

      ! The deposition area is the channel bed area.
      bed_area = BtmWdth * Lst !BK20260518

      !sediment deposition [kg s-1]
      !channelSed%Schdp(ips,ich) = Jd * Area_x
      channelSed%Schdp(ips,ich) = Jd * bed_area !BK20260518
      !channelSed%Schdp(ips,ich) = min(channelSed%Scw(ips,ich), channelSed%Schdp(ips,ich) * dt) / dt
      channelSed%Schdp(ips,ich) = min(max(0.0d0, channelSed%Scw(ips,ich)), &
                                max(0.0d0, channelSed%Schdp(ips,ich) * dt)) / dt !BK20260522

      !BK20260523: net flux approach — erosion pick-up and settling deposition both computed from
      !beginning-of-step Scw (erosion no longer pre-updates Scw), applied together as net flux.
      !Non-negativity is guaranteed: Scher*dt <= Sch and Schdp*dt <= Scw from upstream caps.
      channelSed%Scw(ips,ich) = channelSed%Scw(ips,ich) + &
                                 (channelSed%Scher(ips,ich) - channelSed%Schdp(ips,ich)) * dt
      channelSed%Sch(ips,ich) = channelSed%Sch(ips,ich) - &
                                 (channelSed%Scher(ips,ich) - channelSed%Schdp(ips,ich)) * dt
      channelSed%Scw(ips,ich) = max(0.0d0, channelSed%Scw(ips,ich))  !BK20260523
      channelSed%Sch(ips,ich) = max(0.0d0, channelSed%Sch(ips,ich))  !BK20260523


      !================= Calculation of probability of deposition !BK20260219
      !--- calculation of tau_0
      ! tau_0 = rho_water * g * Rh * Sf
      !
      !--- for non-cohesive particles
      !
      ! # tau_cd_n [N/m2]
      ! tau_cd_n = shields_c * (rho_sed - rho_water) * g * dp
      ! shields_c = 0.047 (Gessler, 1971)
      ! rho_sed = 2650 kg/m3
      ! rho_water = 1000 kg/m3
      ! g = 9.8 m/s2
      ! dp = particle diameter [m]
      !
      ! # P_dep
      ! if (tau_0 <= 0.0) then
      !    ! If there is no flow/shear stress, particles will definitely deposit
      !    P_dep = 1.0  
      ! else
      !    ! Calculate the Y parameter for non-cohesive sediment
      !    Y = (1.0 / sigma) * ((tau_cd_n / tau_0) - 1.0)
      !   
      !    ! Calculate the standard normal CDF using the intrinsic error function
      !    P_dep = 0.5 * (1.0 + erf(Y / sqrt(2.0)))
      ! end if

      !--- for cohesive particles
      !
      ! # tau_cd_c [N/m2]
      ! user-specific parameter 
      ! Standard default values for environmental modeling 
      ! typically fall between 0.01 and 0.1 N/m²
      !
      ! # P_dep
      ! if (tau_0 >= tau_cd_c) then
      !    ! Turbulence is too high. Flocs are torn apart and cannot settle.
      !    P_dep = 0.0
      ! else
      !    ! Flow is calm enough. Calculate Y and the deposition probability.
      !    ! Note: 'sigma' is the standard deviation parameter from your model
         
      !    inner_term = 0.25 * ((tau_cd_c / tau_0) - 1.0) * exp(1.27 * tau_cd_c)
      !    Y = (1.0 / sigma) * log(inner_term)
         
      !    ! Using the highly optimized complementary error function
      !    P_dep = 0.5 * erfc(Y / sqrt(2.0))
      ! end if
      !==================== Calculation of probability of deposition !BK20260219

   end subroutine channelSedDeposition


   subroutine channelSedDischarge(ich,ips,Qdsch,dt)  !BK20240519(ich>chid) !BK20241101(chid>ich)

      implicit none

      integer,    intent(in) :: ich  !channel index
      integer,    intent(in) :: ips
      real(8),    intent(in) :: Qdsch !channel discharge (should be read from Noah-MP outpus) [m3/s]
      !real(8)                :: Wstg_st !stream water storage [m3]
      real(8),    intent(in) :: dt    !time-step [sec] (get this from config file)


      if (SedCNP_hydro%Wstg_st(ich) > 0.0) then
         channelSed%Cchsd(ips,ich) = channelSed%Scw(ips,ich) / SedCNP_hydro%Wstg_st(ich)  !kyh_debug; WRF-Hydro에서 Wstg가 어떻게 구현 및 출력되는지확인필요
      else
         channelSed%Cchsd(ips,ich) = 0.0
      endif
      channelSed%Sdsch(ips,ich) = Qdsch * channelSed%Cchsd(ips,ich)  ! [kg s-1]
      !channelSed%Sdsch(ips,ich) = min(channelSed%Scw(ips,ich), channelSed%Sdsch(ips,ich) * dt) / dt  
      channelSed%Sdsch(ips,ich) = min(max(0.0d0, channelSed%Scw(ips,ich)), &
                                max(0.0d0, channelSed%Sdsch(ips,ich) * dt)) / dt !BK20260522

      !sediment storage update [kg]
      channelSed%Scw(ips,ich) = channelSed%Scw(ips,ich) - channelSed%Sdsch(ips,ich) * dt
      channelSed%Scw(ips,ich) = max(0.0d0, channelSed%Scw(ips,ich))  !BK20260522 flush fp rounding residual

   !     !kyh_debug
   !     if(ich == 50) then
   !        open (2445212, file='log_file_channelSedDischarge.txt', status='unknown')
   !        write(2445212, *) 'SedCNP_hydro%Wstg_st(ich)= ', SedCNP_hydro%Wstg_st(ich)
   !        write(2445212, *) 'channelSed%Cchsd(ips,ich)= ', channelSed%Cchsd(ips,ich)
   !        write(2445212, *) 'Qdsch= ', Qdsch
   !        write(2445212, *) 'dt= ', dt
   !        write(2445212, *) 'channelSed%Sdsch(ips,ich)= ', channelSed%Sdsch(ips,ich)
   !        write(2445212, *) 'channelSed%Scw(ips,ich)= ', channelSed%Scw(ips,ich)
   !     endif

   end subroutine channelSedDischarge


end module module_channelSed
