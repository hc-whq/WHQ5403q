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

module module_CNPupdates

   !use CNPvariables  !Y.Kwon
   use module_SedCNPvariables    !BK20230310
!=====||__WHQ5403q__||=====!
!
   use SedCNP_config,          only: SedCNPmodel
!
!=====||__WHQ5403q__||=====!

   contains

   subroutine UpdateCNP_So (i,j)


      !to update the soil storages/variables for the next time step

      IMPLICIT NONE

      INTEGER :: i, j, isl
      integer :: ips            !BK20240618

      ! === Sediment: overland  !BK20240618 
      do ips = 1, nps
      !    ! update storages
      !    overSed%Sow0(ips,i,j) = overSed%Sow(ips,i,j) !BK20250925 obsolete updates - commented out
      !    overSed%Sol0(ips,i,j) = overSed%Sol(ips,i,j) !BK20250925 obsolete updates - commented out
         ! update transport variables
         overSed%Ssurf0(ips,i,j) = overSed%Ssurf(ips,i,j)
      end do


      ! === C
      ! --- update storages
      So_LPOCLIT0(i,j) = So_LPOCLIT(i,j)
      So_LPOCEXC0(i,j) = So_LPOCEXC(i,j)
      So_LPOCMAN0(i,j) = So_LPOCMAN(i,j)

      DO isl = 1, domain%nsl
         So_MBMC0(i,isl,j) = So_MBMC(i,isl,j)
         So_LPOCRES0(i,isl,j) = So_LPOCRES(i,isl,j)
         So_RPOC0(i,isl,j) = So_RPOC(i,isl,j)
         So_LDOC0(i,isl,j) = So_LDOC(i,isl,j)
         So_RDOC0(i,isl,j) = So_RDOC(i,isl,j)
      END DO

      ! --- update transport variables
      DO isl = 1, domain%nsl
         So_LDOCintf0(i,isl,j) = So_LDOCintf(i,isl,j)
         So_RDOCintf0(i,isl,j) = So_RDOCintf(i,isl,j)
         So_LDOCperc0(i,isl,j) = So_LDOCperc(i,isl,j)
         So_RDOCperc0(i,isl,j) = So_RDOCperc(i,isl,j)
      END DO

      So_LDOCgwso0(i,j) = So_LDOCgwso(i,j)  !BK20260323
      So_RDOCgwso0(i,j) = So_RDOCgwso(i,j)  !BK20260323
      So_LPOCsurf0(i,j) = So_LPOCsurf(i,j)
      So_RPOCsurf0(i,j) = So_RPOCsurf(i,j)
      So_MBMCsurf0(i,j) = So_MBMCsurf(i,j)
      So_LDOCsurf0(i,j) = So_LDOCsurf(i,j)
      So_RDOCsurf0(i,j) = So_RDOCsurf(i,j)
      So_LDOCsogw0(i,j) = So_LDOCsogw(i,j)  !BK20260318
      So_RDOCsogw0(i,j) = So_RDOCsogw(i,j)  !BK20260318

      ! === N
      ! --- update storages
      So_LPONLIT0(i,j) = So_LPONLIT(i,j)
      So_LPONEXC0(i,j) = So_LPONEXC(i,j)
      So_LPONMAN0(i,j) = So_LPONMAN(i,j)
      
      DO isl = 1, domain%nsl
         So_MBMN0(i,isl,j) = So_MBMN(i,isl,j)
         So_LPONRES0(i,isl,j) = So_LPONRES(i,isl,j)
         So_RPON0(i,isl,j) = So_RPON(i,isl,j)
         So_LDON0(i,isl,j) = So_LDON(i,isl,j)
         So_RDON0(i,isl,j) = So_RDON(i,isl,j)
         So_NH40(i,isl,j) = So_NH4(i,isl,j)
         So_NO30(i,isl,j) = So_NO3(i,isl,j)
      END DO

      ! --- update transport variables
      DO isl = 1, domain%nsl
         So_LDONintf0(i,isl,j) = So_LDONintf(i,isl,j)
         So_RDONintf0(i,isl,j) = So_RDONintf(i,isl,j)
         So_NH4intf0(i,isl,j) = So_NH4intf(i,isl,j)
         So_NO3intf0(i,isl,j) = So_NO3intf(i,isl,j)
         So_LDONperc0(i,isl,j) = So_LDONperc(i,isl,j)
         So_RDONperc0(i,isl,j) = So_RDONperc(i,isl,j)
         So_NH4perc0(i,isl,j) = So_NH4perc(i,isl,j)
         So_NO3perc0(i,isl,j) = So_NO3perc(i,isl,j)
      END DO

      So_LDONgwso0(i,j) = So_LDONgwso(i,j)  !BK20260323
      So_RDONgwso0(i,j) = So_RDONgwso(i,j)  !BK20260323
      So_NH4gwso0(i,j) = So_NH4gwso(i,j)  !BK20260323
      So_NO3gwso0(i,j) = So_NO3gwso(i,j)  !BK20260323
      So_LPONsurf0(i,j) = So_LPONsurf(i,j)
      So_RPONsurf0(i,j) = So_RPONsurf(i,j)
      So_MBMNsurf0(i,j) = So_MBMNsurf(i,j)
      So_LDONsurf0(i,j) = So_LDONsurf(i,j)
      So_RDONsurf0(i,j) = So_RDONsurf(i,j)
      So_NH4surf0(i,j) = So_NH4surf(i,j)
      So_NO3surf0(i,j) = So_NO3surf(i,j)
      So_LDONsogw0(i,j) = So_LDONsogw(i,j)  !BK20260318
      So_RDONsogw0(i,j) = So_RDONsogw(i,j)  !BK20260318
      So_NH4sogw0(i,j) = So_NH4sogw(i,j)  !BK20260318
      So_NO3sogw0(i,j) = So_NO3sogw(i,j)  !BK20260318

      ! === P
      ! --- update storages
      So_LPOPLIT0(i,j) = So_LPOPLIT(i,j)
      So_LPOPEXC0(i,j) = So_LPOPEXC(i,j)
      So_LPOPMAN0(i,j) = So_LPOPMAN(i,j)

      DO isl = 1, domain%nsl
         So_MBMP0(i,isl,j) = So_MBMP(i,isl,j)
         So_LPOPRES0(i,isl,j) = So_LPOPRES(i,isl,j)
         So_RPOP0(i,isl,j) = So_RPOP(i,isl,j)
         So_LDOP0(i,isl,j) = So_LDOP(i,isl,j)
         So_RDOP0(i,isl,j) = So_RDOP(i,isl,j)
         So_PO40(i,isl,j) = So_PO4(i,isl,j)
         So_PIPA0(i,isl,j) = So_PIPA(i,isl,j)
         So_PIPS0(i,isl,j) = So_PIPS(i,isl,j)
      END DO

      ! --- update transport variables
      DO isl = 1, domain%nsl
         So_LDOPintf0(i,isl,j) = So_LDOPintf(i,isl,j)
         So_RDOPintf0(i,isl,j) = So_RDOPintf(i,isl,j)
         So_PO4intf0(i,isl,j) = So_PO4intf(i,isl,j)
         So_LDOPperc0(i,isl,j) = So_LDOPperc(i,isl,j)
         So_RDOPperc0(i,isl,j) = So_RDOPperc(i,isl,j)
         So_PO4perc0(i,isl,j) = So_PO4perc(i,isl,j)
      END DO

      So_LDOPgwso0(i,j) = So_LDOPgwso(i,j)  !BK20260323
      So_RDOPgwso0(i,j) = So_RDOPgwso(i,j)  !BK20260323
      So_PO4gwso0(i,j) = So_PO4gwso(i,j)  !BK20260323
      So_LPOPsurf0(i,j) = So_LPOPsurf(i,j)
      So_RPOPsurf0(i,j) = So_RPOPsurf(i,j)
      So_MBMPsurf0(i,j) = So_MBMPsurf(i,j)
      So_LDOPsurf0(i,j) = So_LDOPsurf(i,j)
      So_RDOPsurf0(i,j) = So_RDOPsurf(i,j)
      So_PO4surf0(i,j) = So_PO4surf(i,j)
      So_PIPAsurf0(i,j) = So_PIPAsurf(i,j)
      So_PIPSsurf0(i,j) = So_PIPSsurf(i,j)
      So_LDOPsogw0(i,j) = So_LDOPsogw(i,j)  !BK20260318
      So_RDOPsogw0(i,j) = So_RDOPsogw(i,j)  !BK20260318
      So_PO4sogw0(i,j) = So_PO4sogw(i,j)  !BK20260318

   end subroutine UpdateCNP_So


   subroutine UpdateCNP_Gw (gwid)  !BK20240614

      !to update the subwatershed groundwater storages/variables for the next time step

      IMPLICIT NONE
      !integer :: sbid
      integer :: gwid  !BK20240614

      ! === C
      ! --- update storages
      Gw_LDOC0(gwid) = Gw_LDOC(gwid)
      Gw_RDOC0(gwid) = Gw_RDOC(gwid)

      ! --- update transport variables
      Gw_LDOCgwch0(gwid) = Gw_LDOCgwch(gwid)
      Gw_RDOCgwch0(gwid) = Gw_RDOCgwch(gwid)

      ! === N
      ! --- update storages
      Gw_LDON0(gwid) = Gw_LDON(gwid)
      Gw_RDON0(gwid) = Gw_RDON(gwid)
      Gw_NH40(gwid) = Gw_NH4(gwid)
      Gw_NO30(gwid) = Gw_NO3(gwid)

      ! --- update transport variables
      Gw_LDONgwch0(gwid) = Gw_LDONgwch(gwid)
      Gw_RDONgwch0(gwid) = Gw_RDONgwch(gwid)
      Gw_NH4gwch0(gwid) = Gw_NH4gwch(gwid)
      Gw_NO3gwch0(gwid) = Gw_NO3gwch(gwid)

      ! === P
      ! --- update storages
      Gw_LDOP0(gwid) = Gw_LDOP(gwid)
      Gw_RDOP0(gwid) = Gw_RDOP(gwid)
      Gw_PO40(gwid) = Gw_PO4(gwid)

      ! --- update transport variables
      Gw_LDOPgwch0(gwid) = Gw_LDOPgwch(gwid)
      Gw_RDOPgwch0(gwid) = Gw_RDOPgwch(gwid)
      Gw_PO4gwch0(gwid) = Gw_PO4gwch(gwid)

   end subroutine UpdateCNP_Gw


   !subroutine UpdateCNP_Ch (ich)   !BK20231029 !BK20240513 !BK20240519(ich>chid) !BK20241102(chid>ich)
   subroutine UpdateCNP_Ch (ich)

      !to update the channel storages/variables for the next time step

      IMPLICIT NONE

      integer,intent(in)   :: ich
      INTEGER              :: chid           !channel ID
      integer              :: iich           !channel index for the loop
      integer              :: iupch          !channel index for upstream channels
      integer              :: ia  
      integer              :: ips            !BK20240618

      chid = SedCNP_hydro%linkID(ich)  !BK20231029 !BK20240513

      ! === Sediment: channel !BK20240618
      do ips = 1, nps
         ! update storages
         channelSed%Scw0(ips,ich) = channelSed%Scw(ips,ich)
         channelSed%Sch0(ips,ich) = channelSed%Sch(ips,ich)
      end do

      ! === C
      ! --- update storages
      Ch_ALGC0(ich,1:nalg) = Ch_ALGC(ich,1:nalg)
      Ch_ZOOC0(ich) = Ch_ZOOC(ich)
      Ch_MBMC0(ich) = Ch_MBMC(ich)
      Ch_DIC0(ich)  = Ch_DIC(ich)  !BK20251120
      Ch_LPOC0(ich) = Ch_LPOC(ich)
      Ch_RPOC0(ich) = Ch_RPOC(ich)
      Ch_LDOC0(ich) = Ch_LDOC(ich)
      Ch_RDOC0(ich) = Ch_RDOC(ich)
      Ch_POCDEPOSIT0(ich) = Ch_POCDEPOSIT(ich)

      ! === N
      ! --- update storages
      Ch_ALGN0(ich,1:nalg) = Ch_ALGN(ich,1:nalg)
      Ch_ZOON0(ich) = Ch_ZOON(ich)
      Ch_MBMN0(ich) = Ch_MBMN(ich)
      Ch_LPON0(ich) = Ch_LPON(ich)
      Ch_RPON0(ich) = Ch_RPON(ich)
      Ch_LDON0(ich) = Ch_LDON(ich)
      Ch_RDON0(ich) = Ch_RDON(ich)
      Ch_NH40(ich)  = Ch_NH4(ich)
      Ch_NO30(ich)  = Ch_NO3(ich)
      Ch_PONDEPOSIT0(ich) = Ch_PONDEPOSIT(ich)

      ! === P
      ! --- update storages
      Ch_ALGP0(ich,1:nalg) = Ch_ALGP(ich,1:nalg)
      Ch_ZOOP0(ich) = Ch_ZOOP(ich)
      Ch_MBMP0(ich) = Ch_MBMP(ich)
      Ch_LPOP0(ich) = Ch_LPOP(ich)
      Ch_RPOP0(ich) = Ch_RPOP(ich)
      Ch_LDOP0(ich) = Ch_LDOP(ich)
      Ch_RDOP0(ich) = Ch_RDOP(ich)
      Ch_PO40(ich)  = Ch_PO4(ich)
      Ch_PIPA0(ich) = Ch_PIPA(ich)
      Ch_PIPS0(ich) = Ch_PIPS(ich)
      Ch_POPDEPOSIT0(ich) = Ch_POPDEPOSIT(ich)

      
       !-------------------------------------------------------Y.Kwon 20230604
       ! --- update transport variables (CNP)
       ! *** For each of the channel links, calculate sum of the DSCH variables of the upstream channels ***

       !initialize
       !C
       Ch_ALGCusch0(ich,:) = 0.0
       Ch_ZOOCusch0(ich) = 0.0
       Ch_MBMCusch0(ich) = 0.0
       Ch_DICusch0(ich)  = 0.0  !BK20251121
       Ch_LPOCusch0(ich) = 0.0
       Ch_RPOCusch0(ich) = 0.0
       Ch_LDOCusch0(ich) = 0.0
       Ch_RDOCusch0(ich) = 0.0

       !N
       Ch_ALGNusch0(ich,:) = 0.0
       Ch_ZOONusch0(ich) = 0.0
       Ch_MBMNusch0(ich) = 0.0
       Ch_LPONusch0(ich) = 0.0
       Ch_RPONusch0(ich) = 0.0
       Ch_LDONusch0(ich) = 0.0
       Ch_RDONusch0(ich) = 0.0
       Ch_NH4usch0(ich)  = 0.0
       Ch_NO3usch0(ich)  = 0.0

       !P
       Ch_ALGPusch0(ich,:) = 0.0 
       Ch_ZOOPusch0(ich) = 0.0
       Ch_MBMPusch0(ich) = 0.0
       Ch_LPOPusch0(ich) = 0.0
       Ch_RPOPusch0(ich) = 0.0
       Ch_LDOPusch0(ich) = 0.0
       Ch_RDOPusch0(ich) = 0.0
       Ch_PO4usch0(ich)  = 0.0
       Ch_PIPAusch0(ich) = 0.0
       Ch_PIPSusch0(ich) = 0.0

       !channel water (Y.Kwon)
       SedCNP_hydro%q_usch0(ich) = 0.0
       !SedCNP_hydro%Wstg_st(ich) = SedCNP_hydro%Wstg_st0(ich)
       channelSed%Susch0(:,ich) = 0.0


 !for debugging !BK20240516
 !open (24051618, file='./debug/chid_upchid.csv', status='unknown')


       do iich = 1, domain%nch    !BK20241102
          if (iich == ich) then
            cycle
          elseif (SedCNP_hydro%dwnstrm_linkID(iich) == chid) then  !BK20231030

             !upchid = SedCNP_hydro%linkID(iupch)                !BK20231030
             iupch = iich

             !channel water (Y.Kwon)
             SedCNP_hydro%q_usch0(ich) = SedCNP_hydro%q_usch0(ich) + SedCNP_hydro%q_dsch(iupch)

             !sed
             channelSed%Susch0(:,ich) = channelSed%Susch0(:,ich) + channelSed%Sdsch(:,iupch)

             !C
             Ch_ALGCusch0(ich,1:nalg) = Ch_ALGCusch0(ich,1:nalg) + Ch_ALGCdsch(iupch,1:nalg) 
             Ch_ZOOCusch0(ich) = Ch_ZOOCusch0(ich) + Ch_ZOOCdsch(iupch)
             Ch_MBMCusch0(ich) = Ch_MBMCusch0(ich) + Ch_MBMCdsch(iupch)
             Ch_DICusch0(ich) = Ch_DICusch0(ich) + Ch_DICdsch(iupch)  !BK20251121
             Ch_LPOCusch0(ich) = Ch_LPOCusch0(ich) + Ch_LPOCdsch(iupch)      
             Ch_RPOCusch0(ich) = Ch_RPOCusch0(ich) + Ch_RPOCdsch(iupch)
             Ch_LDOCusch0(ich) = Ch_LDOCusch0(ich) + Ch_LDOCdsch(iupch) 
             Ch_RDOCusch0(ich) = Ch_RDOCusch0(ich) + Ch_RDOCdsch(iupch)

             !N
             Ch_ALGNusch0(ich,1:nalg) = Ch_ALGNusch0(ich,1:nalg) + Ch_ALGNdsch(iupch,1:nalg)
             Ch_ZOONusch0(ich) = Ch_ZOONusch0(ich) + Ch_ZOONdsch(iupch)
             Ch_MBMNusch0(ich) = Ch_MBMNusch0(ich) + Ch_MBMNdsch(iupch)
             Ch_LPONusch0(ich) = Ch_LPONusch0(ich) + Ch_LPONdsch(iupch)
             Ch_RPONusch0(ich) = Ch_RPONusch0(ich) + Ch_RPONdsch(iupch)
             Ch_LDONusch0(ich) = Ch_LDONusch0(ich) + Ch_LDONdsch(iupch)
             Ch_RDONusch0(ich) = Ch_RDONusch0(ich) + Ch_RDONdsch(iupch)
             Ch_NH4usch0(ich)  = Ch_NH4usch0(ich) + Ch_NH4dsch(iupch)
             Ch_NO3usch0(ich)  = Ch_NO3usch0(ich) + Ch_NO3dsch(iupch)

             !P
             Ch_ALGPusch0(ich,1:nalg) = Ch_ALGPusch0(ich,1:nalg) + Ch_ALGPdsch(iupch,1:nalg)
             Ch_ZOOPusch0(ich) = Ch_ZOOPusch0(ich) + Ch_ZOOPdsch(iupch)
             Ch_MBMPusch0(ich) = Ch_MBMPusch0(ich) + Ch_MBMPdsch(iupch)
             Ch_LPOPusch0(ich) = Ch_LPOPusch0(ich) + Ch_LPOPdsch(iupch)
             Ch_RPOPusch0(ich) = Ch_RPOPusch0(ich) + Ch_RPOPdsch(iupch)
             Ch_LDOPusch0(ich) = Ch_LDOPusch0(ich) + Ch_LDOPdsch(iupch)
             Ch_RDOPusch0(ich) = Ch_RDOPusch0(ich) + Ch_RDOPdsch(iupch)
             Ch_PO4usch0(ich)  = Ch_PO4usch0(ich) + Ch_PO4dsch(iupch)
             Ch_PIPAusch0(ich) = Ch_PIPAusch0(ich) + Ch_PIPAdsch(iupch)
             Ch_PIPSusch0(ich) = Ch_PIPSusch0(ich) + Ch_PIPSdsch(iupch)


 !for debugging !BK20240516
 !write (2451618,*) iupch, chid, upchid, SedCNP_hydro%q_usch0(chid), SedCNP_hydro%q_dsch(upchid), &
 !          channelSed%Susch0(1,chid), channelSed%Sdsch(1,upchid), Ch_LDOCusch0(chid), Ch_LDOCdsch(upchid)


          endif
       enddo
       !----------------------------------------------------------------------

 !close (24051618)

    end subroutine UpdateCNP_Ch

end module module_CNPupdates
