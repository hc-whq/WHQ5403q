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

!#BK
! *** numbering rule for I/O unit number (6 digits) ***
! digit 1: read/write (5 = read, 6 = write)
! digit 2: subject (9 = CNP or water quality)
! digit 3: subject (0 = common, 1 = C, 2 = N, 3 = P)
! digit 4-6: any number to identify individual I/O units

module module_SedCNP_out

  use netcdf      
  !use module_hydro_stop, only:HYDRO_stop
!=====||__WHQ5403q__||=====!
!
  use SedCNP_config,         only: SedCNPmodel   !BK20231016
!
!=====||__WHQ5403q__||=====!
  !use module_SedCNPvariables

  contains

  subroutine write_Sed_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                        :: numSedVars = 8
     character (len=64), dimension(numSedVars) :: SedOut_varNames
     character (len=64), dimension(numSedVars) :: SedOut_longName  ! Long names for each variable.
     character (len=64), dimension(numSedVars) :: SedOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm      ! SEDOUT_DOMAIN filename
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnSed      ! NetCDF file handle
     integer                                   :: iret             ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId            ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     SedOut_varNames(:) = [character(len=64) :: &
                          "Sol","Sow","Solcp","Ssurf","Srdet", &  !1-5
                          "Soler","Soldp","Colsd"]                !6-8
     SedOut_longName(:) = [character(len=128) :: &
                          "Overland sediment storage", &                                !1
                          "Overland water sediment storage", &                          !2
                          "Sediment transport capacity", &                              !3
                          "Overland sediment transport", &                              !4
                          "Sediment detachment by rainfall and irrigation water", &     !5
                          "Sediment erosion", &                                         !6
                          "Sediment deposition", &                                      !7
                          "Concentration of sediment particles in the overland flow"]   !8
     SedOut_units(:) = [character(len=64) :: &
                       "kg","kg","kg s-1","kg s-1","kg s-1", &     !1-5
                       "kg s-1","kg s-1","kg m-3"]                 !6-8       

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     !iret = nf90_create(trim(output_flnm),OR(NF90_CLOBBER,NF90_NETCDF4),ncid=ftn)  !BK20250630
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create SEDOUT NetCDF file.')
     ftnSed = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnSed,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnSed,'x',domain%ix,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define x dimension')
     iret = nf90_def_dim(ftnSed,'y',domain%jx,dimId(3))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define y dimension')
     iret = nf90_def_dim(ftnSed,'grain_size',4,dimId(4))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define grain_size dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnSed,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnSed,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnSed,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numSedVars
        iret = nf90_def_var(ftnSed,trim(SedOut_varNames(iTmp)),nf90_double,(/dimId(4),dimId(2),dimId(3),dimId(1)/),varId)
        call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(SedOut_varNames(iTmp)))

        ! Create variable attributes
        iret = nf90_put_att(ftnSed,varId,'long_name',trim(SedOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(SedOut_varNames(iTmp)))
        iret = nf90_put_att(ftnSed,varId,'units',trim(SedOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(SedOut_varNames(iTmp)))
     enddo
   
     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnSed)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take SEDOUT file out of definition mode') 

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnSed,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnSed,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time') 

     ! Write array out to NetCDF file
     !Sol
     iret = nf90_inq_varid(ftnSed,'Sol',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sol')
     iret = nf90_put_var(ftnSed,varId,overSed%Sol,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sol')

     !Sow
     iret = nf90_inq_varid(ftnSed,'Sow',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sow')
     iret = nf90_put_var(ftnSed,varId,overSed%Sow,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sow')

     !Solcp
     iret = nf90_inq_varid(ftnSed,'Solcp',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Solcp')
     iret = nf90_put_var(ftnSed,varId,overSed%Solcp,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Solcp')

     !Ssurf
     where(overSed%Ssurf < 1.0d-14) overSed%Ssurf = 0.0d0  !BK20260510
     iret = nf90_inq_varid(ftnSed,'Ssurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ssurf')
     iret = nf90_put_var(ftnSed,varId,overSed%Ssurf,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ssurf')

     !Srdet
     iret = nf90_inq_varid(ftnSed,'Srdet',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Srdet')
     iret = nf90_put_var(ftnSed,varId,overSed%Srdet,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))

     ! --- for debugging  BK20240213
     !if(iret .ne. 0) then
     !   write(*,*) 'varID = ', varID
     !   write(*,*) 'ftnSed = ', ftnSed
     !   write(*,*) 'Srdet = ', overSed%Srdet
     !endif
     ! ---

     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Srdet')

     !Soler
     iret = nf90_inq_varid(ftnSed,'Soler',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Soler')
     iret = nf90_put_var(ftnSed,varId,overSed%Soler,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Soler')

     !Soldp
     iret = nf90_inq_varid(ftnSed,'Soldp',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Soldp')
     iret = nf90_put_var(ftnSed,varId,overSed%Soldp,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Soldp')

     !Colsd
     iret = nf90_inq_varid(ftnSed,'Colsd',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Colsd')
     iret = nf90_put_var(ftnSed,varId,overSed%Colsd,(/1,1,1,1/),(/4,domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Colsd')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnSed)

  end subroutine write_Sed_output


  subroutine write_Ch_Sed_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                        :: numChSedVars = 13  !14  !BK20260510 removed Sirri
     character (len=64), dimension(numChSedVars) :: ChSedOut_varNames
     character (len=64), dimension(numChSedVars) :: ChSedOut_longName  ! Long names for each variable.
     character (len=64), dimension(numChSedVars) :: ChSedOut_units     ! Units for each variable.
     character(len=256)                          :: output_flnm        ! CHSEDOUT_DOMAIN filename
     integer                                     :: diagFlag
     integer                                     :: ftn, ftnChSed        ! NetCDF file handle
     integer                                     :: iret               ! NetCDF return statuses
     integer                                     :: iTmp
     integer                                     :: varId              ! Variable ID value created as NetCDF variables are created and populated.
     integer                                     :: dimId(4)           ! Dimension ID values created during NetCDF created.
     real(8)                                     :: dt

     diagFlag = 0    !BK20251031

   !   ChSedOut_varNames(:) = [character(len=64) :: &
   !                          "Sch","Scw","Susch0","Sxpnt","Solch0", &       !1-5  !BK20240701  !BK20250925
   !                          "Sirri","Sinput","Sdsch","Scher","Schdp", &   !6-10
   !                          "Schcp","Cchsd","Sxabs","Sxdis"]             !11-14  !BK20240729
     ChSedOut_varNames(:) = [character(len=64) :: &
                            "Sch","Scw","Susch0","Sxpnt","Solch0", &       !1-5  !BK20240701  !BK20250925
                            "Sinput","Sdsch","Scher","Schdp", &   !6-9  !BK20260510 removed Sirri
                            "Schcp","Cchsd","Sxabs","Sxdis"]             !10-13  !BK20240729
     ChSedOut_longName(:) = [character(len=128) :: &
                            "channel sediment storage", &         !1
                            "channel water sediment storage", &   !2
                            "upstream discharge", &               !3
                            "sed. in-flux from point source", &   !4   !BK20240729
                            "overland to channel", &              !5
                            !"irrigation", &                       !6  !BK20260510 removed Sirri
                            "total sediment input", &             !6
                            "downstream discharge", &             !7
                            "channel erosion", &                  !8
                            "channel deposition", &               !9
                            "channel capacity", &                 !10
                            "channel sediment concentration", &   !11
                            "sed. out-flux by abstraction", &     !12   !BK20240729
                            "sed. in-flux by discharge"]          !13   !BK20240729
     ChSedOut_units(:) = [character(len=64) :: &
                         "kg","kg","kg s-1","kg","kg s-1", &            !1-5   !BK20240729
                         "kg","kg s-1","kg s-1","kg s-1",& !6-9 !BK20260510 removed Sirri
                         "kg","kg m-3","kg","kg"]                       !10-13

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create CHSEDOUT NetCDF file.')
     ftnChSed = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnChSed,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnChSed,'nch',domain%nch,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define nch dimension')
     iret = nf90_def_dim(ftnChSed,'grain_size',4,dimId(3))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define grain_size dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  Y.Kwon20241125
     iret = nf90_def_var(ftnChSed,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnChSed,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnChSed,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numChSedVars
        iret = nf90_def_var(ftnChSed,trim(ChSedOut_varNames(iTmp)),nf90_double,(/dimId(3),dimId(2),dimId(1)/),varId)
        call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(ChSedOut_varNames(iTmp)))

        ! Create variable attributes
        iret = nf90_put_att(ftnChSed,varId,'long_name',trim(ChSedOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(ChSedOut_varNames(iTmp)))
        iret = nf90_put_att(ftnChSed,varId,'units',trim(ChSedOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(ChSedOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnChSed)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take SEDOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon 20241125
     iret = nf90_inq_varid(ftnChSed,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnChSed,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !Sch
     iret = nf90_inq_varid(ftnChSed,'Sch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sch')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Sch,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sch')

     !Scw
     where(channelSed%Scw < 1.0d-14) channelSed%Scw = 0.0d0  !BK20260510
     iret = nf90_inq_varid(ftnChSed,'Scw',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Scw')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Scw,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Scw')

     !Susch0
     iret = nf90_inq_varid(ftnChSed,'Susch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Susch0')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Susch0,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Susch0')

     !Spnt
     !iret = nf90_inq_varid(ftnChSed,'Spnt',varId)
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Spnt')
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Spnt,(/1,1,1/),(/4,domain%nch,1/))
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Spnt')
     
     !Sxpnt  !BK20240701
     iret = nf90_inq_varid(ftnChSed,'Sxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sxpnt')
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Sxpnt,(/1,1/),(/4,domain%nch/))
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Sxpnt,(/1,1,1/),(/4,domain%nch,1/))  !BK20250701
     block        !--- BK20250702
        real(8), dimension(4, domain%nch, 1) :: tmp_arr
        tmp_arr(:,:,1) = channelSed%Sxpnt(:,:)
        iret = nf90_put_var(ftnChSed,varId,tmp_arr,(/1,1,1/),(/4,domain%nch,1/))
     end block    !--- BK20250702
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sxpnt')
     
     !Solch0  !BK20250925
     iret = nf90_inq_varid(ftnChSed,'Solch0',varId)  
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Solch0')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Solch0,(/1,1,1/),(/4,domain%nch,1/))  
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Solch0')

     !Sirri  !BK20260510 removed Sirri
     !iret = nf90_inq_varid(ftnChSed,'Sirri',varId)
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sirri')
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Sirri,(/1,1,1/),(/4,domain%nch,1/))
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sirri')

     !Sinput
     where(channelSed%Sinput < 1.0d-14) channelSed%Sinput = 0.0d0  !BK20260510
     iret = nf90_inq_varid(ftnChSed,'Sinput',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sinput')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Sinput,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sinput')

     !Sdsch
     iret = nf90_inq_varid(ftnChSed,'Sdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sdsch')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Sdsch,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sdsch')

     !Scher
     iret = nf90_inq_varid(ftnChSed,'Scher',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Scher')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Scher,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Scher')

     !Schdp
     where(channelSed%Schdp < 1.0d-14) channelSed%Schdp = 0.0d0  !BK20260510
     iret = nf90_inq_varid(ftnChSed,'Schdp',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Schdp')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Schdp,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Schdp')

     !Schcp
     iret = nf90_inq_varid(ftnChSed,'Schcp',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Schcp')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Schcp,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Schcp')

     !Cchsd
     iret = nf90_inq_varid(ftnChSed,'Cchsd',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Cchsd')
     iret = nf90_put_var(ftnChSed,varId,channelSed%Cchsd,(/1,1,1/),(/4,domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Cchsd')

     !dt = real(SedCNPmodel%SedCNP_timestep)   !BK20240729
     !Sxabs   !BK20240729 
     where(channelSed%Sxabs < 1.0d-14) channelSed%Sxabs = 0.0d0  !BK20260510
     iret = nf90_inq_varid(ftnChSed,'Sxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sxabs')
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Sxabs,(/1,1/),(/4,domain%nch/))
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Sxabs,(/1,1,1/),(/4,domain%nch,1/))  !BK20250701
      block       !--- BK20250702
         real(8), dimension(4, domain%nch, 1) :: tmp_arr
         tmp_arr(:,:,1) = channelSed%Sxabs(:,:)
         iret = nf90_put_var(ftnChSed,varId,tmp_arr,(/1,1,1/),(/4,domain%nch,1/))
      end block   !--- BK20250702
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sxabs')

     !Sxdis   !BK20240729
     iret = nf90_inq_varid(ftnChSed,'Sxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Sxdis')
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Sxdis,(/1,1/),(/4,domain%nch/))
     !iret = nf90_put_var(ftnChSed,varId,channelSed%Sxdis,(/1,1,1/),(/4,domain%nch,1/))  !BK20250701
     block        !--- BK20250702
        real(8), dimension(4, domain%nch, 1) :: tmp_arr
        tmp_arr(:,:,1) = channelSed%Sxdis(:,:)
        iret = nf90_put_var(ftnChSed,varId,tmp_arr,(/1,1,1/),(/4,domain%nch,1/))
     end block    !--- BK20250702
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Sxdis')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnChSed)
     
  end subroutine write_Ch_Sed_output


  subroutine write_So_C_output(output_flnm)
  
     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                         :: numSoC3DVars = 13
     integer, parameter                         :: numSoC4DVars = 13
     integer, parameter                         :: numSoCVars   = 26 !(numWSoC3DVars + numWSoC4DVars)
     character (len=64), dimension(numSoCVars)  :: SoCOut_varNames
     character (len=64), dimension(numSoCVars)  :: SoCOut_longName  ! Long names for each variable.
     character (len=64), dimension(numSoCVars)  :: SoCOut_units     ! Units for each variable.
     character(len=256)                         :: output_flnm       ! CNPOUT_DOMAIN filename
     integer                                    :: diagFlag
     integer                                    :: ftn, ftnSoC      ! NetCDF file handle
     integer                                    :: iret              ! NetCDF return statuses
     integer                                    :: iTmp
     integer                                    :: varId             ! Variable ID value created as NetCDF variables are created and populated.
     integer                                    :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     SoCOut_varNames(:) = [character(len=64) :: &
                           "So_LPOCLITin","So_LPOCEXCin","So_LPOCMANin","So_LDOCgwso0","So_RDOCgwso0",& !1-5
                           "So_LPOCsurf","So_RPOCsurf","So_LDOCsurf","So_RDOCsurf","So_MBMCsurf",&      !6-10
                           "So_LPOCLIT","So_LPOCEXC","So_LPOCMAN","So_LPOCRESin","So_LDOCintf",&        !11-15
                           "So_RDOCintf","So_LDOCperc","So_RDOCperc","So_CO2","So_LPOCRES",&            !16-20
                           "So_RPOC","So_LDOC","So_RDOC","So_MBMC","So_CSC",&                           !21-25
                           "So_CMBError"]                                                               !26
     SoCOut_longName(:) = [character(len=128) :: &
                           "So_LPOCLITin","So_LPOCEXCin","So_LPOCMANin","So_LDOCgwso0","So_RDOCgwso0",& !1-5
                           "So_LPOCsurf","So_RPOCsurf","So_LDOCsurf","So_RDOCsurf","So_MBMCsurf",&      !6-10
                           "So_LPOCLIT","So_LPOCEXC","So_LPOCMAN","So_LPOCRESin","So_LDOCintf",&        !11-15
                           "So_RDOCintf","So_LDOCperc","So_RDOCperc","So_CO2","So_LPOCRES",&            !16-20
                           "So_RPOC","So_LDOC","So_RDOC","So_MBMC","So_CSC",&                           !21-25
                           "So_CMBError"]                                                               !26
     SoCOut_units(:) = [character(len=64) :: &
                           "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1", &  !1-5
                           "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1", &  !6-10
                           "kgC","kgC","kgC","kgC dt-1","kgC dt-1", &                 !11-15
                           "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC",&        !16-20
                           "kgC","kgC","kgC","kgC","kgC dt-1",&                       !21-25
                           "kgC dt-1"]                                                !26

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create SoCOUT NetCDF file.')
     ftnSoC = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnSoC,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnSoC,'x',domain%ix,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define x dimension')
     iret = nf90_def_dim(ftnSoC,'y',domain%jx,dimId(3))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define y dimension')
     iret = nf90_def_dim(ftnSoC,'soil_layers',domain%nsl,dimId(4))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define soil_layers dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnSoC,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnSoC,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnSoC,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numSoCVars
        if(iTmp.le.numSoC3DVars) then
           iret = nf90_def_var(ftnSoC,trim(SoCOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(3),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(SoCOut_varNames(iTmp)))
        else   !numWSoC4DVars
           iret = nf90_def_var(ftnSoC,trim(SoCOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(4),dimId(3),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(SoCOut_varNames(iTmp)))
        endif

        ! Create variable attributes
        iret = nf90_put_att(ftnSoC,varId,'long_name',trim(SoCOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(SoCOut_varNames(iTmp)))
        iret = nf90_put_att(ftnSoC,varId,'units',trim(SoCOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(SoCOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnSoC)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take SoCOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon 20241125
     iret = nf90_inq_varid(ftnSoC,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnSoC,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !So_LPOCLITin
     iret = nf90_inq_varid(ftnSoC,'So_LPOCLITin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOCLITin')
     iret = nf90_put_var(ftnSoC,varId,So_LPOCLITin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOCLITin')

     !So_LPOCEXCin
     iret = nf90_inq_varid(ftnSoC,'So_LPOCEXCin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOCEXCin')
     iret = nf90_put_var(ftnSoC,varId,So_LPOCEXCin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOCEXCin')

     !So_LPOCMANin
     iret = nf90_inq_varid(ftnSoC,'So_LPOCMANin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOCMANin')
     iret = nf90_put_var(ftnSoC,varId,So_LPOCMANin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOCMANin')

     !So_LDOCgwso0
     iret = nf90_inq_varid(ftnSoC,'So_LDOCgwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LDOCgwso0')
     iret = nf90_put_var(ftnSoC,varId,So_LDOCgwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LDOCgwso0')

     !So_RDOCgwso0
     iret = nf90_inq_varid(ftnSoC,'So_RDOCgwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RDOCgwso0')
     iret = nf90_put_var(ftnSoC,varId,So_RDOCgwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RDOCgwso0')

     !So_LPOCsurf
     iret = nf90_inq_varid(ftnSoC,'So_LPOCsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOCsurf')
     iret = nf90_put_var(ftnSoC,varId,So_LPOCsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOCsurf')

     !So_RPOCsurf
     iret = nf90_inq_varid(ftnSoC,'So_RPOCsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RPOCsurf')
     iret = nf90_put_var(ftnSoC,varId,So_RPOCsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RPOCsurf')

     !So_LDOCsurf
     iret = nf90_inq_varid(ftnSoC,'So_LDOCsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LDOCsurf')
     iret = nf90_put_var(ftnSoC,varId,So_LDOCsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LDOCsurf')

     !So_RDOCsurf
     iret = nf90_inq_varid(ftnSoC,'So_RDOCsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RDOCsurf')
     iret = nf90_put_var(ftnSoC,varId,So_RDOCsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RDOCsurf')

     !So_MBMCsurf
     iret = nf90_inq_varid(ftnSoC,'So_MBMCsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_MBMCsurf')
     iret = nf90_put_var(ftnSoC,varId,So_MBMCsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_MBMCsurf')

     !So_LPOCLIT
     iret = nf90_inq_varid(ftnSoC,'So_LPOCLIT',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOCLIT')
     iret = nf90_put_var(ftnSoC,varId,So_LPOCLIT,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOCLIT')

     !So_LPOCEXC
     iret = nf90_inq_varid(ftnSoC,'So_LPOCEXC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOCEXC')
     iret = nf90_put_var(ftnSoC,varId,So_LPOCEXC,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOCEXC')

     !So_LPOCMAN
     iret = nf90_inq_varid(ftnSoC,'So_LPOCMAN',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOCMAN')
     iret = nf90_put_var(ftnSoC,varId,So_LPOCMAN,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOCMAN')

     ! Write 4D variables
     do iTmp = numSoC3DVars + 1, numSoCVars
        iret = nf90_inq_varid(ftnSoC, trim(SoCOut_varNames(iTmp)), varId)
        call nwmCheck(diagFlag, iret, 'ERROR: Unable to find variable ID for var: ' // trim(SoCOut_varNames(iTmp)))
        
        SELECT CASE (trim(SoCOut_varNames(iTmp)))
        CASE ('So_LPOCRESin'); iret = nf90_put_var(ftnSoC, varId, So_LPOCRESin, (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDOCintf');  iret = nf90_put_var(ftnSoC, varId, So_LDOCintf,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDOCintf');  iret = nf90_put_var(ftnSoC, varId, So_RDOCintf,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDOCperc');  iret = nf90_put_var(ftnSoC, varId, So_LDOCperc,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDOCperc');  iret = nf90_put_var(ftnSoC, varId, So_RDOCperc,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_CO2');       iret = nf90_put_var(ftnSoC, varId, So_CO2,       (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LPOCRES');   iret = nf90_put_var(ftnSoC, varId, So_LPOCRES,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RPOC');      iret = nf90_put_var(ftnSoC, varId, So_RPOC,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDOC');      iret = nf90_put_var(ftnSoC, varId, So_LDOC,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDOC');      iret = nf90_put_var(ftnSoC, varId, So_RDOC,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_MBMC');      iret = nf90_put_var(ftnSoC, varId, So_MBMC,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_CSC');       iret = nf90_put_var(ftnSoC, varId, So_CSC,       (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_CMBError');  iret = nf90_put_var(ftnSoC, varId, So_CMBError,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        END SELECT
        call nwmCheck(diagFlag, iret, 'ERROR: Unable to place data into output variable: ' // trim(SoCOut_varNames(iTmp)))
     end do

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnSoC)
     
  end subroutine write_So_C_output 


  subroutine write_Gw_C_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                        :: numGwCVars = 10
     character (len=64), dimension(numGwCVars) :: GwCOut_varNames
     character (len=64), dimension(numGwCVars) :: GwCOut_longName  ! Long names for each variable.
     character (len=64), dimension(numGwCVars) :: GwCOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm      ! CHCOUT_DOMAIN filename
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnGwC      ! NetCDF file handle
     integer                                   :: iret             ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId            ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     GwCOut_varNames(:) = [character(len=64) :: &
            "Gw_LDOCsogw0","Gw_RDOCsogw0","Gw_LDOCgwso","Gw_RDOCgwso","Gw_LDOCgwch", & !1-5
            "Gw_RDOCgwch","Gw_LDOC","Gw_RDOC","Gw_CSC","Gw_CMBError"]                  !6-10
     GwCOut_longName(:) = [character(len=128) :: &
            "Gw_LDOCsogw0","Gw_RDOCsogw0","Gw_LDOCgwso","Gw_RDOCgwso","Gw_LDOCgwch", & !1-5
            "Gw_RDOCgwch","Gw_LDOC","Gw_RDOC","Gw_CSC","Gw_CMBError"]                  !6-10
     GwCOut_units(:) = [character(len=64) :: &
            "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1", &                  !1-5
            "kgC dt-1","kgC","kgC","kgC dt-1","kgC dt-1"]                              !6-10

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create GWCOUT NetCDF file.')
     ftnGwC = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnGwC,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnGwC,'nbasin',domain%nbasin,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define nbasin dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnGwC,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnGwC,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnGwC,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')
     
     do iTmp=1,numGwCVars
        iret = nf90_def_var(ftnGwC,trim(GwCOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(1)/),varId)
        call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(GwCOut_varNames(iTmp)))

        ! Create variable attributes
        iret = nf90_put_att(ftnGwC,varId,'long_name',trim(GwCOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(GwCOut_varNames(iTmp)))
        iret = nf90_put_att(ftnGwC,varId,'units',trim(GwCOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(GwCOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnGwC)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take GWCOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnGwC,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnGwC,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !Gw_LDOCsogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwC,'Gw_LDOCsogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOCsogw0')
     iret = nf90_put_var(ftnGwC,varId,Gw_LDOCsogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOCsogw0')

     !Gw_RDOCsogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwC,'Gw_RDOCsogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOCsogw0')
     iret = nf90_put_var(ftnGwC,varId,Gw_RDOCsogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOCsogw0')

     !Gw_LDOCgwso  !BK20260322
     iret = nf90_inq_varid(ftnGwC,'Gw_LDOCgwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOCgwso')
     iret = nf90_put_var(ftnGwC,varId,Gw_LDOCgwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOCgwso')

     !Gw_RDOCgwso  !BK20260322
     iret = nf90_inq_varid(ftnGwC,'Gw_RDOCgwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOCgwso')
     iret = nf90_put_var(ftnGwC,varId,Gw_RDOCgwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOCgwso')

     !Gw_LDOCgwch
     iret = nf90_inq_varid(ftnGwC,'Gw_LDOCgwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOCgwch')
     iret = nf90_put_var(ftnGwC,varId,Gw_LDOCgwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOCgwch')

     !Gw_RDOCgwch
     iret = nf90_inq_varid(ftnGwC,'Gw_RDOCgwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOCgwch')
     iret = nf90_put_var(ftnGwC,varId,Gw_RDOCgwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOCgwch')

     !Gw_LDOC
     iret = nf90_inq_varid(ftnGwC,'Gw_LDOC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOC')
     iret = nf90_put_var(ftnGwC,varId,Gw_LDOC,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOC')

     !Gw_RDOC
     iret = nf90_inq_varid(ftnGwC,'Gw_RDOC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOC')
     iret = nf90_put_var(ftnGwC,varId,Gw_RDOC,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOC')

     !Gw_CSC
     iret = nf90_inq_varid(ftnGwC,'Gw_CSC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_CSC')
     iret = nf90_put_var(ftnGwC,varId,Gw_CSC,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_CSC')

     !Gw_CMBError
     iret = nf90_inq_varid(ftnGwC,'Gw_CMBError',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_CMBError')
     iret = nf90_put_var(ftnGwC,varId,Gw_CMBError,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_CMBError')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnGwC)
     
  end subroutine write_Gw_C_output


  subroutine write_Ch_C_output(output_flnm)
    
     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     !integer, parameter                        :: numChchVars     = 42   !BK20250329_2
     !integer, parameter                        :: numChbasinVars  = 9    !BK20250329_2
     integer, parameter                        :: numChCVars      = 54  !(numChchVars + numChbasinVars) !BK20251120
     character (len=64), dimension(numChCVars) :: ChCOut_varNames
     character (len=64), dimension(numChCVars) :: ChCOut_longName  ! Long names for each variable. 
     character (len=64), dimension(numChCVars) :: ChCOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm      ! CHCOUT_DOMAIN filename 
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnChC      ! NetCDF file handle
     integer                                   :: iret             ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId            ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     ChCOut_varNames(:) = [character(len=64) :: &
                          "Ch_ALGCusch0_1","Ch_ALGCusch0_2","Ch_ALGCusch0_3","Ch_ZOOCusch0","Ch_LPOCusch0",&    !1-5
                          "Ch_RPOCusch0","Ch_LDOCusch0","Ch_RDOCusch0","Ch_MBMCusch0","Ch_DICusch0",&           !6-10  !BK20251120
                          "Ch_LPOCxpnt","Ch_RPOCxpnt","Ch_LDOCxpnt","Ch_RDOCxpnt","Ch_LPOCxabs",&               !11-15   !BK20240729
                          "Ch_RPOCxabs","Ch_LDOCxabs","Ch_RDOCxabs","Ch_LPOCxdis","Ch_RPOCxdis",&               !16-20   !BK20240729
                          "Ch_LDOCxdis","Ch_RDOCxdis","Ch_ALGCdsch_1","Ch_ALGCdsch_2","Ch_ALGCdsch_3",&         !21-25   !BK20240729
                          "Ch_ZOOCdsch","Ch_LPOCdsch","Ch_RPOCdsch","Ch_LDOCdsch","Ch_RDOCdsch",&               !26-30
                          "Ch_MBMCdsch","Ch_DICdsch","Ch_ALGC_1","Ch_ALGC_2","Ch_ALGC_3",&                      !31-35  !BK20251120
                          "Ch_ZOOC","Ch_LPOC","Ch_RPOC","Ch_LDOC","Ch_RDOC",&                                   !36-40
                          "Ch_MBMC","Ch_DIC","Ch_CSC","Ch_CMBError","Ch_LPOCsurf0",&                            !41-45   !BK20250328
                          "Ch_RPOCsurf0","Ch_LDOCsurf0","Ch_RDOCsurf0","Ch_MBMCsurf0","Ch_LDOCintf0",&          !46-50   !BK20250328
                          "Ch_RDOCintf0","Ch_LDOCgwch0","Ch_RDOCgwch0","Twater"]                                !51-54      !BK20250328
     ChCOut_longName(:) = [character(len=128) :: &
                          "Ch_ALGCusch0_1","Ch_ALGCusch0_2","Ch_ALGCusch0_3","Ch_ZOOCusch0","Ch_LPOCusch0",&    !1-5
                          "Ch_RPOCusch0","Ch_LDOCusch0","Ch_RDOCusch0","Ch_MBMCusch0","Ch_DICusch0",&           !6-10  !BK20251120
                          "Ch_LPOCxpnt","Ch_RPOCxpnt","Ch_LDOCxpnt","Ch_RDOCxpnt","Ch_LPOCxabs",&               !11-15   !BK20240729
                          "Ch_RPOCxabs","Ch_LDOCxabs","Ch_RDOCxabs","Ch_LPOCxdis","Ch_RPOCxdis",&               !16-20   !BK20240729
                          "Ch_LDOCxdis","Ch_RDOCxdis","Ch_ALGCdsch_1","Ch_ALGCdsch_2","Ch_ALGCdsch_3",&         !21-25   !BK20240729
                          "Ch_ZOOCdsch","Ch_LPOCdsch","Ch_RPOCdsch","Ch_LDOCdsch","Ch_RDOCdsch",&               !26-30
                          "Ch_MBMCdsch","Ch_DICdsch","Ch_ALGC_1","Ch_ALGC_2","Ch_ALGC_3",&                      !31-35  !BK20251120
                          "Ch_ZOOC","Ch_LPOC","Ch_RPOC","Ch_LDOC","Ch_RDOC",&                                   !36-40
                          "Ch_MBMC","Ch_DIC","Ch_CSC","Ch_CMBError","Ch_LPOCsurf0",&                            !41-45   !BK20250328
                          "Ch_RPOCsurf0","Ch_LDOCsurf0","Ch_RDOCsurf0","Ch_MBMCsurf0","Ch_LDOCintf0",&          !46-50   !BK20250328
                          "Ch_RDOCintf0","Ch_LDOCgwch0","Ch_RDOCgwch0","Twater"]                                !51-54      !BK20250328
     ChCOut_units(:) = [character(len=64) :: &
                          "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1",&   !1-5
                          "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1",&   !6-10
                          "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1",&   !11-15
                          "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1",&   !16-20
                          "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1",&   !21-25
                          "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1",&   !26-30
                          "kgC dt-1","kgC dt-1","kgC","kgC","kgC",&                  !31-35
                          "kgC","kgC","kgC","kgC","kgC",&                            !36-40
                          "kgC","kgC","kgC dt-1","kgC dt-1","kgC dt-1",&             !41-45
                          "kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1","kgC dt-1",&   !46-50
                          "kgC dt-1","kgC dt-1","kgC dt-1","deg-C"]                  !51-54
  
     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create CHCOUT NetCDF file.')
     ftnChC = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnChC,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnChC,'nch',domain%nch,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define nch dimension')
     !iret = nf90_def_dim(ftnChC,'nbasin',domain%nbasin,dimId(3))              !BK20250329
     !iret = nf90_def_dim(ftnChC,'nch',domain%nch,dimId(3))                     !BK20250329 !BK20250329_2
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to define nbasin dimension')   !BK20250329 !BK20250329_2

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnChC,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnChC,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnChC,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numChCVars
        !if(iTmp.le.numChchVars) then  !BK20250329_2
           iret = nf90_def_var(ftnChC,trim(ChCOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(ChCOut_varNames(iTmp)))
        !else  !numChbasinVars
        !   iret = nf90_def_var(ftnChC,trim(ChCOut_varNames(iTmp)),nf90_double,(/dimId(3),dimId(1)/),varId)  !BK20250329_2
        !   call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(ChCOut_varNames(iTmp)))  !BK20250329_2
        !endif

        ! Create variable attributes
        iret = nf90_put_att(ftnChC,varId,'long_name',trim(ChCOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(ChCOut_varNames(iTmp)))
        iret = nf90_put_att(ftnChC,varId,'units',trim(ChCOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(ChCOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnChC)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take CHCOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnChC,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnChC,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !Ch_ALGCusch0_1
     iret = nf90_inq_varid(ftnChC,'Ch_ALGCusch0_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGCusch0_1')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGCusch0(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGCusch0_1')

     !Ch_ALGCusch0_2
     iret = nf90_inq_varid(ftnChC,'Ch_ALGCusch0_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGCusch0_2')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGCusch0(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGCusch0_2')

     !Ch_ALGCusch0_3
     iret = nf90_inq_varid(ftnChC,'Ch_ALGCusch0_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGCusch0_3')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGCusch0(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGCusch0_3')

     !Ch_ZOOCusch0
     iret = nf90_inq_varid(ftnChC,'Ch_ZOOCusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOOCusch0')
     iret = nf90_put_var(ftnChC,varId,Ch_ZOOCusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOOCusch0')

     !Ch_LPOCusch0
     iret = nf90_inq_varid(ftnChC,'Ch_LPOCusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOCusch0')
     iret = nf90_put_var(ftnChC,varId,Ch_LPOCusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOCusch0')

     !Ch_RPOCusch0
     iret = nf90_inq_varid(ftnChC,'Ch_RPOCusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOCusch0')
     iret = nf90_put_var(ftnChC,varId,Ch_RPOCusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOCusch0')

     !Ch_LDOCusch0
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCusch0')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCusch0')

     !Ch_RDOCusch0
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCusch0')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCusch0')

     !Ch_MBMCusch0
     iret = nf90_inq_varid(ftnChC,'Ch_MBMCusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMCusch0')
     iret = nf90_put_var(ftnChC,varId,Ch_MBMCusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMCusch0')

     !Ch_DICusch0  !BK20251120
     iret = nf90_inq_varid(ftnChC,'Ch_DICusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_DICusch0')
     iret = nf90_put_var(ftnChC,varId,Ch_DICusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_DICusch0')

     !Ch_LPOCxpnt
     iret = nf90_inq_varid(ftnChC,'Ch_LPOCxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOCxpnt')
     iret = nf90_put_var(ftnChC,varId,Ch_LPOCxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOCxpnt')

     !Ch_RPOCxpnt
     iret = nf90_inq_varid(ftnChC,'Ch_RPOCxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOCxpnt')
     iret = nf90_put_var(ftnChC,varId,Ch_RPOCxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOCxpnt')

     !Ch_LDOCxpnt
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCxpnt')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCxpnt')

     !Ch_RDOCxpnt   !BK20241208
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCxpnt')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCxpnt')

     !Ch_LPOCxabs   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_LPOCxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOCxabs')
     iret = nf90_put_var(ftnChC,varId,Ch_LPOCxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOCxabs')

     !Ch_RPOCxabs   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_RPOCxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOCxabs')
     iret = nf90_put_var(ftnChC,varId,Ch_RPOCxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOCxabs')

     !Ch_LDOCxabs   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCxabs')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCxabs')

     !Ch_RDOCxabs   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCxabs')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCxabs')

     !Ch_LPOCxdis   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_LPOCxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOCxdis')
     iret = nf90_put_var(ftnChC,varId,Ch_LPOCxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOCxdis')     

     !Ch_RPOCxdis   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_RPOCxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOCxdis')
     iret = nf90_put_var(ftnChC,varId,Ch_RPOCxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOCxdis')
     
     !Ch_LDOCxdis   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCxdis')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCxdis')
    
     !Ch_RDOCxdis   !BK20240729
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCxdis')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCxdis')

     !Ch_RDOCxpnt
     !iret = nf90_inq_varid(ftnChC,'Ch_RDOCxpnt',varId)
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCxpnt')
     !iret = nf90_put_var(ftnChC,varId,Ch_RDOCxpnt,(/1,1/),(/domain%nch,1/))
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCxpnt')

     !Ch_ALGCdsch_1
     iret = nf90_inq_varid(ftnChC,'Ch_ALGCdsch_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGCdsch_1')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGCdsch(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGCdsch_1')

     !Ch_ALGCdsch_2
     iret = nf90_inq_varid(ftnChC,'Ch_ALGCdsch_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGCdsch_2')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGCdsch(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGCdsch_2')

     !Ch_ALGCdsch_3
     iret = nf90_inq_varid(ftnChC,'Ch_ALGCdsch_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGCdsch_3')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGCdsch(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGCdsch_3')

     !Ch_ZOOCdsch
     iret = nf90_inq_varid(ftnChC,'Ch_ZOOCdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOOCdsch')
     iret = nf90_put_var(ftnChC,varId,Ch_ZOOCdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOOCdsch')

     !Ch_LPOCdsch
     iret = nf90_inq_varid(ftnChC,'Ch_LPOCdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOCdsch')
     iret = nf90_put_var(ftnChC,varId,Ch_LPOCdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOCdsch')

     !Ch_RPOCdsch
     iret = nf90_inq_varid(ftnChC,'Ch_RPOCdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOCdsch')
     iret = nf90_put_var(ftnChC,varId,Ch_RPOCdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOCdsch')

     !Ch_LDOCdsch
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCdsch')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCdsch')

     !Ch_RDOCdsch
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCdsch')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCdsch')

     !Ch_MBMCdsch
     iret = nf90_inq_varid(ftnChC,'Ch_MBMCdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMCdsch')
     iret = nf90_put_var(ftnChC,varId,Ch_MBMCdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMCdsch')

     !Ch_DICdsch  !BK20251120
     iret = nf90_inq_varid(ftnChC,'Ch_DICdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_DICdsch')
     iret = nf90_put_var(ftnChC,varId,Ch_DICdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_DICdsch')

     !Ch_ALGC_1
     iret = nf90_inq_varid(ftnChC,'Ch_ALGC_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGC_1')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGC(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGC_1')

     !Ch_ALGC_2
     iret = nf90_inq_varid(ftnChC,'Ch_ALGC_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGC_2')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGC(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGC_2')

     !Ch_ALGC_3
     iret = nf90_inq_varid(ftnChC,'Ch_ALGC_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGC_3')
     iret = nf90_put_var(ftnChC,varId,Ch_ALGC(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGC_3')

     !Ch_ZOOC
     iret = nf90_inq_varid(ftnChC,'Ch_ZOOC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOOC')
     iret = nf90_put_var(ftnChC,varId,Ch_ZOOC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOOC')

     !Ch_LPOC
     iret = nf90_inq_varid(ftnChC,'Ch_LPOC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOC')
     iret = nf90_put_var(ftnChC,varId,Ch_LPOC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOC')

     !Ch_RPOC
     iret = nf90_inq_varid(ftnChC,'Ch_RPOC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOC')
     iret = nf90_put_var(ftnChC,varId,Ch_RPOC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOC')

     !Ch_LDOC
     iret = nf90_inq_varid(ftnChC,'Ch_LDOC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOC')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOC')

     !Ch_RDOC
     iret = nf90_inq_varid(ftnChC,'Ch_RDOC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOC')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOC')

     !Ch_MBMC
     iret = nf90_inq_varid(ftnChC,'Ch_MBMC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMC')
     iret = nf90_put_var(ftnChC,varId,Ch_MBMC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMC')

     !Ch_DIC  !BK20251120
     iret = nf90_inq_varid(ftnChC,'Ch_DIC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_DIC')
     iret = nf90_put_var(ftnChC,varId,Ch_DIC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_DIC')

     !Ch_CSC
     iret = nf90_inq_varid(ftnChC,'Ch_CSC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_CSC')
     iret = nf90_put_var(ftnChC,varId,Ch_CSC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_CSC')

     !Ch_CMBError
     iret = nf90_inq_varid(ftnChC,'Ch_CMBError',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_CMBError')
     iret = nf90_put_var(ftnChC,varId,Ch_CMBError,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_CMBError')

     !Ch_LPOCsurf0
     iret = nf90_inq_varid(ftnChC,'Ch_LPOCsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOCsurf0')
     iret = nf90_put_var(ftnChC,varId,Ch_LPOCsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOCsurf0')

     !Ch_RPOCsurf0
     iret = nf90_inq_varid(ftnChC,'Ch_RPOCsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOCsurf0')
     iret = nf90_put_var(ftnChC,varId,Ch_RPOCsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOCsurf0')

     !Ch_LDOCsurf0
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCsurf0')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCsurf0')

     !Ch_RDOCsurf0
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCsurf0')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCsurf0')

     !Ch_MBMCsurf0
     iret = nf90_inq_varid(ftnChC,'Ch_MBMCsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMCsurf0')
     iret = nf90_put_var(ftnChC,varId,Ch_MBMCsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMCsurf0')

     !Ch_LDOCintf0
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCintf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCintf0')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCintf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCintf0')

     !Ch_RDOCintf0
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCintf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCintf0')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCintf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCintf0')

     !Ch_LDOCgwch0
     iret = nf90_inq_varid(ftnChC,'Ch_LDOCgwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOCgwch0')
     iret = nf90_put_var(ftnChC,varId,Ch_LDOCgwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOCgwch0')

     !Ch_RDOCgwch0
     iret = nf90_inq_varid(ftnChC,'Ch_RDOCgwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOCgwch0')
     iret = nf90_put_var(ftnChC,varId,Ch_RDOCgwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOCgwch0')

     !Twater
     iret = nf90_inq_varid(ftnChC,'Twater',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Twater')
     iret = nf90_put_var(ftnChC,varId,TWATER,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Twater')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnChC)

  end subroutine write_Ch_C_output


  subroutine write_So_N_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                         :: numSoN3DVars = 19
     integer, parameter                         :: numSoN4DVars = 22
     integer, parameter                         :: numSoNVars   = 41 !(numWSoN3DVars + numWSoN4DVars)     
     character (len=64), dimension(numSoNVars)  :: SoNOut_varNames
     character (len=64), dimension(numSoNVars)  :: SoNOut_longName  ! Long names for each variable.
     character (len=64), dimension(numSoNVars)  :: SoNOut_units     ! Units for each variable.
     character(len=256)                         :: output_flnm       ! CNPOUT_DOMAIN filename
     integer                                    :: diagFlag
     integer                                    :: ftn, ftnSoN      ! NetCDF file handle
     integer                                    :: iret              ! NetCDF return statuses
     integer                                    :: iTmp
     integer                                    :: varId             ! Variable ID value created as NetCDF variables are created and populated.
     integer                                    :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     SoNOut_varNames(:) = [character(len=64) :: &
                           "So_LPONLITin","So_LPONEXCin","So_LPONMANin","So_NH4fxin","So_NO3fxin",&       !1-5
                           "So_LDONgwso0","So_RDONgwso0","So_NH4gwso0","So_NO3gwso0","So_LPONsurf",&  !6-10
                           "So_RPONsurf","So_LDONsurf","So_RDONsurf","So_MBMNsurf","So_NH4surf",&     !11-15
                           "So_NO3surf","So_LPONLIT","So_LPONEXC","So_LPONMAN","So_LPONRESin",&       !16-20
                           "So_LDONintf","So_RDONintf","So_NH4intf","So_NO3intf","So_NH4uptk",&       !21-25
                           "So_NO3uptk","So_LDONperc","So_RDONperc","So_NH4perc","So_NO3perc",&       !26-30
                           "So_NH3VOL","So_DENIT","So_LPONRES","So_RPON","So_LDON",&                  !31-35
                           "So_RDON","So_MBMN","So_NH4","So_NO3","So_NSC",&                           !36-40
                           "So_NMBError"]                                                             !41
     SoNOut_longName(:) = [character(len=128) :: &
                           "So_LPONLITin","So_LPONEXCin","So_LPONMANin","So_NH4fxin","So_NO3fxin",&       !1-5
                           "So_LDONgwso0","So_RDONgwso0","So_NH4gwso0","So_NO3gwso0","So_LPONsurf",&  !6-10
                           "So_RPONsurf","So_LDONsurf","So_RDONsurf","So_MBMNsurf","So_NH4surf",&     !11-15
                           "So_NO3surf","So_LPONLIT","So_LPONEXC","So_LPONMAN","So_LPONRESin",&       !16-20
                           "So_LDONintf","So_RDONintf","So_NH4intf","So_NO3intf","So_NH4uptk",&       !21-25
                           "So_NO3uptk","So_LDONperc","So_RDONperc","So_NH4perc","So_NO3perc",&       !26-30
                           "So_NH3VOL","So_DENIT","So_LPONRES","So_RPON","So_LDON",&                  !31-35
                           "So_RDON","So_MBMN","So_NH4","So_NO3","So_NSC",&                           !36-40
                           "So_NMBError"]
     SoNOut_units(:) = [character(len=64) :: &
                           "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !1-5
                           "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !6-10
                           "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !11-15
                           "kgN dt-1","kgN","kgN","kgN","kgN dt-1",&                 !16-20
                           "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !21-25
                           "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !26-30
                           "kgN dt-1","kgN dt-1","kgN","kgN","kgN",&                 !31-35
                           "kgN","kgN","kgN","kgN","kgN dt-1",&                      !36-40
                           "kgN dt-1"]                                               !41

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create SoNOUT NetCDF file.')
     ftnSoN = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnSoN,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnSoN,'x',domain%ix,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define x dimension')
     iret = nf90_def_dim(ftnSoN,'y',domain%jx,dimId(3))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define y dimension')
     iret = nf90_def_dim(ftnSoN,'soil_layers',domain%nsl,dimId(4))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define soil_layers dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnSoN,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnSoN,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnSoN,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numSoNVars
        if(iTmp.le.numSoN3DVars) then
           iret = nf90_def_var(ftnSoN,trim(SoNOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(3),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(SoNOut_varNames(iTmp)))
        else   !numWSoN4DVars
           iret = nf90_def_var(ftnSoN,trim(SoNOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(4),dimId(3),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(SoNOut_varNames(iTmp)))
        endif

        ! Create variable attributes
        iret = nf90_put_att(ftnSoN,varId,'long_name',trim(SoNOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(SoNOut_varNames(iTmp)))
        iret = nf90_put_att(ftnSoN,varId,'units',trim(SoNOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(SoNOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnSoN)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take SoNOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnSoN,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnSoN,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !So_LPONLITin
     iret = nf90_inq_varid(ftnSoN,'So_LPONLITin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPONLITin')
     iret = nf90_put_var(ftnSoN,varId,So_LPONLITin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPONLITin')

     !So_LPONEXCin
     iret = nf90_inq_varid(ftnSoN,'So_LPONEXCin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPONEXCin')
     iret = nf90_put_var(ftnSoN,varId,So_LPONEXCin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPONEXCin')

     !So_LPONMANin
     iret = nf90_inq_varid(ftnSoN,'So_LPONMANin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPONMANin')
     iret = nf90_put_var(ftnSoN,varId,So_LPONMANin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPONMANin')

     !So_NH4fxin
     iret = nf90_inq_varid(ftnSoN,'So_NH4fxin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_NH4fxin')
     iret = nf90_put_var(ftnSoN,varId,So_NH4fxin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_NH4fxin')

     !So_NO3in
     iret = nf90_inq_varid(ftnSoN,'So_NO3fxin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_NO3fxin')
     iret = nf90_put_var(ftnSoN,varId,So_NO3fxin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_NO3fxin')

     !So_LDONgwso0
     iret = nf90_inq_varid(ftnSoN,'So_LDONgwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LDONgwso0')
     iret = nf90_put_var(ftnSoN,varId,So_LDONgwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LDONgwso0')

     !So_RDONgwso0
     iret = nf90_inq_varid(ftnSoN,'So_RDONgwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RDONgwso0')
     iret = nf90_put_var(ftnSoN,varId,So_RDONgwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RDONgwso0')

     !So_NH4gwso0
     iret = nf90_inq_varid(ftnSoN,'So_NH4gwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_NH4gwso0')
     iret = nf90_put_var(ftnSoN,varId,So_NH4gwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_NH4gwso0')

     !So_NO3gwso0
     iret = nf90_inq_varid(ftnSoN,'So_NO3gwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_NO3gwso0')
     iret = nf90_put_var(ftnSoN,varId,So_NO3gwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_NO3gwso0')

     !So_LPONsurf
     iret = nf90_inq_varid(ftnSoN,'So_LPONsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPONsurf')
     iret = nf90_put_var(ftnSoN,varId,So_LPONsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPONsurf')

     !So_RPONsurf
     iret = nf90_inq_varid(ftnSoN,'So_RPONsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RPONsurf')
     iret = nf90_put_var(ftnSoN,varId,So_RPONsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RPONsurf')

     !So_LDONsurf
     iret = nf90_inq_varid(ftnSoN,'So_LDONsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LDONsurf')
     iret = nf90_put_var(ftnSoN,varId,So_LDONsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LDONsurf')

     !So_RDONsurf
     iret = nf90_inq_varid(ftnSoN,'So_RDONsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RDONsurf')
     iret = nf90_put_var(ftnSoN,varId,So_RDONsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RDONsurf')

     !So_MBMNsurf
     iret = nf90_inq_varid(ftnSoN,'So_MBMNsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_MBMNsurf')
     iret = nf90_put_var(ftnSoN,varId,So_MBMNsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_MBMNsurf')

     !So_NH4surf
     iret = nf90_inq_varid(ftnSoN,'So_NH4surf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_NH4surf')
     iret = nf90_put_var(ftnSoN,varId,So_NH4surf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_NH4surf')

     !So_NO3surf
     iret = nf90_inq_varid(ftnSoN,'So_NO3surf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_NO3surf')
     iret = nf90_put_var(ftnSoN,varId,So_NO3surf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_NO3surf')

     !So_LPONLIT
     iret = nf90_inq_varid(ftnSoN,'So_LPONLIT',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPONLIT')
     iret = nf90_put_var(ftnSoN,varId,So_LPONLIT,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPONLIT')

     !So_LPONEXC
     iret = nf90_inq_varid(ftnSoN,'So_LPONEXC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPONEXC')
     iret = nf90_put_var(ftnSoN,varId,So_LPONEXC,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPONEXC')

     !So_LPONMAN
     iret = nf90_inq_varid(ftnSoN,'So_LPONMAN',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPONMAN')
     iret = nf90_put_var(ftnSoN,varId,So_LPONMAN,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPONMAN')

     do iTmp = numSoN3DVars + 1, numSoNVars
        iret = nf90_inq_varid(ftnSoN, trim(SoNOut_varNames(iTmp)), varId)
        call nwmCheck(diagFlag, iret, 'ERROR: Unable to find variable ID for var: ' // trim(SoNOut_varNames(iTmp)))
        
        SELECT CASE (trim(SoNOut_varNames(iTmp)))
        CASE ('So_LPONRESin'); iret = nf90_put_var(ftnSoN, varId, So_LPONRESin, (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDONintf');  iret = nf90_put_var(ftnSoN, varId, So_LDONintf,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDONintf');  iret = nf90_put_var(ftnSoN, varId, So_RDONintf,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NH4intf');   iret = nf90_put_var(ftnSoN, varId, So_NH4intf,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NO3intf');   iret = nf90_put_var(ftnSoN, varId, So_NO3intf,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NH4uptk');   iret = nf90_put_var(ftnSoN, varId, So_NH4uptk,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NO3uptk');   iret = nf90_put_var(ftnSoN, varId, So_NO3uptk,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDONperc');  iret = nf90_put_var(ftnSoN, varId, So_LDONperc,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDONperc');  iret = nf90_put_var(ftnSoN, varId, So_RDONperc,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NH4perc');   iret = nf90_put_var(ftnSoN, varId, So_NH4perc,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NO3perc');   iret = nf90_put_var(ftnSoN, varId, So_NO3perc,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NH3VOL');    iret = nf90_put_var(ftnSoN, varId, So_NH3VOL,    (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_DENIT');     iret = nf90_put_var(ftnSoN, varId, So_DENIT,     (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LPONRES');   iret = nf90_put_var(ftnSoN, varId, So_LPONRES,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RPON');      iret = nf90_put_var(ftnSoN, varId, So_RPON,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDON');      iret = nf90_put_var(ftnSoN, varId, So_LDON,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDON');      iret = nf90_put_var(ftnSoN, varId, So_RDON,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_MBMN');      iret = nf90_put_var(ftnSoN, varId, So_MBMN,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NH4');       iret = nf90_put_var(ftnSoN, varId, So_NH4,       (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NO3');       iret = nf90_put_var(ftnSoN, varId, So_NO3,       (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NSC');       iret = nf90_put_var(ftnSoN, varId, So_NSC,       (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_NMBError');  iret = nf90_put_var(ftnSoN, varId, So_NMBError,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        END SELECT
        call nwmCheck(diagFlag, iret, 'ERROR: Unable to place data into output variable: ' // trim(SoNOut_varNames(iTmp)))
     end do

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnSoN)
     
  end subroutine write_So_N_output


  subroutine write_Gw_N_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                        :: numGwNVars = 18
     character (len=64), dimension(numGwNVars) :: GwNOut_varNames
     character (len=64), dimension(numGwNVars) :: GwNOut_longName  ! Long names for each variable.
     character (len=64), dimension(numGwNVars) :: GwNOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm      ! CHNOUT_DOMAIN filename
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnGwN      ! NetCDF file handle
     integer                                   :: iret             ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId            ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     GwNOut_varNames(:) = [character(len=64) :: &
            "Gw_LDONsogw0","Gw_RDONsogw0","Gw_NH4sogw0","Gw_NO3sogw0","Gw_LDONgwso", & !1-5
            "Gw_RDONgwso","Gw_NH4gwso","Gw_NO3gwso","Gw_LDONgwch","Gw_RDONgwch", &     !6-10
            "Gw_NH4gwch","Gw_NO3gwch","Gw_LDON","Gw_RDON","Gw_NH4", &                  !11-15
            "Gw_NO3","Gw_NSC","Gw_NMBError"]                                           !16-18
     GwNOut_longName(:) = [character(len=128) :: &
            "Gw_LDONsogw0","Gw_RDONsogw0","Gw_NH4sogw0","Gw_NO3sogw0","Gw_LDONgwso", & !1-5
            "Gw_RDONgwso","Gw_NH4gwso","Gw_NO3gwso","Gw_LDONgwch","Gw_RDONgwch", &     !6-10
            "Gw_NH4gwch","Gw_NO3gwch","Gw_LDON","Gw_RDON","Gw_NH4", &                  !11-15
            "Gw_NO3","Gw_NSC","Gw_NMBError"]                                           !16-18
     GwNOut_units(:) = [character(len=64) :: &
            "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1", &                  !1-5
            "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1", &                  !6-10
            "kgN dt-1","kgN dt-1","kgN","kgN","kgN","kgN", &                           !11-15
            "kgN dt-1","kgN dt-1"]                                                     !16-18

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create GWNOUT NetCDF file.')
     ftnGwN = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnGwN,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnGwN,'nbasin',domain%nbasin,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define nbasin dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnGwN,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnGwN,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnGwN,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numGwNVars
        iret = nf90_def_var(ftnGwN,trim(GwNOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(1)/),varId)
        call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(GwNOut_varNames(iTmp)))

        ! Create variable attributes
        iret = nf90_put_att(ftnGwN,varId,'long_name',trim(GwNOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(GwNOut_varNames(iTmp)))
        iret = nf90_put_att(ftnGwN,varId,'units',trim(GwNOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(GwNOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnGwN)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take GWNOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnGwN,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnGwN,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !Gw_LDONsogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwN,'Gw_LDONsogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDONsogw0')
     iret = nf90_put_var(ftnGwN,varId,Gw_LDONsogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDONsogw0')

     !Gw_RDONsogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwN,'Gw_RDONsogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDONsogw0')
     iret = nf90_put_var(ftnGwN,varId,Gw_RDONsogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDONsogw0')

     !Gw_NH4sogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwN,'Gw_NH4sogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NH4sogw0')
     iret = nf90_put_var(ftnGwN,varId,Gw_NH4sogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NH4sogw0')

     !Gw_NO3sogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwN,'Gw_NO3sogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NO3sogw0')
     iret = nf90_put_var(ftnGwN,varId,Gw_NO3sogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NO3sogw0')

     !Gw_LDONgwso  !BK20260323
     iret = nf90_inq_varid(ftnGwN,'Gw_LDONgwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDONgwso')
     iret = nf90_put_var(ftnGwN,varId,Gw_LDONgwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDONgwso')

     !Gw_RDONgwso  !BK20260323
     iret = nf90_inq_varid(ftnGwN,'Gw_RDONgwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDONgwso')
     iret = nf90_put_var(ftnGwN,varId,Gw_RDONgwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDONgwso')

     !Gw_NH4gwso  !BK20260323
     iret = nf90_inq_varid(ftnGwN,'Gw_NH4gwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NH4gwso')
     iret = nf90_put_var(ftnGwN,varId,Gw_NH4gwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NH4gwso')

     !Gw_NO3gwso  !BK20260323
     iret = nf90_inq_varid(ftnGwN,'Gw_NO3gwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NO3gwso')
     iret = nf90_put_var(ftnGwN,varId,Gw_NO3gwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NO3gwso')

     !Gw_LDONgwch
     iret = nf90_inq_varid(ftnGwN,'Gw_LDONgwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDONgwch')
     iret = nf90_put_var(ftnGwN,varId,Gw_LDONgwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDONgwch')

     !Gw_RDONgwch
     iret = nf90_inq_varid(ftnGwN,'Gw_RDONgwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDONgwch')
     iret = nf90_put_var(ftnGwN,varId,Gw_RDONgwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDONgwch')

     !Gw_NH4gwch
     iret = nf90_inq_varid(ftnGwN,'Gw_NH4gwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NH4gwch')
     iret = nf90_put_var(ftnGwN,varId,Gw_NH4gwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NH4gwch')

     !Gw_NO3gwch
     iret = nf90_inq_varid(ftnGwN,'Gw_NO3gwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NO3gwch')
     iret = nf90_put_var(ftnGwN,varId,Gw_NO3gwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NO3gwch')

     !Gw_LDON
     iret = nf90_inq_varid(ftnGwN,'Gw_LDON',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDON')
     iret = nf90_put_var(ftnGwN,varId,Gw_LDON,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDON')

     !Gw_RDON
     iret = nf90_inq_varid(ftnGwN,'Gw_RDON',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDON')
     iret = nf90_put_var(ftnGwN,varId,Gw_RDON,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDON')

     !Gw_NH4
     iret = nf90_inq_varid(ftnGwN,'Gw_NH4',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NH4')
     iret = nf90_put_var(ftnGwN,varId,Gw_NH4,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NH4')

     !Gw_NO3
     iret = nf90_inq_varid(ftnGwN,'Gw_NO3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NO3')
     iret = nf90_put_var(ftnGwN,varId,Gw_NO3,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NO3')

     !Gw_NSC
     iret = nf90_inq_varid(ftnGwN,'Gw_NSC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NSC')
     iret = nf90_put_var(ftnGwN,varId,Gw_NSC,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NSC')

     !Gw_NMBError
     iret = nf90_inq_varid(ftnGwN,'Gw_NMBError',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_NMBError')
     iret = nf90_put_var(ftnGwN,varId,Gw_NMBError,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_NMBError')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnGwN)

  end subroutine write_Gw_N_output


  subroutine write_Ch_N_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     !integer, parameter                        :: numChchVars     = 55   !BK20250329_2
     !integer, parameter                        :: numChbasinVars  = 15   !BK20250329_2
     integer, parameter                        :: numChNVars      = 70  !(numChchVars + numChbasinVars)
     character (len=64), dimension(numChNVars) :: ChNOut_varNames
     character (len=64), dimension(numChNVars) :: ChNOut_longName  ! Long names for each variable.
     character (len=64), dimension(numChNVars) :: ChNOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm      ! CHNOUT_DOMAIN filename
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnChN      ! NetCDF file handle
     integer                                   :: iret             ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId            ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     ChNOut_varNames(:) = [character(len=64) :: &
                          "Ch_ALGNusch0_1","Ch_ALGNusch0_2","Ch_ALGNusch0_3","Ch_ZOONusch0","Ch_LPONusch0",&  !1-5
                          "Ch_RPONusch0","Ch_LDONusch0","Ch_RDONusch0","Ch_MBMNusch0","Ch_NH4usch0",&         !6-10
                          "Ch_NO3usch0","Ch_LPONxpnt","Ch_RPONxpnt","Ch_LDONxpnt","Ch_RDONxpnt",&             !11-15
                          "Ch_NH4xpnt","Ch_NO3xpnt","Ch_LPONxabs","Ch_RPONxabs","Ch_LDONxabs",&               !16-20   !BK20240729
                          "Ch_RDONxabs","Ch_NH4xabs","Ch_NO3xabs","Ch_LPONxdis","Ch_RPONxdis",&               !21-25   !BK20240729
                          "Ch_LDONxdis","Ch_RDONxdis","Ch_NH4xdis","Ch_NO3xdis","Ch_ALGNdsch_1",&             !26-30   !BK20240729
                          "Ch_ALGNdsch_2","Ch_ALGNdsch_3","Ch_ZOONdsch","Ch_LPONdsch","Ch_RPONdsch",&         !31-35
                          "Ch_LDONdsch","Ch_RDONdsch","Ch_MBMNdsch","Ch_NH4dsch","Ch_NO3dsch",&               !36-40
                          "Ch_NH3VOL","Ch_DENIT","Ch_ALGN_1","Ch_ALGN_2","Ch_ALGN_3",&                        !41-45
                          "Ch_ZOON","Ch_LPON","Ch_RPON","Ch_LDON","Ch_RDON",&                                 !46-50
                          "Ch_MBMN","Ch_NH4","Ch_NO3","Ch_NSC","Ch_NMBError",&                                !51-55
                          "Ch_LPONsurf0","Ch_RPONsurf0","Ch_LDONsurf0","Ch_RDONsurf0","Ch_MBMNsurf0",&        !56-60   !BK20250329
                          "Ch_NH4surf0","Ch_NO3surf0","Ch_LDONintf0","Ch_RDONintf0","Ch_NH4intf0",&           !61-65   !BK20250329
                          "Ch_NO3intf0","Ch_LDONgwch0","Ch_RDONgwch0","Ch_NH4gwch0","Ch_NO3gwch0"]            !66-70   !BK20250329
     ChNOut_longName(:) = [character(len=128) :: &
                          "Ch_ALGNusch0_1","Ch_ALGNusch0_2","Ch_ALGNusch0_3","Ch_ZOONusch0","Ch_LPONusch0",&  !1-5
                          "Ch_RPONusch0","Ch_LDONusch0","Ch_RDONusch0","Ch_MBMNusch0","Ch_NH4usch0",&         !6-10
                          "Ch_NO3usch0","Ch_LPONxpnt","Ch_RPONxpnt","Ch_LDONxpnt","Ch_RDONxpnt",&             !11-15
                          "Ch_NH4xpnt","Ch_NO3xpnt","Ch_LPONxabs","Ch_RPONxabs","Ch_LDONxabs",&               !16-20   !BK20240729
                          "Ch_RDONxabs","Ch_NH4xabs","Ch_NO3xabs","Ch_LPONxdis","Ch_RPONxdis",&               !21-25   !BK20240729
                          "Ch_LDONxdis","Ch_RDONxdis","Ch_NH4xdis","Ch_NO3xdis","Ch_ALGNdsch_1",&             !26-30   !BK20240729
                          "Ch_ALGNdsch_2","Ch_ALGNdsch_3","Ch_ZOONdsch","Ch_LPONdsch","Ch_RPONdsch",&         !31-35
                          "Ch_LDONdsch","Ch_RDONdsch","Ch_MBMNdsch","Ch_NH4dsch","Ch_NO3dsch",&               !36-40
                          "Ch_NH3VOL","Ch_DENIT","Ch_ALGN_1","Ch_ALGN_2","Ch_ALGN_3",&                        !41-45
                          "Ch_ZOON","Ch_LPON","Ch_RPON","Ch_LDON","Ch_RDON",&                                 !46-50
                          "Ch_MBMN","Ch_NH4","Ch_NO3","Ch_NSC","Ch_NMBError",&                                !51-55
                          "Ch_LPONsurf0","Ch_RPONsurf0","Ch_LDONsurf0","Ch_RDONsurf0","Ch_MBMNsurf0",&        !56-60   !BK20250329
                          "Ch_NH4surf0","Ch_NO3surf0","Ch_LDONintf0","Ch_RDONintf0","Ch_NH4intf0",&           !61-65   !BK20250329
                          "Ch_NO3intf0","Ch_LDONgwch0","Ch_RDONgwch0","Ch_NH4gwch0","Ch_NO3gwch0"]            !66-70   !BK20250329
     ChNOut_units(:) = [character(len=64) :: &
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !1-5
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !6-10
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !11-15
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !16-20
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !21-25
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !26-30
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !31-35
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !36-40
                          "kgN dt-1","kgN dt-1","kgN","kgN","kgN",&                 !41-45
                          "kgN","kgN","kgN","kgN","kgN",&                           !46-50
                          "kgN","kgN","kgN","kgN dt-1","kgN dt-1",&                 !51-55
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !56-60
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1",&  !61-65
                          "kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1","kgN dt-1"]   !66-70

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create CHNOUT NetCDF file.')
     ftnChN = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnChN,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnChN,'nch',domain%nch,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define nch dimension')
     !iret = nf90_def_dim(ftnChN,'nbasin',domain%nbasin,dimId(3))              !BK20250329
     !iret = nf90_def_dim(ftnChN,'nch',domain%nch,dimId(3))                     !BK20250329  !BK20250329_2
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to define nbasin dimension')   !BK20250329  !BK20250329_2

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnChN,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnChN,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnChN,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numChNVars
        !if(iTmp.le.numChchVars) then   !BK20250329_2
           iret = nf90_def_var(ftnChN,trim(ChNOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(ChNOut_varNames(iTmp)))
        !else !numChbasinVars
        !   iret = nf90_def_var(ftnChN,trim(ChNOut_varNames(iTmp)),nf90_double,(/dimId(3),dimId(1)/),varId)  !BK20250329_2
        !   call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(ChNOut_varNames(iTmp)))  !BK20250329_2
        !endif

        ! Create variable attributes
        iret = nf90_put_att(ftnChN,varId,'long_name',trim(ChNOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(ChNOut_varNames(iTmp)))
        iret = nf90_put_att(ftnChN,varId,'units',trim(ChNOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(ChNOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnChN)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take CHNOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnChN,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnChN,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !Ch_ALGNusch0_1
     iret = nf90_inq_varid(ftnChN,'Ch_ALGNusch0_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGNusch0_1')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGNusch0(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGNusch0_1')

     !Ch_ALGNusch0_2
     iret = nf90_inq_varid(ftnChN,'Ch_ALGNusch0_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGNusch0_2')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGNusch0(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGNusch0_2')

     !Ch_ALGNusch0_3
     iret = nf90_inq_varid(ftnChN,'Ch_ALGNusch0_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGNusch0_3')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGNusch0(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGNusch0_3')

     !Ch_ZOONusch0
     iret = nf90_inq_varid(ftnChN,'Ch_ZOONusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOONusch0')
     iret = nf90_put_var(ftnChN,varId,Ch_ZOONusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOONusch0')

     !Ch_LPONusch0
     iret = nf90_inq_varid(ftnChN,'Ch_LPONusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPONusch0')
     iret = nf90_put_var(ftnChN,varId,Ch_LPONusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPONusch0')

     !Ch_RPONusch0
     iret = nf90_inq_varid(ftnChN,'Ch_RPONusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPONusch0')
     iret = nf90_put_var(ftnChN,varId,Ch_RPONusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPONusch0')

     !Ch_LDONusch0
     iret = nf90_inq_varid(ftnChN,'Ch_LDONusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONusch0')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONusch0')

     !Ch_RDONusch0
     iret = nf90_inq_varid(ftnChN,'Ch_RDONusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONusch0')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONusch0')

     !Ch_MBMNusch0
     iret = nf90_inq_varid(ftnChN,'Ch_MBMNusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMNusch0')
     iret = nf90_put_var(ftnChN,varId,Ch_MBMNusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMNusch0')

     !Ch_NH4usch0
     iret = nf90_inq_varid(ftnChN,'Ch_NH4usch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4usch0')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4usch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4usch0')

     !Ch_NO3usch0
     iret = nf90_inq_varid(ftnChN,'Ch_NO3usch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3usch0')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3usch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3usch0')

     !Ch_LPONxpnt
     iret = nf90_inq_varid(ftnChN,'Ch_LPONxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPONxpnt')
     iret = nf90_put_var(ftnChN,varId,Ch_LPONxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPONxpnt')

     !Ch_RPONxpnt
     iret = nf90_inq_varid(ftnChN,'Ch_RPONxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPONxpnt')
     iret = nf90_put_var(ftnChN,varId,Ch_RPONxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPONxpnt')

     !Ch_LDONxpnt
     iret = nf90_inq_varid(ftnChN,'Ch_LDONxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONxpnt')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONxpnt')

     !Ch_RDONxpnt   !BK20240729 commented out (duplicated)
     !iret = nf90_inq_varid(ftnChN,'Ch_LDONxpnt',varId)
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONxpnt')
     !iret = nf90_put_var(ftnChN,varId,Ch_LDONxpnt,(/1,1/),(/domain%nch,1/))
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONxpnt')

     !Ch_RDONxpnt
     iret = nf90_inq_varid(ftnChN,'Ch_RDONxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONxpnt')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONxpnt')

     !Ch_NH4xpnt
     iret = nf90_inq_varid(ftnChN,'Ch_NH4xpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4xpnt')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4xpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4xpnt')

     !Ch_NO3xpnt
     iret = nf90_inq_varid(ftnChN,'Ch_NO3xpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3xpnt')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3xpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3xpnt')

     !Ch_LPONxabs   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_LPONxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPONxabs')
     iret = nf90_put_var(ftnChN,varId,Ch_LPONxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPONxabs')

     !Ch_RPONxabs   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_RPONxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPONxabs')
     iret = nf90_put_var(ftnChN,varId,Ch_RPONxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPONxabs')

     !Ch_LDONxabs   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_LDONxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONxabs')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONxabs')

     !Ch_RDONxabs   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_RDONxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONxabs')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONxabs')

     !Ch_NH4xabs   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_NH4xabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4xabs')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4xabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4xabs')

     !Ch_NO3xabs   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_NO3xabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3xabs')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3xabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3xabs')

     !Ch_LPONxdis   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_LPONxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPONxdis')
     iret = nf90_put_var(ftnChN,varId,Ch_LPONxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPONxdis')

     !Ch_RPONxdis   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_RPONxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPONxdis')
     iret = nf90_put_var(ftnChN,varId,Ch_RPONxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPONxdis')

     !Ch_LDONxdis   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_LDONxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONxdis')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONxdis')

     !Ch_RDONxdis   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_RDONxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONxdis')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONxdis')

     !Ch_NH4xdis   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_NH4xdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4xdis')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4xdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4xdis')

     !Ch_NO3xdis   !BK20240729
     iret = nf90_inq_varid(ftnChN,'Ch_NO3xdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3xdis')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3xdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3xdis')

     !Ch_ALGNdsch_1
     iret = nf90_inq_varid(ftnChN,'Ch_ALGNdsch_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGNdsch_1')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGNdsch(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGNdsch_1')

     !Ch_ALGNdsch_2
     iret = nf90_inq_varid(ftnChN,'Ch_ALGNdsch_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGNdsch_2')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGNdsch(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGNdsch_2')

     !Ch_ALGNdsch_3
     iret = nf90_inq_varid(ftnChN,'Ch_ALGNdsch_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGNdsch_3')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGNdsch(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGNdsch_3')

     !Ch_ZOONdsch
     iret = nf90_inq_varid(ftnChN,'Ch_ZOONdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOONdsch')
     iret = nf90_put_var(ftnChN,varId,Ch_ZOONdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOONdsch')

     !Ch_LPONdsch
     iret = nf90_inq_varid(ftnChN,'Ch_LPONdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPONdsch')
     iret = nf90_put_var(ftnChN,varId,Ch_LPONdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPONdsch')

     !Ch_RPONdsch
     iret = nf90_inq_varid(ftnChN,'Ch_RPONdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPONdsch')
     iret = nf90_put_var(ftnChN,varId,Ch_RPONdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPONdsch')

     !Ch_LDONdsch
     iret = nf90_inq_varid(ftnChN,'Ch_LDONdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONdsch')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONdsch')

     !Ch_RDONdsch
     iret = nf90_inq_varid(ftnChN,'Ch_RDONdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONdsch')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONdsch')

     !Ch_MBMNdsch
     iret = nf90_inq_varid(ftnChN,'Ch_MBMNdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMNdsch')
     iret = nf90_put_var(ftnChN,varId,Ch_MBMNdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMNdsch')

     !Ch_NH4dsch
     iret = nf90_inq_varid(ftnChN,'Ch_NH4dsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4dsch')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4dsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4dsch')

     !Ch_NO3dsch
     iret = nf90_inq_varid(ftnChN,'Ch_NO3dsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3dsch')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3dsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3dsch')

     !Ch_NH3VOL
     iret = nf90_inq_varid(ftnChN,'Ch_NH3VOL',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH3VOL')
     iret = nf90_put_var(ftnChN,varId,Ch_NH3VOL,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH3VOL')

     !Ch_DENIT
     iret = nf90_inq_varid(ftnChN,'Ch_DENIT',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_DENIT')
     iret = nf90_put_var(ftnChN,varId,Ch_DENIT,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_DENIT')

     !Ch_ALGN_1
     iret = nf90_inq_varid(ftnChN,'Ch_ALGN_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGN_1')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGN(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGN_1')

     !Ch_ALGN_2
     iret = nf90_inq_varid(ftnChN,'Ch_ALGN_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGN_2')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGN(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGN_2')

     !Ch_ALGN_3
     iret = nf90_inq_varid(ftnChN,'Ch_ALGN_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGN_3')
     iret = nf90_put_var(ftnChN,varId,Ch_ALGN(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGN_3')

     !Ch_ZOON
     iret = nf90_inq_varid(ftnChN,'Ch_ZOON',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOON')
     iret = nf90_put_var(ftnChN,varId,Ch_ZOON,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOON')

     !Ch_LPON
     iret = nf90_inq_varid(ftnChN,'Ch_LPON',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPON')
     iret = nf90_put_var(ftnChN,varId,Ch_LPON,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPON')

     !Ch_RPON
     iret = nf90_inq_varid(ftnChN,'Ch_RPON',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPON')
     iret = nf90_put_var(ftnChN,varId,Ch_RPON,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPON')

     !Ch_LDON
     iret = nf90_inq_varid(ftnChN,'Ch_LDON',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDON')
     iret = nf90_put_var(ftnChN,varId,Ch_LDON,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDON')

     !Ch_RDON
     iret = nf90_inq_varid(ftnChN,'Ch_RDON',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDON')
     iret = nf90_put_var(ftnChN,varId,Ch_RDON,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDON')

     !Ch_MBMN
     iret = nf90_inq_varid(ftnChN,'Ch_MBMN',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMN')
     iret = nf90_put_var(ftnChN,varId,Ch_MBMN,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMN')

     !Ch_NH4
     iret = nf90_inq_varid(ftnChN,'Ch_NH4',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4')

     !Ch_NO3
     iret = nf90_inq_varid(ftnChN,'Ch_NO3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3')

     !Ch_NSC
     iret = nf90_inq_varid(ftnChN,'Ch_NSC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NSC')
     iret = nf90_put_var(ftnChN,varId,Ch_NSC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NSC')

     !Ch_NMBError
     iret = nf90_inq_varid(ftnChN,'Ch_NMBError',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NMBError')
     iret = nf90_put_var(ftnChN,varId,Ch_NMBError,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NMBError')

     !Ch_LPONsurf0
     iret = nf90_inq_varid(ftnChN,'Ch_LPONsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPONsurf0')
     iret = nf90_put_var(ftnChN,varId,Ch_LPONsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPONsurf0')

     !Ch_RPONsurf0
     iret = nf90_inq_varid(ftnChN,'Ch_RPONsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPONsurf0')
     iret = nf90_put_var(ftnChN,varId,Ch_RPONsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPONsurf0')

     !Ch_LDONsurf0
     iret = nf90_inq_varid(ftnChN,'Ch_LDONsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONsurf0')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONsurf0')

     !Ch_RDONsurf0
     iret = nf90_inq_varid(ftnChN,'Ch_RDONsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONsurf0')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONsurf0')

     !Ch_MBMNsurf0
     iret = nf90_inq_varid(ftnChN,'Ch_MBMNsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMNsurf0')
     iret = nf90_put_var(ftnChN,varId,Ch_MBMNsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMNsurf0')

     !Ch_NH4surf0
     iret = nf90_inq_varid(ftnChN,'Ch_NH4surf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4surf0')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4surf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4surf0')

     !Ch_NO3surf0
     iret = nf90_inq_varid(ftnChN,'Ch_NO3surf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3surf0')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3surf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3surf0')

     !Ch_LDONintf0
     iret = nf90_inq_varid(ftnChN,'Ch_LDONintf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONintf0')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONintf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONintf0')

     !Ch_RDONintf0
     iret = nf90_inq_varid(ftnChN,'Ch_RDONintf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONintf0')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONintf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONintf0')

     !Ch_NH4intf0
     iret = nf90_inq_varid(ftnChN,'Ch_NH4intf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4intf0')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4intf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4intf0')

     !Ch_NO3intf0
     iret = nf90_inq_varid(ftnChN,'Ch_NO3intf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3intf0')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3intf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3intf0')

     !Ch_LDONgwch0
     iret = nf90_inq_varid(ftnChN,'Ch_LDONgwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDONgwch0')
     iret = nf90_put_var(ftnChN,varId,Ch_LDONgwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDONgwch0')

     !Ch_RDONgwch0
     iret = nf90_inq_varid(ftnChN,'Ch_RDONgwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDONgwch0')
     iret = nf90_put_var(ftnChN,varId,Ch_RDONgwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDONgwch0')

     !Ch_NH4gwch0
     iret = nf90_inq_varid(ftnChN,'Ch_NH4gwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NH4gwch0')
     iret = nf90_put_var(ftnChN,varId,Ch_NH4gwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NH4gwch0')

     !Ch_NO3gwch0
     iret = nf90_inq_varid(ftnChN,'Ch_NO3gwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_NO3gwch0')
     iret = nf90_put_var(ftnChN,varId,Ch_NO3gwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_NO3gwch0')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnChN)
     
  end subroutine write_Ch_N_output


  subroutine write_So_P_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                        :: numSoP3DVars = 18
     integer, parameter                        :: numSoP4DVars = 18
     integer, parameter                        :: numSoPVars   = 36 !(numWSoP3DVars + numWSoP4DVars)
     character(len=64), dimension(numSoPVars)  :: SoPOut_varNames
     character(len=64), dimension(numSoPVars)  :: SoPOut_longName  ! Long names for each variable.
     character(len=64), dimension(numSoPVars)  :: SoPOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm       ! CNPOUT_DOMAIN filename
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnSoP      ! NetCDF file handle
     integer                                   :: iret              ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId             ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     SoPOut_varNames(:) = [character(len=64) :: &
                           "So_LPOPLITin","So_LPOPEXCin","So_LPOPMANin","So_PO4fxin","So_LDOPgwso0",&  !1-5
                           "So_RDOPgwso0","So_PO4gwso0","So_LPOPsurf","So_RPOPsurf","So_LDOPsurf",&  !6-10
                           "So_RDOPsurf","So_MBMPsurf","So_PO4surf","So_PIPAsurf","So_PIPSsurf",&    !11-15
                           "So_LPOPLIT","So_LPOPEXC","So_LPOPMAN","So_LPOPRESin","So_LDOPintf",&     !16-20
                           "So_RDOPintf","So_PO4intf","So_PO4uptk","So_LDOPperc","So_RDOPperc",&     !21-25
                           "So_PO4perc","So_LPOPRES","So_RPOP","So_LDOP","So_RDOP",&                 !26-30
                           "So_MBMP","So_PO4","So_PIPA","So_PIPS","So_PSC",&                         !31-35
                           "So_PMBError"]                                                            !36
     SoPOut_longName(:) = [character(len=128) :: &
                           "So_LPOPLITin","So_LPOPEXCin","So_LPOPMANin","So_PO4fxin","So_LDOPgwso0",&  !1-5
                           "So_RDOPgwso0","So_PO4gwso0","So_LPOPsurf","So_RPOPsurf","So_LDOPsurf",&  !6-10
                           "So_RDOPsurf","So_MBMPsurf","So_PO4surf","So_PIPAsurf","So_PIPSsurf",&    !11-15
                           "So_LPOPLIT","So_LPOPEXC","So_LPOPMAN","So_LPOPRESin","So_LDOPintf",&     !16-20
                           "So_RDOPintf","So_PO4intf","So_PO4uptk","So_LDOPperc","So_RDOPperc",&     !21-25
                           "So_PO4perc","So_LPOPRES","So_RPOP","So_LDOP","So_RDOP",&                 !26-30
                           "So_MBMP","So_PO4","So_PIPA","So_PIPS","So_PSC",&                         !31-35
                           "So_PMBError"]                                                            !36
     SoPOut_units(:) = [character(len=64) :: &
                           "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !1-5
                           "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !6-10
                           "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !11-15
                           "kgP","kgP","kgP","kgP dt-1","kgP dt-1",&                 !16-20
                           "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !21-25
                           "kgP dt-1","kgP","kgP","kgP","kgP",&                      !26-30
                           "kgP","kgP","kgP","kgP","kgP dt-1",&                      !31-35
                           "kgP dt-1"]                                               !36

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create SoPOUT NetCDF file.')
     ftnSoP = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnSoP,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnSoP,'x',domain%ix,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define x dimension')
     iret = nf90_def_dim(ftnSoP,'y',domain%jx,dimId(3))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define y dimension')
     iret = nf90_def_dim(ftnSoP,'soil_layers',domain%nsl,dimId(4))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define soil_layers dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnSoP,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnSoP,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnSoP,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numSoPVars
        if(iTmp.le.numSoP3DVars) then
           iret = nf90_def_var(ftnSoP,trim(SoPOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(3),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(SoPOut_varNames(iTmp)))
        else   !numWSoP4DVars
           iret = nf90_def_var(ftnSoP,trim(SoPOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(4),dimId(3),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(SoPOut_varNames(iTmp)))
        endif

        ! Create variable attributes
        iret = nf90_put_att(ftnSoP,varId,'long_name',trim(SoPOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(SoPOut_varNames(iTmp)))
        iret = nf90_put_att(ftnSoP,varId,'units',trim(SoPOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(SoPOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnSoP)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take SoPOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnSoP,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnSoP,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !So_LPOPLITin
     iret = nf90_inq_varid(ftnSoP,'So_LPOPLITin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOPLITin')
     iret = nf90_put_var(ftnSoP,varId,So_LPOPLITin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOPLITin')

     !So_LPOPEXCin
     iret = nf90_inq_varid(ftnSoP,'So_LPOPEXCin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOPEXCin')
     iret = nf90_put_var(ftnSoP,varId,So_LPOPEXCin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOPEXCin')

     !So_LPOPMANin
     iret = nf90_inq_varid(ftnSoP,'So_LPOPMANin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOPMANin')
     iret = nf90_put_var(ftnSoP,varId,So_LPOPMANin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOPMANin')

     !So_PO4in
     iret = nf90_inq_varid(ftnSoP,'So_PO4fxin',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_PO4fxin')
     iret = nf90_put_var(ftnSoP,varId,So_PO4fxin,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_PO4fxin')

     !So_LDOPgwso0
     iret = nf90_inq_varid(ftnSoP,'So_LDOPgwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LDOPgwso0')
     iret = nf90_put_var(ftnSoP,varId,So_LDOPgwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LDOPgwso0')

     !So_RDOPgwso0
     iret = nf90_inq_varid(ftnSoP,'So_RDOPgwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RDOPgwso0')
     iret = nf90_put_var(ftnSoP,varId,So_RDOPgwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RDOPgwso0')

     !So_PO4gwso0
     iret = nf90_inq_varid(ftnSoP,'So_PO4gwso0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_PO4gwso0')
     iret = nf90_put_var(ftnSoP,varId,So_PO4gwso0,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_PO4gwso0')

     !So_LPOPsurf
     iret = nf90_inq_varid(ftnSoP,'So_LPOPsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOPsurf')
     iret = nf90_put_var(ftnSoP,varId,So_LPOPsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOPsurf')

     !So_RPOPsurf
     iret = nf90_inq_varid(ftnSoP,'So_RPOPsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RPOPsurf')
     iret = nf90_put_var(ftnSoP,varId,So_RPOPsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RPOPsurf')

     !So_LDOPsurf
     iret = nf90_inq_varid(ftnSoP,'So_LDOPsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LDOPsurf')
     iret = nf90_put_var(ftnSoP,varId,So_LDOPsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LDOPsurf')

     !So_RDOPsurf
     iret = nf90_inq_varid(ftnSoP,'So_RDOPsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_RDOPsurf')
     iret = nf90_put_var(ftnSoP,varId,So_RDOPsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_RDOPsurf')

     !So_MBMPsurf
     iret = nf90_inq_varid(ftnSoP,'So_MBMPsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_MBMPsurf')
     iret = nf90_put_var(ftnSoP,varId,So_MBMPsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_MBMPsurf')

     !So_PO4surf
     iret = nf90_inq_varid(ftnSoP,'So_PO4surf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_PO4surf')
     iret = nf90_put_var(ftnSoP,varId,So_PO4surf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_PO4surf')

     !So_PIPAsurf
     iret = nf90_inq_varid(ftnSoP,'So_PIPAsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_PIPAsurf')
     iret = nf90_put_var(ftnSoP,varId,So_PIPAsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_PIPAsurf')

     !So_PIPSsurf
     iret = nf90_inq_varid(ftnSoP,'So_PIPSsurf',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_PIPSsurf')
     iret = nf90_put_var(ftnSoP,varId,So_PIPSsurf,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_PIPSsurf')

     !So_LPOPLIT
     iret = nf90_inq_varid(ftnSoP,'So_LPOPLIT',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOPLIT')
     iret = nf90_put_var(ftnSoP,varId,So_LPOPLIT,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOPLIT')

     !So_LPOPEXC
     iret = nf90_inq_varid(ftnSoP,'So_LPOPEXC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOPEXC')
     iret = nf90_put_var(ftnSoP,varId,So_LPOPEXC,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOPEXC')

     !So_LPOPMAN
     iret = nf90_inq_varid(ftnSoP,'So_LPOPMAN',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: So_LPOPMAN')
     iret = nf90_put_var(ftnSoP,varId,So_LPOPMAN,(/1,1,1/),(/domain%ix,domain%jx,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: So_LPOPMAN')

     do iTmp = numSoP3DVars + 1, numSoPVars
        iret = nf90_inq_varid(ftnSoP, trim(SoPOut_varNames(iTmp)), varId)
        call nwmCheck(diagFlag, iret, 'ERROR: Unable to find variable ID for var: ' // trim(SoPOut_varNames(iTmp)))
        
        SELECT CASE (trim(SoPOut_varNames(iTmp)))  !BK20251103 changed the order of the soil layers index
        CASE ('So_LPOPRESin'); iret = nf90_put_var(ftnSoP, varId, So_LPOPRESin, (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDOPintf');  iret = nf90_put_var(ftnSoP, varId, So_LDOPintf,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDOPintf');  iret = nf90_put_var(ftnSoP, varId, So_RDOPintf,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PO4intf');   iret = nf90_put_var(ftnSoP, varId, So_PO4intf,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PO4uptk');   iret = nf90_put_var(ftnSoP, varId, So_PO4uptk,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDOPperc');  iret = nf90_put_var(ftnSoP, varId, So_LDOPperc,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDOPperc');  iret = nf90_put_var(ftnSoP, varId, So_RDOPperc,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PO4perc');   iret = nf90_put_var(ftnSoP, varId, So_PO4perc,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LPOPRES');   iret = nf90_put_var(ftnSoP, varId, So_LPOPRES,   (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RPOP');      iret = nf90_put_var(ftnSoP, varId, So_RPOP,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_LDOP');      iret = nf90_put_var(ftnSoP, varId, So_LDOP,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_RDOP');      iret = nf90_put_var(ftnSoP, varId, So_RDOP,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_MBMP');      iret = nf90_put_var(ftnSoP, varId, So_MBMP,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PO4');       iret = nf90_put_var(ftnSoP, varId, So_PO4,       (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PIPA');      iret = nf90_put_var(ftnSoP, varId, So_PIPA,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PIPS');      iret = nf90_put_var(ftnSoP, varId, So_PIPS,      (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PSC');       iret = nf90_put_var(ftnSoP, varId, So_PSC,       (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        CASE ('So_PMBError');  iret = nf90_put_var(ftnSoP, varId, So_PMBError,  (/1,1,1,1/), (/domain%ix,domain%nsl,domain%jx,1/))
        END SELECT
        call nwmCheck(diagFlag, iret, 'ERROR: Unable to place data into output variable: ' // trim(SoPOut_varNames(iTmp)))
     end do

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnSoP)

  end subroutine write_So_P_output


  subroutine write_Gw_P_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     integer, parameter                        :: numGwPVars = 14  !BK20250809
     character (len=64), dimension(numGwPVars) :: GwPOut_varNames
     character (len=64), dimension(numGwPVars) :: GwPOut_longName  ! Long names for each variable.
     character (len=64), dimension(numGwPVars) :: GwPOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm      ! CHPOUT_DOMAIN filename
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnGwP      ! NetCDF file handle
     integer                                   :: iret             ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId            ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     GwPOut_varNames(:) = [character(len=64) :: &
            "Gw_LDOPsogw0","Gw_RDOPsogw0","Gw_PO4sogw0","Gw_LDOPgwso","Gw_RDOPgwso", & !1-5
            "Gw_PO4gwso","Gw_LDOPgwch","Gw_RDOPgwch","Gw_PO4gwch","Gw_LDOP", &         !6-10
            "Gw_RDOP","Gw_PO4","Gw_PSC","Gw_PMBError"]                                 !11-14
     GwPOut_longName(:) = [character(len=128) :: &
            "Gw_LDOPsogw0","Gw_RDOPsogw0","Gw_PO4sogw0","Gw_LDOPgwso","Gw_RDOPgwso", & !1-5
            "Gw_PO4gwso","Gw_LDOPgwch","Gw_RDOPgwch","Gw_PO4gwch","Gw_LDOP", &         !6-10
            "Gw_RDOP","Gw_PO4","Gw_PSC","Gw_PMBError"]                                 !11-14
     GwPOut_units(:) = [character(len=64) :: &
            "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&                   !1-5
            "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP", &                       !6-10
            "kgP","kgP","kgP dt-1","kgP dt-1"]                                         !11-14

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create GWPOUT NetCDF file.')
     ftnGwP = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnGwP,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnGwP,'nbasin',domain%nbasin,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define nbasin dimension')

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnGwP,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnGwP,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnGwP,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numGwPVars
        iret = nf90_def_var(ftnGwP,trim(GwPOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(1)/),varId)
        call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(GwPOut_varNames(iTmp)))

        ! Create variable attributes
        iret = nf90_put_att(ftnGwP,varId,'long_name',trim(GwPOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(GwPOut_varNames(iTmp)))
        iret = nf90_put_att(ftnGwP,varId,'units',trim(GwPOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(GwPOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnGwP)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take GWPOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnGwP,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnGwP,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !Gw_LDOPsogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwP,'Gw_LDOPsogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOPsogw0')
     iret = nf90_put_var(ftnGwP,varId,Gw_LDOPsogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOPsogw0')

     !Gw_RDOPsogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwP,'Gw_RDOPsogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOPsogw0')
     iret = nf90_put_var(ftnGwP,varId,Gw_RDOPsogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOPsogw0')

     !Gw_PO4sogw0  !BK20250310
     iret = nf90_inq_varid(ftnGwP,'Gw_PO4sogw0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_PO4sogw0')
     iret = nf90_put_var(ftnGwP,varId,Gw_PO4sogw0,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_PO4sogw0')

     !Gw_LDOPgwso  !BK20260323
     iret = nf90_inq_varid(ftnGwP,'Gw_LDOPgwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOPgwso')
     iret = nf90_put_var(ftnGwP,varId,Gw_LDOPgwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOPgwso')

     !Gw_RDOPgwso  !BK20260323
     iret = nf90_inq_varid(ftnGwP,'Gw_RDOPgwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOPgwso')
     iret = nf90_put_var(ftnGwP,varId,Gw_RDOPgwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOPgwso')

     !Gw_PO4gwso  !BK20260323
     iret = nf90_inq_varid(ftnGwP,'Gw_PO4gwso',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_PO4gwso')
     iret = nf90_put_var(ftnGwP,varId,Gw_PO4gwso,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_PO4gwso')

     !Gw_LDOPgwch
     iret = nf90_inq_varid(ftnGwP,'Gw_LDOPgwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOPgwch')
     iret = nf90_put_var(ftnGwP,varId,Gw_LDOPgwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOPgwch')

     !Gw_RDOPgwch
     iret = nf90_inq_varid(ftnGwP,'Gw_RDOPgwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOPgwch')
     iret = nf90_put_var(ftnGwP,varId,Gw_RDOPgwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOPgwch')

     !Gw_PO4gwch
     iret = nf90_inq_varid(ftnGwP,'Gw_PO4gwch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_PO4gwch')
     iret = nf90_put_var(ftnGwP,varId,Gw_PO4gwch,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_PO4gwch')

     !Gw_LDOP
     iret = nf90_inq_varid(ftnGwP,'Gw_LDOP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_LDOP')
     iret = nf90_put_var(ftnGwP,varId,Gw_LDOP,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_LDOP')

     !Gw_RDOP
     iret = nf90_inq_varid(ftnGwP,'Gw_RDOP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_RDOP')
     iret = nf90_put_var(ftnGwP,varId,Gw_RDOP,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_RDOP')

     !Gw_PO4
     iret = nf90_inq_varid(ftnGwP,'Gw_PO4',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_PO4')
     iret = nf90_put_var(ftnGwP,varId,Gw_PO4,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_PO4')

     !Gw_PSC
     iret = nf90_inq_varid(ftnGwP,'Gw_PSC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_PSC')
     iret = nf90_put_var(ftnGwP,varId,Gw_PSC,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_PSC')

     !Gw_PMBError
     iret = nf90_inq_varid(ftnGwP,'Gw_PMBError',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Gw_PMBError')
     iret = nf90_put_var(ftnGwP,varId,Gw_PMBError,(/1,1/),(/domain%nbasin,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Gw_PMBError')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnGwP)

  end subroutine write_Gw_P_output
  

  subroutine write_Ch_P_output(output_flnm)

     use module_SedCNPvariables
     use module_NWM_io_dict

     implicit none

     !integer, parameter                        :: numChchVars     = 53  !BK20250329_2
     !integer, parameter                        :: numChbasinVars  = 14  !BK20250329_2
     integer, parameter                        :: numChPVars      = 67 !(numChchVars + numChbasinVars)
     character (len=64), dimension(numChPVars) :: ChPOut_varNames
     character (len=64), dimension(numChPVars) :: ChPOut_longName  ! Long names for each variable.
     character (len=64), dimension(numChPVars) :: ChPOut_units     ! Units for each variable.
     character(len=256)                        :: output_flnm      ! CHPOUT_DOMAIN filename
     integer                                   :: diagFlag
     integer                                   :: ftn, ftnChP      ! NetCDF file handle
     integer                                   :: iret             ! NetCDF return statuses
     integer                                   :: iTmp
     integer                                   :: varId            ! Variable ID value created as NetCDF variables are created and populated.
     integer                                   :: dimId(4)         ! Dimension ID values created during NetCDF created.

     diagFlag = 0    !BK20251031

     ChPOut_varNames(:) = [character(len=64) :: &
                         "Ch_ALGPusch0_1","Ch_ALGPusch0_2","Ch_ALGPusch0_3","Ch_ZOOPusch0","Ch_LPOPusch0",&   !1-5
                         "Ch_RPOPusch0","Ch_LDOPusch0","Ch_RDOPusch0","Ch_MBMPusch0","Ch_PO4usch0",&          !6-10
                         "Ch_PIPAusch0","Ch_PIPSusch0","Ch_LPOPxpnt","Ch_RPOPxpnt","Ch_LDOPxpnt",&            !11-15
                         "Ch_RDOPxpnt","Ch_PO4xpnt","Ch_LPOPxabs","Ch_RPOPxabs","Ch_LDOPxabs",&               !16-20   !BK20240729
                         "Ch_RDOPxabs","Ch_PO4xabs","Ch_LPOPxdis","Ch_RPOPxdis","Ch_LDOPxdis",&               !21-25   !BK20240729
                         "Ch_RDOPxdis","Ch_PO4xdis","Ch_ALGPdsch_1","Ch_ALGPdsch_2","Ch_ALGPdsch_3",&         !26-30   !BK20240729
                         "Ch_ZOOPdsch","Ch_LPOPdsch","Ch_RPOPdsch","Ch_LDOPdsch","Ch_RDOPdsch",&              !31-35
                         "Ch_MBMPdsch","Ch_PO4dsch","Ch_PIPAdsch","Ch_PIPSdsch","Ch_ALGP_1",&                 !36-40
                         "Ch_ALGP_2","Ch_ALGP_3","Ch_ZOOP","Ch_LPOP","Ch_RPOP",&                              !41-45
                         "Ch_LDOP","Ch_RDOP","Ch_MBMP","Ch_PO4","Ch_PIPA",&                                   !46-50
                         "Ch_PIPS","Ch_PSC","Ch_PMBError","Ch_LPOPsurf0","Ch_RPOPsurf0",&                     !51-55   !BK20250329
                         "Ch_LDOPsurf0","Ch_RDOPsurf0","Ch_MBMPsurf0","Ch_PO4surf0","Ch_PIPAsurf0",&          !56-60   !BK20250329
                         "Ch_PIPSsurf0","Ch_LDOPintf0","Ch_RDOPintf0","Ch_PO4intf0","Ch_LDOPgwch0",&          !61-65   !BK20250329
                         "Ch_RDOPgwch0","Ch_PO4gwch0"]                                                        !66-67   !BK20250329
     ChPOut_longName(:) = [character(len=128) :: &
                         "Ch_ALGPusch0_1","Ch_ALGPusch0_2","Ch_ALGPusch0_3","Ch_ZOOPusch0","Ch_LPOPusch0",&   !1-5
                         "Ch_RPOPusch0","Ch_LDOPusch0","Ch_RDOPusch0","Ch_MBMPusch0","Ch_PO4usch0",&          !6-10
                         "Ch_PIPAusch0","Ch_PIPSusch0","Ch_LPOPxpnt","Ch_RPOPxpnt","Ch_LDOPxpnt",&            !11-15
                         "Ch_RDOPxpnt","Ch_PO4xpnt","Ch_LPOPxabs","Ch_RPOPxabs","Ch_LDOPxabs",&               !16-20   !BK20240729
                         "Ch_RDOPxabs","Ch_PO4xabs","Ch_LPOPxdis","Ch_RPOPxdis","Ch_LDOPxdis",&               !21-25   !BK20240729
                         "Ch_RDOPxdis","Ch_PO4xdis","Ch_ALGPdsch_1","Ch_ALGPdsch_2","Ch_ALGPdsch_3",&         !26-30   !BK20240729
                         "Ch_ZOOPdsch","Ch_LPOPdsch","Ch_RPOPdsch","Ch_LDOPdsch","Ch_RDOPdsch",&              !31-35
                         "Ch_MBMPdsch","Ch_PO4dsch","Ch_PIPAdsch","Ch_PIPSdsch","Ch_ALGP_1",&                 !36-40
                         "Ch_ALGP_2","Ch_ALGP_3","Ch_ZOOP","Ch_LPOP","Ch_RPOP",&                              !41-45
                         "Ch_LDOP","Ch_RDOP","Ch_MBMP","Ch_PO4","Ch_PIPA",&                                   !46-50
                         "Ch_PIPS","Ch_PSC","Ch_PMBError","Ch_LPOPsurf0","Ch_RPOPsurf0",&                     !51-55   !BK20250329
                         "Ch_LDOPsurf0","Ch_RDOPsurf0","Ch_MBMPsurf0","Ch_PO4surf0","Ch_PIPAsurf0",&          !56-60   !BK20250329
                         "Ch_PIPSsurf0","Ch_LDOPintf0","Ch_RDOPintf0","Ch_PO4intf0","Ch_LDOPgwch0",&          !61-65   !BK20250329
                         "Ch_RDOPgwch0","Ch_PO4gwch0"]                                                        !66-67   !BK20250329                         
     ChPOut_units(:) = [character(len=64) :: &
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !1-5
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !6-10
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !11-15
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !16-20
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !21-25
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !26-30
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !31-35
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP",&       !36-40
                         "kgP","kgP","kgP","kgP","kgP",&                           !41-45
                         "kgP","kgP","kgP","kgP","kgP",&                           !46-50
                         "kgP","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&       !51-55
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !56-60
                         "kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1","kgP dt-1",&  !61-65
                         "kgP dt-1","kgP dt-1"]                                    !66-67

     ! Create output file
     iret = nf90_create(trim(output_flnm),cmode=NF90_NETCDF4,ncid = ftn)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to create CHPOUT NetCDF file.')
     ftnChP = ftn

     ! Create dimensions
     iret = nf90_def_dim(ftnChP,'time',NF90_UNLIMITED,dimId(1))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define time dimension')
     iret = nf90_def_dim(ftnChP,'nch',domain%nch,dimId(2))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to define nch dimension')
     !iret = nf90_def_dim(ftnChP,'nbasin',domain%nbasin,dimId(3))              !BK20250329
     !iret = nf90_def_dim(ftnChP,'nch',domain%nch,dimId(3))                     !BK20250329 !BK20250329_2
     !call nwmCheck(diagFlag,iret,'ERROR: Unable to define nbasin dimension')   !BK20250329 !BK20250329_2

     ! Loop through all possible variables and create them, along with their
     ! metadata attributes.

     !create time variable attributes [minutes since 1970-01-01 00:00:00 UTC]  !Y.Kwon20241125
     iret = nf90_def_var(ftnChP,"time",nf90_float,(/dimId(1)/),varId)
     call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//"time")
     iret = nf90_put_att(ftnChP,varId,'long_name','valid output time')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//'time')
     iret = nf90_put_att(ftnChP,varId,'units','minutes since 1970-01-01 00:00:00 UTC')
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//'time')

     do iTmp=1,numChPVars
        !if(iTmp.le.numChchVars) then !BK20250329_2
           iret = nf90_def_var(ftnChP,trim(ChPOut_varNames(iTmp)),nf90_double,(/dimId(2),dimId(1)/),varId)
           call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(ChPOut_varNames(iTmp)))
        !else !numChbasinVars
        !   iret = nf90_def_var(ftnChP,trim(ChPOut_varNames(iTmp)),nf90_double,(/dimId(3),dimId(1)/),varId) !BK20250329_2
        !   call nwmCheck(diagFlag,iret,"ERROR: Unable to create variable: "//trim(ChPOut_varNames(iTmp))) !BK20250329_2
        !endif

        ! Create variable attributes
        iret = nf90_put_att(ftnChP,varId,'long_name',trim(ChPOut_longName(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place long_name attribute into variable '//trim(ChPOut_varNames(iTmp)))
        iret = nf90_put_att(ftnChP,varId,'units',trim(ChPOut_units(iTmp)))
        call nwmCheck(diagFlag,iret,'ERROR: Unable to place units attribute into variable '//trim(ChPOut_varNames(iTmp)))
     enddo

     ! Remove NetCDF file from definition mode
     iret = nf90_enddef(ftnChP)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to take CHPOUT file out of definition mode')

     ! Write time variable to NetCDF file  !Y.Kwon20241125
     iret = nf90_inq_varid(ftnChP,'time',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: time')
     iret = nf90_put_var(ftnChP,varId,SedCNP_hydro%time,(/1/),(/1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: time')

     ! Write array out to NetCDF file
     !Ch_ALGPusch0_1
     iret = nf90_inq_varid(ftnChP,'Ch_ALGPusch0_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGPusch0_1')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGPusch0(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGPusch0_1')

     !Ch_ALGPusch0_2
     iret = nf90_inq_varid(ftnChP,'Ch_ALGPusch0_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGPusch0_2')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGPusch0(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGPusch0_2')

     !Ch_ALGPusch0_3
     iret = nf90_inq_varid(ftnChP,'Ch_ALGPusch0_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGPusch0_3')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGPusch0(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGPusch0_3')

     !Ch_ZOOPusch0
     iret = nf90_inq_varid(ftnChP,'Ch_ZOOPusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOOPusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_ZOOPusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOOPusch0')

     !Ch_LPOPusch0
     iret = nf90_inq_varid(ftnChP,'Ch_LPOPusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOPusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_LPOPusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOPusch0')

     !Ch_RPOPusch0
     iret = nf90_inq_varid(ftnChP,'Ch_RPOPusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOPusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_RPOPusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOPusch0')

     !Ch_LDOPusch0
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPusch0')

     !Ch_RDOPusch0
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPusch0')

     !Ch_MBMPusch0
     iret = nf90_inq_varid(ftnChP,'Ch_MBMPusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMPusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_MBMPusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMPusch0')

     !Ch_PO4usch0
     iret = nf90_inq_varid(ftnChP,'Ch_PO4usch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4usch0')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4usch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4usch0')

     !Ch_PIPAusch0
     iret = nf90_inq_varid(ftnChP,'Ch_PIPAusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPAusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPAusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPAusch0')

     !Ch_PIPSusch0
     iret = nf90_inq_varid(ftnChP,'Ch_PIPSusch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPSusch0')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPSusch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPSusch0')

     !Ch_LPOPxpnt
     iret = nf90_inq_varid(ftnChP,'Ch_LPOPxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOPxpnt')
     iret = nf90_put_var(ftnChP,varId,Ch_LPOPxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOPxpnt')

     !Ch_RPOPxpnt
     iret = nf90_inq_varid(ftnChP,'Ch_RPOPxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOPxpnt')
     iret = nf90_put_var(ftnChP,varId,Ch_RPOPxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOPxpnt')

     !Ch_LDOPxpnt
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPxpnt')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPxpnt')

     !Ch_RDOPxpnt
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPxpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPxpnt')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPxpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPxpnt')

     !Ch_PO4xpnt
     iret = nf90_inq_varid(ftnChP,'Ch_PO4xpnt',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4xpnt')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4xpnt,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4xpnt')

     !Ch_LPOPxabs  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_LPOPxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOPxabs')
     iret = nf90_put_var(ftnChP,varId,Ch_LPOPxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOPxabs')

     !Ch_RPOPxabs  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_RPOPxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOPxabs')
     iret = nf90_put_var(ftnChP,varId,Ch_RPOPxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOPxabs')

     !Ch_LDOPxabs  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPxabs')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPxabs')

     !Ch_RDOPxabs  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPxabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPxabs')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPxabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPxabs')
     
     !Ch_PO4xabs  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_PO4xabs',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4xabs')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4xabs,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4xabs')

     !Ch_LPOPxdis  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_LPOPxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOPxdis')
     iret = nf90_put_var(ftnChP,varId,Ch_LPOPxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOPxdis')

     !Ch_RPOPxdis  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_RPOPxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOPxdis')
     iret = nf90_put_var(ftnChP,varId,Ch_RPOPxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOPxdis')
     
     !Ch_LDOPxdis  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPxdis')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPxdis')

     !Ch_RDOPxdis  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPxdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPxdis')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPxdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPxdis')

     !Ch_PO4xdis  !BK20240729
     iret = nf90_inq_varid(ftnChP,'Ch_PO4xdis',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4xdis')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4xdis,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4xdis')

     !Ch_ALGPdsch_1
     iret = nf90_inq_varid(ftnChP,'Ch_ALGPdsch_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGPdsch_1')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGPdsch(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGPdsch_1')

     !Ch_ALGPdsch_2
     iret = nf90_inq_varid(ftnChP,'Ch_ALGPdsch_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGPdsch_2')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGPdsch(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGPdsch_2')

     !Ch_ALGPdsch_3
     iret = nf90_inq_varid(ftnChP,'Ch_ALGPdsch_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGPdsch_3')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGPdsch(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGPdsch_3')

     !Ch_ZOOPdsch
     iret = nf90_inq_varid(ftnChP,'Ch_ZOOPdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOOPdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_ZOOPdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOOPdsch')

     !Ch_LPOPdsch  !BK20240912 added
     iret = nf90_inq_varid(ftnChP,'Ch_LPOPdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOPdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_LPOPdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOPdsch')

     !Ch_RPOPdsch
     iret = nf90_inq_varid(ftnChP,'Ch_RPOPdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOPdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_RPOPdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOPdsch')

     !Ch_LDOPdsch
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPdsch')

     !Ch_RDOPdsch
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPdsch')

     !Ch_MBMPdsch
     iret = nf90_inq_varid(ftnChP,'Ch_MBMPdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMPdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_MBMPdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMPdsch')

     !Ch_PO4dsch
     iret = nf90_inq_varid(ftnChP,'Ch_PO4dsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4dsch')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4dsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4dsch')

     !Ch_PIPAdsch
     iret = nf90_inq_varid(ftnChP,'Ch_PIPAdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPAdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPAdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPAdsch')

     !Ch_PIPSdsch
     iret = nf90_inq_varid(ftnChP,'Ch_PIPSdsch',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPSdsch')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPSdsch,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPSdsch')

     !Ch_ALGP_1
     iret = nf90_inq_varid(ftnChP,'Ch_ALGP_1',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGP_1')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGP(:,1),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGP_1')

     !Ch_ALGP_2
     iret = nf90_inq_varid(ftnChP,'Ch_ALGP_2',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGP_2')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGP(:,2),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGP_2')

     !Ch_ALGP_3
     iret = nf90_inq_varid(ftnChP,'Ch_ALGP_3',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ALGP_3')
     iret = nf90_put_var(ftnChP,varId,Ch_ALGP(:,3),(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ALGP_3')

     !Ch_ZOOP
     iret = nf90_inq_varid(ftnChP,'Ch_ZOOP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_ZOOP')
     iret = nf90_put_var(ftnChP,varId,Ch_ZOOP,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_ZOOP')

     !Ch_LPOP
     iret = nf90_inq_varid(ftnChP,'Ch_LPOP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOP')
     iret = nf90_put_var(ftnChP,varId,Ch_LPOP,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOP')

     !Ch_RPOP
     iret = nf90_inq_varid(ftnChP,'Ch_RPOP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOP')
     iret = nf90_put_var(ftnChP,varId,Ch_RPOP,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOP')

     !Ch_LDOP
     iret = nf90_inq_varid(ftnChP,'Ch_LDOP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOP')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOP,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOP')

     !Ch_RDOP
     iret = nf90_inq_varid(ftnChP,'Ch_RDOP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOP')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOP,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOP')

     !Ch_MBMP
     iret = nf90_inq_varid(ftnChP,'Ch_MBMP',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMP')
     iret = nf90_put_var(ftnChP,varId,Ch_MBMP,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMP')

     !Ch_PO4
     iret = nf90_inq_varid(ftnChP,'Ch_PO4',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4')

     !Ch_PIPA
     iret = nf90_inq_varid(ftnChP,'Ch_PIPA',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPA')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPA,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPA')

     !Ch_PIPS
     iret = nf90_inq_varid(ftnChP,'Ch_PIPS',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPS')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPS,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPS')

     !Ch_PSC
     iret = nf90_inq_varid(ftnChP,'Ch_PSC',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PSC')
     iret = nf90_put_var(ftnChP,varId,Ch_PSC,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PSC')

     !Ch_PMBError
     iret = nf90_inq_varid(ftnChP,'Ch_PMBError',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PMBError')
     iret = nf90_put_var(ftnChP,varId,Ch_PMBError,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PMBError')

     !Ch_LPOPsurf0
     iret = nf90_inq_varid(ftnChP,'Ch_LPOPsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LPOPsurf0')
     iret = nf90_put_var(ftnChP,varId,Ch_LPOPsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LPOPsurf0')

     !Ch_RPOPsurf0
     iret = nf90_inq_varid(ftnChP,'Ch_RPOPsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RPOPsurf0')
     iret = nf90_put_var(ftnChP,varId,Ch_RPOPsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RPOPsurf0')

     !Ch_LDOPsurf0
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPsurf0')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPsurf0')

     !Ch_RDOPsurf0
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPsurf0')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPsurf0')

     !Ch_MBMPsurf0
     iret = nf90_inq_varid(ftnChP,'Ch_MBMPsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_MBMPsurf0')
     iret = nf90_put_var(ftnChP,varId,Ch_MBMPsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_MBMPsurf0')

     !Ch_PO4surf0
     iret = nf90_inq_varid(ftnChP,'Ch_PO4surf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4surf0')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4surf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4surf0')

     !Ch_PIPAsurf0
     iret = nf90_inq_varid(ftnChP,'Ch_PIPAsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPAsurf0')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPAsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPAsurf0')

     !Ch_PIPSsurf0
     iret = nf90_inq_varid(ftnChP,'Ch_PIPSsurf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PIPSsurf0')
     iret = nf90_put_var(ftnChP,varId,Ch_PIPSsurf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PIPSsurf0')

     !Ch_LDOPintf0
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPintf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPintf0')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPintf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPintf0')

     !Ch_RDOPintf0
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPintf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPintf0')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPintf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPintf0')

     !Ch_PO4intf0
     iret = nf90_inq_varid(ftnChP,'Ch_PO4intf0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4intf0')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4intf0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4intf0')

     !Ch_LDOPgwch0
     iret = nf90_inq_varid(ftnChP,'Ch_LDOPgwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_LDOPgwch0')
     iret = nf90_put_var(ftnChP,varId,Ch_LDOPgwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_LDOPgwch0')

     !Ch_RDOPgwch0
     iret = nf90_inq_varid(ftnChP,'Ch_RDOPgwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_RDOPgwch0')
     iret = nf90_put_var(ftnChP,varId,Ch_RDOPgwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_RDOPgwch0')

     !Ch_PO4gwch0
     iret = nf90_inq_varid(ftnChP,'Ch_PO4gwch0',varId)
     call nwmCheck(diagFlag,iret,'ERROR: Unable to find variable ID for var: Ch_PO4gwch0')
     iret = nf90_put_var(ftnChP,varId,Ch_PO4gwch0,(/1,1/),(/domain%nch,1/))
     call nwmCheck(diagFlag,iret,'ERROR: Unable to place data into output variable: Ch_PO4gwch0')

     ! close the netCDF file JES 20231212
     iret = nf90_close(ftnChP)
     
  end subroutine write_Ch_P_output
 
end module module_SedCNP_out
