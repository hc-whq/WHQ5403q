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

module module_SedCNP_in

  use netcdf
  use CNPfunctions
  use module_SedCNPvariables  
  !use module_hydro_stop, only:HYDRO_stop
!=====||__WHQ5403q__||=====!
!
  use SedCNP_config,         only: SedCNPmodel   !BK20231016
!
!=====||__WHQ5403q__||=====!
  
  contains

  integer function get_lsm_time(var_name,out_value,fileName)  !Y.Kwon20241125

  implicit none
  
  character(len=*),                      intent(in)  :: var_name,fileName
  integer                                            :: out_value(1)
  integer                                            :: ntime, timeDimid
  integer                                            :: out_buff(1)
  integer                                            :: iret,ivar,varid,nid
  character(len=256)                                 :: errMsg

  get_lsm_time = -1

  iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
  if (iret .ne. 0) then
     errMsg = "get_lsm_time: failed to open the netcdf file: " // trim(fileName)
     print*, trim(errMsg)
     !if(fatalErr_local) call hydro_stop(trim(errMsg))
     out_buff = -9999.
     return
  endif

  iret = nf90_inq_dimid(nid, "time", timeDimid)
  if(iret .ne. 0) then
     errMsg = "WARNING: get_lsm_time: failed to get the time dimension id: " //      &
              ' in ' // trim(fileName)
     write(6,*) errMsg
  endif

  iret = nf90_Inquire_Dimension(nid, timeDimid, len = ntime)
  if(iret .ne. 0) then
     errMsg = "WARNING: get_lsm_time: failed to get the time dimension: " //      &
              ' in ' // trim(fileName)
     write(6,*) errMsg
  endif

  ivar = nf90_inq_varid(nid,trim(var_name),  varid)
  if(ivar .ne. 0) then
     ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
     if(ivar .ne. 0) then
        errMsg = "WARNING: get_lsm_time: failed to find the variables: " //      &
                  trim(var_name) // ' and ' // trim(var_name//"_M") // &
                  ' in ' // trim(fileName)
        write(6,*) errMsg
        !if(fatalErr_local) call hydro_stop(errMsg)
        return
     endif
  end if

  iret = nf90_get_var(nid, varid, out_buff(1))
  if(iret .ne. 0) then
     errMsg = "WARNING: get_lsm_time: failed to read the variable: " // &
              trim(var_name) // ' or ' // trim(var_name//"_M") // &
              ' in ' // trim(fileName)
     print*,trim(errMsg)
     !if(fatalErr_local) call hydro_stop(trim(errMsg))
     return
  endif

  iret = nf90_close(nid)
  if(iret .ne. 0) then
     errMsg = "WARNING: get_lsm_time: failed to close the file: " // &
              trim(fileName)
     print*,trim(errMsg)
     !if(fatalErr_local) call hydro_stop(trim(errMsg))
  endif

  out_value(1) = out_buff(1)

  get_lsm_time = ivar

  end function get_lsm_time  !Y.Kwon20241125

!--------------------------------------------------------------------Y.Kwon(20250621)
  subroutine get_feature_id(out_value,fileName)  !Y.Kwon20241125

     use netcdf

     implicit none

     character(len=*),                      intent(in)  :: fileName
     integer,                               intent(out) :: out_value
     integer                                            :: feature_len,Dimid
     integer                                            :: iret,nid
     character(len=1024)                                 :: errMsg

     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get_feature_id: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        return
     endif

     iret = nf90_inq_dimid(nid, "feature_id", Dimid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_feature_id: failed to get the feature id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
     endif

     iret = nf90_Inquire_Dimension(nid, Dimid, len = feature_len)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_feature_id: failed to get the feature id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
     endif

     iret = nf90_close(nid)

     out_value = feature_len

  end subroutine get_feature_id

!--------------------------------------------------------------------Y.Kwon(20250621)

!----------------------------------------------------------BK20251208
  subroutine get_global_iswater(fileName, out_value)
    use netcdf
    implicit none
    character(len=*), intent(in) :: fileName
    integer, intent(out) :: out_value
    integer :: ncid, iret, val
    
    out_value = -999 ! Default error value
    
    iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=ncid)
    if (iret /= 0) return 

    iret = nf90_get_att(ncid, NF90_GLOBAL, "ISWATER", val)
    if (iret == 0) then
       out_value = val
    else
        ! Try lowercase
        iret = nf90_get_att(ncid, NF90_GLOBAL, "iswater", val)
        if (iret == 0) out_value = val
    endif
    
    iret = nf90_close(ncid)
  end subroutine get_global_iswater
!----------------------------------------------------------BK20251208

!--------------------------------------------------------------------Y.Kwon(20250621)
  subroutine get_ixjxnsl(out_value_i,out_value_j,out_value_nsl,fileName)

     use netcdf

     implicit none

     character(len=*),                      intent(in)  :: fileName
     integer,                               intent(out) :: out_value_i !longitude (LSM grid)
     integer,                               intent(out) :: out_value_j !latitude (LSM grid)
     integer,                               intent(out) :: out_value_nsl !number of soil layers (default nsl = 4)
     integer                                            :: ix_len, jx_len, nsl_len
     integer                                            :: Dimid_ix, Dimid_jx, Dimid_nsl
     integer                                            :: iret,nid
     character(len=1024)                                 :: errMsg

     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get_ixjxnsl: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        return
     endif

     iret = nf90_inq_dimid(nid, "west_east", Dimid_ix)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_ixjxnsl: failed to get west_east: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
     endif

     iret = nf90_inq_dimid(nid, "south_north", Dimid_jx)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_ixjxnsl: failed to get south_north: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
     endif

     iret = nf90_inq_dimid(nid, "soil_layers_stag", Dimid_nsl)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_ixjxnsl: failed to get soil_layers_stag: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
     endif

     iret = nf90_Inquire_Dimension(nid, Dimid_ix, len = ix_len)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_ixjxnsl: failed to get west_east: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
     endif

     iret = nf90_Inquire_Dimension(nid, Dimid_jx, len = jx_len)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_ixjxnsl: failed to get south_north: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
          endif

     iret = nf90_Inquire_Dimension(nid, Dimid_nsl, len = nsl_len)
     if(iret .ne. 0) then
        errMsg = "WARNING: get_ixjxnsl: failed to get soil_layers_stag: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
        iret = nf90_close(nid)
        return
     endif

     iret = nf90_close(nid)

     out_value_i   = ix_len
     out_value_j   = jx_len
     out_value_nsl = nsl_len

  end subroutine get_ixjxnsl

!--------------------------------------------------------------------Y.Kwon(20250621)

  !--------------------------------------------------------------------Y.Kwon(20250621) !BK20251210
  real(8) function get_lsm_dx(var_name,fileName)
  ! Read grid spacing (dx) from geo_static file specified in hydro.namelist
  ! This is the most reliable method since geo_static files always contain DX attribute
  ! The value is cached after first read to avoid repeated file access

  implicit none

  character(len=*), intent(in) :: var_name, fileName
  real :: dx_value
  real(8), save :: cached_dx = -9999.0_8  ! Static variable to cache dx value
  logical, save :: first_call = .true.    ! Flag to track first call

  ! If already read, return cached value
  if (.not. first_call .and. cached_dx > 0.0_8) then
     get_lsm_dx = cached_dx
     return
  endif
  
  get_lsm_dx = -9999.0_8
  
  ! Read DX from geo_static file (only on first call)
  call read_dx_from_geo_static(dx_value)
  
  if (dx_value > 0.0) then
     cached_dx = real(dx_value, 8)
     get_lsm_dx = cached_dx
     first_call = .false.
     !write(6,'(A,F10.2,A)') 'INFO: get_lsm_dx: Successfully read DX = ', get_lsm_dx, ' from the geo_static file'  !BK20251210
  else
     write(6,*) 'ERROR: get_lsm_dx: Could not read DX from geo_static file!'  !BK20251210
     !write(6,*) '       Please ensure GEO_STATIC_FLNM is set in hydro.namelist'  !BK20251210
     !write(6,*) '       and that the file contains a DX global attribute'  !BK20251210
  endif

  end function get_lsm_dx
!----------------------------------------------------------Y.Kwon(20250621)

  subroutine read_dx_from_geo_static(dx_out)
  ! Helper subroutine to read DX from geo_static file specified in hydro.namelist
  implicit none
  real, intent(out) :: dx_out
  integer :: iunit, ierr, iret, nid
  character(len=1024) :: line, geo_static_file
  real :: dx_temp
  
  dx_out = -9999.0
  geo_static_file = ''
  
  ! Read GEO_STATIC_FLNM from hydro.namelist
  open(newunit=iunit, file='hydro.namelist', status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
     open(newunit=iunit, file='./hydro.namelist', status='old', action='read', iostat=ierr)
  endif
  
  if (ierr == 0) then
     ! Read file line by line looking for GEO_STATIC_FLNM
     do
        read(iunit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit
        
        ! Look for GEO_STATIC_FLNM = "filename"
        if (index(line, 'GEO_STATIC_FLNM') > 0 .and. index(line, '=') > 0) then
           ! Extract the filename between quotes
           call extract_quoted_string(line, geo_static_file)
           if (len_trim(geo_static_file) > 0) exit
        endif
     enddo
     close(iunit)
  endif
  
  ! If we found the geo_static file, try to read DX from it
  if (len_trim(geo_static_file) > 0) then
     !write(6,*) '         Found GEO_STATIC_FLNM: ', trim(geo_static_file)  !BK20251210
     iret = nf90_open(path=trim(geo_static_file), mode=NF90_NOWRITE, ncid=nid)
     if (iret == 0) then
        ! Try to read DX attribute
        iret = nf90_get_att(nid, NF90_GLOBAL, "DX", dx_temp)
        if (iret == 0 .and. dx_temp > 0.0) then
           dx_out = dx_temp
        else
           ! Try lowercase dx
           iret = nf90_get_att(nid, NF90_GLOBAL, "dx", dx_temp)
           if (iret == 0 .and. dx_temp > 0.0) then
              dx_out = dx_temp
           endif
        endif
        iret = nf90_close(nid)
     else
        write(6,*) '         WARNING: Could not open geo_static file: ', trim(geo_static_file)  !BK20251210
     endif
  else
     write(6,*) '         WARNING: Could not find GEO_STATIC_FLNM in hydro.namelist'  !BK20251210
  endif
  
  end subroutine read_dx_from_geo_static
!----------------------------------------------------------BK20251210

  subroutine extract_quoted_string(line, output)
  ! Extract string between quotes from a line
  implicit none
  character(len=*), intent(in) :: line
  character(len=*), intent(out) :: output
  integer :: start_pos, end_pos
  
  output = ''
  
  ! Find first quote
  start_pos = index(line, '"')
  if (start_pos == 0) then
     start_pos = index(line, "'")
  endif
  
  if (start_pos > 0) then
     ! Find second quote
     end_pos = index(line(start_pos+1:), '"')
     if (end_pos == 0) then
        end_pos = index(line(start_pos+1:), "'")
     endif
     
     if (end_pos > 0) then
        output = line(start_pos+1:start_pos+end_pos-1)
     endif
  endif
  
  end subroutine extract_quoted_string
!----------------------------------------------------------BK20251210

  integer function get3d_lsm_real(var_name,out_value,ix,jx,fileName)

     implicit none

     integer,                               intent(in)  :: ix,jx
     character(len=*),                      intent(in)  :: var_name,fileName
     real(8),         dimension(ix,jx),     intent(out) :: out_value
     integer                                            :: ntime, timeDimid
     real,allocatable,dimension(:,:,:)                  :: out_buff
     integer                                            :: iret,ivar,varid,nid
     !logical                                            :: fatalErr_local
     character(len=256)                                 :: errMsg

     !fatalErr_local = .false.

     get3d_lsm_real = -1


     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get3d_lsm_real: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        out_buff = -9999.
        return
     endif

     iret = nf90_inq_dimid(nid, "time", timeDimid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_real: failed to get the time dimension id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     iret = nf90_Inquire_Dimension(nid, timeDimid, len = ntime)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_real: failed to get the time dimension: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     !write(8765, *) 'ntime= ', ntime

     allocate(out_buff(ix,jx,ntime))

     ivar = nf90_inq_varid(nid,trim(var_name),  varid)
     if(ivar .ne. 0) then
        ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
        if(ivar .ne. 0) then
           errMsg = "WARNING: get3d_lsm_real: failed to find the variables: " //      &
                     trim(var_name) // ' and ' // trim(var_name//"_M") // &
                     ' in ' // trim(fileName)
           write(6,*) errMsg
           !if(fatalErr_local) call hydro_stop(errMsg)
           return
        endif
     end if
    
     iret = nf90_get_var(nid, varid, out_buff)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_real: failed to read the variable: " // &
                 trim(var_name) // ' or ' // trim(var_name//"_M") // &
                 ' in ' // trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        return
     endif 
     
     iret = nf90_close(nid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_real: failed to close the file: " // &
                 trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
     endif

     out_value(:,:) = out_buff(:,:,ntime)

     !kyh_debug
     !write(8765, *) 'out_buff(3,16,ntime)=', out_buff(3,16,ntime)

     deallocate(out_buff)

     get3d_lsm_real = ivar

  end function get3d_lsm_real


  integer function get4d_lsm_real(var_name,out_value,ix,jx,fileName)

     implicit none

     integer,                               intent(in)  :: ix,jx
     character(len=*),                      intent(in)  :: var_name,fileName
     real(8),         dimension(ix,4,jx),   intent(out) :: out_value
     integer                                            :: ntime, timeDimid
     integer                                            :: nsoillayer, soillayerDimid
     real,allocatable,dimension(:,:,:,:)                :: out_buff
     integer                                            :: iret,ivar,varid,nid
     !logical                                            :: fatalErr_local
     character(len=256)                                 :: errMsg

     !fatalErr_local = .false.

     get4d_lsm_real = -1

     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get4d_lsm_real: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        out_buff = -9999.
        return
     endif

     !time dimension
     iret = nf90_inq_dimid(nid, "time", timeDimid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_lsm_real: failed to get the time dimension id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     iret = nf90_Inquire_Dimension(nid, timeDimid, len = ntime)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_lsm_real: failed to get the time dimension: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     !soil layer dimension
     iret = nf90_inq_dimid(nid, "soil_layers_stag", soillayerDimid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_lsm_real: failed to get the time dimension id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     iret = nf90_Inquire_Dimension(nid, soillayerDimid, len = nsoillayer)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_lsm_real: failed to get the soil layer dimension: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     allocate(out_buff(ix,nsoillayer,jx,ntime))

     ivar = nf90_inq_varid(nid,trim(var_name),  varid)
     if(ivar .ne. 0) then
        ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
        if(ivar .ne. 0) then
           errMsg = "WARNING: get4d_lsm_real: failed to find the variables: " //      &
                     trim(var_name) // ' and ' // trim(var_name//"_M") // &
                     ' in ' // trim(fileName)
           write(6,*) errMsg
           !if(fatalErr_local) call hydro_stop(errMsg)
           return
        endif
     end if

     iret = nf90_get_var(nid, varid, out_buff)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_lsm_real: failed to read the variable: " // &
                 trim(var_name) // ' or ' // trim(var_name//"_M") // &
                 ' in ' // trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        return
     endif

     iret = nf90_close(nid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_lsm_real: failed to close the file: " // &
                 trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
     endif

     out_value(:,:,:) = out_buff(:,:,:,ntime)

     deallocate(out_buff)

     get4d_lsm_real = ivar

  end function get4d_lsm_real


  integer function get4d_soil_real(var_name,out_value,ix,jx,fileName)

     implicit none

     integer,                               intent(in)  :: ix,jx
     character(len=*),                      intent(in)  :: var_name,fileName
     real(8),         dimension(ix,jx,4),   intent(out) :: out_value
     integer                                            :: ntime, timeDimid
     integer                                            :: nsoillayer, soillayerDimid
     real,allocatable,dimension(:,:,:,:)                :: out_buff
     integer                                            :: iret,ivar,varid,nid
     !logical                                            :: fatalErr_local
     character(len=256)                                 :: errMsg

     !fatalErr_local = .false.

     get4d_soil_real = -1

     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get4d_soil_real: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        out_buff = -9999.
        return
     endif

     !time dimension
     iret = nf90_inq_dimid(nid, "Time", timeDimid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_soil_real: failed to get the time dimension id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     iret = nf90_Inquire_Dimension(nid, timeDimid, len = ntime)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_soil_real: failed to get the time dimension: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     !soil layer dimension
     iret = nf90_inq_dimid(nid, "soil_layers_stag", soillayerDimid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_soil_real: failed to get the time dimension id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     iret = nf90_Inquire_Dimension(nid, soillayerDimid, len = nsoillayer)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_soil_real: failed to get the soil layer dimension: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     allocate(out_buff(ix,jx,nsoillayer,ntime))

     ivar = nf90_inq_varid(nid,trim(var_name),  varid)
     if(ivar .ne. 0) then
        ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
        if(ivar .ne. 0) then
           errMsg = "WARNING: get4d_soil_real: failed to find the variables: " //      &
                     trim(var_name) // ' and ' // trim(var_name//"_M") // &
                     ' in ' // trim(fileName)
           write(6,*) errMsg
           !if(fatalErr_local) call hydro_stop(errMsg)
           return
        endif
     end if

     iret = nf90_get_var(nid, varid, out_buff)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_soil_real: failed to read the variable: " // &
                 trim(var_name) // ' or ' // trim(var_name//"_M") // &
                 ' in ' // trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        return
     endif

     iret = nf90_close(nid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get4d_soil_real: failed to close the file: " // &
                 trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
     endif

     out_value(:,:,:) = out_buff(:,:,:,ntime)

     !kyh_debug
     !write(87651122, *) 'out_buff(3,4,1,ntime)=', out_buff(3,4,1,ntime)
     !write(87651122, *) 'out_buff(3,4,2,ntime)=', out_buff(3,4,2,ntime)
     !write(87651122, *) 'out_buff(3,4,3,ntime)=', out_buff(3,4,3,ntime)
     !write(87651122, *) 'out_buff(3,4,4,ntime)=', out_buff(3,4,4,ntime)

     deallocate(out_buff)

     get4d_soil_real = ivar

  end function get4d_soil_real


  integer function get1d_ch_real(var_name,out_value,nch,fileName)

     implicit none

     integer,                              intent(in)  :: nch  !number of channel reaches
     character(len=*),                     intent(in)  :: var_name,fileName
     real(8),         dimension(nch),      intent(out) :: out_value
     integer                                           :: ntime, timeDimid
     real,allocatable,dimension(:)                     :: out_buff
     integer                                           :: iret,ivar,varid,nid
     !logical                                           :: fatalErr_local
     character(len=256)                                :: errMsg

     !fatalErr_local = .false.

     get1d_ch_real = -1


     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get1d_ch_real: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        out_buff = -9999.
        return
     endif

     !iret = nf90_inq_dimid(nid, "time", timeDimid)
     !if(iret .ne. 0) then
     !   errMsg = "WARNING: get1d_ch_real: failed to get the time dimension id: " //      &
     !            ' in ' // trim(fileName)
     !   write(6,*) errMsg
     !endif

     !iret = nf90_Inquire_Dimension(nid, timeDimid, len = ntime)
     !if(iret .ne. 0) then
     !   errMsg = "WARNING: get1d_ch_real: failed to get the time dimension: " //      &
     !            ' in ' // trim(fileName)
     !   write(6,*) errMsg
     !endif

     allocate(out_buff(nch))

     ivar = nf90_inq_varid(nid,trim(var_name),  varid)
     if(ivar .ne. 0) then
        ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
        if(ivar .ne. 0) then
           errMsg = "WARNING: get1d_ch_real: failed to find the variables: " //      &
                     trim(var_name) // ' and ' // trim(var_name//"_M") // &
                     ' in ' // trim(fileName)
           write(6,*) errMsg
           !if(fatalErr_local) call hydro_stop(errMsg)
           return
        endif
     end if

     iret = nf90_get_var(nid, varid, out_buff)
     if(iret .ne. 0) then
        errMsg = "WARNING: get1d_ch_real: failed to read the variable: " // &
                 trim(var_name) // ' or ' // trim(var_name//"_M") // &
                 ' in ' // trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        return
     endif

     iret = nf90_close(nid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get1d_lsm_real: failed to close the file: " // &
                 trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
     endif

     out_value(:) = out_buff(:)

     !kyh_debug
     !write(87651, *) 'out_buff(10)=', out_buff(10)

     deallocate(out_buff)

     get1d_ch_real = ivar

  end function get1d_ch_real  


  integer function get1d_ch_int(var_name,out_value,nch,fileName)

     implicit none

     integer,                              intent(in)  :: nch  !number of channel reaches
     character(len=*),                     intent(in)  :: var_name,fileName
     integer,            dimension(nch),   intent(out) :: out_value
     integer                                           :: ntime, timeDimid
     integer,allocatable,dimension(:)                  :: out_buff
     integer                                           :: iret,ivar,varid,nid
     character(len=256)                                :: errMsg 

     get1d_ch_int = -1

     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get1d_ch_int: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        out_buff = -9999
        return
     endif 

     allocate(out_buff(nch))

     ivar = nf90_inq_varid(nid,trim(var_name),  varid)
     if(ivar .ne. 0) then
        ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
        if(ivar .ne. 0) then
           errMsg = "WARNING: get1d_ch_int: failed to find the variables: " //      &
                     trim(var_name) // ' and ' // trim(var_name//"_M") // &
                     ' in ' // trim(fileName)
           write(6,*) errMsg
           !if(fatalErr_local) call hydro_stop(errMsg)
           return
        endif
     end if

     iret = nf90_get_var(nid, varid, out_buff)
     if(iret .ne. 0) then
        errMsg = "WARNING: get1d_ch_int: failed to read the variable: " // &
                 trim(var_name) // ' or ' // trim(var_name//"_M") // &
                 ' in ' // trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        return
     endif

     iret = nf90_close(nid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get1d_ch_int: failed to close the file: " // &
                 trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
     endif

     out_value(:) = out_buff(:)

     deallocate(out_buff)

     get1d_ch_int = ivar

  end function get1d_ch_int

  
  integer function get2d_int(var_name,out_value,ix,jx,fileName)

     implicit none

     integer,                               intent(in)  :: ix,jx
     character(len=*),                      intent(in)  :: var_name,fileName
     integer,         dimension(ix,jx),     intent(out) :: out_value
     integer                                            :: ntime, timeDimid
     integer,allocatable,dimension(:,:)                 :: out_buff
     integer                                            :: iret,ivar,varid,nid
     !logical                                            :: fatalErr_local
     character(len=256)                                 :: errMsg

     !fatalErr_local = .false.

     get2d_int = -1

     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get2d_int: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        out_buff = -9999
        return
     endif

     allocate(out_buff(ix,jx))

     ivar = nf90_inq_varid(nid,trim(var_name),  varid)
     if(ivar .ne. 0) then
        ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
        if(ivar .ne. 0) then
           errMsg = "WARNING: get2d_int: failed to find the variables: " //      &
                     trim(var_name) // ' and ' // trim(var_name//"_M") // &
                     ' in ' // trim(fileName)
           write(6,*) errMsg
           !if(fatalErr_local) call hydro_stop(errMsg)
           return
        endif
     end if

     iret = nf90_get_var(nid, varid, out_buff)
     if(iret .ne. 0) then
        errMsg = "WARNING: get2d_int: failed to read the variable: " // &
                 trim(var_name) // ' or ' // trim(var_name//"_M") // &
                 ' in ' // trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        return
     endif

     iret = nf90_close(nid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get2d_int: failed to close the file: " // &
                 trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
     endif

     out_value(:,:) = out_buff(:,:)

     !kyh_debug
     !write(8765, *) 'out_buff(3,16,ntime)=', out_buff(3,16,ntime)

     deallocate(out_buff)

     get2d_int = ivar

  end function get2d_int


  integer function get3d_lsm_int(var_name,out_value,ix,jx,fileName)    !BK20231005

     implicit none

     integer,                               intent(in)  :: ix,jx
     character(len=*),                      intent(in)  :: var_name,fileName
     integer,            dimension(ix,jx),  intent(out) :: out_value
     integer                                            :: ntime, timeDimid
     integer,allocatable,dimension(:,:,:)               :: out_buff
     integer                                            :: iret,ivar,varid,nid
     !logical                                            :: fatalErr_local
     character(len=256)                                 :: errMsg

     !fatalErr_local = .false.

     get3d_lsm_int = -1

     iret = nf90_open(path=trim(fileName), mode=NF90_NOWRITE, ncid=nid)
     if (iret .ne. 0) then
        errMsg = "get3d_lsm_int: failed to open the netcdf file: " // trim(fileName)
        print*, trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        out_buff = -9999
        return
     endif

     iret = nf90_inq_dimid(nid, "time", timeDimid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_int: failed to get the time dimension id: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     iret = nf90_Inquire_Dimension(nid, timeDimid, len = ntime)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_int: failed to get the time dimension: " //      &
                 ' in ' // trim(fileName)
        write(6,*) errMsg
     endif

     allocate(out_buff(ix,jx,ntime))

     ivar = nf90_inq_varid(nid,trim(var_name),  varid)
     if(ivar .ne. 0) then
        ivar = nf90_inq_varid(nid,trim(var_name//"_M"),  varid)
        if(ivar .ne. 0) then
           errMsg = "WARNING: get3d_lsm_int: failed to find the variables: " //      &
                     trim(var_name) // ' and ' // trim(var_name//"_M") // &
                     ' in ' // trim(fileName)
           write(6,*) errMsg
           !if(fatalErr_local) call hydro_stop(errMsg)
           return
        endif
     end if
    
     iret = nf90_get_var(nid, varid, out_buff)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_int: failed to read the variable: " // &
                 trim(var_name) // ' or ' // trim(var_name//"_M") // &
                 ' in ' // trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
        return
     endif 
     
     iret = nf90_close(nid)
     if(iret .ne. 0) then
        errMsg = "WARNING: get3d_lsm_int: failed to close the file: " // &
                 trim(fileName)
        print*,trim(errMsg)
        !if(fatalErr_local) call hydro_stop(trim(errMsg))
     endif

     out_value(:,:) = out_buff(:,:,ntime)

     deallocate(out_buff)

     get3d_lsm_int = ivar

  end function get3d_lsm_int   !BK20231005


  subroutine ReadCNPinputs()  !BK20231024

     !to read in monthly CNP inputs from each of the land covers (CNPinputs.dat)

     !use CNPvariables
     !use CNPfunctions
     !use module_SedCNPvariables
     use module_hydro_stop, only:HYDRO_stop  !BK20251130

     implicit none

     !integer, intent(in) :: itime  !Y.Kwon  !BK20231024
     real, dimension(nlc,12) :: dLPOCLIT, dLPOCRES, dLPOCEXC, dLPOCMAN  !BK20231024
     real, dimension(nlc,12) :: dLPONLIT, dLPONRES, dLPONEXC, dLPONMAN, dNH4fert, dNO3fert  !BK20231024
     real, dimension(nlc,12) :: dLPOPLIT, dLPOPRES, dLPOPEXC, dLPOPMAN, dPO4fert  !BK20231024
     integer :: ilc   !land cover index   !BK20231005
     ! real :: ilc   !land code index
     character(len=3) :: dum
     !character(len=2) :: smonth  !BK20240117
     !integer :: itime
     integer :: imonth
     !real(8) :: dt           ! time-step (sec)  !Y.Kwon
     
     !BK20251128 
     logical :: file_exists
     integer :: ierr

     !DT = real(SedCNPmodel%SedCNP_timestep)  !get model time-step  !BK20240117

     !initialize (Y.Kwon)
     So_LPOCLITinso = 0.
     So_LPOCRESinso = 0.
     So_LPOCEXCinso = 0.
     So_LPOCMANinso = 0.
     So_LPONLITinso = 0.
     So_LPONRESinso = 0.
     So_LPONEXCinso = 0.
     So_LPONMANinso = 0.
     So_NH4fert     = 0.
     So_NO3fert     = 0.
     So_LPOPLITinso = 0.
     So_LPOPRESinso = 0.
     So_LPOPEXCinso = 0.
     So_LPOPMANinso = 0.
     So_PO4fert     = 0.


     if (len_trim(SedCNPmodel%CNPinputs_file) > 0) then
        inquire(file=trim(SedCNPmodel%CNPinputs_file), exist=file_exists)
        if (file_exists) then
           open (590005, file=trim(SedCNPmodel%CNPinputs_file), status='unknown', iostat=ierr)
        else
            ierr = 1 ! Flag to try defaults
        endif
     else
        ierr = 1 ! Flag to try defaults
     endif

     if (ierr /= 0) then
         !BK20251128 Check for default files
         inquire(file='./WHQIN/CNPinputs.dat', exist=file_exists)
         if (file_exists) then
            open (590005, file='./WHQIN/CNPinputs.dat', status='unknown', iostat=ierr)
         else
            open (590005, file='./WHQIN/CNPinputs_TEST.dat', status='unknown', iostat=ierr)
         endif
     endif
     
     if (ierr /= 0) then
        call hydro_stop("Failed to open CNPinputs data file")
     endif
     !read (590005,*)

     do ilc = 1, nlc
        read(590005, *) !table header
        do imonth = 1, 12
           ! read in monthly-average inputs of organic matter (kg ha-1 day-1)
           read(590005, *) dum, dLPOCLIT(ilc,imonth), dLPOCRES(ilc,imonth), dLPOCEXC(ilc,imonth), dLPOCMAN(ilc,imonth), &
              dLPONLIT(ilc,imonth), dLPONRES(ilc,imonth), dLPONEXC(ilc,imonth), dLPONMAN(ilc,imonth), dNH4fert(ilc,imonth), dNO3fert(ilc,imonth), &
              dLPOPLIT(ilc,imonth), dLPOPRES(ilc,imonth), dLPOPEXC(ilc,imonth), dLPOPMAN(ilc,imonth), dPO4fert(ilc,imonth)   ! (kg ha-1 day-1)

           !BK20240117
           ! - LDOMMAN and RDOMMAN to be added later to reflect inputs of liquid manure
           ! - read-in format should be changed later as follows:
           !       LPOMLIT, LPOMRES, LPOMEXC - read monthly-average inputs and convert into daily inputs
           !       by using a mean-preserving interpolation method
           !       LPOMMAN, NH4fert, NO3fert, PO4fert - read irregular time-series inputs
           So_LPOCLITinso(ilc,imonth) = dLPOCLIT(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPOCRESinso(ilc,imonth) = dLPOCRES(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPOCEXCinso(ilc,imonth) = dLPOCEXC(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPOCMANinso(ilc,imonth) = dLPOCMAN(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPONLITinso(ilc,imonth) = dLPONLIT(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPONRESinso(ilc,imonth) = dLPONRES(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPONEXCinso(ilc,imonth) = dLPONEXC(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPONMANinso(ilc,imonth) = dLPONMAN(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_NH4fert(ilc,imonth) = dNH4fert(ilc,imonth) * (domain%areaxy / 10000.) / 86400.       !kg sec-1
           So_NO3fert(ilc,imonth) = dNO3fert(ilc,imonth) * (domain%areaxy / 10000.) / 86400.       !kg sec-1
           So_LPOPLITinso(ilc,imonth) = dLPOPLIT(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPOPRESinso(ilc,imonth) = dLPOPRES(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPOPEXCinso(ilc,imonth) = dLPOPEXC(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_LPOPMANinso(ilc,imonth) = dLPOPMAN(ilc,imonth) * (domain%areaxy / 10000.) / 86400.   !kg sec-1
           So_PO4fert(ilc,imonth) = dPO4fert(ilc,imonth) * (domain%areaxy / 10000.) / 86400.       !kg sec-1
           ! *** reserved for future development ***
           ! So_RDS (road-deposited sediment), So_NH4DryDeposit, So_NO3DryDeposit (kg ha-1 day-1)
           ! So_NH4WetDeposit, So_NO3WetDeposit (kg m-3, in precipitation)   

        end do
     end do
     close (590005)

     !do itime = 1, domain%ntime_sedcnp   !BK20231024
     !do itime = 1, domain%ntime_sedcnp   !BK20231024
     ! calculate loading rates for each of the timesteps (kg)

         !smonth = trim(dateSedCNP%olddate(6:7))
         !read(smonth,*) imonth

     !  do ilc = 1, nlc
     !     So_LPOCLITinso(ilc,itime) = dLPOCLIT(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPOCRESinso(ilc,itime) = dLPOCRES(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPOCEXCinso(ilc,itime) = dLPOCEXC(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPOCMANinso(ilc,itime) = dLPOCMAN(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPONLITinso(ilc,itime) = dLPONLIT(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPONRESinso(ilc,itime) = dLPONRES(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPONEXCinso(ilc,itime) = dLPONEXC(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPONMANinso(ilc,itime) = dLPONMAN(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_NH4fert(ilc,itime) = dNH4fert(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)     ! !kg
     !     So_NO3fert(ilc,itime) = dNO3fert(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)     ! !kg
     !     So_LPOPLITinso(ilc,itime) = dLPOPLIT(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPOPRESinso(ilc,itime) = dLPOPRES(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPOPEXCinso(ilc,itime) = dLPOPEXC(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_LPOPMANinso(ilc,itime) = dLPOPMAN(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)   !kg
     !     So_PO4fert(ilc,itime) = dPO4fert(ilc,imonth) * (domain%areaxy / 10000.) * (DT / 86400.)     ! !kg

     !     ! *** reserved for future development ***
     !     ! So_RDS (road-deposited sediment), So_NH4DryDeposit, So_NO3DryDeposit (kg ha-1 day-1)
     !     ! So_NH4WetDeposit, So_NO3WetDeposit (kg m-3, in precipitation)
     !  end do  
     !end do 

  end subroutine ReadCNPinputs


   subroutine ReadCH_PNTSRC()   !BK20231016

      ! (to be called if itime = 1)
      ! to read in channel inputs from point sources (CH_PNTSRC.dat)
      ! loadings at t0 are put into the channel at t0 + 1 (one timestep later)

      !use CNPmain
      use module_SedCNPvariables
      use module_hydro_stop, only:HYDRO_stop  !BK20241023
! !--- BK20250707
! #ifdef MPP_LAND
!          use module_mpp_land, &
!                only: my_id, IO_id, mpp_land_bcast_int1, mpp_land_bcast_char
! #endif
! !--- BK20250707

      implicit none

      integer :: ich, ich_pnt, chid  !BK20241016
      integer :: year0, month0, day0, hour0, minute0, julian_day0, jh0, t0 
      integer :: year1, month1, day1, hour1, minute1, julian_day1, jh1, t1
      integer :: itime, ierr
      real(8) :: dQ0, dQ1
      real(8) :: dSED0(4), dLPOC0, dRPOC0, dLDOC0, dRDOC0, dLPON0, dRPON0, dLDON0, dRDON0, dNH40, dNO30, &
                 dLPOP0, dRPOP0, dLDOP0, dRDOP0, dPO40
      real(8) :: dSED1(4), dLPOC1, dRPOC1, dLDOC1, dRDOC1, dLPON1, dRPON1, dLDON1, dRDON1, dNH41, dNO31, &
                 dLPOP1, dRPOP1, dLDOP1, dRDOP1, dPO41
      character(len=630) :: str_line               !BK20250525
      character(len=30)  :: str_date, str_var(20)   !BK20250525
      character(len=1)   :: first_char   !BK20241015
      character(len=630), allocatable :: file_content(:)    !BK20250714
      integer            :: i, num_lines, line_idx          !BK20250714
      integer            :: khour                           !BK20250714

      logical :: file_exists  !BK20251128

      khour = SedCNPmodel%SedCNP_khour    !BK20250714

      if ((khour < 0) .and. (SedCNPmodel%SedCNP_kday < 0)) then
         write(*, '("FATAL ERROR: In ReadCH_PNTSRC() - "// &
               "Namelist error: Either KHOUR or KDAY must be defined.")')
         stop
      else if (( khour < 0 ) .and. (SedCNPmodel%SedCNP_kday > 0)) then
         khour = SedCNPmodel%SedCNP_kday * 24
      else if ((khour > 0) .and. (SedCNPmodel%SedCNP_kday > 0)) then
         write(*, '("WARNING: In ReadCH_PNTSRC() - Check Namelist: KHOUR and KDAY both defined.")')
         stop
      endif      

      !--- BK20250714
! #ifdef MPP_LAND
!     if (my_id == IO_id) then
! #endif
      !BK20251128 Check for user-specified file
      if (len_trim(SedCNPmodel%CH_PNTSRC_file) > 0) then
         inquire(file=trim(SedCNPmodel%CH_PNTSRC_file), exist=file_exists)
         if (file_exists) then
            open (unit=10, file=trim(SedCNPmodel%CH_PNTSRC_file), status='unknown', iostat=ierr)
         else
             ierr = 1 ! Flag to try defaults
         endif
      else
         ierr = 1 ! Flag to try defaults
      endif

      if (ierr /= 0) then
          !BK20251128 Check for default files
          inquire(file='./WHQIN/CH_PNTSRC.dat', exist=file_exists)
          if (file_exists) then
             open (unit=10, file='./WHQIN/CH_PNTSRC.dat', status='unknown', iostat=ierr)
          endif
      endif
      
      if (ierr /= 0) then
         call hydro_stop("Failed to open CH_PNTSRC data file")
      endif

      ! Count lines to allocate file_content
      num_lines = 0
      DO
         read(10, '(A)', iostat=ierr)
         if (ierr /= 0) exit
         num_lines = num_lines + 1
      ENDDO
      rewind(10)

      allocate(file_content(num_lines))
      do i = 1, num_lines
         read(10, '(A)') file_content(i)
      enddo
      close(10)
! #ifdef MPP_LAND
!     endif
!     call mpp_land_bcast_int1(num_lines)
!     if (my_id /= IO_id) then
!        allocate(file_content(num_lines))
!     endif
!     do i = 1, num_lines
!        call mpp_land_bcast_char(len(file_content(i)), file_content(i))
!     enddo
! #endif
      !--- BK20250714

      !--- initialisation      
      !water inflow (m3 s-1)
      !SedCNP_hydro%Qpnt = 0.
      !loadings (kg s-1)
      channelSed%Spnt = 0.
      Ch_LPOCpnt = 0.
      Ch_RPOCpnt = 0.
      Ch_LDOCpnt = 0.
      Ch_RDOCpnt = 0.
      Ch_LPONpnt = 0.
      Ch_RPONpnt = 0.
      Ch_LDONpnt = 0.
      Ch_RDONpnt = 0.
      Ch_NH4pnt = 0.
      Ch_NO3pnt = 0.
      Ch_LPOPpnt = 0.
      Ch_RPOPpnt = 0.
      Ch_LDOPpnt = 0.
      Ch_RDOPpnt = 0.
      Ch_PO4pnt = 0.
      str_date = ' '    !BK20250622
      str_var = ' '     !BK20250622
      ich_pnt = 0       !BK20250714

      !--- BK20250714
      !--- do loop
      ich_pnt = 0 ! 0 indicates not currently in a REACH block owned by this process
      DO line_idx = 1, num_lines    !BK20250720
         str_line = trim(adjustl(file_content(line_idx)))
         ! Skip blank lines and comments
         if (str_line == '') cycle
         if (str_line(1:1) == '!') cycle

         ! REACH block starts
         if (str_line(1:5) == 'REACH') then
            ! A new reach is starting. Finalize the previous one if it exists  !--- BK20250720
            if (ich_pnt > 0) then
               !do itime = t0, domain%ntime_sedcnp
               do itime = max(1, t0), domain%ntime_sedcnp  !BK20251204
                  channelSed%Spnt(:,ich_pnt,itime) =  dSED0(:) / 3600.
                  Ch_LPOCpnt(ich_pnt,itime) = dLPOC0 / 3600.
                  Ch_RPOCpnt(ich_pnt,itime) = dRPOC0 / 3600.
                  Ch_LDOCpnt(ich_pnt,itime) = dLDOC0 / 3600.
                  Ch_RDOCpnt(ich_pnt,itime) = dRDOC0 / 3600.
                  Ch_LPONpnt(ich_pnt,itime) = dLPON0 / 3600.
                  Ch_RPONpnt(ich_pnt,itime) = dRPON0 / 3600.
                  Ch_LDONpnt(ich_pnt,itime) = dLDON0 / 3600.
                  Ch_RDONpnt(ich_pnt,itime) = dRDON0 / 3600.
                  Ch_NH4pnt(ich_pnt,itime) = dNH40 / 3600.
                  Ch_NO3pnt(ich_pnt,itime) = dNO30 / 3600.
                  Ch_LPOPpnt(ich_pnt,itime) = dLPOP0 / 3600.
                  Ch_RPOPpnt(ich_pnt,itime) = dRPOP0 / 3600.
                  Ch_LDOPpnt(ich_pnt,itime) = dLDOP0 / 3600.
                  Ch_RDOPpnt(ich_pnt,itime) = dRDOP0 / 3600.
                  Ch_PO4pnt(ich_pnt,itime) = dPO40 / 3600.
               enddo
            endif    !--- BK20250720

            read(str_line(7:), *, iostat=ierr) chid   !BK20250720
            if (ierr /= 0) then
               call hydro_stop("Failed to read REACH ID (CH_PNTSRC)")
            endif

            ich_pnt = 0
            do ich = 1,domain%nch
               if(SedCNP_hydro%linkID(ich)==chid) then
                  ich_pnt = ich
                  exit
               endif
            enddo
       
            ! All processes must initialize for interpolation to reset state
            year0 = SedCNPmodel%SedCNP_start_year
            month0 = SedCNPmodel%SedCNP_start_month
            day0 = SedCNPmodel%SedCNP_start_day
            hour0 = SedCNPmodel%SedCNP_start_hour
            minute0 = SedCNPmodel%SedCNP_start_min
            t0 = 1
            dQ0 = 0.0
            dSED0 = 0.0
            dLPOC0 = 0.0
            dRPOC0 = 0.0
            dLDOC0 = 0.0
            dRDOC0 = 0.0
            dLPON0 = 0.0
            dRPON0 = 0.0 
            dLDON0 = 0.0
            dRDON0 = 0.0
            dNH40 = 0.0
            dNO30 = 0.0
            dLPOP0 = 0.0
            dRPOP0 = 0.0
            dLDOP0 = 0.0
            dRDOP0 = 0.0
            dPO40 = 0.0

         elseif (str_line(1:3) == 'END') then
            ! End of all data. Finalize the last reach block   !--- BK20250720
            if (ich_pnt > 0) then
               !do itime = t0, domain%ntime_sedcnp
               do itime = max(1, t0), domain%ntime_sedcnp !BK20251204               
                  channelSed%Spnt(:,ich_pnt,itime) =  dSED0(:) / 3600.
                  Ch_LPOCpnt(ich_pnt,itime) = dLPOC0 / 3600.
                  Ch_RPOCpnt(ich_pnt,itime) = dRPOC0 / 3600.
                  Ch_LDOCpnt(ich_pnt,itime) = dLDOC0 / 3600.
                  Ch_RDOCpnt(ich_pnt,itime) = dRDOC0 / 3600.
                  Ch_LPONpnt(ich_pnt,itime) = dLPON0 / 3600.
                  Ch_RPONpnt(ich_pnt,itime) = dRPON0 / 3600.
                  Ch_LDONpnt(ich_pnt,itime) = dLDON0 / 3600.
                  Ch_RDONpnt(ich_pnt,itime) = dRDON0 / 3600.
                  Ch_NH4pnt(ich_pnt,itime) = dNH40 / 3600.
                  Ch_NO3pnt(ich_pnt,itime) = dNO30 / 3600.
                  Ch_LPOPpnt(ich_pnt,itime) = dLPOP0 / 3600.
                  Ch_RPOPpnt(ich_pnt,itime) = dRPOP0 / 3600.
                  Ch_LDOPpnt(ich_pnt,itime) = dLDOP0 / 3600.
                  Ch_RDOPpnt(ich_pnt,itime) = dRDOP0 / 3600.
                  Ch_PO4pnt(ich_pnt,itime) = dPO40 / 3600.
               enddo
            endif       !--- BK20250720
            ich_pnt = 0 ! Mark as finalized to prevent double-counting
            exit

         else
            ! This is a data line. All processes parse the line and update state
            ! to stay synchronized. Only the owning process stores the data.
            read(str_line, *, iostat=ierr) str_date, str_var
            if (ierr > 0) then
               write(*,*) 'Format error in CH_PNTSRC file: ', trim(str_line)
               call hydro_stop("Failed to read input data (CH_PNTSRC)")
            endif

            str_date = adjustl(str_date)     !BK20241015
            read(str_date(1:4), *) year1
            read(str_date(6:7), *) month1
            read(str_date(9:10), *) day1
            read(str_date(12:13), *) hour1
            read(str_date(15:16), *) minute1
            read(str_var(1), *) dQ1          ! (m3 hr-1)
            if (dQ1 > 0) then
               read(str_var(2), *) dSED1(1)     ! (kg hr-1)
               read(str_var(3), *) dSED1(2)     ! (kg hr-1)
               read(str_var(4), *) dSED1(3)     ! (kg hr-1)
               read(str_var(5), *) dSED1(4)     ! (kg hr-1)
               read(str_var(6), *) dLPOC1       ! (kg hr-1)
               read(str_var(7), *) dRPOC1       ! (kg hr-1)
               read(str_var(8), *) dLDOC1       ! (kg hr-1)
               read(str_var(9), *) dRDOC1       ! (kg hr-1)
               read(str_var(10), *) dLPON1       ! (kg hr-1)
               read(str_var(11), *) dRPON1       ! (kg hr-1)
               read(str_var(12), *) dLDON1       ! (kg hr-1)
               read(str_var(13), *) dRDON1       ! (kg hr-1)
               read(str_var(14), *) dNH41         ! (kg hr-1)
               read(str_var(15), *) dNO31         ! (kg hr-1)
               read(str_var(16), *) dLPOP1       ! (kg hr-1)
               read(str_var(17), *) dRPOP1       ! (kg hr-1)
               read(str_var(18), *) dLDOP1       ! (kg hr-1)
               read(str_var(19), *) dRDOP1       ! (kg hr-1)
               read(str_var(20), *) dPO41         ! (kg hr-1)
            else !dQ1 <= 0.
               dSED1 = 0.0
               dLPOC1 = 0.0
               dRPOC1 = 0.0
               dLDOC1 = 0.0
               dRDOC1 = 0.0
               dLPON1 = 0.0
               dRPON1 = 0.0
               dLDON1 = 0.0
               dRDON1 = 0.0
               dNH41 = 0.0
               dNO31 = 0.0
               dLPOP1 = 0.0
               dRPOP1 = 0.0
               dLDOP1 = 0.0
               dRDOP1 = 0.0
               dPO41 = 0.0
            endif

            julian_day0 = calcJulianDay(year0,month0,day0)
            julian_day1 = calcJulianDay(year1,month1,day1)
            jh0 = (julian_day0 - 1) * 24 + hour0
            jh1 = (julian_day1 - 1) * 24 + hour1
            !t1 = t0 + nint((jh1 - jh0) * 3600. / real(domain%dt))
            if (t0 <= domain%ntime_sedcnp) then !---BK20251204
               t1 = t0 + nint((jh1 - jh0) * 3600. / real(domain%dt))
            else
               t1 = domain%ntime_sedcnp + 1
            endif  !---BK20251204

            if (ich_pnt > 0) then
               !do itime = t0, t1 - 1
               do itime = max(1, t0), t1 - 1 !BK20251204
                  if (itime > domain%ntime_sedcnp) exit
                  !water inflow (m3 s-1)
                  !  Qpnt is estimated in the subroutine drive_CHANNEL and
                  !  read-in in the subroutine SedCNPmodel_hydro_input
                  !loadings (kg s-1)
                  channelSed%Spnt(:,ich_pnt,itime) =  dSED0(:) / 3600.
                  Ch_LPOCpnt(ich_pnt,itime) = dLPOC0 / 3600.
                  Ch_RPOCpnt(ich_pnt,itime) = dRPOC0 / 3600.
                  Ch_LDOCpnt(ich_pnt,itime) = dLDOC0 / 3600.
                  Ch_RDOCpnt(ich_pnt,itime) = dRDOC0 / 3600.
                  Ch_LPONpnt(ich_pnt,itime) = dLPON0 / 3600.
                  Ch_RPONpnt(ich_pnt,itime) = dRPON0 / 3600.
                  Ch_LDONpnt(ich_pnt,itime) = dLDON0 / 3600.
                  Ch_RDONpnt(ich_pnt,itime) = dRDON0 / 3600.
                  Ch_NH4pnt(ich_pnt,itime) = dNH40 / 3600.
                  Ch_NO3pnt(ich_pnt,itime) = dNO30 / 3600.
                  Ch_LPOPpnt(ich_pnt,itime) = dLPOP0 / 3600.
                  Ch_RPOPpnt(ich_pnt,itime) = dRPOP0 / 3600.
                  Ch_LDOPpnt(ich_pnt,itime) = dLDOP0 / 3600.
                  Ch_RDOPpnt(ich_pnt,itime) = dRDOP0 / 3600.
                  Ch_PO4pnt(ich_pnt,itime) = dPO40 / 3600.
               enddo
            endif

            ! All processes must update state for next interpolation step
            year0 = year1
            month0 = month1
            day0 = day1
            hour0 = hour1
            t0 = t1
            dQ0 = dQ1
            dSED0 = dSED1
            dLPOC0 = dLPOC1
            dRPOC0 = dRPOC1
            dLDOC0 = dLDOC1
            dRDOC0 = dRDOC1
            dLPON0 = dLPON1
            dRPON0 = dRPON1
            dLDON0 = dLDON1
            dRDON0 = dRDON1
            dNH40 = dNH41
            dNO30 = dNO31
            dLPOP0 = dLPOP1
            dRPOP0 = dRPOP1
            dLDOP0 = dLDOP1
            dRDOP0 = dRDOP1
            dPO40 = dPO41
         endif
      ENDDO
      !--- BK20250714

      ! In case the file does not have an 'END' tag, finalize the last reach  !--- BK20250720
      if (ich_pnt > 0) then
         !do itime = t0, domain%ntime_sedcnp
         do itime = max(1, t0), domain%ntime_sedcnp  !BK20251204
            channelSed%Spnt(:,ich_pnt,itime) =  dSED0(:) / 3600.
            Ch_LPOCpnt(ich_pnt,itime) = dLPOC0 / 3600.
            Ch_RPOCpnt(ich_pnt,itime) = dRPOC0 / 3600.
            Ch_LDOCpnt(ich_pnt,itime) = dLDOC0 / 3600.
            Ch_RDOCpnt(ich_pnt,itime) = dRDOC0 / 3600.
            Ch_LPONpnt(ich_pnt,itime) = dLPON0 / 3600.
            Ch_RPONpnt(ich_pnt,itime) = dRPON0 / 3600.
            Ch_LDONpnt(ich_pnt,itime) = dLDON0 / 3600.
            Ch_RDONpnt(ich_pnt,itime) = dRDON0 / 3600.
            Ch_NH4pnt(ich_pnt,itime) = dNH40 / 3600.
            Ch_NO3pnt(ich_pnt,itime) = dNO30 / 3600.
            Ch_LPOPpnt(ich_pnt,itime) = dLPOP0 / 3600.
            Ch_RPOPpnt(ich_pnt,itime) = dRPOP0 / 3600.
            Ch_LDOPpnt(ich_pnt,itime) = dLDOP0 / 3600.
            Ch_RDOPpnt(ich_pnt,itime) = dRDOP0 / 3600.
            Ch_PO4pnt(ich_pnt,itime) = dPO40 / 3600.
         enddo
      endif    !--- BK20250720

      if (allocated(file_content)) deallocate(file_content)  !BK20250707

   end subroutine ReadCH_PNTSRC


   subroutine ReadCH_ABSDIS()   !BK20240703 !BK20241018

      ! (to be called if itime = 1)
      ! to read in channel abstraction/discharge (CH_ABSDIS.dat)
      ! water abstraction in negative (-) value
      ! water discharge in positive (+) value
      ! loadings at t0 are put into the channel at t0 + 1 (one timestep later)

      !use CNPmain
      use module_SedCNPvariables
      use module_hydro_stop, only:HYDRO_stop  !BK20241023
! !--- BK20250707
! #ifdef MPP_LAND
!          use module_mpp_land, &
!                only: my_id, IO_id, mpp_land_bcast_int1, mpp_land_bcast_char
! #endif
! !--- BK20250707

      implicit none

      integer :: ich, ich_absdis, chid
      integer :: year0, month0, day0, hour0, minute0, julian_day0, jh0, t0 
      integer :: year1, month1, day1, hour1, minute1, julian_day1, jh1, t1
      integer :: itime, ierr
      !integer :: lag_ts    ! lag-time for discharge (number of time-steps)   !BK20250622
      !real(8) :: lag_hr0, lag_hr1     ! lag-time for discharge (hr)          !BK20250622
      real(8) :: dQ0, dQ1  
      real(8) :: dSED0(4), dLPOC0, dRPOC0, dLDOC0, dRDOC0, dLPON0, dRPON0, dLDON0, dRDON0, dNH40, dNO30, &
                 dLPOP0, dRPOP0, dLDOP0, dRDOP0, dPO40
      real(8) :: dSED1(4), dLPOC1, dRPOC1, dLDOC1, dRDOC1, dLPON1, dRPON1, dLDON1, dRDON1, dNH41, dNO31, &
                 dLPOP1, dRPOP1, dLDOP1, dRDOP1, dPO41
      character(len=630) :: str_line
      character(len=30)  :: str_date, str_var(20)
      character(len=1)   :: first_char   !BK20241015
      character(len=630), allocatable :: file_content(:)    !BK20250714
      integer            :: i, num_lines, line_idx          !BK20250714
      integer            :: khour                           !BK20250714

      logical :: file_exists  !BK20251128

      khour = SedCNPmodel%SedCNP_khour    !BK20250714

      if ((khour < 0) .and. (SedCNPmodel%SedCNP_kday < 0)) then
         write(*, '("FATAL ERROR: In ReadCH_ABSDIS() - "// &
               "Namelist error: Either KHOUR or KDAY must be defined.")')
         stop
      else if (( khour < 0 ) .and. (SedCNPmodel%SedCNP_kday > 0)) then
         khour = SedCNPmodel%SedCNP_kday * 24
      else if ((khour > 0) .and. (SedCNPmodel%SedCNP_kday > 0)) then
         write(*, '("WARNING: In ReadCH_ABSDIS() - Check Namelist: KHOUR and KDAY both defined.")')
         stop
      endif

      !--- BK20250714
! #ifdef MPP_LAND
!     if (my_id == IO_id) then
! #endif
      !BK20251128 Check for user-specified file
      if (len_trim(SedCNPmodel%CH_ABSDIS_file) > 0) then
         inquire(file=trim(SedCNPmodel%CH_ABSDIS_file), exist=file_exists)
         if (file_exists) then
            open (unit=10, file=trim(SedCNPmodel%CH_ABSDIS_file), status='unknown', iostat=ierr)
         else
             ierr = 1 ! Flag to try defaults
         endif
      else
         ierr = 1 ! Flag to try defaults
      endif

      if (ierr /= 0) then
          !BK20251128 Check for default files
          inquire(file='./WHQIN/CH_ABSDIS.dat', exist=file_exists)
          if (file_exists) then
             open (unit=10, file='./WHQIN/CH_ABSDIS.dat', status='unknown', iostat=ierr)
          endif
      endif

      if (ierr /= 0) then
         call hydro_stop("Failed to open CH_ABSDIS data file")
      endif

      ! Count lines to allocate file_content
      num_lines = 0
      DO
         read(10, '(A)', iostat=ierr)
         if (ierr /= 0) exit
         num_lines = num_lines + 1
      ENDDO
      rewind(10)

      allocate(file_content(num_lines))
      do i = 1, num_lines
         read(10, '(A)') file_content(i)
      enddo
      close(10)
! #ifdef MPP_LAND
!     endif
!     call mpp_land_bcast_int1(num_lines)
!     if (my_id /= IO_id) then
!        allocate(file_content(num_lines))
!     endif
!     do i = 1, num_lines
!        call mpp_land_bcast_char(len(file_content(i)), file_content(i))
!     enddo
! #endif
      !--- BK20250714

      !--- initialisation
      !look-up table
      !abs2dis_ich = 0  !reach_index to which abstrated water is discharged  !BK20250622
      !abs2dis_lag = 0  !lag time for discharge, in number of timesteps      !BK20250622

      !water abstraction, discharge (m3 s-1)
      !SedCNP_hydro%Qabs = 0.
      !SedCNP_hydro%Qdis = 0.

      !abstraction loadings (kg s-1)
      channelSed%Sabs = 0.
      Ch_LPOCabs = 0.
      Ch_RPOCabs = 0.
      Ch_LDOCabs = 0.
      Ch_RDOCabs = 0.
      Ch_LPONabs = 0.
      Ch_RPONabs = 0.
      Ch_LDONabs = 0.
      Ch_RDONabs = 0.
      Ch_NH4abs = 0.
      Ch_NO3abs = 0.
      Ch_LPOPabs = 0.
      Ch_RPOPabs = 0.
      Ch_LDOPabs = 0.
      Ch_RDOPabs = 0.
      Ch_PO4abs = 0.
      
      !discharge loadings (kg s-1)
      channelSed%Sdis = 0.
      Ch_LPOCdis = 0.
      Ch_RPOCdis = 0.
      Ch_LDOCdis = 0.
      Ch_RDOCdis = 0.
      Ch_LPONdis = 0.
      Ch_RPONdis = 0.
      Ch_LDONdis = 0.
      Ch_RDONdis = 0.
      Ch_NH4dis = 0.
      Ch_NO3dis = 0.
      Ch_LPOPdis = 0.
      Ch_RPOPdis = 0.
      Ch_LDOPdis = 0.
      Ch_RDOPdis = 0.
      Ch_PO4dis = 0.

      str_date = ' '    !BK20250622
      str_date = ' '    !BK20250622
      str_var = ' '     !BK20250622

      !--- BK20250714
      !--- do loop
      ich_absdis = 0 ! 0 indicates not currently in a REACH block owned by this process
      DO line_idx = 1, num_lines    !BK20250720
         str_line = trim(adjustl(file_content(line_idx)))
         ! Skip blank lines and comments
         if (str_line == '') cycle
         if (str_line(1:1) == '!') cycle

         ! REACH block starts
         if (str_line(1:5) == 'REACH') then
            ! A new reach is starting. Finalize the previous one if it exists  !--- BK20250720
            if (ich_absdis > 0) then
               if (dQ0 <= 0.) then ! Abstraction 
               !    do itime = max(1, t0), domain%ntime_sedcnp
               !    enddo
               else ! (dQ0 > 0.) ! Discharge
                  !do itime = t0, domain%ntime_sedcnp
                  do itime = max(1, t0), domain%ntime_sedcnp  !BK20251204
                     channelSed%Sdis(:,ich_absdis,itime) = channelSed%Sdis(:,ich_absdis,itime) + dSED0(:) / 3600.
                     Ch_LPOCdis(ich_absdis,itime) = Ch_LPOCdis(ich_absdis,itime) + dLPOC0 / 3600.
                     Ch_RPOCdis(ich_absdis,itime) = Ch_RPOCdis(ich_absdis,itime) + dRPOC0 / 3600.
                     Ch_LDOCdis(ich_absdis,itime) = Ch_LDOCdis(ich_absdis,itime) + dLDOC0 / 3600.
                     Ch_RDOCdis(ich_absdis,itime) = Ch_RDOCdis(ich_absdis,itime) + dRDOC0 / 3600.
                     Ch_LPONdis(ich_absdis,itime) = Ch_LPONdis(ich_absdis,itime) + dLPON0 / 3600.
                     Ch_RPONdis(ich_absdis,itime) = Ch_RPONdis(ich_absdis,itime) + dRPON0 / 3600.
                     Ch_LDONdis(ich_absdis,itime) = Ch_LDONdis(ich_absdis,itime) + dLDON0 / 3600.
                     Ch_RDONdis(ich_absdis,itime) = Ch_RDONdis(ich_absdis,itime) + dRDON0 / 3600.
                     Ch_NH4dis(ich_absdis,itime) = Ch_NH4dis(ich_absdis,itime) + dNH40 / 3600.
                     Ch_NO3dis(ich_absdis,itime) = Ch_NO3dis(ich_absdis,itime) + dNO30 / 3600.
                     Ch_LPOPdis(ich_absdis,itime) = Ch_LPOPdis(ich_absdis,itime) + dLPOP0 / 3600.
                     Ch_RPOPdis(ich_absdis,itime) = Ch_RPOPdis(ich_absdis,itime) + dRPOP0 / 3600.
                     Ch_LDOPdis(ich_absdis,itime) = Ch_LDOPdis(ich_absdis,itime) + dLDOP0 / 3600.
                     Ch_RDOPdis(ich_absdis,itime) = Ch_RDOPdis(ich_absdis,itime) + dRDOP0 / 3600.
                     Ch_PO4dis(ich_absdis,itime) = Ch_PO4dis(ich_absdis,itime) + dPO40 / 3600.
                  enddo
               endif
            endif       !--- BK20250720

            read(str_line(7:), *, iostat=ierr) chid      !BK20250720
            if (ierr /= 0) then
               call hydro_stop("Failed to read REACH ID (CH_ABSDIS)")
            endif

            ich_absdis = 0
            do ich = 1,domain%nch
               if(SedCNP_hydro%linkID(ich)==chid) then
                  ich_absdis = ich
                  exit
               endif
            enddo

            ! All processes must initialize for interpolation to reset state
            year0 = SedCNPmodel%SedCNP_start_year
            month0 = SedCNPmodel%SedCNP_start_month
            day0 = SedCNPmodel%SedCNP_start_day
            hour0 = SedCNPmodel%SedCNP_start_hour
            minute0 = SedCNPmodel%SedCNP_start_min
            t0 = 1
            dQ0 = 0.0
            dSED0 = 0.0
            dLPOC0 = 0.0
            dRPOC0 = 0.0
            dLDOC0 = 0.0
            dRDOC0 = 0.0
            dLPON0 = 0.0
            dRPON0 = 0.0
            dLDON0 = 0.0
            dRDON0 = 0.0
            dNH40 = 0.0
            dNO30 = 0.0
            dLPOP0 = 0.0
            dRPOP0 = 0.0
            dLDOP0 = 0.0
            dRDOP0 = 0.0
            dPO40 = 0.0
            !chid_dis0 = 0    !BK20250622
            !lag_hr0 = 0.     !BK20250622

         elseif (str_line(1:3) == 'END') then
            ! End of all data. Finalize the last reach block.  !--- BK20250720
            if (ich_absdis > 0) then
               if (dQ0 <= 0.) then ! Abstraction 
               !    do itime = max(1, t0), domain%ntime_sedcnp
               !    enddo
               else ! (dQ0 > 0.) ! Discharge
                  !do itime = t0, domain%ntime_sedcnp
                  do itime = max(1, t0), domain%ntime_sedcnp  !BK20251204
                     channelSed%Sdis(:,ich_absdis,itime) = channelSed%Sdis(:,ich_absdis,itime) + dSED0(:) / 3600.
                     Ch_LPOCdis(ich_absdis,itime) = Ch_LPOCdis(ich_absdis,itime) + dLPOC0 / 3600.
                     Ch_RPOCdis(ich_absdis,itime) = Ch_RPOCdis(ich_absdis,itime) + dRPOC0 / 3600.
                     Ch_LDOCdis(ich_absdis,itime) = Ch_LDOCdis(ich_absdis,itime) + dLDOC0 / 3600.
                     Ch_RDOCdis(ich_absdis,itime) = Ch_RDOCdis(ich_absdis,itime) + dRDOC0 / 3600.
                     Ch_LPONdis(ich_absdis,itime) = Ch_LPONdis(ich_absdis,itime) + dLPON0 / 3600.
                     Ch_RPONdis(ich_absdis,itime) = Ch_RPONdis(ich_absdis,itime) + dRPON0 / 3600.
                     Ch_LDONdis(ich_absdis,itime) = Ch_LDONdis(ich_absdis,itime) + dLDON0 / 3600.
                     Ch_RDONdis(ich_absdis,itime) = Ch_RDONdis(ich_absdis,itime) + dRDON0 / 3600.
                     Ch_NH4dis(ich_absdis,itime) = Ch_NH4dis(ich_absdis,itime) + dNH40 / 3600.
                     Ch_NO3dis(ich_absdis,itime) = Ch_NO3dis(ich_absdis,itime) + dNO30 / 3600.
                     Ch_LPOPdis(ich_absdis,itime) = Ch_LPOPdis(ich_absdis,itime) + dLPOP0 / 3600.
                     Ch_RPOPdis(ich_absdis,itime) = Ch_RPOPdis(ich_absdis,itime) + dRPOP0 / 3600.
                     Ch_LDOPdis(ich_absdis,itime) = Ch_LDOPdis(ich_absdis,itime) + dLDOP0 / 3600.
                     Ch_RDOPdis(ich_absdis,itime) = Ch_RDOPdis(ich_absdis,itime) + dRDOP0 / 3600.
                     Ch_PO4dis(ich_absdis,itime) = Ch_PO4dis(ich_absdis,itime) + dPO40 / 3600.
                  enddo
               endif
            endif    !--- BK20250720
            ich_absdis = 0 ! Mark as finalized to prevent double-counting
            exit     

         else
            ! This is a data line. All processes parse the line and update state
            ! to stay synchronized. Only the owning process stores the data.
            read(str_line, *, iostat=ierr) str_date, str_var
            if (ierr > 0) then
               write(*,*) 'Format error in CH_ABSDIS file: ', trim(str_line)
               call hydro_stop("Failed to read input data (CH_ABSDIS)")
            endif

            str_date = adjustl(str_date)     !BK20241015
            read(str_date(1:4), *) year1
            read(str_date(6:7), *) month1
            read(str_date(9:10), *) day1
            read(str_date(12:13), *) hour1
            read(str_date(15:16), *) minute1
            read(str_var(1), *) dQ1          ! (m3 hr-1)
            if (dQ1 > 0.) then
               read(str_var(2), *) dSED1(1)     ! (kg hr-1)
               read(str_var(3), *) dSED1(2)     ! (kg hr-1)
               read(str_var(4), *) dSED1(3)     ! (kg hr-1)
               read(str_var(5), *) dSED1(4)     ! (kg hr-1)
               read(str_var(6), *) dLPOC1       ! (kg hr-1)
               read(str_var(7), *) dRPOC1       ! (kg hr-1)
               read(str_var(8), *) dLDOC1       ! (kg hr-1)
               read(str_var(9), *) dRDOC1       ! (kg hr-1)
               read(str_var(10), *) dLPON1       ! (kg hr-1)
               read(str_var(11), *) dRPON1       ! (kg hr-1)
               read(str_var(12), *) dLDON1       ! (kg hr-1)
               read(str_var(13), *) dRDON1       ! (kg hr-1)
               read(str_var(14), *) dNH41         ! (kg hr-1)
               read(str_var(15), *) dNO31         ! (kg hr-1)
               read(str_var(16), *) dLPOP1       ! (kg hr-1)
               read(str_var(17), *) dRPOP1       ! (kg hr-1)
               read(str_var(18), *) dLDOP1       ! (kg hr-1)
               read(str_var(19), *) dRDOP1       ! (kg hr-1)
               read(str_var(20), *) dPO41         ! (kg hr-1)
            else !dQ1 <= 0.
               dSED1 = 0.0
               dLPOC1 = 0.0
               dRPOC1 = 0.0
               dLDOC1 = 0.0
               dRDOC1 = 0.0
               dLPON1 = 0.0
               dRPON1 = 0.0
               dLDON1 = 0.0
               dRDON1 = 0.0
               dNH41 = 0.0
               dNO31 = 0.0
               dLPOP1 = 0.0
               dRPOP1 = 0.0
               dLDOP1 = 0.0
               dRDOP1 = 0.0
               dPO41 = 0.0
            endif

            julian_day0 = calcJulianDay(year0,month0,day0)
            julian_day1 = calcJulianDay(year1,month1,day1)
            jh0 = (julian_day0 - 1) * 24 + hour0
            jh1 = (julian_day1 - 1) * 24 + hour1
            !t1 = t0 + nint((jh1 - jh0) * 3600. / real(domain%dt))
            if (t0 <= domain%ntime_sedcnp) then  !---BK20251204
               t1 = t0 + nint((jh1 - jh0) * 3600. / real(domain%dt))
            else
               t1 = domain%ntime_sedcnp + 1
            endif  !---BK20251204

            if (ich_absdis > 0) then
               if (dQ0 <= 0.) then  ! Abstraction
                  ! do itime = t0, t1 - 1
                  !    !if (itime > domain%ntime_sedcnp) exit
                  !    !actual amount of Qabs is estimated in the subroutine drive_CHANNEL and
                  !    !  read-in in the subroutine SedCNPmodel_hydro_input
                  !    !actual amount of abstracted loadings is to be estimated in the CNP subroutines
                  !    !  taking account of Qabs and availability of loadings
                  ! enddo
                  !
                  !--- BK20250622 commented out (to be restored later if mpi can be applied to time-lag simulations)
                  ! if (chid_dis0 .gt. 0) then !--- discharge (abstraction from internal source)
                  !    ich_absdis = 0
                  !    do ich = 1,domain%nch
                  !       if(SedCNP_hydro%linkID(ich).eq.chid_dis0) then
                  !          ich_absdis = ich
                  !          exit
                  !       endif
                  !    enddo
                  !    lag_ts = max(1, int(lag_hr0 * 3600. / domain%DT))

                  !    do itime = t0 + 1, t1
                  !       abs2dis_ich(ich_absdis,itime) = ich_dis
                  !       abs2dis_lag(ich_absdis,itime) = lag_ts
                  !    enddo
                  !    !
                  !    !actual abstraction water/loadings to be estimated
                  !    !in CNP subroutines considering the availability of water/loadings
                  !    !
                  ! endif
                  !--- BK20250622 commented out
               else ! (dQ0 > 0.) ! Discharge
                  !do itime = t0, t1 - 1
                  do itime = max(1, t0), t1 - 1  !BK20251204
                     if (itime > domain%ntime_sedcnp) exit
                     !water inflow (m3 s-1)
                     !  Qdis is estimated in the subroutine drive_CHANNEL and
                     !  read-in in the subroutine SedCNPmodel_hydro_input
                     !loadings (kg s-1)
                     channelSed%Sdis(:,ich_absdis,itime) = channelSed%Sdis(:,ich_absdis,itime) + dSED0(:) / 3600.
                     Ch_LPOCdis(ich_absdis,itime) = Ch_LPOCdis(ich_absdis,itime) + dLPOC0 / 3600.
                     Ch_RPOCdis(ich_absdis,itime) = Ch_RPOCdis(ich_absdis,itime) + dRPOC0 / 3600.
                     Ch_LDOCdis(ich_absdis,itime) = Ch_LDOCdis(ich_absdis,itime) + dLDOC0 / 3600.
                     Ch_RDOCdis(ich_absdis,itime) = Ch_RDOCdis(ich_absdis,itime) + dRDOC0 / 3600.
                     Ch_LPONdis(ich_absdis,itime) = Ch_LPONdis(ich_absdis,itime) + dLPON0 / 3600.
                     Ch_RPONdis(ich_absdis,itime) = Ch_RPONdis(ich_absdis,itime) + dRPON0 / 3600.
                     Ch_LDONdis(ich_absdis,itime) = Ch_LDONdis(ich_absdis,itime) + dLDON0 / 3600.
                     Ch_RDONdis(ich_absdis,itime) = Ch_RDONdis(ich_absdis,itime) + dRDON0 / 3600.
                     Ch_NH4dis(ich_absdis,itime) = Ch_NH4dis(ich_absdis,itime) + dNH40 / 3600.
                     Ch_NO3dis(ich_absdis,itime) = Ch_NO3dis(ich_absdis,itime) + dNO30 / 3600.
                     Ch_LPOPdis(ich_absdis,itime) = Ch_LPOPdis(ich_absdis,itime) + dLPOP0 / 3600.
                     Ch_RPOPdis(ich_absdis,itime) = Ch_RPOPdis(ich_absdis,itime) + dRPOP0 / 3600.
                     Ch_LDOPdis(ich_absdis,itime) = Ch_LDOPdis(ich_absdis,itime) + dLDOP0 / 3600.
                     Ch_RDOPdis(ich_absdis,itime) = Ch_RDOPdis(ich_absdis,itime) + dRDOP0 / 3600.
                     Ch_PO4dis(ich_absdis,itime) = Ch_PO4dis(ich_absdis,itime) + dPO40 / 3600.
                  enddo
               endif
            endif

            ! All processes must update state for next interpolation step
            year0 = year1
            month0 = month1
            day0 = day1
            hour0 = hour1
            jh0 = jh1
            t0 = t1
            dQ0 = dQ1
            dSED0 = dSED1
            dLPOC0 = dLPOC1
            dRPOC0 = dRPOC1
            dLDOC0 = dLDOC1
            dRDOC0 = dRDOC1
            dLPON0 = dLPON1
            dRPON0 = dRPON1
            dLDON0 = dLDON1
            dRDON0 = dRDON1
            dNH40 = dNH41
            dNO30 = dNO31
            dLPOP0 = dLPOP1
            dRPOP0 = dRPOP1
            dLDOP0 = dLDOP1
            dRDOP0 = dRDOP1
            dPO40 = dPO41
            !--- BK20250622 commented out (to be restored later if mpi can be applied to time-lag simulations)
            !chid_dis0 = chid_dis1
            !lag_hr0 = lag_hr1 
            !--- BK20250622 commented out
         endif  
      ENDDO
            
      ! In case the file does not have an 'END' tag, finalize the last reach     !--- BK20250720
      if (ich_absdis > 0) then
         if (dQ0 <= 0.) then ! Abstraction 
         !    do itime = max(1, t0), domain%ntime_sedcnp
         !    enddo
         else ! (dQ0 > 0.) ! Discharge
            !do itime = t0, domain%ntime_sedcnp
            do itime = max(1, t0), domain%ntime_sedcnp  !BK20251204
               channelSed%Sdis(:,ich_absdis,itime) = channelSed%Sdis(:,ich_absdis,itime) + dSED0(:) / 3600.
               Ch_LPOCdis(ich_absdis,itime) = Ch_LPOCdis(ich_absdis,itime) + dLPOC0 / 3600.
               Ch_RPOCdis(ich_absdis,itime) = Ch_RPOCdis(ich_absdis,itime) + dRPOC0 / 3600.
               Ch_LDOCdis(ich_absdis,itime) = Ch_LDOCdis(ich_absdis,itime) + dLDOC0 / 3600.
               Ch_RDOCdis(ich_absdis,itime) = Ch_RDOCdis(ich_absdis,itime) + dRDOC0 / 3600.
               Ch_LPONdis(ich_absdis,itime) = Ch_LPONdis(ich_absdis,itime) + dLPON0 / 3600.
               Ch_RPONdis(ich_absdis,itime) = Ch_RPONdis(ich_absdis,itime) + dRPON0 / 3600.
               Ch_LDONdis(ich_absdis,itime) = Ch_LDONdis(ich_absdis,itime) + dLDON0 / 3600.
               Ch_RDONdis(ich_absdis,itime) = Ch_RDONdis(ich_absdis,itime) + dRDON0 / 3600.
               Ch_NH4dis(ich_absdis,itime) = Ch_NH4dis(ich_absdis,itime) + dNH40 / 3600.
               Ch_NO3dis(ich_absdis,itime) = Ch_NO3dis(ich_absdis,itime) + dNO30 / 3600.
               Ch_LPOPdis(ich_absdis,itime) = Ch_LPOPdis(ich_absdis,itime) + dLPOP0 / 3600.
               Ch_RPOPdis(ich_absdis,itime) = Ch_RPOPdis(ich_absdis,itime) + dRPOP0 / 3600.
               Ch_LDOPdis(ich_absdis,itime) = Ch_LDOPdis(ich_absdis,itime) + dLDOP0 / 3600.
               Ch_RDOPdis(ich_absdis,itime) = Ch_RDOPdis(ich_absdis,itime) + dRDOP0 / 3600.
               Ch_PO4dis(ich_absdis,itime) = Ch_PO4dis(ich_absdis,itime) + dPO40 / 3600.
            enddo
         endif
      endif       !--- BK20250720

      if (allocated(file_content)) deallocate(file_content)   !BK20250707

   end subroutine ReadCH_ABSDIS


!   subroutine ReadSBID_CH()    !BK20231011 !BK20240509 !BK20240515 commented out
!
!      !to read in subbasin ID for channel IDs from a text file (SBIDCH_TEST.dat)
!
!      !use config_base, only: SedCNPmodel
!
!      implicit none
! 
!      integer             :: chid   !channel id
!      integer             :: sbid   !subbasin id
!
!      open (51011, file='./WHQIN/SBIDCH_TEST.dat', status='unknown')
!      read (51011,*)
!      
!      DO
!         read(51011, *, end=999) chid, sbid
!         if (chid .gt. 0) then
!            if (sbid .gt. 0) then
!               domain%sbid_ch(chid) = sbid
!            else
!               domain%sbid_ch(chid) = 0  ! sbid set to zero
!               write(*, '(A, I5, A)') 'WARNING: Subbasin ID set to 0 for Channel-ID ', chid, '(subroutine ReadSBID_CH)'  !BK20240513
!            endif
!         else 
!            write(*, '("FATAL ERROR: Channel ID must be greater than zero (subroutine ReadSBID_CH)")')
!            stop
!         endif
!      END DO
! 999  continue
!      close (51011)
!
!   end subroutine ReadSBID_CH   !BK20231011


   subroutine ReadSOILPSF()    !BK20240205

     !to read in soil particle size fractions for soil texture classes from a text file (SOILPSF_TEST.dat)
     
     use module_hydro_stop, only:HYDRO_stop  !BK20251130

     implicit none

     integer  :: st           !soil texture code
     real     :: psf(nps)     !soil particle size fractions     
     integer  :: ips          !particle size index
     character(len=20) :: description
     
     !BK20251128 
     logical :: file_exists
     integer :: ierr

      !initialisation
      overSed%SOILPSF = 0.0  !BK20240517

      if (len_trim(SedCNPmodel%SOILPSF_file) > 0) then
         inquire(file=trim(SedCNPmodel%SOILPSF_file), exist=file_exists)
         if (file_exists) then
            open (2402051, file=trim(SedCNPmodel%SOILPSF_file), status='unknown', iostat=ierr)
         else
             ierr = 1 ! Flag to try defaults
         endif
      else
         ierr = 1 ! Flag to try defaults
      endif

      if (ierr /= 0) then
          !BK20251128 Check for default files
          inquire(file='./WHQIN/SOILPSF.dat', exist=file_exists)
          if (file_exists) then
             open (2402051, file='./WHQIN/SOILPSF.dat', status='unknown', iostat=ierr)
         endif
      endif

      if (ierr /= 0) then
         call hydro_stop("Failed to open SOILPSF data file")
      endif
      read (2402051,*)
      DO
         !read(2402051, *, end=999) st, psf(1:nps), description  !BK20240517
         ! do ips = 1,nps
         !    if (st .gt. 0) then
         !       overSed%SOILPSF(st,:) = psf(ips)
         !    else
         !       write(*, '("FATAL ERROR: Soil texture code must be greater than zero")')
         !       stop
         !    endif
         ! enddo
         !--- BK20250720
         read(2402051, *, end=999) st, psf(1:nps), description
         if (st > 0 .and. st <= size(overSed%SOILPSF, dim=1)) then
            overSed%SOILPSF(st, 1:nps) = psf(1:nps)
         else if (st <= 0) then
            write(*, '("FATAL ERROR: Soil texture code must be greater than zero")')
            stop
         else ! st is out of bounds
            write(*, '("WARNING: Soil texture code ", I0, " is out of bounds, skipping.")') st
         endif
         !--- BK20250720
      ENDDO
 999  continue
      close (2402051)

   end subroutine ReadSOILPSF   !BK20240205


  function calcJulianDay(year, month, day) result(julian_day)

    integer, intent(in) :: year, month, day
    integer :: julian_day
    integer :: a, y, m

    a = (14 - month) / 12
    y = year + 4800 - a
    m = month + 12 * a - 3

    julian_day = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045

  end function calcJulianDay

end module module_SedCNP_in
