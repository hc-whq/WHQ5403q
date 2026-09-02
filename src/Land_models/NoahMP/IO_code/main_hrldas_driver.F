program noah_hrldas_driver
! this is the main program to drive HRLDAS-Noah, HRLDAS-NoahMP, and other Land models.

#ifdef Noah1d
! this is used to drive Noah1d
  use module_noah1d_hrldas_driver, only: land_driver_ini, land_driver_exe
#else
  ! this is used to drive NoahMP
  use module_noahmp_hrldas_driver, only: land_driver_ini, land_driver_exe
#endif

!=====||___WHQ___||=====!
!
#ifdef MPP_LAND
  use module_mpp_land, only: mpp_land_sync
#endif
!=====||__WHQ5403q__||=====!
!
  use module_SedCNPmodel_driver
  ! SedCNPmodel is already visible via module_SedCNPmodel_driver/orchestrator_base below;
  ! a redundant direct "use config_base" here triggers a gfortran diamond-import
  ! derived-type mismatch when combined with those.
!
!=====||__WHQ5403q__||=====!
!
!=====||___WHQ___||=====!

  use module_HYDRO_drv, only: HYDRO_finish

  ! NEW MODULE WITH DT INCLUDING STATE, PARAMETER,FORCING AND GEOMETRY. Different modules
  use state_module, only: state_type
  use parameters
  use forcing
  use geometry
  use orchestrator_base

  implicit none
  integer :: ITIME, NTIME
  type(state_type) :: state

  call orchestrator%init()

!=====||___WHQ___||=====!
!
   if(SedCNPmodel%SedCNP_option == 0) then
!
!=====||___WHQ___||=====!

      call land_driver_ini(NTIME, state)

#ifdef WRF_HYDRO
      do ITIME = 1, NTIME
#else
      do ITIME = 0, NTIME
#endif
         ! DTs must be passed to the exe.
         call land_driver_exe(ITIME, state)
      end do

!=====||___WHQ___||=====!
!
   else
!=====||__WHQ5403q__||=====!
!
      call SedCNP_driver_ini(NTIME)
      call SedCNPmodel_surface_input()

#ifdef WRF_HYDRO
      do ITIME = 1, NTIME
#else
      do ITIME = 0, NTIME
#endif
         call SedCNP_driver_exe(ITIME)
      end do
!
!=====||__WHQ5403q__||=====!
   endif
!
!=====||___WHQ___||=====!

#ifdef WRF_HYDRO
  !Pass state for destruction
  call hydro_finish()
#endif

end program noah_hrldas_driver
