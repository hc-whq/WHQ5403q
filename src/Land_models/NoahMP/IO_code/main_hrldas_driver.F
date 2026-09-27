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
  !use module_mpp_land, only: my_id, io_id   !WHQ5403 progress bar on the I/O rank
#endif
  !use ISO_FORTRAN_ENV, only: ERROR_UNIT, OUTPUT_UNIT   !WHQ5403
  use module_progress, only: progress_show   !WHQ5403 progress bar (module_progress)
!=====||__WHQ5403__||=====!
!
  use module_SedCNPmodel_driver
  ! SedCNPmodel is already visible via module_SedCNPmodel_driver/orchestrator_base below;
  ! a redundant direct "use config_base" here triggers a gfortran diamond-import
  ! derived-type mismatch when combined with those.
  !use config_base, only: SedCNPmodel
!
!=====||__WHQ5403__||=====!
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
         !call progress_clear()   !WHQ5403 messages of the time-step first, progress bar last
         call land_driver_exe(ITIME, state)
         call progress_show(ITIME, NTIME, "water-only")   !WHQ5403
      end do

!=====||___WHQ___||=====!
!
   else
!=====||__WHQ5403__||=====!
!
      call SedCNP_driver_ini(NTIME)
      call SedCNPmodel_surface_input()

#ifdef WRF_HYDRO
      do ITIME = 1, NTIME
#else
      do ITIME = 0, NTIME
#endif
         !call progress_clear()   !WHQ5403 messages of the time-step first, progress bar last
         call SedCNP_driver_exe(ITIME)
         call progress_show(ITIME, NTIME, "water quality")   !WHQ5403
      end do
!
!=====||__WHQ5403__||=====!
   endif
!
!=====||___WHQ___||=====!

#ifdef WRF_HYDRO
  !Pass state for destruction
  call hydro_finish()
#endif

!=====||__WHQ5403__||=====!
!
! progress bar moved to module_progress (SedCNPmodel), which also erases it before messages
!contains

!  ! Progress bar on stderr, printed at the end of each time-step after the messages of the
!  ! time-step (INFO, WARNING, ...), and erased at the start of the next time-step so that
!  ! those messages start on a clean line.
!  subroutine progress_clear()
!    if (my_id /= io_id) return
!    flush(OUTPUT_UNIT)
!    write(ERROR_UNIT, '(A)', advance='no') char(13)//repeat(' ', 80)//char(13)
!    flush(ERROR_UNIT)
!  end subroutine progress_clear

!  subroutine progress_show(it, nt, label)
!    integer,          intent(in) :: it, nt
!    character(len=*), intent(in) :: label
!    if (my_id /= io_id) return
!    flush(OUTPUT_UNIT)
!    write(ERROR_UNIT, '(A,F6.2,A,I0,A,I0,A)', advance='no') " WRF-HydroQual running ("//label//") ... ", &
!         real(it) / real(max(nt,1)) * 100.0, "% (", it, " / ", nt, ")"
!    if (it == nt) then
!       ! On the final timestep, print a newline, then the success message.
!       write(ERROR_UNIT,*)
!       write(ERROR_UNIT,*) "The model finished successfully"
!    endif
!    flush(ERROR_UNIT)
!  end subroutine progress_show
!
!=====||__WHQ5403__||=====!

end program noah_hrldas_driver
