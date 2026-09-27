!=====||__WHQ5403__||=====!
!
!> @brief Progress bar on stderr. It is printed at the end of each time-step (after the messages
!! of the time-step) and stays on the screen while the next time-step runs; messages printed
!! during a time-step call progress_clear first, which erases the bar only if it is shown.
!
module module_progress

#ifdef MPP_LAND
   use module_mpp_land, only: my_id, io_id
#endif
   use ISO_FORTRAN_ENV, only: ERROR_UNIT, OUTPUT_UNIT

   implicit none

   private
   public :: progress_clear, progress_show

   logical, save :: bar_shown = .false.

   contains

   !> erases the progress bar (if shown) so that the following message starts on a clean line
   subroutine progress_clear()
      if (.not. bar_shown) return
#ifdef MPP_LAND
      if (my_id /= io_id) return
#endif
      write(ERROR_UNIT, '(A)', advance='no') char(13)//repeat(' ', 80)//char(13)
      flush(ERROR_UNIT)
      bar_shown = .false.
   end subroutine progress_clear

   !> prints the progress bar of time-step it of nt (at the end of the time-step)
   subroutine progress_show(it, nt, label)
      integer,          intent(in) :: it, nt
      character(len=*), intent(in) :: label
#ifdef MPP_LAND
      if (my_id /= io_id) return
#endif
      flush(OUTPUT_UNIT)
      write(ERROR_UNIT, '(A,A,F6.2,A,I0,A,I0,A)', advance='no') char(13), &
           " WRF-HydroQual running ("//label//") ... ", real(it) / real(max(nt,1)) * 100.0, &
           "% (", it, " / ", nt, ")"
      bar_shown = .true.
      if (it == nt) then
         ! On the final timestep, print a newline, then the success message.
         write(ERROR_UNIT,*)
         write(ERROR_UNIT,*) "The model finished successfully"
         bar_shown = .false.
      endif
      flush(ERROR_UNIT)
   end subroutine progress_show

end module module_progress
!
!=====||__WHQ5403__||=====!
