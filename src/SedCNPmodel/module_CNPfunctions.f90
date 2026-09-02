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

module CNPfunctions

   !use CNPvariables
   use module_SedCNPvariables !Y.Kwon

   implicit none

   contains 
   
   real function calc_KTB(dTemp, dTC, dKC)
      ! to estimate temperature adjustment parameter KTB for biological processes
     
      real(8), intent(in)                 :: dTemp  ! temperature, deg_C
      real(8), intent(in), dimension(ntm) :: dTC    ! array of temperatures for the rising/falling limb curve, deg_C
      real(8), intent(in), dimension(ntm) :: dKC    ! array of multiplier factors for the rising/falling limb curve
      real(8)                             :: gamma_r, gamma_f  ! temperature multiplier coefficient for the rising/falling limb of curve
      real(8) :: dKTB

      dKTB = 0.0 ! Initialize to a safe value

      IF ((dTemp .LE. dTC(1)) .or. (dTemp .GE. dTC(4))) THEN
         dKTB = 1.0e-9
      ELSE IF (dTemp .LE. dTC(2)) THEN
         ! Rising limb
         if (dTC(2) - dTC(1) > 1.0e-9 .and. dKC(1) > 0. .and. dKC(2) > 0. .and. (1. - dKC(1)) > 0. .and. (1. - dKC(2)) > 0.) then
            gamma_r = 1. / (dTC(2) - dTC(1)) * log((dKC(2) * (1. - dKC(1))) / (dKC(1) * (1. - dKC(2))))
            dKTB = dKC(1) / (dKC(1) + (1. - dKC(1)) * exp(-gamma_r * (dTemp - dTC(1))))
         else
            dKTB = 1.0e-9
         endif
      ELSE IF (dTemp .LE. dTC(3)) THEN
         ! Optimal range
         dKTB = 1.0
      ELSE ! dTemp < dTC(4)
         ! Falling limb
         if (dTC(4) - dTC(3) > 1.0e-9 .and. dKC(3) > 0. .and. dKC(4) > 0. .and. (1. - dKC(3)) > 0. .and. (1. - dKC(4)) > 0.) then
            gamma_f = 1. / (dTC(4) - dTC(3)) * log((dKC(3) * (1. - dKC(4))) / (dKC(4) * (1. - dKC(3))))
            dKTB = dKC(4) / (dKC(4) + (1. - dKC(4)) * exp(-gamma_f * (dTC(4) - dTemp)))
         else
            dKTB = 1.0e-9
         endif
      END IF

      calc_KTB = dKTB

   end function calc_KTB


   real function calc_KTHETA(dTheta, dTheta_s, dTheta_f)
      ! to estimate soil water content adjustment parameters KTHETA

      real(8), intent(in) :: dTheta  ! soil water content, m3 m-3
      real(8), intent(in) :: dTheta_s, dTheta_f  ! soil water content at saturation and field capacity, m3 m-3
      real(8) :: dKTHETA  !BK20251031

      !IF (dTheta .GT. (dTheta_s + dTheta_f) / 2.) THEN
      IF (dTheta .GT. (dTheta_s + dTheta_f) / 2. .AND. dTheta_s - dTheta_f > 1.0e-9) THEN !BK20251029
         dKTHETA = 0.6 + 0.8 * (dTheta_s - dTheta) / (dTheta_s - dTheta_f)
      ELSE IF (dTheta .GE. dTheta_f) THEN
         dKTHETA = 1.
      ELSE IF (dTheta_f > 1.0e-9) THEN  !BK20251029
         dKTHETA = max (dTheta / dTheta_f, 0.05)  !BK20230311
      END IF

      calc_KTHETA = dKTHETA  !BK20251031

   end function calc_KTHETA

  
   real function calc_RK4SOL(rate0, y0)
      ! to solve simultaneous nonlinear equations using the Runge-Kutta 4-th order methods
      
      real(8), intent(in) :: rate0, y0  ! initial values at t0
      real(8) :: dRK1, dRK2, dRK3, dRK4  ! Runge-Kutta coefficients
      real(8) :: y0_local  !BK20251031

      y0_local = y0  
      if (y0_local <= 0.) then
         y0_local = 0.        !BK20231120
      endif

      ! dRK1 = rate0 * y0
      ! dRK2 = rate0 * (y0 + dRK1 * domain%DT / 2.)
      ! dRK3 = rate0 * (y0 + dRK2 * domain%DT / 2.)
      ! dRK4 = rate0 * (y0 + dRK3 * domain%DT) 
      ! RK4SOL = y0 + (dRK1 + 2. * dRK2 + 2. * dRK3 + dRK4 ) * domain%DT / 6.

      !BK20251031
      dRK1 = rate0 * y0_local
      dRK2 = rate0 * (y0_local + dRK1 * domain%DT / 2.)
      dRK3 = rate0 * (y0_local + dRK2 * domain%DT / 2.)
      dRK4 = rate0 * (y0_local + dRK3 * domain%DT) 

      calc_RK4SOL = y0_local + (dRK1 + 2. * dRK2 + 2. * dRK3 + dRK4 ) * domain%DT / 6.

   end function calc_RK4SOL


   real function calc_dPI(eta, Qout, Wstorage)    !BK20250705 !BK20260301
      implicit none
      real(8), intent(in)  :: eta         ! user-specified parameter
      real(8), intent(in)  :: Qout        ! out-flow at the end of the time-step (m3)
      real(8), intent(in)  :: Wstorage    ! water storage at the end of the time-step (m3)
      
      if (Qout > 0.0) then
         if (Qout + Wstorage > 0.0) then
            calc_dPI = 1.0 - EXP(-eta * Qout / (Qout + Wstorage))
         else
            calc_dPI = 0.0
         endif
      else
         calc_dPI = 0.0
      endif

   end function calc_dPI


   !--- BK20250705 commented out: not being used. use calcJulianDay instead
   ! integer function get_JDAY(year, month, day)   !BK20231015
   !    ! convert Gregorian year, month, day to Julian day
   
   !    implicit none
   !    integer, intent(in) :: year, month, day
   
   !    JDAY = FLOOR(2. - REAL(year)/100. + REAL(year)/400. + REAL(day) + 365.25 &
   !           * (REAL(year) + 4716.) + 30.6001* (REAL(month) + 1.) - 1524.5)
   
   ! end function get_JDAY


   subroutine GREGORIAN(jday, year, month, day)  !BK20231015
      ! convert Julian day to Gregorian year, month, day
      implicit none
      integer, intent(in) :: jday
      integer, intent(out) :: year, month, day
      integer :: ii,jj,kk,ll,nn
   
      ll = jday + 68569
      nn = 4 * ll / 146097
      ll = ll - (146097 * nn + 3) / 4
      ii = 4000 * (ll + 1) / 1461001
      ll = ll - 1461 * ii / 4 + 31
      jj = 80 * ll / 2447
      kk = ll - 2447 * jj / 80
      ll = jj / 11
      jj = jj + 2 - 12 * ll
      ii = 100 * (nn - 49) + ii + ll
   
      year = ii
      month = jj
      day = kk
   
   end subroutine GREGORIAN

end module CNPfunctions 
