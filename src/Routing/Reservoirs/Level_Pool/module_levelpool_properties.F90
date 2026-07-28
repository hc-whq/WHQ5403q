! This module defines and instantiates objects
! for a level pool type reservoir's
! parameters/properties. Properties holds
! static/unchanging variables that are
! set when the given reservoir object is
! initialized/instantiated.
module module_levelpool_properties

    use module_reservoir, only: reservoir_properties
    use iso_fortran_env, only: int64
    implicit none

    ! Extend/derive level pool properties from the abstract base
    ! type for reservoir properties.
    type, extends(reservoir_properties) :: levelpool_properties_interface
        real :: lake_area                ! area of lake (km^2)
        real :: weir_elevation           ! bottom of weir elevation (meters AMSL)
        real :: weir_coeffecient         ! weir coefficient
        real :: weir_length              ! weir length (meters)
        real :: dam_length               ! dam length (meters)
        real :: orifice_elevation        ! orifice elevation (meters AMSL)
        real :: orifice_coefficient      ! orifice coefficient
        real :: orifice_area             ! orifice area (meters^2)
        real :: max_depth                ! max depth of reservoir before overtop (meters)
        integer(kind=int64) :: lake_number           ! lake number
        integer :: lake_opt              ! reservoir physics options (1: levelpool, 2: passthrough)

        !=====||___WHQ___||=====!
        !
        real :: slc_elev         ! sluice gate spillway elevation (m)
        real :: slc_length       ! sluice gate length (m)
        real :: slc_coeff        ! sluice gate discharge coefficient (-)
        real :: slc_open         ! sluice gate opening height (m)
        real :: flp_elev         ! flap gate overflow elevation, when closed (m)
        real :: flp_length       ! flap gate length (m)
        real :: flp_coeff        ! flap gate discharge coefficient (-)
        real :: flp_open         ! flap gate opening height (m)
        real :: rad_elev         ! radial gate spillway elevation (m)
        real :: rad_length       ! radial gate length (m)
        real :: rad_coeff        ! radial gate discharge coefficient (-)
        real :: rad_open         ! radial gate opening height (m)
        real :: rad_htrn         ! radial gate trunnion height, from crest to pivot axis (m)
        real :: rad_tau          ! radial gate exponent for trunnion height (-)
        real :: rad_beta         ! radial gate exponent for opening height (-)
        real :: rad_eta          ! radial gate exponent for upstream water head (-)
        !
        !=====||___WHQ___||=====!

    contains

        procedure :: init => levelpool_properties_init
        procedure :: destroy => levelpool_properties_destroy

    end type levelpool_properties_interface

contains

    !Level Pool Properties Constructor
    subroutine levelpool_properties_init(this, lake_area, &
        weir_elevation, weir_coeffecient, weir_length, dam_length, orifice_elevation, &
        orifice_coefficient, orifice_area, max_depth, lake_number, lake_opt, &

        !=====||___WHQ___||=====!
        !
        slc_elev, slc_length, slc_coeff, slc_open, &
        flp_elev, flp_length, flp_coeff, flp_open, &
        rad_elev, rad_length, rad_coeff, rad_open, &
        rad_htrn, rad_tau, rad_beta, rad_eta &
        !
        !=====||___WHQ___||=====!

        )
        implicit none
        class(levelpool_properties_interface), intent(inout) :: this ! the type object being initialized
        real, intent(in)    :: lake_area      	        ! area of lake (km^2)
        real, intent(in)    :: weir_elevation           ! bottom of weir elevation (meters AMSL)
        real, intent(in)    :: weir_coeffecient         ! weir coefficient
        real, intent(in)    :: weir_length              ! weir length (meters)
        real, intent(in)    :: dam_length               ! dam length (meters)
        real, intent(in)    :: orifice_elevation        ! orifice elevation (meters AMSL)
        real, intent(in)    :: orifice_coefficient      ! orifice coefficient
        real, intent(in)    :: orifice_area             ! orifice area (meters^2)
        real, intent(in)    :: max_depth                ! max depth of reservoir before overtop (meters)
        integer(kind=int64), intent(in) :: lake_number              ! lake number
        integer             :: lake_opt                 ! reservoir physics options (1: levelpool, 2: passthrough)

        !=====||___WHQ___||=====!
        !
        real, intent(in)    :: slc_elev         ! sluice gate spillway elevation (m)
        real, intent(in)    :: slc_length       ! sluice gate length (m)
        real, intent(in)    :: slc_coeff        ! sluice gate discharge coefficient (-)
        real, intent(in)    :: slc_open         ! sluice gate opening height (m)
        real, intent(in)    :: flp_elev         ! flap gate overflow elevation, when closed (m)
        real, intent(in)    :: flp_length       ! flap gate length (m)
        real, intent(in)    :: flp_coeff        ! flap gate discharge coefficient (-)
        real, intent(in)    :: flp_open         ! flap gate opening height (m)
        real, intent(in)    :: rad_elev         ! radial gate spillway elevation (m)
        real, intent(in)    :: rad_length       ! radial gate length (m)
        real, intent(in)    :: rad_coeff        ! radial gate discharge coefficient (-)
        real, intent(in)    :: rad_open         ! radial gate opening height (m)
        real, intent(in)    :: rad_htrn         ! radial gate trunnion height, from crest to pivot axis (m)
        real, intent(in)    :: rad_tau          ! radial gate exponent for trunnion height (-)
        real, intent(in)    :: rad_beta         ! radial gate exponent for opening height (-)
        real, intent(in)    :: rad_eta          ! radial gate exponent for upstream water head (-)
        !
        !=====||___WHQ___||=====!

        ! Assign the values passed in to a particular level pool reservoir
        ! properties object's variables.
        this%lake_area = lake_area
        this%weir_elevation = weir_elevation
        this%weir_coeffecient = weir_coeffecient
        this%weir_length = weir_length
        this%orifice_elevation = orifice_elevation
        this%orifice_coefficient = orifice_coefficient
        this%orifice_area = orifice_area
        this%max_depth = max_depth
        this%lake_number = lake_number
        this%dam_length = dam_length
        this%lake_opt = lake_opt

        !=====||___WHQ___||=====!
        !
        this%slc_elev = slc_elev
        this%slc_length = slc_length
        this%slc_coeff = slc_coeff
        this%slc_open = slc_open
        this%flp_elev = flp_elev
        this%flp_length = flp_length
        this%flp_coeff = flp_coeff
        this%flp_open = flp_open
        this%rad_elev = rad_elev
        this%rad_length = rad_length
        this%rad_coeff = rad_coeff
        this%rad_open = rad_open
        this%rad_htrn = rad_htrn
        this%rad_tau = rad_tau
        this%rad_beta = rad_beta
        this%rad_eta = rad_eta
        !
        !=====||___WHQ___||=====!

    end subroutine levelpool_properties_init

    !Level Pool Properties Destructor
    subroutine levelpool_properties_destroy(this)
        implicit none
        class(levelpool_properties_interface), intent(inout) :: this ! the type object being destroyed
    end subroutine levelpool_properties_destroy

end module module_levelpool_properties
