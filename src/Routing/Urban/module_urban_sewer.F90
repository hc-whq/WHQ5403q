module module_urban_sewer
    implicit none
    
    ! Parameters for Solutes
    integer, parameter :: num_solute_classes = 4
    
    ! Type definition for a single urban sewer cell
    type urban_sewer_cell
        logical :: is_active          ! Is this an urban cell with a sewer?
        
        ! Geometry
        real :: width                 ! Width of rectangular sewer [m]
        real :: depth_max             ! Maximum depth of sewer [m]
        real :: area_upstream         ! Upstream contributing area [m2]
        real :: slope                 ! Slope of the sewer channel [m/m]
        real :: length                ! Length of sewer segment [m]
        real :: manning_n             ! Manning's roughness for sewer
        
        ! Configuration
        integer :: sewer_type         ! 1: Combined, 2: Separated
        logical :: has_pump           ! Does this cell have a pumping station?
        real :: pump_capacity         ! Pumping capacity [m3/s]
        
        ! Topology & Elevation
        integer :: i_down, j_down     ! Indices of downstream cell
        real :: z_bed                 ! Bed elevation of sewer [m]
        real :: dist_down             ! Distance to downstream cell [m]

        ! State Variables - Hydraulic
        real :: water_depth           ! Current water depth in sewer [m]
        real :: discharge             ! Current Outflow [m3/s]
        real :: inflow                ! Total Inflow (Upstream + Surface + Point) [m3/s]
        real :: overflow              ! Overflow rate (CSO or SSO) [m3/s]
        
        ! State Variables - Solutes
        ! Surface accumulation (on impervious surface)
        real, dimension(num_solute_classes) :: surface_solute_mass ! [kg]
        
        ! In-sewer transport
        real, dimension(num_solute_classes) :: sewer_solute_mass   ! Mass in sewer water [kg]
        real, dimension(num_solute_classes) :: solute_outflow      ! Mass flux out [kg/s]
        
    end type urban_sewer_cell

    ! Global Grid of Sewer Cells
    type(urban_sewer_cell), allocatable, dimension(:,:) :: urban_sewer_grid
    
    ! Module Constants
    real, parameter :: vacuum_efficiency = 0.5  ! 50% removal efficiency
    real, parameter :: accumulation_rate(num_solute_classes) = [0.1, 0.2, 0.1, 0.05] ! kg/day (Example)
    integer, parameter :: TYPE_COMBINED = 1
    integer, parameter :: TYPE_SEPARATED = 2

contains

    !---------------------------------------------------------------------
    ! Subroutine: urban_sewer_init
    ! Purpose: Initialize the urban sewer module, allocate grid, and setup geometries
    !---------------------------------------------------------------------
    subroutine urban_sewer_init(ix, jx, is_urban, flow_dir, terrain_z, dx, dy, dt)
        implicit none
        integer, intent(in) :: ix, jx
        integer, dimension(ix,dx), intent(in) :: is_urban ! 1 if urban
        integer, dimension(ix,dx), intent(in) :: flow_dir ! 8-direction flow
        real, dimension(ix,dx), intent(in) :: terrain_z ! Surface Elevation
        real, intent(in) :: dx, dy, dt
        
        integer :: i, j
        
        allocate(urban_sewer_grid(ix, jx))
        
        ! Initialize cells
        do j = 1, jx
            do i = 1, ix
                if (is_urban(i,j) == 1) then
                    urban_sewer_grid(i,j)%is_active = .true.
                    urban_sewer_grid(i,j)%length = dx ! Approximation
                    urban_sewer_grid(i,j)%manning_n = 0.013 ! Concrete pipe
                    
                    ! Elevation setup (Assume sewer is 2m below terrain)
                    urban_sewer_grid(i,j)%z_bed = terrain_z(i,j) - 2.0
                    
                    ! Topology Mapping (Simplified D8 Mapping for illustration)
                    ! Assumes standard: 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
                    call get_downstream_index(i, j, ix, jx, flow_dir(i,j), &
                                            urban_sewer_grid(i,j)%i_down, &
                                            urban_sewer_grid(i,j)%j_down)
                    
                    ! Calculate distance to downstream
                    if (urban_sewer_grid(i,j)%i_down > 0) then
                        real :: d_x, d_y
                        d_x = (urban_sewer_grid(i,j)%i_down - i) * dx
                        d_y = (urban_sewer_grid(i,j)%j_down - j) * dy
                        urban_sewer_grid(i,j)%dist_down = sqrt(d_x**2 + d_y**2)
                    else
                        urban_sewer_grid(i,j)%dist_down = dx ! Outlet
                    endif
                    
                    ! Default Configuration (User settable logic would go here)
                    urban_sewer_grid(i,j)%sewer_type = TYPE_COMBINED 
                    urban_sewer_grid(i,j)%has_pump = .false.
                    urban_sewer_grid(i,j)%pump_capacity = 0.0
                    
                    ! Calculate Upstream Area (Simplified placeholder for full routing traversal)
                    ! In full implementation, we would traverse flow_dir to sum up areas.
                    ! Here we assume a conceptual area based on grid size for demo.
                    urban_sewer_grid(i,j)%area_upstream = dx * dy 
                    
                    ! Set Geometry based on Upstream Area (User Logic)
                    call set_sewer_geometry(urban_sewer_grid(i,j))
                    
                    ! Initialize States
                    urban_sewer_grid(i,j)%water_depth = 0.0
                    urban_sewer_grid(i,j)%discharge = 0.0
                    urban_sewer_grid(i,j)%surface_solute_mass = 0.0
                    urban_sewer_grid(i,j)%sewer_solute_mass = 0.0
                else
                    urban_sewer_grid(i,j)%is_active = .false.
                endif
            end do
        end do
        
    end subroutine urban_sewer_init

    ! Helper to map D8 flow direction
    subroutine get_downstream_index(i, j, nx, ny, fdir, i_next, j_next)
        integer, intent(in) :: i, j, nx, ny, fdir
        integer, intent(out) :: i_next, j_next
        
        i_next = 0
        j_next = 0
        
        ! Simple D8 Logic (Example)
        if (fdir == 1) then      ! East
            i_next = i + 1; j_next = j
        elseif (fdir == 4) then  ! South
            i_next = i; j_next = j - 1
        elseif (fdir == 16) then ! West
            i_next = i - 1; j_next = j
        elseif (fdir == 64) then ! North
            i_next = i; j_next = j + 1
        endif
        ! Include diagonals in full implementation...
        
        ! Boundary check
        if (i_next < 1 .or. i_next > nx .or. j_next < 1 .or. j_next > ny) then
            i_next = 0; j_next = 0
        endif
    end subroutine get_downstream_index

    !---------------------------------------------------------------------
    ! Subroutine: urban_sewer_geometry
    ! Purpose: Define width/depth based on upstream area
    !---------------------------------------------------------------------
    subroutine set_sewer_geometry(cell)
        implicit none
        class(urban_sewer_cell), intent(inout) :: cell
        
        ! Conceptual sizing rule (Example)
        ! Width scales with sqrt of area, Depth scales similarly
        if (cell%area_upstream < 10000.0) then
            cell%width = 0.5
            cell%depth_max = 0.5
        elseif (cell%area_upstream < 100000.0) then
            cell%width = 1.0
            cell%depth_max = 1.0
        else
            cell%width = 2.0
            cell%depth_max = 1.5
        endif
    end subroutine set_sewer_geometry

    !---------------------------------------------------------------------
    ! Subroutine: urban_sewer_run
    ! Purpose: Main driver for the module, called every time step
    !---------------------------------------------------------------------
    subroutine urban_sewer_drive(ix, jx, dt, rainfall, surface_runoff, point_sources, vacuum_flag)
        implicit none
        integer, intent(in) :: ix, jx
        real, intent(in) :: dt
        real, dimension(ix,jx), intent(in) :: rainfall       ! mm/s
        real, dimension(ix,jx), intent(in) :: surface_runoff ! m3/s (Runoff entering sewer)
        real, dimension(ix,jx), intent(in) :: point_sources  ! m3/s (Pollution loads/industrial)
        logical, intent(in) :: vacuum_flag                   ! True if vacuum cleaning happens this step
        
        integer :: i, j, k
        real :: q_in, q_out, q_pump, flow_area, r_hyd, velocity
        
        do j = 1, jx
            do i = 1, ix
                if (.not. urban_sewer_grid(i,j)%is_active) cycle
                
                ! 1. Solute Processes (Impervious Surface)
                ! Accumulation
                urban_sewer_grid(i,j)%surface_solute_mass = urban_sewer_grid(i,j)%surface_solute_mass + &
                                                            (accumulation_rate * dt / 86400.0)
                
                ! Vacuum Removal
                if (vacuum_flag) then
                   urban_sewer_grid(i,j)%surface_solute_mass = urban_sewer_grid(i,j)%surface_solute_mass * (1.0 - vacuum_efficiency)
                endif
                
                ! Washoff (Transport to sewer via surface runoff)
                ! Exponential washoff function: M_wash = M_accum * (1 - exp(-k * runoff))
                do k = 1, num_solute_classes
                     ! Simplified washoff
                     if (surface_runoff(i,j) > 0.0) then
                        urban_sewer_grid(i,j)%surface_solute_mass(k) = urban_sewer_grid(i,j)%surface_solute_mass(k) * 0.9 ! 10% washoff per step example
                        ! Add washed mass to sewer
                        urban_sewer_grid(i,j)%sewer_solute_mass(k) = urban_sewer_grid(i,j)%sewer_solute_mass(k) + &
                                                                     (urban_sewer_grid(i,j)%surface_solute_mass(k) * 0.1)
                     endif
                end do
                
                ! 2. Hydraulic Routing (Diffusive Wave - Backwater Support)
                
                ! Inputs
                q_in = 0.0
                
                ! Surface Runoff entry
                if (urban_sewer_grid(i,j)%sewer_type == TYPE_COMBINED .or. &
                    urban_sewer_grid(i,j)%sewer_type == TYPE_SEPARATED) then
                    q_in = q_in + surface_runoff(i,j)
                endif
                
                ! Point Sources
                if (urban_sewer_grid(i,j)%sewer_type == TYPE_COMBINED) then
                    q_in = q_in + point_sources(i,j)
                endif
                
                urban_sewer_grid(i,j)%inflow = q_in
                
                ! Discharge Calculation (Diffusive Wave)
                if (urban_sewer_grid(i,j)%water_depth > 0.001) then
                    real :: wse_curr, wse_next, friction_slope, dist
                    real :: A_flow, R_h
                    
                    ! Current Water Surface Elevation
                    wse_curr = urban_sewer_grid(i,j)%z_bed + urban_sewer_grid(i,j)%water_depth
                    
                    ! Downstream WSE
                    if (urban_sewer_grid(i,j)%i_down > 0) then
                        integer :: id, jd
                        id = urban_sewer_grid(i,j)%i_down
                        jd = urban_sewer_grid(i,j)%j_down
                        wse_next = urban_sewer_grid(id,jd)%z_bed + urban_sewer_grid(id,jd)%water_depth
                        dist = urban_sewer_grid(i,j)%dist_down
                    else
                         ! Outlet condition (Normal Depth assumption or Critical)
                         ! Here we assume the slope is just the bed slope
                         wse_next = wse_curr - (urban_sewer_grid(i,j)%slope * urban_sewer_grid(i,j)%length)
                         dist = urban_sewer_grid(i,j)%length
                    endif
                    
                    ! Friction Slope (S_f = -dH/dx)
                    friction_slope = (wse_curr - wse_next) / dist
                    
                    ! Geometric properties
                    A_flow = urban_sewer_grid(i,j)%water_depth * urban_sewer_grid(i,j)%width
                    R_h = A_flow / (urban_sewer_grid(i,j)%width + 2.0 * urban_sewer_grid(i,j)%water_depth)
                    
                    ! Manning's Eq for Diffusive Wave: Q = (1/n) * A * R^(2/3) * sqrt(|Sf|) * sgn(Sf)
                    ! Note: We typically handle reverse flow (negative Sf). 
                    ! For simplicity in this conceptual module, we calculate magnitude and direction.
                    
                    velocity = (1.0 / urban_sewer_grid(i,j)%manning_n) * (R_h**(2.0/3.0)) * sqrt(abs(friction_slope))
                    
                    if (friction_slope < 0.0) then
                        q_out = -1.0 * A_flow * velocity ! Reverse flow (Backwater pushing back)
                    else
                        q_out = A_flow * velocity
                    endif
                    
                else
                    q_out = 0.0
                endif
                
                ! Pumping Station Logic
                if (urban_sewer_grid(i,j)%has_pump) then
                    if (q_out > urban_sewer_grid(i,j)%pump_capacity) then
                        ! Overflow occurs
                        urban_sewer_grid(i,j)%overflow = q_out - urban_sewer_grid(i,j)%pump_capacity
                        q_out = urban_sewer_grid(i,j)%pump_capacity ! Pump is maxed out
                    else
                        urban_sewer_grid(i,j)%overflow = 0.0
                    endif
                else
                     ! Check for pipe capacity overflow (surcharge)
                     real :: max_capacity
                     max_capacity = (1.0/urban_sewer_grid(i,j)%manning_n) * &
                                    (urban_sewer_grid(i,j)%width * urban_sewer_grid(i,j)%depth_max) * &
                                    (( (urban_sewer_grid(i,j)%width * urban_sewer_grid(i,j)%depth_max) / &
                                       (urban_sewer_grid(i,j)%width + 2*urban_sewer_grid(i,j)%depth_max) )**(2.0/3.0)) * &
                                    urban_sewer_grid(i,j)%slope**0.5
                     
                     if (q_out > max_capacity) then
                         urban_sewer_grid(i,j)%overflow = q_out - max_capacity
                         q_out = max_capacity
                     endif
                endif
                
                urban_sewer_grid(i,j)%discharge = q_out
                
                ! Update Depth
                ! V_new = V_old + (Qin - Qout - Qoverflow) * dt
                real :: vol_new
                vol_new = (urban_sewer_grid(i,j)%water_depth * urban_sewer_grid(i,j)%width * urban_sewer_grid(i,j)%length) + &
                          (q_in - q_out - urban_sewer_grid(i,j)%overflow) * dt
                
                if (vol_new < 0.0) vol_new = 0.0
                urban_sewer_grid(i,j)%water_depth = vol_new / (urban_sewer_grid(i,j)%width * urban_sewer_grid(i,j)%length)
                
                ! 3. Solute Transport
                ! Simple mixing and advection
                ! C = M / V
                ! M_out = C * Q_out * dt
                if (vol_new > 0.0) then
                    do k = 1, num_solute_classes
                        real :: conc, mass_out
                        conc = urban_sewer_grid(i,j)%sewer_solute_mass(k) / vol_new
                        mass_out = conc * q_out * dt
                        urban_sewer_grid(i,j)%solute_outflow(k) = mass_out
                        urban_sewer_grid(i,j)%sewer_solute_mass(k) = urban_sewer_grid(i,j)%sewer_solute_mass(k) - mass_out
                    end do
                else
                    urban_sewer_grid(i,j)%solute_outflow = 0.0
                endif
                
            end do
        end do
        
    end subroutine urban_sewer_drive

end module module_urban_sewer
