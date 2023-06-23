! Configuration file for the Mobility Model
! Change any value corresponding to your needs, but keep the line limitation of 132 symbols in mind!
! Some compilers (intel) allow more symbols per line, but to keep compatibility with others make sure your lines are short enough
! You can break a line (even within a string) with the symbol & and continue in the next line (for a string with & again)
! The second & in a string is required, otherwise all spaces before the continued string are included in the string

module mmf_config
  implicit none
  save

  ! Declaration of variables and constants

  ! Path for loading file, which includes all required data (subgroups are not working with this script!)
  character (len = *), parameter :: path_EHEP = "/data/sfb806/human_mobility_model/christian/Applic_Auri_Grav/Exp2/&
&Preprocessed_auri_grav_45k_25k_gaus_init_5hp100km2_glac_extended.nc"
  
  ! Name definitions for loading data. These should be the dimension/variable names in your netcdf file.
  character (len = *), parameter :: name_EHEP_x = "x"                         ! x dimension
  character (len = *), parameter :: name_EHEP_y = "y"                         ! y dimension
  character (len = *), parameter :: name_EHEP_lat = "lat"                     ! latitude
  character (len = *), parameter :: name_EHEP_lon = "lon"                     ! longitude
  character (len = *), parameter :: name_EHEP_pot = "EHEP"                    ! Accessible Potential
  character (len = *), parameter :: name_EHEP_water = "watermask"             ! watermask
  character (len = *), parameter :: name_EHEP_dx = "dx"                       ! dx (gridspacing in x direction, equidistant)
  character (len = *), parameter :: name_EHEP_dy = "dy"                       ! dy (gridspacing in y direction, equidistant)
  character (len = *), parameter :: name_EHEP_start = "dens"                  ! Start values for density
  character (len = *), parameter :: name_EHEP_x_1 = "x_1"                     ! x dimension -1 for the velocity
  character (len = *), parameter :: name_EHEP_y_1 = "y_1"                     ! y dimension -1 for the velocity
  character (len = *), parameter :: name_EHEP_dx_grid = "dx_grid"             ! dx (gridspacing in x direction, distorted)
  character (len = *), parameter :: name_EHEP_dy_grid = "dy_grid"             ! dy (gridspacing in y direction, distorted)
  character (len = *), parameter :: name_EHEP_posx_lat = "posx_lat"           ! latitude of border between two points, x-direc
  character (len = *), parameter :: name_EHEP_posx_lon = "posx_lon"           ! longitude of border between two points, x-direc
  character (len = *), parameter :: name_EHEP_posy_lat = "posy_lat"           ! latitude of border between two points, y-direc
  character (len = *), parameter :: name_EHEP_posy_lon = "posy_lon"           ! longitude of border between two points, y-direc
  character (len = *), parameter :: name_EHEP_time_steps = "time_s"           ! time step dimension
  character (len = *), parameter :: name_EHEP_area = "area"                   ! area of every grid cell
  character (len = *), parameter :: name_EHEP_pot_numb = "EHEPnumb"           ! Potential converted to number of Humans
  character (len = *), parameter :: name_EHEP_start_numb = "hnumb"            ! Start values (number of humans)

  ! Model specific configuration. Make sure that the model is stable with the chosen parameters (Check the documentation)!
  integer, parameter :: fraction_time = 4        ! number of timesteps per year
  integer, parameter :: years = 20000                ! number of years to simulate
  integer, parameter :: save_interval = 10         ! interval of saving intermediate results (1 = yearly, 10 = decennial, etc.)
  integer, parameter :: time_interval = 100          ! interval of time steps in input
  ! Make sure that years divided by save_interval is a natural number (integer) 
  ! The same helds for years divided by time_interval
  double precision, parameter :: vmax = 5.d0     ! in km/year, maximal possible movement speed
  double precision, parameter :: vcrt = 0.50d0*vmax     ! critical velocity, minimum for diffusion coefficient
  double precision, parameter :: deltatau = 1.d0 ! response time of humans to change velocity
  double precision, parameter :: growth = 0.01d0   ! growth rate every year, percentile of the current population, [%/yr]
  double precision, parameter :: alpha = 2.5d0    ! exponent for weibull distribution
  double precision, parameter :: beta = 0.4d0       ! demoninator for weibull distribution
  double precision, parameter :: diff_frac = 0.4d0  ! standard deviation of v, required for the diffusion calculation
  double precision, parameter :: theta = 1.5d0      ! Flux-limiter parameter between 1 and 2, 2 for least dissipative
  double precision, parameter :: repro = 0.0025d0       ! Reproduction onset, in humans per km^2
  double precision, parameter :: eta = 0.1d0      ! "friction" of the velocity term, decays the velocity over time
  double precision, parameter :: maxpot = 0.25d0   ! Fixed maximal potential gradient (humans per grid cell per km) 
  ! maxpot calculated by: 0.5 humans per 100km^2, divided by 100, multiplied with the average area of a grid cell (km),
  ! divided by the average distance of grid cells
  integer, parameter :: switch_adv = 1            ! 0 = turn advection term off, 1 = turn advection term on
  integer, parameter :: switch_diff = 1           ! 0 = turn diffusion term off, 1 = turn diffusion term on
  integer, parameter :: switch_birth = 1          ! 0 = turn birth term off, 1 = turn birth term on
  integer, parameter :: switch_equid = 0          ! 0 = non-equidistant, uses dx/dy_grid, 1 = equidistant, uses dx/dy
  integer, parameter :: switch_slices = 1         ! 0 = only one timeslice of potential in input, 1 = several slices in input
  integer, parameter :: switch_border = 0         ! 0 = closed borders, outer three gridpoints potential is zero
                                                  ! 1 = open borders, outer two gridpoints potential is zero, human can leave
                                                  ! the model, summed up as borderloss in output

  double precision, parameter :: Eps2 = 1d-5      ! required to avoid floating point errors due to divisions in weibull calcs
  double precision, parameter :: Eps = 1d-10   ! required to avoid floating point errors in the angle calculations
  double precision, parameter :: PI = 4 * atan(1.0_8) ! definition of pi

  ! Output path and name definitions. These define the location and description of the output file. No further documentation here
  ! or elsewhere for each parameter.
  character (len = *), parameter :: path_output = "/data/sfb806/human_mobility_model/christian/reworkHDM/output/&
&Dispersal_45k_Rework_run_5hm_oldAdv_fixDiff_lwm_corr_02_05_long_clbd_teststability.nc"
  character (len = *), parameter :: output_name_description = "Description"
  character (len = *), parameter :: output_name_longname = "longname"
  character (len = *), parameter :: output_name_unit = "unit"
  character (len = *), parameter :: output_description = "Reworked Run of OWtE from the 2nd of May 2022"
  character (len = *), parameter :: output_group_name = "Parameters"
  character (len = *), parameter :: name_output_dim_param = "Value"
  character (len = *), parameter :: name_output_dim_x = "x"
  character (len = *), parameter :: name_output_dim_y = "y"
  character (len = *), parameter :: name_output_dim_t = "t"
  character (len = *), parameter :: name_output_dim_x_1 = "x_1"
  character (len = *), parameter :: name_output_dim_y_1 = "y_1"
  character (len = *), parameter :: name_output_var_lat = "lat"
  character (len = *), parameter :: name_output_var_lon = "lon"
  character (len = *), parameter :: name_output_var_posx_lat = "posx_lat"
  character (len = *), parameter :: name_output_var_posx_lon = "posx_lon"
  character (len = *), parameter :: name_output_var_posy_lat = "posy_lat"
  character (len = *), parameter :: name_output_var_posy_lon = "posy_lon"
  character (len = *), parameter :: name_output_var_time = "time"
  character (len = *), parameter :: name_output_var_density = "Density"
  character (len = *), parameter :: name_output_var_AvailPot = "AvailPot"
  character (len = *), parameter :: name_output_var_Border = "Borderloss"
  character (len = *), parameter :: name_output_var_DiffPara = "Diffusion_Parameter"
  character (len = *), parameter :: name_output_var_Velocity_x = "Velocity_x"
  character (len = *), parameter :: name_output_var_Velocity_y = "Velocity_y"
  character (len = *), parameter :: name_output_var_Advection = "Advection_term"
  character (len = *), parameter :: name_output_var_Diffusion = "Diffusion_term"
  character (len = *), parameter :: name_output_var_Birth = "Birth_term"
  character (len = *), parameter :: name_output_var_Death = "Death_term"
  character (len = *), parameter :: name_output_var_watermask = "watermask"
  character (len = *), parameter :: name_output_var_Area = "area"
  character (len = *), parameter :: name_output_var_NumbHuman = "hnumb"
  character (len = *), parameter :: name_output_var_EHEP = "Ehep"
  character (len = *), parameter :: name_output_var_EHEP_dens = "Ehep_dens"
  character (len = *), parameter :: name_output_param_adv = "Advection_Switch"
  character (len = *), parameter :: name_output_param_diff = "Diffusion_Switch"
  character (len = *), parameter :: name_output_param_birth = "Birth_Switch"
  character (len = *), parameter :: name_output_param_dynpot = "Dynamic_Potential"
  character (len = *), parameter :: name_output_param_equidis = "Grid_Structure"
  character (len = *), parameter :: name_output_param_numtime = "Numerical_Timestep"
  character (len = *), parameter :: name_output_param_time = "Calculation_Years"
  character (len = *), parameter :: name_output_param_saveint = "Saving_Intervall"
  character (len = *), parameter :: name_output_param_potint = "Potential_Intervall"
  character (len = *), parameter :: name_output_param_vmax = "Maximum_Speed"
  character (len = *), parameter :: name_output_param_vcrt = "Critical_Velocity"
  character (len = *), parameter :: name_output_param_dtau = "Response_Time"
  character (len = *), parameter :: name_output_param_diffc = "Adv_Diff_Fraction"
  character (len = *), parameter :: name_output_param_theta = "Flux_Limiter"
  character (len = *), parameter :: name_output_param_repro = "Reproduction_Onset"
  character (len = *), parameter :: name_output_param_growth = "Growth_Rate"
  character (len = *), parameter :: name_output_param_alpha = "Weibull_Exponent"
  character (len = *), parameter :: name_output_param_beta = "Weibull_Demoninator"
  character (len = *), parameter :: name_output_param_eta = "Friction_Time_scale"
  character (len = *), parameter :: name_output_param_maxpot = "Maximal potential gradient"
  character (len = *), parameter :: name_output_var_xi = "Xi"
  character (len = *), parameter :: name_output_var_gradx = "gradx"
  character (len = *), parameter :: name_output_var_grady = "grady"
  character (len = *), parameter :: name_output_var_pop_los = "Pop_loss"
  character (len = *), parameter :: name_output_var_advflux_x_p = "Adv_Flux_x_pos"
  character (len = *), parameter :: name_output_var_advflux_x_m = "Adv_Flux_x_neg"
  character (len = *), parameter :: name_output_var_advflux_y_p = "Adv_Flux_y_pos"
  character (len = *), parameter :: name_output_var_advflux_y_m = "Adv_Flux_y_neg"
  character (len = *), parameter :: name_output_var_diffflux_x_p = "Diff_Flux_x_pos"
  character (len = *), parameter :: name_output_var_diffflux_x_m = "Diff_Flux_x_neg"
  character (len = *), parameter :: name_output_var_diffflux_y_p = "Diff_Flux_y_pos"
  character (len = *), parameter :: name_output_var_diffflux_y_m = "Diff_Flux_y_neg"
  character (len = *), parameter :: name_output_var_dx = "dx_equi"
  character (len = *), parameter :: name_output_var_dy = "dy_equi"
  character (len = *), parameter :: name_output_var_dx_grid = "dx_grid"
  character (len = *), parameter :: name_output_var_dy_grid = "dy_grid"
  character (len = *), parameter :: output_long_lat = "Latitude"
  character (len = *), parameter :: output_unit_lat = "degrees_north"
  character (len = *), parameter :: output_long_lon = "Longitude"
  character (len = *), parameter :: output_unit_lon = "degrees_east"
  character (len = *), parameter :: output_long_time = "Timesteps"
  character (len = *), parameter :: output_unit_time = "years"
  character (len = *), parameter :: output_long_Density = "Density of modern humans"
  character (len = *), parameter :: output_unit_Density = "humans / km^2"
  character (len = *), parameter :: output_long_NumbHuman = "Total number of Humans within the grid cell"
  character (len = *), parameter :: output_unit_NumbHuman = "Humans per grid cell"
  character (len = *), parameter :: output_long_Area = "Area of the grid cell"
  character (len = *), parameter :: output_unit_Area = "km^2"
  character (len = *), parameter :: output_long_AvailPot = "Available environmental human existance potential"
  character (len = *), parameter :: output_unit_AvailPot = "humans / km^2"
  character (len = *), parameter :: output_long_Borderloss = "Amount of humans moving out of the system during the timestep"
  character (len = *), parameter :: output_unit_Borderloss = "number of humans"!"humans / 100 km^2"
  character (len = *), parameter :: output_long_DiffPara = "Diffusion Parameter"
  character (len = *), parameter :: output_unit_DiffPara = "km^2/yr"
  character (len = *), parameter :: output_long_Velocity_x = "Current dispersal velocity per gridcell in x-direction"
  character (len = *), parameter :: output_unit_Velocity_x = "km/yr"
  character (len = *), parameter :: output_long_Velocity_y = "Current dispersal velocity per gridcell in y-direction"
  character (len = *), parameter :: output_unit_Velocity_y = "km/yr"
  character (len = *), parameter :: output_long_posx_lat = "latitude of border between two points, x-direction"
  character (len = *), parameter :: output_unit_posx_lat = "degrees_north"
  character (len = *), parameter :: output_long_posx_lon = "longitude of border between two points, x-direction"
  character (len = *), parameter :: output_unit_posx_lon = "degrees_east"
  character (len = *), parameter :: output_long_posy_lat = "latitude of border between two points, y-direction"
  character (len = *), parameter :: output_unit_posy_lat = "degrees_north"
  character (len = *), parameter :: output_long_posy_lon = "longitude of border between two points, y-direction"
  character (len = *), parameter :: output_unit_posy_lon = "degrees_east"
  character (len = *), parameter :: output_long_Advection = "Impact of advection term on the change of the human density"
  character (len = *), parameter :: output_unit_Advection = "humans / km^2"
  character (len = *), parameter :: output_long_Diffusion = "Impact of diffusion term on the change of the human density"
  character (len = *), parameter :: output_unit_Diffusion = "humans / km^2"
  character (len = *), parameter :: output_long_Birth = "Impact of the birth term on the change of the human density"
  character (len = *), parameter :: output_unit_Birth = "humans / km^2"
  character (len = *), parameter :: output_long_Death = "Impact of the death term on the change of the human density"
  character (len = *), parameter :: output_unit_Death = "humans / km^2"
  character (len = *), parameter :: output_long_watermask = "Watermask that determines valid (=1) and invalid (=0) gridpoints"
  character (len = *), parameter :: output_unit_watermask = "none" 
  character (len = *), parameter :: output_long_EHEP = "Human Existence Potential"
  character (len = *), parameter :: output_unit_EHEP = "number of humans"
  character (len = *), parameter :: output_long_EHEP_dens = "Human Existence Potential"
  character (len = *), parameter :: output_unit_EHEP_dens = "humans/km^2"
  character (len = *), parameter :: output_long_param_adv = "Indicator showing whether the Advection Term is activated or not"
  character (len = *), parameter :: output_unit_param_adv = "1 = enabled, 0 = disabled"
  character (len = *), parameter :: output_long_param_diff = "Indicator showing whether the Diffusion Term is activated or not"
  character (len = *), parameter :: output_unit_param_diff = "1 = enabled, 0 = disabled"
  character (len = *), parameter :: output_long_param_birth = "Indicator showing wether the Birth Term is activatedor not"
  character (len = *), parameter :: output_unit_param_birth = "1 = enabled, 0 = disabled"
  character (len = *), parameter :: output_long_param_dynpot = "Indicator showing if the potential is changing over time"
  character (len = *), parameter :: output_unit_param_dynpot = "1 = dynamic Potential, 0 = static Potential"
  character (len = *), parameter :: output_long_param_equidis = "Indicator showing which grid distance type is used"
  character (len = *), parameter :: output_unit_param_equidis = "1 = equidistant, 0 = variable distances"
  character (len = *), parameter :: output_long_param_numtime = "Numerical timesteps per year of simulation"
  character (len = *), parameter :: output_unit_param_numtime = "Number of timesteps per year"
  character (len = *), parameter :: output_long_param_time = "Total years of calculation"
  character (len = *), parameter :: output_unit_param_time = "years"
  character (len = *), parameter :: output_long_param_saveint = "Intervall of intermediate year steps to be saved"
  character (len = *), parameter :: output_unit_param_saveint = "years"
  character (len = *), parameter :: output_long_param_potint = "Intervall of changing potential in the input"
  character (len = *), parameter :: output_unit_param_potint = "years"
  character (len = *), parameter :: output_long_param_vmax = "Maximum allowed travel speed of humans per year. Important for &
&the stability of the model"
  character (len = *), parameter :: output_unit_param_vmax = "km/yr"
  character (len = *), parameter :: output_long_param_growth = "Growth rate of the population, evaluated every year"
  character (len = *), parameter :: output_unit_param_growth = "None"
  character (len = *), parameter :: output_long_param_vcrt = "Minimum velocity that is considered for the diffusion term"
  character (len = *), parameter :: output_unit_param_vcrt = "km/yr"
  character (len = *), parameter :: output_long_param_dtau = "Response time of humans to react to changing potentials on their &
&travel path"
  character (len = *), parameter :: output_unit_param_dtau = "yr"
  character (len = *), parameter :: output_long_param_diffc = "Fraction between advection and diffusion term, controls strength &
&of the diffusion term"
  character (len = *), parameter :: output_unit_param_diffc = "None"
  character (len = *), parameter :: output_long_param_theta = "Control parameter for the flux-limiter, 1 = most dissipative, 2 = &
&least dissipative"
  character (len = *), parameter :: output_unit_param_theta = "None"
  character (len = *), parameter :: output_long_param_repro = "Reproduction Onset, minimum of humans required to trigger birth"
  character (len = *), parameter :: output_unit_param_repro = "humans/km^2"
  character (len = *), parameter :: output_long_param_alpha = "Weibull distribution Exponent value"
  character (len = *), parameter :: output_unit_param_alpha = "None"
  character (len = *), parameter :: output_long_param_beta = "Weibull distribution Demoninator value"
  character (len = *), parameter :: output_unit_param_beta = "None"
  character (len = *), parameter :: output_long_param_eta = "Time-scale of the Velocity Friction term"
  character (len = *), parameter :: output_unit_param_eta = "1/yr"
  character (len = *), parameter :: output_long_param_maxpot = "Maximal potential gradient required to achieve vmax"
  character (len = *), parameter :: output_unit_param_maxpot = "humans/Average gridcell/km"
  character (len = *), parameter :: output_long_xi = "Xi angle value"
  character (len = *), parameter :: output_unit_xi = "radiant"
  character (len = *), parameter :: output_long_advflux_x_p = "Advection Flux, positive x-direction"
  character (len = *), parameter :: output_unit_advflux_x_p = "humans/km/10yr"
  character (len = *), parameter :: output_long_advflux_x_m = "Advection Flux, negative x-direction"
  character (len = *), parameter :: output_unit_advflux_x_m = "humans/km/10yr"
  character (len = *), parameter :: output_long_advflux_y_p = "Advection Flux, positive y-direction"
  character (len = *), parameter :: output_unit_advflux_y_p = "humans/km/10yr"
  character (len = *), parameter :: output_long_advflux_y_m = "Advection Flux, negative y-direction"
  character (len = *), parameter :: output_unit_advflux_y_m = "humans/km/10yr"
  character (len = *), parameter :: output_long_diffflux_x_p = "Diffusion Flux, positive x-direction"
  character (len = *), parameter :: output_unit_diffflux_x_p = "humans/km/10yr"
  character (len = *), parameter :: output_long_diffflux_x_m = "Diffusion Flux, negative x-direction"
  character (len = *), parameter :: output_unit_diffflux_x_m = "humans/km/10yr"
  character (len = *), parameter :: output_long_diffflux_y_p = "Diffusion Flux, positive y-direction"
  character (len = *), parameter :: output_unit_diffflux_y_p = "humans/km/10yr"
  character (len = *), parameter :: output_long_diffflux_y_m = "Diffusion Flux, negative y-direction"
  character (len = *), parameter :: output_unit_diffflux_y_m = "humans/km/10yr"
  character (len = *), parameter :: output_long_dx = "Equidistant grid spacing x-direction"
  character (len = *), parameter :: output_unit_dx = "km"
  character (len = *), parameter :: output_long_dy = "Equidistant grid spacing y-direction"
  character (len = *), parameter :: output_unit_dy = "km"
  character (len = *), parameter :: output_long_dx_grid = "Actual grid spacing x-direction"
  character (len = *), parameter :: output_unit_dx_grid = "km"
  character (len = *), parameter :: output_long_dy_grid = "Actual grid spacing y-direction"
  character (len = *), parameter :: output_unit_dy_grid = "km"
  character (len = *), parameter :: output_long_pop_los = "Numerical population loss"
  character (len = *), parameter :: output_unit_pop_los = "number of humans"
end module mmf_config
