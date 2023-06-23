module mmf_setup
  !use mpi
  use netcdf
  use mmf_config
  use omp_lib
  implicit none
  save

  ! Declaration of variables and constants
  integer :: nc_id_EHEP, nc_id_out, lat_id, lon_id, pot_id, water_id, water_id_out, lat_id_out, lon_id_out
  integer :: x, y, x_1, y_1, x_1_id, y_1_id, x_id, y_id, dx_id, dy_id, x_id_out, y_id_out, x_1_id_out, y_1_id_out
  integer :: posx_lat_id, posx_lon_id, posy_lat_id, posy_lon_id, dx_grid_id, dy_grid_id
  integer :: dx_grid_id_out, dy_grid_id_out, posx_lat_id_out, posx_lon_id_out, posy_lat_id_out, posy_lon_id_out
  integer :: time_slices_id, potnum_id, area_id, numdens_id, area_id_out, numdens_id_out, ehep_id_out, diffpara_id_out
  integer :: paramgrp_id_out, param_dim_id, param_adv_id, param_birth_id, param_death_id, param_eta_id
  integer :: param_dynpot_id, param_equidis_id, param_numtime_id, param_time_id, param_saveint_id, param_potint_id
  integer :: param_vmax_id, param_growth_id, xi_id_out, gradx_id_out, grady_id_out, pop_los_id_out
  integer :: Fx_id_out, Fy_id_out, Px_id_out,Py_id_out
  integer :: param_vcrt_id, param_dtau_id, param_diffc_id, param_theta_id, param_repro_id, param_alpha_id, param_beta_id
  integer :: Fxp_id_out, Fxm_id_out, Fyp_id_out, Fym_id_out, Pxp_id_out, Pxm_id_out, Pyp_id_out, Pym_id_out
  integer :: dx_id_out, dy_id_out, param_maxpot_id
  integer :: j,k
  double precision,dimension(:,:),allocatable :: lat, lon, potential, watermask, dx_grid, dy_grid
  double precision,dimension(:,:),allocatable :: density, avail_pot, diffusion, potnum, area, numdens, minpot, minnum
  double precision,dimension(:,:),allocatable :: density_s, avail_pot_s, wei_scale, Theta_tmp, reprod
  double precision,dimension(:,:),allocatable :: vxfield, vyfield, vxfield_old, vyfield_old, gradxfield, gradyfield
  double precision,dimension(:,:),allocatable :: k1, k2, k3, k4
  double precision,dimension(:,:),allocatable :: posx_lat, posx_lon, posy_lat, posy_lon
  double precision,dimension(:,:),allocatable :: advec_t, diffu_t, birth_t, death_t
  double precision,dimension(:,:),allocatable :: slope, intercept, xi, local_grid, vlength
  double precision,dimension(:,:),allocatable :: Fxp_s, Fxm_s, Fyp_s, Fym_s, Pxp_s, Pxm_s, Pyp_s, Pym_s
  double precision,dimension(:,:,:),allocatable :: advec_temp, diffu_temp, potential_sl, watermask_sl, potnum_sl
  double precision,dimension(:,:,:),allocatable :: birth_temp, death_temp
  double precision,dimension(:,:,:), allocatable :: Fxp_tmp, Fxm_tmp, Fyp_tmp, Fym_tmp, Pxp_tmp, Pxm_tmp, Pyp_tmp, Pym_tmp
  double precision,dimension(:),allocatable :: time_array
  double precision :: dt, dx, dy, wei_maxpos, wei_max, gamma_para, local_speed, borderloss, pop_los
  integer :: currentyear, timesteps, t_id_out, time_id, density_id, availpot_id, borderloss_id
  integer :: velocity_x_id, velocity_y_id, advection_id, diffusion_id, birth_id, death_id
  integer :: start_id, i, savetimes, t, time_slices, currentslice, uppboundbord,lowboundbord
  integer,dimension(3) :: dim_ids_all, dim_ids_vel
  integer,dimension(2) :: dim_latlon, dim_latlon_1

  ! Definitions of routines used in the mmf_run file
  contains
    ! Improved Error output, if anything goes wrong you get the number of the error call and the original message
    subroutine check(status,state,operation,variable,timestep,year)
      integer, intent(in) :: status
      integer, intent(in),optional :: state, timestep, year
      character (len=*),optional :: operation,variable
      if (status /= nf90_noerr) then
        print *,"Error detected, aborting program..."
        print *,"Information:"
        if (present(state)) then
          select case (state)
          case (1)
            print *,"Programm stopped at mmf_setup_load() while loading the input file."
          case (2)
            print *,"Programm stopped at mmf_setup_load() while creating the output file."
          case (3)
            print *,"Programm stopped at mmf_setup_load() while saving the initial state."
          case (4)
            print *,"Programm stopped within mmf_run.f95 while saving the results."
          end select
        end if
        if (present(operation)) then
          print *,"The operation that failed was nf90_",operation,"."
        end if
        if (present(variable)) then
          print *,"The variable ",variable," was involved in the operation."
        end if
        if (present(timestep) .AND. present(year)) then
          print *,"The exception occured during the timestep ",timestep," within the modelled year ",year,"."
        end if
        print *,"The orignal Error statement follows:"
        print *,trim(nf90_strerror(status))
        
        stop 2
      end if
    end subroutine check

    elemental function Theta_func(rho, rho_max, rho_min)
      double precision, intent(in) :: rho, rho_max, rho_min
      double precision :: Theta_func
      if (abs(rho_max - rho_min) > Eps2) then
        Theta_func = (rho - rho_min) / (rho_max - rho_min)
        if ((Theta_func /= Theta_func) .OR. (Theta_func >= 10.d0))then
          Theta_func = 10.d0
        end if
      else
        Theta_func = 10.d0
      end if
    end function Theta_func

    elemental function weibull_func(ThetaIn, Cons)
      double precision, intent(in):: ThetaIn, Cons
      double precision :: weibull_func
      if (abs(ThetaIn) > Eps2) then
        weibull_func = Cons * (alpha/beta) * (ThetaIn/beta) ** (alpha-1.0d0) * exp(-(ThetaIn/beta)**alpha)
        if (weibull_func /= weibull_func) then
          weibull_func = 0.d0
        end if
      else
        weibull_func = 0.d0
      end if
    end function weibull_func

!    elemental function calc_avail_func(pot_in, dens_in)
!      double precision, intent(in) :: pot_in, dens_in
!      double precision :: Theta
!      double precision :: calc_avail_func

!      Theta = Theta_func(dens_in, pot_in, minnum)
!      calc_avail_func = weibull_func(Theta, wei_scale)
!    end function calc_avail_func

    subroutine mmf_setup_load()
      ! Main routine to initialise the required variables
      ! Open the input file and save the length of the dimensions
      call check(nf90_open(path_EHEP,nf90_nowrite,nc_id_EHEP),1,"open")
      call check(nf90_inq_dimid(nc_id_EHEP,name_EHEP_x,x_id),1,"inq_dimid",name_EHEP_x)
      call check(nf90_inquire_dimension(nc_id_EHEP,x_id,len=x),1,"inquire_dimension",name_EHEP_x)
      call check(nf90_inq_dimid(nc_id_EHEP,name_EHEP_y,y_id),1,"inq_dimid",name_EHEP_y)
      call check(nf90_inquire_dimension(nc_id_EHEP,y_id,len=y),1,"inquire_dimension",name_EHEP_y)

      call check(nf90_inq_dimid(nc_id_EHEP,name_EHEP_x_1,x_1_id),1,"inq_dimid",name_EHEP_x_1)
      call check(nf90_inquire_dimension(nc_id_EHEP,x_1_id,len=x_1),1,"inquire_dimension",name_EHEP_x_1)
      call check(nf90_inq_dimid(nc_id_EHEP,name_EHEP_y_1,y_1_id),1,"inq_dimid",name_EHEP_y_1)
      call check(nf90_inquire_dimension(nc_id_EHEP,y_1_id,len=y_1),1,"inquire_dimension",name_EHEP_y_1)

      if (switch_slices == 1) then
        call check(nf90_inq_dimid(nc_id_EHEP,name_EHEP_time_steps,time_slices_id),1,"inq_dimid",name_EHEP_time_steps)
        call check(nf90_inquire_dimension(nc_id_EHEP,time_slices_id,len=time_slices),1,"inquire_dimension",name_EHEP_time_steps)
      end if
      ! Get the variable IDs, allocate memory and load the data
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_lat,lat_id),1,"inq_varid",name_EHEP_lat)
      allocate(lat(x,y))
      call check(nf90_get_var(nc_id_EHEP,lat_id,lat),1,"get_var",name_EHEP_lat)
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_lon,lon_id),1,"inq_varid",name_EHEP_lon)
      allocate(lon(x,y))
      call check(nf90_get_var(nc_id_EHEP,lon_id,lon),1,"get_var",name_EHEP_lon)
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_area,area_id),1,"inq_varid",name_EHEP_area)
      allocate(area(x,y))
      call check(nf90_get_var(nc_id_EHEP,area_id,area),1,"get_var",name_EHEP_area)
      
      if (switch_slices == 0) then
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_pot,pot_id),1,"inq_varid",name_EHEP_pot)
        allocate(potential(x,y))
        call check(nf90_get_var(nc_id_EHEP,pot_id,potential),1,"get_var",name_EHEP_pot)
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_water,water_id),1,"inq_varid",name_EHEP_water)
        allocate(watermask(x,y))
        call check(nf90_get_var(nc_id_EHEP,water_id,watermask),1,"get_var",name_EHEP_water)
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_pot_numb,potnum_id),1,"inq_varid",name_EHEP_pot_numb)
        allocate(potnum(x,y))
        call check(nf90_get_var(nc_id_EHEP,potnum_id,potnum),1,"get_var",name_EHEP_pot_numb)
      else if (switch_slices == 1) then
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_pot,pot_id),1,"inq_varid",name_EHEP_pot)
        allocate(potential_sl(x,y,time_slices))
        allocate(potential(x,y))
        call check(nf90_get_var(nc_id_EHEP,pot_id,potential_sl),1,"get_var",name_EHEP_pot)
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_water,water_id),1,"inq_varid",name_EHEP_water)
        allocate(watermask(x,y))
        allocate(watermask_sl(x,y,time_slices))
        call check(nf90_get_var(nc_id_EHEP,water_id,watermask_sl),1,"get_var",name_EHEP_water)
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_pot_numb,potnum_id),1,"inq_varid",name_EHEP_pot_numb)
        allocate(potnum(x,y))
        allocate(potnum_sl(x,y,time_slices))
        call check(nf90_get_var(nc_id_EHEP,potnum_id,potnum_sl),1,"get_var",name_EHEP_pot_numb)
      else
        print *,"Wrong entry for switch_slices: not binary (0 or 1)! Aborting now..."
        stop 3
      end if
      
      
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_start,start_id),1,"inq_varid",name_EHEP_start)
      allocate(density(x,y))
      call check(nf90_get_var(nc_id_EHEP,start_id,density),1,"get_var",name_EHEP_start)
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_start_numb,numdens_id),1,"inq_varid",name_EHEP_start_numb)
      allocate(numdens(x,y))
      call check(nf90_get_var(nc_id_EHEP,numdens_id,numdens),1,"get_var",name_EHEP_start_numb)

      if (switch_equid == 1) then      
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_dx,dx_id),1,"inq_varid",name_EHEP_dx)
        call check(nf90_get_var(nc_id_EHEP,dx_id,dx),1,"get_var",name_EHEP_dx)
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_dy,dy_id),1,"inq_varid",name_EHEP_dy)
        call check(nf90_get_var(nc_id_EHEP,dy_id,dy),1,"get_var",name_EHEP_dy)
      else if (switch_equid == 0) then
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_dx_grid,dx_grid_id),1,"inq_varid",name_EHEP_dx_grid)
        allocate(dx_grid(x_1,y_1))
        call check(nf90_get_var(nc_id_EHEP,dx_grid_id,dx_grid),1,"get_var",name_EHEP_dx_grid)
        call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_dy_grid,dy_grid_id),1,"inq_varid",name_EHEP_dy_grid)
        allocate(dy_grid(x_1,y_1))
        call check(nf90_get_var(nc_id_EHEP,dy_grid_id,dy_grid),1,"get_var",name_EHEP_dy_grid)
      else
        print *,"Wrong entry for switch_equid: not binary (0 or 1)! Aborting now..."
        stop 3
      end if
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_posx_lat,posx_lat_id),1,"inq_varid",name_EHEP_posx_lat)
      allocate(posx_lat(x_1,y_1))
      call check(nf90_get_var(nc_id_EHEP,posx_lat_id,posx_lat),1,"get_var",name_EHEP_posx_lat)
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_posx_lon,posx_lon_id),1,"inq_varid",name_EHEP_posx_lon)
      allocate(posx_lon(x_1,y_1))
      call check(nf90_get_var(nc_id_EHEP,posx_lon_id,posx_lon),1,"get_var",name_EHEP_posy_lon)
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_posy_lat,posy_lat_id),1,"inq_varid",name_EHEP_posy_lat)
      allocate(posy_lat(x_1,y_1))
      call check(nf90_get_var(nc_id_EHEP,posy_lat_id,posy_lat),1,"get_var",name_EHEP_posy_lat)
      call check(nf90_inq_varid(nc_id_EHEP,name_EHEP_posy_lon,posy_lon_id),1,"inq_varid",name_EHEP_posy_lon)
      allocate(posy_lon(x_1,y_1))
      call check(nf90_get_var(nc_id_EHEP,posy_lon_id,posy_lon),1,"get_var",name_EHEP_posy_lon)
      ! Allocations and first calculations to set up the model

      if (switch_border == 1) then
        uppboundbord = 1
        lowboundbord = 2
      else
        uppboundbord = 2
        lowboundbord = 3
      end if
      
      allocate(avail_pot(x,y))
      dt = 1.d0/fraction_time
      timesteps = int(years*fraction_time)
      savetimes = years/save_interval
      currentyear = 0
      currentslice = 1 
      allocate(diffusion(x,y))
      allocate(minpot(x,y))
      allocate(minnum(x,y))
      allocate(wei_scale(x,y))
      allocate(Theta_tmp(x,y))
      allocate(reprod(x,y))
      Theta_tmp(:,:) = 0.d0
      minpot(:,:) = 0.d0
      minnum(:,:) = 0.d0
      wei_scale(:,:) = 1.d0
      reprod(:,:) = repro * area
      allocate(local_grid(x,y))
      local_grid(:,:) = 1.d0
      
      
      allocate(density_s(x,y))
      allocate(avail_pot_s(x,y))
      allocate(vxfield(x-1,y-1))
      allocate(vyfield(x-1,y-1))
      allocate(vxfield_old(x-1,y-1))
      allocate(vyfield_old(x-1,y-1))
      allocate(vlength(x-1,y-1))
      allocate(gradxfield(x-1,y-1))
      allocate(gradyfield(x-1,y-1))
      allocate(xi(x-1,y-1))
      allocate(time_array(savetimes+1))
      time_array = (/(i,i=0,years+1,save_interval)/)


      if (switch_slices == 0) then
        wei_maxpos = ((alpha-1.0d0)/alpha * beta **alpha)**(1./alpha)
        wei_max = (alpha/beta) * (wei_maxpos/beta) ** (alpha-1.0d0) * exp(-(wei_maxpos/beta)**alpha)
        wei_scale = potnum / wei_max
        Theta_tmp = Theta_func(numdens,potnum,minnum)
        avail_pot = weibull_func(Theta_tmp,wei_scale)
        diffusion = 0.
        !maxpot = maxval(potnum)
        !maxpot = (maxpot / 100.d0) * 2500.d0
        watermask(:,1:lowboundbord) = 0.d0
        watermask(:,y-uppboundbord:) = 0.d0
        watermask(1:lowboundbord,:) = 0.d0
        watermask(x-uppboundbord:,:) = 0.d0
      else if (switch_slices == 1) then
        watermask_sl(:,1:lowboundbord,:) = 0.d0
        watermask_sl(:,y-uppboundbord:,:) = 0.d0
        watermask_sl(1:lowboundbord,:,:) = 0.d0
        watermask_sl(x-uppboundbord:,:,:) = 0.d0
        potnum = potnum_sl(:,:,1)
        watermask = watermask_sl(:,:,1)
        wei_maxpos = ((alpha-1.0d0)/alpha * beta **alpha)**(1./alpha)
        wei_max = (alpha/beta) * (wei_maxpos/beta) ** (alpha-1.0d0) * exp(-(wei_maxpos/beta)**alpha)
        wei_scale = potnum / wei_max
        Theta_tmp = Theta_func(numdens,potnum,minnum)
        avail_pot = weibull_func(Theta_tmp,wei_scale)
        diffusion = 0.
        !maxpot = (maxpot / 100.d0) * 2500.d0
        !maxpot = maxval(potnum_sl)
      end if


      allocate(advec_t(x,y))
      allocate(diffu_t(x,y))
      allocate(birth_t(x,y))
      allocate(death_t(x,y))

      allocate(Fxp_s(x,y))
      allocate(Fxm_s(x,y))
      allocate(Fyp_s(x,y))
      allocate(Fym_s(x,y))
      allocate(Pxp_s(x,y))
      allocate(Pxm_s(x,y))
      allocate(Pyp_s(x,y))
      allocate(Pym_s(x,y))

      allocate(advec_temp(x,y,4))
      allocate(diffu_temp(x,y,4))
      allocate(birth_temp(x,y,4))
      allocate(death_temp(x,y,4))
      allocate(Fxp_tmp(x,y,4))
      allocate(Fyp_tmp(x,y,4))
      allocate(Pxp_tmp(x,y,4))
      allocate(Pyp_tmp(x,y,4))
      allocate(Fxm_tmp(x,y,4))
      allocate(Fym_tmp(x,y,4))
      allocate(Pxm_tmp(x,y,4))
      allocate(Pym_tmp(x,y,4))

      allocate(slope(x,y))
      allocate(intercept(x,y))
      
      allocate(k1(x,y))
      k1(:,:) = 0.d0
      allocate(k2(x,y))
      k2(:,:) = 0.d0
      allocate(k3(x,y))
      k3(:,:) = 0.d0
      allocate(k4(x,y))
      k4(:,:) = 0.d0

      ! Close the input file
      call check(nf90_close(nc_id_EHEP))

      

      ! Clearing the borders
      !density(:,1:2) = 0.d0
      !density(:,ubound(density,dim=2)-1:) = 0.d0
      !density(1:2,:) = 0.d0
      !density(ubound(density,dim=1)-1:,:) = 0.d0
      numdens(:,1:lowboundbord) = 0.d0
      numdens(:,y-uppboundbord:) = 0.d0
      numdens(1:lowboundbord,:) = 0.d0
      numdens(x-uppboundbord:,:) = 0.d0
      
      if (switch_slices == 0) then
        !potential(:,1:2) = 0.d0
        !potential(:,ubound(potential,dim=2)-1:) = 0.d0
        !potential(1:2,:) = 0.d0
        !potential(ubound(potential,dim=1)-1:,:) = 0.d0
        potnum(:,1:lowboundbord) = 0.d0
        potnum(:,y-uppboundbord:) = 0.d0
        potnum(1:lowboundbord,:) = 0.d0
        potnum(x-uppboundbord:,:) = 0.d0
      else if (switch_slices == 1) then
        !potential_sl(:,1:2,:) = 0.d0
        !potential_sl(:,ubound(potential,dim=2)-1:,:) = 0.d0
        !potential_sl(1:2,:,:) = 0.d0
        !potential_sl(ubound(potential,dim=1)-1:,:,:) = 0.d0
        !potential(:,1:2) = 0.d0
        !potential(:,ubound(potential,dim=2)-1:) = 0.d0
        !potential(1:2,:) = 0.d0
        !potential(ubound(potential,dim=1)-1:,:) = 0.d0
        potnum_sl(:,1:lowboundbord,:) = 0.d0
        potnum_sl(:,y-uppboundbord:,:) = 0.d0
        potnum_sl(1:lowboundbord,:,:) = 0.d0
        potnum_sl(x-uppboundbord:,:,:) = 0.d0
        potnum(:,1:lowboundbord) = 0.d0
        potnum(:,y-uppboundbord:) = 0.d0
        potnum(1:lowboundbord,:) = 0.d0
        potnum(x-uppboundbord:,:) = 0.d0
        
      ! calculation the slope and intercept between the first two slices
      !$omp parallel
      !$omp workshare
        slope(:,:) = (potnum_sl(:,:,currentslice+1) - potnum_sl(:,:,currentslice)) / (time_interval*fraction_time) ! -0
        intercept(:,:) = potnum_sl(:,:,1)
      !$omp end workshare
      !$omp end parallel
      end if
      ! Creating the initial saving values
      density_s = numdens
      
      avail_pot_s = avail_pot
      borderloss = 0.d0
      pop_los = 0.d0
      vxfield(:,:) = 0.d0
      vyfield(:,:) = 0.d0
      advec_t(:,:) = 0.d0
      diffu_t(:,:) = 0.d0
      birth_t(:,:) = 0.d0
      death_t(:,:) = 0.d0
      Fxp_s(:,:) = 0.d0
      Fxm_s(:,:) = 0.d0
      Fyp_s(:,:) = 0.d0
      Fxm_s(:,:) = 0.d0
      Pxp_s(:,:) = 0.d0
      Pxm_s(:,:) = 0.d0
      Pyp_s(:,:) = 0.d0
      Pym_s(:,:) = 0.d0

      ! Calculate the first advection estimations:

      y_do_gra: do k = 1, y-1
        x_do_gra: do j = 1, x-1
        if (switch_equid == 1) then
          gradxfield(j,k) = (avail_pot(j+1,k) - avail_pot(j,k))/dx
        else
          gradxfield(j,k) = (avail_pot(j+1,k) - avail_pot(j,k))/dx_grid(j,k)
        end if
        if (switch_equid == 1) then
          gradyfield(j,k) = (avail_pot(j,k+1) - avail_pot(j,k))/dy
        else
          gradyfield(j,k) = (avail_pot(j,k+1) - avail_pot(j,k))/dy_grid(j,k)
        end if

        if (abs(gradxfield(j,k)) <= Eps .AND. abs(gradyfield(j,k)) > Eps) then
          xi(j,k) = (PI/2.d0)
        else if (abs(gradyfield(j,k)) <= Eps .AND. abs(gradxfield(j,k)) > Eps) then
          xi(j,k) = 0.d0
        else if (abs(gradyfield(j,k)) <= Eps .AND. abs(gradxfield(j,k)) <= Eps) then
          xi(j,k) = (PI/4.d0)
        else
          xi(j,k) = atan(-(gradyfield(j,k)/gradxfield(j,k)))
        end if
        end do x_do_gra
      end do y_do_gra

      !$omp do
      y_do_dif: do k = 2, y-1
        x_do_dif: do j = 2, x-1
          if (switch_equid == 1) then
            local_grid(j,k) = (dx+dy) / 2.d0
          else
            local_grid(j,k) = (dx_grid(j,k)+dx_grid(j-1,k)+dy_grid(j,k)+dy_grid(j,k-1)) / 4.d0
          end if
        end do x_do_dif
      end do y_do_dif
      !$omp end do

      !maxpot_grad = max(maxval(gradxfield),maxval(gradyfield))
      !if (switch_equid == 1) then
      !  maxpot = (maxpot / 100.d0) * (sum(area)/size(area))   ! changing the dimension to number of humans 
      !  maxpot = maxpot / ((dx+dy)/2.d0)
      !else
      !  maxpot = (maxpot / 100.d0) * (sum(area)/size(area))   ! changing the dimension to number of humans 
      !  maxpot = maxpot / (((sum(dx_grid)/size(dx_grid))+(sum(dy_grid)/size(dy_grid)))/2.d0)
      !end if
      print *,"maxpot: ",maxpot
      gamma_para = vmax / (maxpot * deltatau)
      print *,"gamma: ",gamma_para
      print *,"dt: ",dt
      vxfield_old = vxfield
      vyfield_old = vyfield
      
      !$omp do
      y_do_adv: do k = 1, y-1
      x_do_adv: do j = 1, x-1
        if (((k <= 2) .OR. (j<=2)) .OR. ((k >= y-2) .OR. (j >= x-2))) then !Boundary conditions
          vxfield(j,k) = 0.d0
          vyfield(j,k) = 0.d0
        else
          vxfield(j,k) = gamma_para*gradxfield(j,k) !+ vxfield_old(j,k)
          vyfield(j,k) = gamma_para*gradyfield(j,k) !+ vyfield_old(j,k)
          !if (vxfield(j,k) /= 0.d0) then
          !  print *,"Reached ",j,k
          !  print *,"X:",vxfield(j,k)
          !else if (vyfield(j,k) /= 0.d0) then
          !  print *,"Reached ",j,k
          !  print *,"Y:",vyfield(j,k)
          !end if
        end if

        vlength(j,k) = sqrt(vxfield(j,k)**2 + vyfield(j,k)**2)
        if (vlength(j,k) > vmax) then
          vxfield(j,k) = (vxfield(j,k) / vlength(j,k)) * vmax
          vyfield(j,k) = (vyfield(j,k) / vlength(j,k)) * vmax
        end if
          
        !if (vxfield(j,k) > vmax) then
        !  vxfield(j,k) = vmax
        !else if (vxfield(j,k) < -vmax) then
        !  vxfield(j,k) = -vmax
        !end if
        !if (vyfield(j,k) > vmax) then
        !  vyfield(j,k) = vmax
        !else if (vyfield(j,k) < -vmax) then
        !  vyfield(j,k) = -vmax

        !if (vxfield(j,k) > 0.d0) then
        !  print *,"vxfield found! ",vxfield(j,k)
        !end if
        !end if
      end do x_do_adv
      end do y_do_adv
      !$omp end do
      
      ! Create the savefile, dimensions, variables and attributes
      call check(nf90_create(path_output,NF90_NETCDF4,nc_id_out),2,"create","path_output")
      call check(nf90_def_grp(nc_id_out,output_group_name,paramgrp_id_out),2,"def_grp",output_group_name)
      call check(nf90_def_dim(nc_id_out,name_output_dim_x,x,x_id_out),2,"def_dim",name_output_dim_x)
      call check(nf90_def_dim(nc_id_out,name_output_dim_y,y,y_id_out),2,"def_dim",name_output_dim_y)
      call check(nf90_def_dim(nc_id_out,name_output_dim_t,savetimes+1,t_id_out),2,"def_dim",name_output_dim_t)
      call check(nf90_def_dim(nc_id_out,name_output_dim_x_1,x_1,x_1_id_out),2,"def_dim",name_output_dim_x_1)
      call check(nf90_def_dim(nc_id_out,name_output_dim_y_1,y_1,y_1_id_out),2,"def_dim",name_output_dim_y_1)
      call check(nf90_def_dim(paramgrp_id_out,name_output_dim_param,1,param_dim_id),2,"def_dim",name_output_dim_param)
      dim_ids_all = (/ x_id_out, y_id_out, t_id_out /)
      dim_ids_vel = (/ x_1_id_out, y_1_id_out,t_id_out/)
      dim_latlon = (/ x_id_out, y_id_out /)
      dim_latlon_1 = (/x_1_id_out, y_1_id_out /)
      call check(nf90_def_var(nc_id_out,name_output_var_lat,NF90_DOUBLE,dim_latlon,lat_id_out),2,"def_var",name_output_var_lat)
      call check(nf90_def_var(nc_id_out,name_output_var_lon,NF90_DOUBLE,dim_latlon,lon_id_out),2,"def_var",name_output_var_lon)
      if (switch_equid == 1) then
        call check(nf90_def_var(paramgrp_id_out,name_output_var_dx,NF90_DOUBLE,param_dim_id,dx_id_out)&
,2,"def_var",name_output_var_dx)
        call check(nf90_def_var(paramgrp_id_out,name_output_var_dy,NF90_DOUBLE,param_dim_id,dy_id_out)&
,2,"def_var",name_output_var_dy)
      else
        call check(nf90_def_var(nc_id_out,name_output_var_dx_grid,NF90_DOUBLE,dim_latlon_1,dx_grid_id_out)&
,2,"def_var",name_output_var_dx_grid)
        call check(nf90_def_var(nc_id_out,name_output_var_dy_grid,NF90_DOUBLE,dim_latlon_1,dy_grid_id_out)&
,2,"def_var",name_output_var_dy_grid)
      end if
      call check(nf90_def_var(nc_id_out,name_output_var_posx_lat,NF90_DOUBLE,dim_latlon_1,posx_lat_id_out)&
,2,"def_var",name_output_var_posx_lat)
      call check(nf90_def_var(nc_id_out,name_output_var_posx_lon,NF90_DOUBLE,dim_latlon_1,posx_lon_id_out)&
,2,"def_var",name_output_var_posx_lon)
      call check(nf90_def_var(nc_id_out,name_output_var_posy_lat,NF90_DOUBLE,dim_latlon_1,posy_lat_id_out)&
,2,"def_var",name_output_var_posy_lat)
      call check(nf90_def_var(nc_id_out,name_output_var_posy_lon,NF90_DOUBLE,dim_latlon_1,posy_lon_id_out)&
,2,"def_var",name_output_var_posy_lon)
      call check(nf90_def_var(nc_id_out,name_output_var_time,NF90_DOUBLE,t_id_out,time_id),2,"def_var",name_output_var_time)

      call check(nf90_def_var(nc_id_out,name_output_var_density,NF90_DOUBLE,dim_ids_all,density_id)&
,2,"def_var",name_output_var_density)
      call check(nf90_def_var(nc_id_out,name_output_var_AvailPot,NF90_DOUBLE,dim_ids_all,availpot_id)&
,2,"def_var",name_output_var_AvailPot)
      call check(nf90_def_var(nc_id_out,name_output_var_Border,NF90_DOUBLE,t_id_out,borderloss_id)&
,2,"def_var",name_output_var_Border)
      call check(nf90_def_var(nc_id_out,name_output_var_DiffPara,NF90_DOUBLE,dim_ids_all,diffpara_id_out)&
,2,"def_var",name_output_var_DiffPara)
      call check(nf90_def_var(nc_id_out,name_output_var_Velocity_x,NF90_DOUBLE,dim_ids_vel,velocity_x_id)&
,2,"def_var",name_output_var_Velocity_x)
      call check(nf90_def_var(nc_id_out,name_output_var_Velocity_y,NF90_DOUBLE,dim_ids_vel,velocity_y_id)&
,2,"def_var",name_output_var_Velocity_y)
      call check(nf90_def_var(nc_id_out,name_output_var_Advection,NF90_DOUBLE,dim_ids_all,advection_id)&
,2,"def_var",name_output_var_Advection)
      call check(nf90_def_var(nc_id_out,name_output_var_Diffusion,NF90_DOUBLE,dim_ids_all,diffusion_id)&
,2,"def_var",name_output_var_Diffusion)
      call check(nf90_def_var(nc_id_out,name_output_var_Birth,NF90_DOUBLE,dim_ids_all,birth_id)&
,2,"def_var",name_output_var_Birth)
      call check(nf90_def_var(nc_id_out,name_output_var_Death,NF90_DOUBLE,dim_ids_all,death_id)&
,2,"def_var",name_output_var_Death) 
      call check(nf90_def_var(nc_id_out,name_output_var_Area,NF90_DOUBLE,dim_latlon,area_id_out)&
,2,"def_var",name_output_var_Area)
      call check(nf90_def_var(nc_id_out,name_output_var_NumbHuman,NF90_DOUBLE,dim_ids_all,numdens_id_out)&
,2,"def_var",name_output_var_NumbHuman)
      if (switch_slices == 0) then
        call check(nf90_def_var(nc_id_out,name_output_var_watermask,NF90_INT,dim_latlon,water_id_out)&
,2,"def_var",name_output_var_watermask)
        call check(nf90_def_var(nc_id_out,name_output_var_EHEP,NF90_DOUBLE,dim_latlon,ehep_id_out)&
,2,"def_var",name_output_var_EHEP)
      else if (switch_slices == 1) then
        call check(nf90_def_var(nc_id_out,name_output_var_watermask,NF90_INT,dim_ids_all,water_id_out)&
,2,"def_var",name_output_var_watermask)
        call check(nf90_def_var(nc_id_out,name_output_var_EHEP,NF90_DOUBLE,dim_ids_all,ehep_id_out)&
,2,"def_var",name_output_var_EHEP)
      end if
      call check(nf90_def_var(nc_id_out,name_output_var_xi,NF90_DOUBLE,dim_ids_vel,xi_id_out)&
,2,"def_var",name_output_var_xi)
      call check(nf90_def_var(nc_id_out,name_output_var_gradx,NF90_DOUBLE,dim_ids_vel,gradx_id_out)&
,2,"def_var",name_output_var_gradx)
      call check(nf90_def_var(nc_id_out,name_output_var_grady,NF90_DOUBLE,dim_ids_vel,grady_id_out)&
,2,"def_var",name_output_var_grady)

      call check(nf90_def_var(nc_id_out,name_output_var_advflux_x_p,NF90_DOUBLE,dim_ids_all,Fxp_id_out)&
,2,"def_var",name_output_var_advflux_x_p)
      call check(nf90_def_var(nc_id_out,name_output_var_advflux_x_m,NF90_DOUBLE,dim_ids_all,Fxm_id_out)&
,2,"def_var",name_output_var_advflux_x_m)
      call check(nf90_def_var(nc_id_out,name_output_var_advflux_y_p,NF90_DOUBLE,dim_ids_all,Fyp_id_out)&
,2,"def_var",name_output_var_advflux_y_p)
      call check(nf90_def_var(nc_id_out,name_output_var_advflux_y_m,NF90_DOUBLE,dim_ids_all,Fym_id_out)&
,2,"dev_var",name_output_var_advflux_y_m)
      call check(nf90_def_var(nc_id_out,name_output_var_diffflux_x_p,NF90_DOUBLE,dim_ids_all,Pxp_id_out)&
,2,"def_var",name_output_var_diffflux_x_p)
      call check(nf90_def_var(nc_id_out,name_output_var_diffflux_x_m,NF90_DOUBLE,dim_ids_all,Pxm_id_out)&
,2,"def_var",name_output_var_diffflux_x_m)
      call check(nf90_def_var(nc_id_out,name_output_var_diffflux_y_p,NF90_DOUBLE,dim_ids_all,Pyp_id_out)&
,2,"def_var",name_output_var_diffflux_y_p)
      call check(nf90_def_var(nc_id_out,name_output_var_diffflux_y_m,NF90_DOUBLE,dim_ids_all,Pym_id_out)&
,2,"dev_var",name_output_var_diffflux_y_m)

      call check(nf90_def_var(nc_id_out,name_output_var_pop_los,NF90_DOUBLE,t_id_out,pop_los_id_out)&
,2,"def_var",name_output_var_pop_los)
      
      call check(nf90_def_var(paramgrp_id_out,name_output_param_adv,NF90_INT,param_dim_id,param_adv_id)&
,2,"def_var",name_output_param_adv)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_birth,NF90_INT,param_dim_id,param_birth_id)&
,2,"def_var",name_output_param_birth)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_dynpot,NF90_INT,param_dim_id,param_dynpot_id)&
,2,"def_var",name_output_param_dynpot)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_equidis,NF90_INT,param_dim_id,param_equidis_id)&
,2,"def_var",name_output_param_equidis)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_numtime,NF90_INT,param_dim_id,param_numtime_id)&
,2,"def_var",name_output_param_numtime)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_time,NF90_INT,param_dim_id,param_time_id)&
,2,"def_var",name_output_param_time)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_saveint,NF90_INT,param_dim_id,param_saveint_id)&
,2,"def_var",name_output_param_saveint)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_potint,NF90_INT,param_dim_id,param_potint_id)&
,2,"def_var",name_output_param_potint)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_vmax,NF90_DOUBLE,param_dim_id,param_vmax_id)&
,2,"def_var",name_output_param_vmax)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_growth,NF90_DOUBLE,param_dim_id,param_growth_id)&
,2,"dev_var",name_output_param_growth)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_vcrt,NF90_DOUBLE,param_dim_id,param_vcrt_id)&
,2,"dev_var",name_output_param_vcrt)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_dtau,NF90_DOUBLE,param_dim_id,param_dtau_id)&
,2,"dev_var",name_output_param_dtau)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_diffc,NF90_DOUBLE,param_dim_id,param_diffc_id)&
,2,"dev_var",name_output_param_diffc)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_theta,NF90_DOUBLE,param_dim_id,param_theta_id)&
,2,"dev_var",name_output_param_theta)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_repro,NF90_DOUBLE,param_dim_id,param_repro_id)&
,2,"dev_var",name_output_param_repro)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_alpha,NF90_DOUBLE,param_dim_id,param_alpha_id)&
,2,"dev_var",name_output_param_alpha)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_beta,NF90_DOUBLE,param_dim_id,param_beta_id)&
,2,"dev_var",name_output_param_beta)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_eta,NF90_DOUBLE,param_dim_id,param_eta_id)&
,2,"dev_var",name_output_param_eta)
      call check(nf90_def_var(paramgrp_id_out,name_output_param_maxpot,NF90_DOUBLE,param_dim_id,param_maxpot_id)&
,2,"dev_var",name_output_param_maxpot)


      call check(nf90_put_att(nc_id_out,NF90_GLOBAL,output_name_description,output_description)&
,2,"put_att","output_description")
      call check(nf90_put_att(nc_id_out,lat_id_out,output_name_longname,output_long_lat)&
,2,"put_att","output_long_lat")
      call check(nf90_put_att(nc_id_out,lat_id_out,output_name_unit,output_unit_lat)&
,2,"put_att","output_unit_lat")
      call check(nf90_put_att(nc_id_out,lon_id_out,output_name_longname,output_long_lon)&
,2,"put_att","output_long_lon")
      call check(nf90_put_att(nc_id_out,lon_id_out,output_name_unit,output_unit_lon)&
,2,"put_att","output_unit_lon")
      call check(nf90_put_att(nc_id_out,posx_lat_id_out,output_name_longname,output_long_posx_lat)&
,2,"put_att","output_long_posx_lat")
      call check(nf90_put_att(nc_id_out,posx_lat_id_out,output_name_unit,output_unit_posx_lat)&
,2,"put_att","output_unit_posx_lat")
      call check(nf90_put_att(nc_id_out,posx_lon_id_out,output_name_longname,output_long_posx_lon)&
,2,"put_att","output_long_posx_lon")
      call check(nf90_put_att(nc_id_out,posx_lon_id_out,output_name_unit,output_unit_posx_lon)&
,2,"put_att","output_unit_posx_lon")
      call check(nf90_put_att(nc_id_out,posy_lat_id_out,output_name_longname,output_long_posy_lat)&
,2,"put_att","output_long_posy_lat")
      call check(nf90_put_att(nc_id_out,posy_lat_id_out,output_name_unit,output_unit_posy_lat)&
,2,"put_att","output_unit_posy_lat")
      call check(nf90_put_att(nc_id_out,posy_lon_id_out,output_name_longname,output_long_posy_lon)&
,2,"put_att","output_long_posy_lon")
      call check(nf90_put_att(nc_id_out,posy_lon_id_out,output_name_unit,output_unit_posy_lon)&
,2,"put_att","output_unit_posy_lon")
      if (switch_equid == 1) then
        call check(nf90_put_att(nc_id_out,dx_id_out,output_name_longname,output_long_dx),2,"put_att","output_long_dx")
        call check(nf90_put_att(nc_id_out,dy_id_out,output_name_longname,output_long_dy),2,"put_att","output_long_dy")
        call check(nf90_put_att(nc_id_out,dx_id_out,output_name_unit,output_unit_dx),2,"put_att","output_unit_dx")
        call check(nf90_put_att(nc_id_out,dy_id_out,output_name_unit,output_unit_dy),2,"put_att","output_unit_dy")
      else
        call check(nf90_put_att(nc_id_out,dx_grid_id_out,output_name_longname,output_long_dx_grid)&
,2,"put_att","output_long_dx_grid")
        call check(nf90_put_att(nc_id_out,dx_grid_id_out,output_name_unit,output_unit_dx_grid)&
,2,"put_att","output_unit_dx_grid")
        call check(nf90_put_att(nc_id_out,dy_grid_id_out,output_name_longname,output_long_dy_grid)&
,2,"put_att","output_long_dy_grid")
        call check(nf90_put_att(nc_id_out,dy_grid_id_out,output_name_unit,output_unit_dy_grid)&
,2,"put_att","output_unit_dy_grid")
      end if

      call check(nf90_put_att(nc_id_out,time_id,output_name_longname,output_long_time)&
,2,"put_att","output_long_time")
      call check(nf90_put_att(nc_id_out,time_id,output_name_unit,output_unit_time)&
,2,"put_att","output_unit_time")
      call check(nf90_put_att(nc_id_out,density_id,output_name_longname,output_long_density)&
,2,"put_att","output_long_density")
      call check(nf90_put_att(nc_id_out,density_id,output_name_unit,output_unit_density)&
,2,"put_att","output_unit_density")
      call check(nf90_put_att(nc_id_out,availpot_id,output_name_longname,output_long_availpot)&
,2,"put_att","output_long_availpot")
      call check(nf90_put_att(nc_id_out,availpot_id,output_name_unit,output_unit_availpot)&
,2,"put_att","output_unit_availpot")
      call check(nf90_put_att(nc_id_out,borderloss_id,output_name_longname,output_long_borderloss)&
,2,"put_att","output_long_borderloss")
      call check(nf90_put_att(nc_id_out,borderloss_id,output_name_unit,output_unit_borderloss)&
,2,"put_att","output_unit_borderloss")
      call check(nf90_put_att(nc_id_out,diffpara_id_out,output_name_longname,output_long_DiffPara)&
,2,"put_att","output_long_DiffPara")
      call check(nf90_put_att(nc_id_out,diffpara_id_out,output_name_unit,output_unit_DiffPara)&
,2,"put_att","output_unit_DiffPara") 
      call check(nf90_put_att(nc_id_out,area_id_out,output_name_longname,output_long_Area)&
,2,"put_att","output_long_Area")
      call check(nf90_put_att(nc_id_out,area_id_out,output_name_unit,output_unit_Area)&
,2,"put_att","output_unit_Area")
      call check(nf90_put_att(nc_id_out,numdens_id_out,output_name_longname,output_long_NumbHuman)&
,2,"put_att","output_long_NumbHuman")
      call check(nf90_put_att(nc_id_out,numdens_id_out,output_name_unit,output_unit_NumbHuman)&
,2,"put_att","output_unit_NumbHuman")

      call check(nf90_put_att(nc_id_out,velocity_x_id,output_name_longname,output_long_velocity_x)&
,2,"put_att","output_long_velocity_x")
      call check(nf90_put_att(nc_id_out,velocity_x_id,output_name_unit,output_unit_velocity_x)&
,2,"put_att","output_unit_velocity_x")
      call check(nf90_put_att(nc_id_out,velocity_y_id,output_name_longname,output_long_velocity_y)&
,2,"put_att","output_long_velocity_y")
      call check(nf90_put_att(nc_id_out,velocity_y_id,output_name_unit,output_unit_velocity_y)&
,2,"put_att","output_unit_velocity_y")
      call check(nf90_put_att(nc_id_out,advection_id,output_name_longname,output_long_advection)&
,2,"put_att","output_long_advection")
      call check(nf90_put_att(nc_id_out,advection_id,output_name_unit,output_unit_advection)&
,2,"put_att","output_unit_advection")
      call check(nf90_put_att(nc_id_out,diffusion_id,output_name_longname,output_long_diffusion)&
,2,"put_att","output_long_diffusion")
      call check(nf90_put_att(nc_id_out,diffusion_id,output_name_unit,output_unit_advection)&
,2,"put_att","output_unit_advection")
      call check(nf90_put_att(nc_id_out,birth_id,output_name_longname,output_long_birth)&
,2,"put_att","output_long_birth")
      call check(nf90_put_att(nc_id_out,birth_id,output_name_unit,output_unit_birth)&
,2,"put_att","output_unit_birth")
      call check(nf90_put_att(nc_id_out,death_id,output_name_longname,output_long_death)&
,2,"put_att","output_long_death")
      call check(nf90_put_att(nc_id_out,death_id,output_name_unit,output_unit_death)&
,2,"put_att","output_unit_death")
      call check(nf90_put_att(nc_id_out,water_id_out,output_name_longname,output_long_watermask)&
,2,"put_att","output_long_watermask")
      call check(nf90_put_att(nc_id_out,water_id_out,output_name_unit,output_unit_watermask)&
,2,"put_att","output_unit_watermask")
      call check(nf90_put_att(nc_id_out,ehep_id_out,output_name_longname,output_long_EHEP)&
,2,"put_att","output_long_EHEP")
      call check(nf90_put_att(nc_id_out,ehep_id_out,output_name_unit,output_unit_EHEP)&
,2,"put_att","output_unit_EHEP")
      call check(nf90_put_att(nc_id_out,xi_id_out,output_name_longname,output_long_xi)&
,2,"put_att","output_long_xi")
      call check(nf90_put_att(nc_id_out,xi_id_out,output_name_unit,output_unit_xi)&
,2,"put_att","output_unit_xi")

      call check(nf90_put_att(nc_id_out,Fxp_id_out,output_name_longname,output_long_advflux_x_p)&
,2,"put_att","output_long_advflux_x_p")
      call check(nf90_put_att(nc_id_out,Fxp_id_out,output_name_unit,output_unit_advflux_x_p)&
,2,"put_att","output_unit_advflux_x_p")
      call check(nf90_put_att(nc_id_out,Fxm_id_out,output_name_longname,output_long_advflux_x_m)&
,2,"put_att","output_long_advflux_x_m")
      call check(nf90_put_att(nc_id_out,Fxm_id_out,output_name_unit,output_unit_advflux_x_m)&
,2,"put_att","output_unit_advflux_x_m")
      call check(nf90_put_att(nc_id_out,Fyp_id_out,output_name_longname,output_long_advflux_y_p)&
,2,"put_att","output_long_advflux_y_p")
      call check(nf90_put_att(nc_id_out,Fyp_id_out,output_name_unit,output_unit_advflux_y_p)&
,2,"put_att","output_unit_advflux_y_p")
      call check(nf90_put_att(nc_id_out,Fym_id_out,output_name_longname,output_long_advflux_y_m)&
,2,"put_att","output_long_advflux_y_m")
      call check(nf90_put_att(nc_id_out,Fym_id_out,output_name_unit,output_unit_advflux_y_m)&
,2,"put_att","output_unit_advflux_y_m")
      call check(nf90_put_att(nc_id_out,Pxp_id_out,output_name_longname,output_long_diffflux_x_p)&
,2,"put_att","output_long_diffflux_x_p")
      call check(nf90_put_att(nc_id_out,Pxp_id_out,output_name_unit,output_unit_diffflux_x_p)&
,2,"put_att","output_unit_diffflux_x_p")
      call check(nf90_put_att(nc_id_out,Pxm_id_out,output_name_longname,output_long_diffflux_x_m)&
,2,"put_att","output_long_diffflux_x_m")
      call check(nf90_put_att(nc_id_out,Pxm_id_out,output_name_unit,output_unit_diffflux_x_m)&
,2,"put_att","output_unit_diffflux_x_m")
      call check(nf90_put_att(nc_id_out,Pyp_id_out,output_name_longname,output_long_diffflux_y_p)&
,2,"put_att","output_long_diffflux_y_p")
      call check(nf90_put_att(nc_id_out,Pyp_id_out,output_name_unit,output_unit_diffflux_y_p)&
,2,"put_att","output_unit_diffflux_y_p")
      call check(nf90_put_att(nc_id_out,Pym_id_out,output_name_longname,output_long_diffflux_y_m)&
,2,"put_att","output_long_diffflux_y_m")
      call check(nf90_put_att(nc_id_out,Pym_id_out,output_name_unit,output_unit_diffflux_y_m)&
,2,"put_att","output_unit_diffflux_y_m")

      call check(nf90_put_att(nc_id_out,pop_los_id_out,output_name_longname,output_long_pop_los)&
,2,"put_att","output_long_pop_los")
      call check(nf90_put_att(nc_id_out,pop_los_id_out,output_name_unit,output_unit_pop_los)&
,2,"put_att","output_unit_pop_los")

      call check(nf90_put_att(paramgrp_id_out,param_adv_id,output_name_longname,output_long_param_adv)&
,2,"put_att","output_long_param_adv")
      call check(nf90_put_att(paramgrp_id_out,param_adv_id,output_name_unit,output_unit_param_adv)&
,2,"put_att","output_unit_param_adv")
      call check(nf90_put_att(paramgrp_id_out,param_birth_id,output_name_longname,output_long_param_birth)&
,2,"put_att","output_long_param_birth")
      call check(nf90_put_att(paramgrp_id_out,param_birth_id,output_name_unit,output_unit_param_birth)&
,2,"put_att","output_unit_param_birth")
      call check(nf90_put_att(paramgrp_id_out,param_dynpot_id,output_name_longname,output_long_param_dynpot)&
,2,"put_att","output_long_param_dynpot")
      call check(nf90_put_att(paramgrp_id_out,param_dynpot_id,output_name_unit,output_unit_param_dynpot)&
,2,"put_att","output_unit_param_dynpot")
      call check(nf90_put_att(paramgrp_id_out,param_equidis_id,output_name_longname,output_long_param_equidis)&
,2,"put_att","output_long_param_equidis")
      call check(nf90_put_att(paramgrp_id_out,param_equidis_id,output_name_unit,output_unit_param_equidis)&
,2,"put_att","output_unit_param_equidis")
      call check(nf90_put_att(paramgrp_id_out,param_numtime_id,output_name_longname,output_long_param_numtime)&
,2,"put_att","output_long_param_numtime")
      call check(nf90_put_att(paramgrp_id_out,param_numtime_id,output_name_unit,output_unit_param_numtime)&
,2,"put_att","output_unit_param_numtime")
      call check(nf90_put_att(paramgrp_id_out,param_time_id,output_name_longname,output_long_param_time)&
,2,"put_att","output_long_param_time")
      call check(nf90_put_att(paramgrp_id_out,param_time_id,output_name_unit,output_unit_param_time)&
,2,"put_att","output_unit_param_time")
      call check(nf90_put_att(paramgrp_id_out,param_saveint_id,output_name_longname,output_long_param_saveint)&
,2,"put_att","output_long_param_saveint")
      call check(nf90_put_att(paramgrp_id_out,param_saveint_id,output_name_unit,output_unit_param_saveint)&
,2,"put_att","output_unit_param_saveint")
      call check(nf90_put_att(paramgrp_id_out,param_potint_id,output_name_longname,output_long_param_potint)&
,2,"put_att","output_long_param_potint")
      call check(nf90_put_att(paramgrp_id_out,param_potint_id,output_name_unit,output_unit_param_potint)&
,2,"put_att","output_unit_param_potint")
      call check(nf90_put_att(paramgrp_id_out,param_vmax_id,output_name_longname,output_long_param_vmax)&
,2,"put_att","output_long_param_vmax")
      call check(nf90_put_att(paramgrp_id_out,param_vmax_id,output_name_unit,output_unit_param_vmax)&
,2,"put_att","output_unit_param_vmax")
      call check(nf90_put_att(paramgrp_id_out,param_growth_id,output_name_longname,output_long_param_growth)&
,2,"put_att","output_long_param_growth")
      call check(nf90_put_att(paramgrp_id_out,param_growth_id,output_name_unit,output_unit_param_growth)&
,2,"put_att","output_unit_param_growth")
      call check(nf90_put_att(paramgrp_id_out,param_vcrt_id,output_name_longname,output_long_param_vcrt)&
,2,"put_att","output_long_param_vcrt")
      call check(nf90_put_att(paramgrp_id_out,param_vcrt_id,output_name_unit,output_unit_param_vcrt)&
,2,"put_att","output_unit_param_vcrt")
      call check(nf90_put_att(paramgrp_id_out,param_dtau_id,output_name_longname,output_long_param_dtau)&
,2,"put_att","output_long_param_dtau")
      call check(nf90_put_att(paramgrp_id_out,param_dtau_id,output_name_unit,output_unit_param_dtau)&
,2,"put_att","output_unit_param_dtau")
      call check(nf90_put_att(paramgrp_id_out,param_diffc_id,output_name_longname,output_long_param_diffc)&
,2,"put_att","output_long_param_diffc")
      call check(nf90_put_att(paramgrp_id_out,param_diffc_id,output_name_unit,output_unit_param_diffc)&
,2,"put_att","output_unit_param_diffc")
      call check(nf90_put_att(paramgrp_id_out,param_theta_id,output_name_longname,output_long_param_theta)&
,2,"put_att","output_long_param_theta")
      call check(nf90_put_att(paramgrp_id_out,param_theta_id,output_name_unit,output_unit_param_theta)&
,2,"put_att","output_unit_param_theta")
      call check(nf90_put_att(paramgrp_id_out,param_repro_id,output_name_longname,output_long_param_repro)&
,2,"put_att","output_long_param_repro")
      call check(nf90_put_att(paramgrp_id_out,param_repro_id,output_name_unit,output_unit_param_repro)&
,2,"put_att","output_unit_param_repro")
      call check(nf90_put_att(paramgrp_id_out,param_alpha_id,output_name_longname,output_long_param_alpha)&
,2,"put_att","output_long_param_alpha")
      call check(nf90_put_att(paramgrp_id_out,param_alpha_id,output_name_unit,output_unit_param_alpha)&
,2,"put_att","output_unit_param_alpha")
      call check(nf90_put_att(paramgrp_id_out,param_beta_id,output_name_longname,output_long_param_beta)&
,2,"put_att","output_long_param_beta")
      call check(nf90_put_att(paramgrp_id_out,param_beta_id,output_name_unit,output_unit_param_beta)&
,2,"put_att","output_unit_param_beta")
      call check(nf90_put_att(paramgrp_id_out,param_eta_id,output_name_longname,output_long_param_eta)&
,2,"put_att","output_long_param_eta")
      call check(nf90_put_att(paramgrp_id_out,param_eta_id,output_name_unit,output_unit_param_eta)&
,2,"put_att","output_unit_param_eta")
      call check(nf90_put_att(paramgrp_id_out,param_maxpot_id,output_name_longname,output_long_param_maxpot)&
,2,"put_att","output_long_param_maxpot")
      call check(nf90_put_att(paramgrp_id_out,param_maxpot_id,output_name_unit,output_unit_param_maxpot)&
,2,"put_att","output_unit_param_maxpot")

      ! Declare that all definitions of the file are finished, only data is filled in the future
      call check(nf90_enddef(nc_id_out),2,"enddef")
      
      ! Save the first (initial) values
      call check(nf90_put_var(nc_id_out,density_id,density,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","density")
      call check(nf90_put_var(nc_id_out,availpot_id,avail_pot,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","avail_pot")
      call check(nf90_put_var(nc_id_out,borderloss_id,borderloss,start=(/1/))&
,3,"put_var","borderloss")
      call check(nf90_put_var(nc_id_out,diffpara_id_out,diffusion,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","diffpara") 
      call check(nf90_put_var(nc_id_out,pop_los_id_out,pop_los,start=(/1/))&
,3,"put_var","pop_los")
      call check(nf90_put_var(nc_id_out,lat_id_out,lat,start=(/1,1/),count=(/x,y/))&
,3,"put_var","lat")
      call check(nf90_put_var(nc_id_out,lon_id_out,lon,start=(/1,1/),count=(/x,y/))&
,3,"put_var","lon")
      call check(nf90_put_var(nc_id_out,time_id,time_array,start=(/1/),count=(/savetimes+1/))&
,3,"put_var","time_array")
      call check(nf90_put_var(nc_id_out,area_id_out,area,start=(/1,1/),count=(/x,y/))&
,3,"put_var","area")
      call check(nf90_put_var(nc_id_out,numdens_id_out,density_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","density_s")
      if (switch_equid == 1) then
        call check(nf90_put_var(nc_id_out,dx_id_out,dx)&
,3,"put_var","dx")
        call check(nf90_put_var(nc_id_out,dy_id_out,dy)&
,3,"put_var","dy")
      else
        call check(nf90_put_var(nc_id_out,dx_grid_id_out,dx_grid,start=(/1,1/),count=(/x_1,y_1/))&
,3,"put_var","dx_grid")
        call check(nf90_put_var(nc_id_out,dy_grid_id_out,dy_grid,start=(/1,1/),count=(/x_1,y_1/))&
,3,"put_var","dy_grid")
      end if
      call check(nf90_put_var(nc_id_out,posx_lat_id_out,posx_lat,start=(/1,1/),count=(/x_1,y_1/))&
,3,"put_var","posx_lat")
      call check(nf90_put_var(nc_id_out,posx_lon_id_out,posx_lon,start=(/1,1/),count=(/x_1,y_1/))&
,3,"put_var","posx_lon")
      call check(nf90_put_var(nc_id_out,posy_lat_id_out,posy_lat,start=(/1,1/),count=(/x_1,y_1/))&
,3,"put_var","posy_lat")
      call check(nf90_put_var(nc_id_out,posy_lon_id_out,posy_lon,start=(/1,1/),count=(/x_1,y_1/))&
,3,"put_var","posy_lon")
      call check(nf90_put_var(nc_id_out,velocity_x_id,vxfield,start=(/1,1,1/),count=(/x_1,y_1,1/))&
,3,"put_var","vxfield")
      call check(nf90_put_var(nc_id_out,velocity_y_id,vyfield,start=(/1,1,1/),count=(/x_1,y_1,1/))&
,3,"put_var","vyfield")
      call check(nf90_put_var(nc_id_out,advection_id,advec_t,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","advec_t")
      call check(nf90_put_var(nc_id_out,diffusion_id,diffu_t,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","diffu_t")
      call check(nf90_put_var(nc_id_out,birth_id,birth_t,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","birth_t")
      call check(nf90_put_var(nc_id_out,death_id,death_t,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","death_t")
      call check(nf90_put_var(nc_id_out,xi_id_out,cos(xi),start=(/1,1,1/),count=(/x_1,y_1,1/))&
,3,"put_var","xi")
      call check(nf90_put_var(nc_id_out,gradx_id_out,gradxfield,start=(/1,1,1/),count=(/x_1,y_1,1/))&
,3,"put_var","gradxfield")
      call check(nf90_put_var(nc_id_out,grady_id_out,gradyfield,start=(/1,1,1/),count=(/x_1,y_1,1/))&
,3,"put_var","gradyfield")

      call check(nf90_put_var(nc_id_out,Fxp_id_out,Fxp_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Fxp_s")
      call check(nf90_put_var(nc_id_out,Fyp_id_out,Fyp_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Fyp_s")
      call check(nf90_put_var(nc_id_out,Pxp_id_out,Pxp_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Pxp_s")
      call check(nf90_put_var(nc_id_out,Pyp_id_out,Pyp_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Pyp_s")
      call check(nf90_put_var(nc_id_out,Fxm_id_out,Fxm_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Fxm_s")
      call check(nf90_put_var(nc_id_out,Fym_id_out,Fym_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Fym_s")
      call check(nf90_put_var(nc_id_out,Pxm_id_out,Pxm_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Pxm_s")
      call check(nf90_put_var(nc_id_out,Pym_id_out,Pym_s,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","Pym_s")

      call check(nf90_put_var(paramgrp_id_out,param_adv_id,switch_adv)&
,3,"put_var","switch_adv")
      call check(nf90_put_var(paramgrp_id_out,param_birth_id,switch_birth)&
,3,"put_var","switch_birth")
      call check(nf90_put_var(paramgrp_id_out,param_dynpot_id,switch_slices)&
,3,"put_var","switch_slices")
      call check(nf90_put_var(paramgrp_id_out,param_equidis_id,switch_equid)&
,3,"put_var","switch_equid")
      call check(nf90_put_var(paramgrp_id_out,param_numtime_id,fraction_time)&
,3,"put_var","fraction_time")
      call check(nf90_put_var(paramgrp_id_out,param_time_id,years)&
,3,"put_var","years")
      call check(nf90_put_var(paramgrp_id_out,param_saveint_id,save_interval)&
,3,"put_var","save_interval")
      call check(nf90_put_var(paramgrp_id_out,param_potint_id,time_interval)&
,3,"put_var","time_interval")
      call check(nf90_put_var(paramgrp_id_out,param_vmax_id,vmax)&
,3,"put_var","vmax")
      call check(nf90_put_var(paramgrp_id_out,param_growth_id,growth)&
,3,"put_var","growth")
      call check(nf90_put_var(paramgrp_id_out,param_vcrt_id,vcrt)&
,3,"put_var","vcrt")
      call check(nf90_put_var(paramgrp_id_out,param_dtau_id,deltatau)&
,3,"put_var","deltatau")
      call check(nf90_put_var(paramgrp_id_out,param_diffc_id,diff_frac)&
,3,"put_var","diff_frac")
      call check(nf90_put_var(paramgrp_id_out,param_theta_id,theta)&
,3,"put_var","theta")
      call check(nf90_put_var(paramgrp_id_out,param_repro_id,repro)&
,3,"put_var","repro")
      call check(nf90_put_var(paramgrp_id_out,param_alpha_id,alpha)&
,3,"put_var","alpha")
      call check(nf90_put_var(paramgrp_id_out,param_beta_id,beta)&
,3,"put_var","beta")
      call check(nf90_put_var(paramgrp_id_out,param_eta_id,eta)&
,3,"put_var","eta")
      call check(nf90_put_var(paramgrp_id_out,param_maxpot_id,maxpot)&
,3,"put_var","maxpot")

      if (switch_slices == 0) then
        call check(nf90_put_var(nc_id_out,water_id_out,watermask,start=(/1,1/),count=(/x,y/))&
,3,"put_var","watermask")
        call check(nf90_put_var(nc_id_out,ehep_id_out,potnum/area,start=(/1,1/),count=(/x,y/))&
,3,"put_var","potnum")
      else if (switch_slices == 1) then
        call check(nf90_put_var(nc_id_out,water_id_out,watermask_sl(:,:,1),start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","watermask_sl")
        call check(nf90_put_var(nc_id_out,ehep_id_out,potnum_sl(:,:,1)/area,start=(/1,1,1/),count=(/x,y,1/))&
,3,"put_var","potnum_sl")
      end if
    print *,"Inputfile: ",path_EHEP
    print *,"Outputfile: ",path_output
    end subroutine mmf_setup_load
    
    double precision function minmod(a,b,c)
    ! MinMod function that is required for the limiter to sustain stability
      double precision,intent(in) :: a,b,c
      
      if ((a > 0.d0) .AND. (b > 0.d0) .AND. (c > 0.d0)) then
        minmod = min(a,b,c)
      elseif ((a < 0.d0) .AND. (b < 0.d0) .AND. (c < 0.d0)) then
        minmod = max(a,b,c)
      else
        minmod = 0.d0
      end if
    end function minmod

    double precision function solve_kolmog(dens,vx,vy,dif,watermask,j,k,dx,dy,dx_grid,&
                                           &dy_grid,adv_t,diff_t,bir_t,dea_t,pot,an,num,&
                                           &Fxp,Fxm,Fyp,Fym,Pxp,Pxm,Pyp,Pym)
    ! Main solver for the Advection-Diffusion equation with sources and sinks covered seperatly in the run script
    ! Solver is based on the Central Schemes for nonlinear Convection-Diffusion Equation by Kurganov and Tadmor (2000)
      ! Variable declarations
      integer,intent(in) :: j,k,num
      double precision, dimension(:,:),intent(in) :: dens, dif, watermask, pot, an
      double precision, dimension(:,:,:),intent(inout) :: adv_t, diff_t, bir_t, dea_t, Fxp, Fxm, Fyp, Fym, Pxp, Pxm, Pyp, Pym
      double precision, dimension(:,:),intent(in) :: vx, vy, dx_grid, dy_grid
      double precision, intent(in) :: dx, dy
      double precision :: ux_r_p, ux_l_p, ux_r_m, ux_l_m, uy_r_p, uy_l_p, uy_r_m, uy_l_m
      double precision :: ax_p, ax_m, ay_p, ay_m
      double precision :: Fx_p, Fx_m, Px_p, Px_m, Fy_p, Fy_m, Py_p, Py_m, Fxs, Fys, Pxs, Pys
      double precision :: Advterm, Diffterm, Birthterm, Deathterm, popul
      
      ! Only calculate the new values on a grid point that is not masked as ocean or water
      if (watermask(j,k) == 1.d0) then
        if (switch_equid == 1) then ! EQUIDISTANT CALCULATIONS
        if (switch_adv == 1) then
        ! intermediate values between two gridpoints using the minmod limiter
        ux_r_p = dens(j+1,k) - 0.5d0*dx*minmod(theta*(dens(j+1,k) - dens(j,k))/dx,(dens(j+2,k) - dens(j,k))/(2.d0*dx), &
                 theta*(dens(j+2,k) - dens(j+1,k))/dx)
        ux_l_p = dens(j,k) + 0.5d0*dx*minmod(theta*(dens(j,k) - dens(j-1,k))/dx,(dens(j+1,k) - dens(j-1,k))/(2.d0*dx), &
                 theta*(dens(j+1,k) - dens(j,k))/dx)
        ux_r_m = dens(j,k) - 0.5d0*dx*minmod(theta*(dens(j,k) - dens(j-1,k))/dx,(dens(j+1,k) - dens(j-1,k))/(2.d0*dx), &
                 theta*(dens(j+1,k) - dens(j,k))/dx)
        ux_l_m = dens(j-1,k) + 0.5d0*dx*minmod(theta*(dens(j-1,k) - dens(j-2,k))/dx,(dens(j,k) - dens(j-2,k))/(2.d0*dx), &
                 theta*(dens(j,k) - dens(j-1,k))/dx)
        uy_r_p = dens(j,k+1) - 0.5d0*dy*minmod(theta*(dens(j,k+1) - dens(j,k))/dy,(dens(j,k+2) - dens(j,k))/(2.d0*dy), &
                 theta*(dens(j,k+2) - dens(j,k+1))/dy)
        uy_l_p = dens(j,k) + 0.5d0*dy*minmod(theta*(dens(j,k) - dens(j,k-1))/dy,(dens(j,k+1) - dens(j,k-1))/(2.d0*dy), &
                 theta*(dens(j,k+1) - dens(j,k))/dy)
        uy_r_m = dens(j,k) - 0.5d0*dy*minmod(theta*(dens(j,k) - dens(j,k-1))/dy,(dens(j,k+1) - dens(j,k-1))/(2.d0*dy), &
                 theta*(dens(j,k+1) - dens(j,k))/dy)
        uy_l_m = dens(j,k-1) + 0.5d0*dy*minmod(theta*(dens(j,k-1) - dens(j,k-2))/dy,(dens(j,k) - dens(j,k-2))/(2.d0*dy), &
                 theta*(dens(j,k) - dens(j,k-1))/dy)

        ! local speeds
        ax_p = abs(vx(j,k))
        ax_m = abs(vx(j-1,k))
        ay_p = abs(vy(j,k))
        ay_m = abs(vy(j,k-1))

        ! numercial convection fluxes on the borders between two gridpoints (advection)
        Fx_p = 0.5d0 * ((vx(j,k)*ux_r_p + vx(j,k)*ux_l_p) - ax_p * (ux_r_p - ux_l_p))
        Fx_m = 0.5d0 * ((vx(j-1,k)*ux_r_m + vx(j-1,k)*ux_l_m) - ax_m * (ux_r_m - ux_l_m))
        
        Fy_p = 0.5d0 * ((vy(j,k)*uy_r_p + vy(j,k)*uy_l_p) - ay_p * (uy_r_p - uy_l_p))
        Fy_m = 0.5d0 * ((vy(j,k-1)*uy_r_m + vy(j,k-1)*uy_l_m) - ay_m * (uy_r_m - uy_l_m))


        ! Combining the calculations above to the advection term
        Fxs = - 1.d0/dx * (Fx_p - Fx_m)
        Fys = - 1.d0/dy * (Fy_p - Fy_m)
        Advterm = Fxs + Fys
        else
          Fx_p = 0.d0
          Fx_m = 0.d0
          Fy_p = 0.d0
          Fy_m = 0.d0
          Advterm = 0.0d0
        end if
        if (switch_diff == 1) then
        ! numercial diffusion fluxes on the borders between two gridpoints (diffusion)
        Px_p = watermask(j+1,k) * 0.5d0 * (dens(j+1,k) - dens(j,k))/dx * (abs(dif(j,k)*sin(an(j,k))) &
               + abs(dif(j+1,k)*sin(an(j+1,k))))
        Px_m = watermask(j-1,k) * 0.5d0 * (dens(j,k) - dens(j-1,k))/dx * (abs(dif(j-1,k)*sin(an(j-1,k)))  &
               + abs(dif(j,k)*sin(an(j,k))))
        Py_p = watermask(j,k+1) * 0.5d0 * (dens(j,k+1) - dens(j,k))/dy * (abs(dif(j,k)*cos(an(j,k)))  &
               + abs(dif(j,k+1)*cos(an(j,k+1))))
        Py_m = watermask(j,k-1) * 0.5d0 * (dens(j,k) - dens(j,k-1))/dy * (abs(dif(j,k-1)*cos(an(j,k-1)))  &
               + abs(dif(j,k)*cos(an(j,k))))
        ! Combining the calculations above to the diffusion term
        Pxs = 1.d0/dx * (Px_p - Px_m)
        Pys = 1.d0/dy * (Py_p - Py_m)
        Diffterm = Pxs + Pys
        else
          Diffterm = 0.d0
          Px_p = 0.d0
          Px_m = 0.d0
          Py_p = 0.d0
          Py_m = 0.d0
        end if
        else ! NON-EQUIDISTANT CALCULATIONS
        if (switch_adv == 1) then
        ! intermediate values between two gridpoints using the minmod limiter
        ux_r_p = dens(j+1,k) - 0.5d0*dx_grid(j,k)*minmod(theta*(dens(j+1,k) - dens(j,k))/dx_grid(j,k), &
                 (dens(j+2,k) - dens(j,k))/(dx_grid(j,k)+dx_grid(j+1,k)), &
                 theta*(dens(j+2,k) - dens(j+1,k))/dx_grid(j+1,k))
        ux_l_p = dens(j,k) + 0.5d0*dx_grid(j,k)*minmod(theta*(dens(j,k) - dens(j-1,k))/dx_grid(j-1,k), &
                 (dens(j+1,k) - dens(j-1,k))/(dx_grid(j,k)+dx_grid(j-1,k)), &
                 theta*(dens(j+1,k) - dens(j,k))/dx_grid(j,k))
        ux_r_m = dens(j,k) - 0.5d0*dx_grid(j-1,k)*minmod(theta*(dens(j,k) - dens(j-1,k))/dx_grid(j-1,k), &
                 (dens(j+1,k) - dens(j-1,k))/(dx_grid(j,k)+dx_grid(j-1,k)), &
                 theta*(dens(j+1,k) - dens(j,k))/dx_grid(j,k))
        ux_l_m = dens(j-1,k) + 0.5d0*dx_grid(j-1,k)*minmod(theta*(dens(j-1,k) - dens(j-2,k))/dx_grid(j-2,k), &
                 (dens(j,k) - dens(j-2,k))/(dx_grid(j-1,k)+dx_grid(j-2,k)), &
                 theta*(dens(j,k) - dens(j-1,k))/dx_grid(j-1,k))
        uy_r_p = dens(j,k+1) - 0.5d0*dy_grid(j,k)*minmod(theta*(dens(j,k+1) - dens(j,k))/dy_grid(j,k), &
                 (dens(j,k+2) - dens(j,k))/(dy_grid(j,k)+dy_grid(j,k+1)), &
                 theta*(dens(j,k+2) - dens(j,k+1))/dy_grid(k,k+1))
        uy_l_p = dens(j,k) + 0.5d0*dy_grid(j,k)*minmod(theta*(dens(j,k) - dens(j,k-1))/dy_grid(j,k-1), &
                 (dens(j,k+1) - dens(j,k-1))/(dy_grid(j,k)+dy_grid(j,k+1)), &
                 theta*(dens(j,k+1) - dens(j,k))/dy_grid(j,k))
        uy_r_m = dens(j,k) - 0.5d0*dy_grid(j,k-1)*minmod(theta*(dens(j,k) - dens(j,k-1))/dy_grid(j,k-1), &
                 (dens(j,k+1) - dens(j,k-1))/(dy_grid(j,k)+dy_grid(j,k-1)), &
                 theta*(dens(j,k+1) - dens(j,k))/dy_grid(j,k))
        uy_l_m = dens(j,k-1) + 0.5d0*dy_grid(j,k-1)*minmod(theta*(dens(j,k-1) - dens(j,k-2))/dy_grid(j,k-2), &
                 (dens(j,k) - dens(j,k-2))/(dy_grid(j,k-1)+dy_grid(j,k-2)), &
                 theta*(dens(j,k) - dens(j,k-1))/dy_grid(j,k-1))

        ! local speeds
        ax_p = abs(vx(j,k))
        ax_m = abs(vx(j-1,k))
        ay_p = abs(vy(j,k))
        ay_m = abs(vy(j,k-1))

        ! numercial convection fluxes on the borders between two gridpoints (advection)
        Fx_p = 0.5d0 * ((vx(j,k)*ux_r_p + vx(j,k)*ux_l_p) - ax_p * (ux_r_p - ux_l_p))
        Fx_m = 0.5d0 * ((vx(j-1,k)*ux_r_m + vx(j-1,k)*ux_l_m) - ax_m * (ux_r_m - ux_l_m))
        
        Fy_p = 0.5d0 * ((vy(j,k)*uy_r_p + vy(j,k)*uy_l_p) - ay_p * (uy_r_p - uy_l_p))
        Fy_m = 0.5d0 * ((vy(j,k-1)*uy_r_m + vy(j,k-1)*uy_l_m) - ay_m * (uy_r_m - uy_l_m))
        

        ! Combining the calculations above to the advection term
        Fxs = - 1.d0/(0.5d0*(dx_grid(j-1,k)+dx_grid(j,k))) * (Fx_p - Fx_m)
        Fys = - 1.d0/(0.5d0*(dy_grid(j,k-1)+dy_grid(j,k))) * (Fy_p - Fy_m)
        Advterm = Fxs + Fys
        else
          Advterm = 0.d0
          Fx_p = 0.d0
          Fx_m = 0.d0
          Fy_p = 0.d0
          Fy_m = 0.d0
        end if
        if (switch_diff == 1) then
        ! numercial diffusion fluxes on the borders between two gridpoints (diffusion)
        Px_p = watermask(j+1,k) * 0.5d0 * (dens(j+1,k) - dens(j,k))/dx_grid(j,k) * (abs(dif(j,k)*sin(an(j,k))) &
               + abs(dif(j+1,k)*sin(an(j+1,k))))
        Px_m = watermask(j-1,k) * 0.5d0 * (dens(j,k) - dens(j-1,k))/dx_grid(j-1,k) * (abs(dif(j-1,k)*sin(an(j-1,k))) &
               + abs(dif(j,k)*sin(an(j,k))))
        Py_p = watermask(j,k+1) * 0.5d0 * (dens(j,k+1) - dens(j,k))/dy_grid(j,k) * (abs(dif(j,k)*cos(an(j,k))) &
               + abs(dif(j,k+1)*cos(an(j,k+1))))
        Py_m = watermask(j,k-1) * 0.5d0 * (dens(j,k) - dens(j,k-1))/dy_grid(j,k-1) * (abs(dif(j,k-1)*cos(an(j,k-1))) &
               + abs(dif(j,k)*cos(an(j,k))))

        ! Combining the calculations above to the diffusion term

        Pxs = 1.d0/(0.5d0*(dx_grid(j,k)+dx_grid(j-1,k))) * (Px_p - Px_m)
        Pys = 1.d0/(0.5d0*(dy_grid(j,k)+dy_grid(j,k-1))) * (Py_p - Py_m)
        Diffterm = Pxs + Pys
        else
          Diffterm = 0.d0
          Px_p = 0.d0
          Px_m = 0.d0
          Py_p = 0.d0
          Py_m = 0.d0
        end if
        end if
        
        if ((switch_birth == 1) .AND. (pot(j,k) > reprod(j,k))) then
          popul = growth * (3.d0 * dexp(-dens(j,k)/(2.d0*pot(j,k))) - 2.d0)
          if (popul >= 0) then
            Birthterm = popul * dens(j,k)
            Deathterm = 0.d0
          else
            Deathterm = popul * dens(j,k)
            if (-Deathterm*dt > dens(j,k)) then
              Deathterm = -dens(j,k) / dt
            end if
            Birthterm = 0.d0
          end if
        else
          Birthterm = 0.d0
          if ((switch_birth == 1) .AND. (dens(j,k) > 0.d0)) then
            Deathterm = -0.1d0
            if (-Deathterm*dt > dens(j,k)) then
              Deathterm = -dens(j,k) / dt
            end if
          else
            Deathterm = 0.d0
          end if
        end if
        ! return the local change at gridpoint j,k of the right hand side
        solve_kolmog = Advterm + Diffterm + Birthterm + Deathterm
        adv_t(j,k,num) = Advterm
        diff_t(j,k,num) = Diffterm
        bir_t(j,k,num) = Birthterm
        dea_t(j,k,num) = Deathterm
        Fxp(j,k,num) = Fx_p
        Fxm(j,k,num) = Fx_m
        Fyp(j,k,num) = Fy_p
        Fym(j,k,num) = Fy_m
        Pxp(j,k,num) = -Px_p
        Pxm(j,k,num) = -Px_m
        Pyp(j,k,num) = -Py_p
        Pym(j,k,num) = -Py_m
      else
        solve_kolmog = 0.d0 ! No change on gridpoints that are masked !
        adv_t(j,k,num) = 0.d0
        diff_t(j,k,num) = 0.d0
        bir_t(j,k,num) = 0.d0
        dea_t(j,k,num) = 0.d0
        Fxp(j,k,num) = 0.d0
        Fxm(j,k,num) = 0.d0
        Fyp(j,k,num) = 0.d0
        Fym(j,k,num) = 0.d0
        Pxp(j,k,num) = 0.d0
        Pxm(j,k,num) = 0.d0
        Pyp(j,k,num) = 0.d0
        Pym(j,k,num) = 0.d0
      end if
    end function solve_kolmog
    
end module mmf_setup
