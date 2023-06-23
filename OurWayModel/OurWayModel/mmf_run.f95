program mmf_run
  use mmf_config
  use mmf_setup
  use netcdf
  use omp_lib
  implicit none
  ! Declaration of variables

  ! Setup all required variabels and predefinitions
  print *,"Preparing the calculations..."
  call mmf_setup_load()
  print *,"Preperations complete, starting calculations..."
  ! Start the calculations with the time loop over all timesteps
  time_do: do t = 1, timesteps
    !$omp parallel
    ! Overwrite the old with the new density
    !$omp workshare
    !density = density_s
    numdens = density_s
    !$omp end workshare
    !Calculate the gradient fields of the available potential
    !$omp do
    y_do_gra: do k = 1, y-1
      x_do_gra: do j = 1, x-1
        if (switch_equid == 1) then
          gradxfield(j,k) = (avail_pot(j+1,k) - avail_pot(j,k))/dx
          gradyfield(j,k) = (avail_pot(j,k+1) - avail_pot(j,k))/dy
        else
          gradxfield(j,k) = (avail_pot(j+1,k) - avail_pot(j,k))/dx_grid(j,k)
          gradyfield(j,k) = (avail_pot(j,k+1) - avail_pot(j,k))/dy_grid(j,k)
        end if
       
        if (abs(gradxfield(j,k)) <= Eps .AND. abs(gradyfield(j,k)) > Eps) then
          if (gradyfield(j,k) > 0) then
            xi(j,k) = (PI/2.d0)
          else
            xi(j,k) = -(PI/2.d0)
          end if
        else if (abs(gradyfield(j,k)) <= Eps .AND. abs(gradxfield(j,k)) > Eps) then
          xi(j,k) = 0.d0
        else if (abs(gradyfield(j,k)) <= Eps .AND. abs(gradxfield(j,k)) <= Eps) then
          xi(j,k) = (PI/4.d0)
        else
          xi(j,k) = atan(gradyfield(j,k)/gradxfield(j,k))
        end if
      end do x_do_gra
    end do y_do_gra
    !$omp end do
    !maxpot_grad = max(maxval(gradxfield),maxval(gradyfield))
    gamma_para = vmax / (maxpot * deltatau)
    vxfield_old = vxfield
    vyfield_old = vyfield
    !Calculate the advection fields
    !$omp do
    y_do_adv: do k = 1, y-1
      x_do_adv: do j = 1, x-1
        if (((k <= 2) .OR. (j<=2)) .OR. ((k >= y-2) .OR. (j >= x-2))) then !Boundary conditions
          vxfield(j,k) = 0.d0
          vyfield(j,k) = 0.d0
        else
          vxfield(j,k) = dt*gamma_para*gradxfield(j,k) + vxfield_old(j,k) - dt*eta*vxfield_old(j,k)
          vyfield(j,k) = dt*gamma_para*gradyfield(j,k) + vyfield_old(j,k) - dt*eta*vyfield_old(j,k)
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
        !end if
      end do x_do_adv
    end do y_do_adv
    !$omp end do

    ! Calculate the diffusion fields
    !$omp do
    y_do_dif: do k = 1, y
      x_do_dif: do j = 1, x
        if (((k <= 2) .OR. (j <= 2)) .OR. ((k >= y-1) .OR. (j >= x-1))) then ! boundary
          diffusion(j,k) = 0.d0
        else
        !local_speed = (abs(vxfield(j,k) + vxfield(j-1,k))/ 2.d0 + abs(vyfield(j,k) + vyfield(j,k-1))/2.d0)/2.d0
        local_speed = sqrt(vxfield(j,k)**2 + vyfield(j,k)**2)/2.d0
        if (local_speed > vcrt) then
          diffusion(j,k) = diff_frac * local_grid(j,k) * local_speed
        else
          diffusion(j,k) = diff_frac * local_grid(j,k) * vcrt
        end if
        end if
      end do x_do_dif
    end do y_do_dif
    !$omp end do
    diffusion = diffusion * watermask

    ! Runge Kutta 4th Order time stepping and calculations of the RHS
    ! First Step
    !$omp do
    y_do_k1: do k = 1, y
      x_do_k1: do j = 1, x
        if (((k <= 2) .OR. (j <= 2)) .OR. ((k >= y-1) .OR. (j >= x-1))) then
          k1(j,k) = 0.d0
        else
          k1(j,k) = solve_kolmog(density_s,vxfield,vyfield,diffusion,watermask, &
          & j,k,dx,dy,dx_grid,dy_grid,advec_temp,diffu_temp,birth_temp,death_temp,potnum,xi,1, &
          & Fxp_tmp,Fxm_tmp,Fyp_tmp,Fym_tmp,Pxp_tmp,Pxm_tmp,Pyp_tmp,Pym_tmp)
        end if
      end do x_do_k1
    end do y_do_k1
    !$omp end do
    !$omp workshare
    !density_s = density + dt * k1 * 0.5d0
    density_s = numdens + dt * k1 * 0.5d0
    !$omp end workshare
    ! Second Step
    !$omp do
    y_do_k2: do k = 1, y
      x_do_k2: do j = 1, x
        if (((k <= 2) .OR. (j <= 2)) .OR. ((k >= y-1) .OR. (j >= x-1))) then
          k2(j,k) = 0.d0
        else
          k2(j,k) = solve_kolmog(density_s,vxfield,vyfield,diffusion,watermask, & 
          & j,k,dx,dy,dx_grid,dy_grid,advec_temp,diffu_temp,birth_temp,death_temp,potnum,xi,2, &
          & Fxp_tmp,Fxm_tmp,Fyp_tmp,Fym_tmp,Pxp_tmp,Pxm_tmp,Pyp_tmp,Pym_tmp)
        end if
      end do x_do_k2
    end do y_do_k2
    !$omp end do
    !$omp workshare
    !density_s = density + dt * k2 * 0.5d0
    density_s = numdens + dt * k2 * 0.5d0
    !$omp end workshare
    ! Third Step
    !$omp do
    y_do_k3: do k = 1, y
      x_do_k3: do j = 1, x
        if (((k <= 2) .OR. (j <= 2)) .OR. ((k >= y-1) .OR. (j >= x-1))) then
          k3(j,k) = 0.d0
        else
          k3(j,k) = solve_kolmog(density_s,vxfield,vyfield,diffusion,watermask, &
          & j,k,dx,dy,dx_grid,dy_grid,advec_temp,diffu_temp,birth_temp,death_temp,potnum,xi,3, &
          & Fxp_tmp,Fxm_tmp,Fyp_tmp,Fym_tmp,Pxp_tmp,Pxm_tmp,Pyp_tmp,Pym_tmp)
        end if
      end do x_do_k3
    end do y_do_k3
    !$omp end do
    !$omp workshare
    !density_s = density + dt * k3
    density_s = numdens + dt * k3
    !$omp end workshare
    ! Fourth Step
    !$omp do
    y_do_k4: do k = 1, y
      x_do_k4: do j = 1, x
        if (((k <= 2) .OR. (j <= 2)) .OR. ((k >= y-1) .OR. (j >= x-1))) then
          k4(j,k) = 0.d0
        else
          k4(j,k) = solve_kolmog(density_s,vxfield,vyfield,diffusion,watermask, &
          & j,k,dx,dy,dx_grid,dy_grid,advec_temp,diffu_temp,birth_temp,death_temp,potnum,xi,4, &
          & Fxp_tmp,Fxm_tmp,Fyp_tmp,Fym_tmp,Pxp_tmp,Pxm_tmp,Pyp_tmp,Pym_tmp)
        end if
      end do x_do_k4
    end do y_do_k4
    !$omp end do
    ! Combine the steps and calculate the new density
    !$omp workshare
    !density_s = density + dt*(k1 + 2.d0*k2 + 2.d0*k3 + k4)/6.d0
    density_s = numdens + dt*(k1 + 2.d0*k2 + 2.d0*k3 + k4)/6.d0
    advec_t = advec_t + dt*(advec_temp(:,:,1) + advec_temp(:,:,2)*2.d0 + advec_temp(:,:,3)*2.d0 + advec_temp(:,:,4))/6.d0
    diffu_t = diffu_t + dt*(diffu_temp(:,:,1) + diffu_temp(:,:,2)*2.d0 + diffu_temp(:,:,3)*2.d0 + diffu_temp(:,:,4))/6.d0
    birth_t = birth_t + dt*(birth_temp(:,:,1) + birth_temp(:,:,2)*2.d0 + birth_temp(:,:,3)*2.d0 + birth_temp(:,:,4))/6.d0
    death_t = death_t + dt*(death_temp(:,:,1) + death_temp(:,:,2)*2.d0 + death_temp(:,:,3)*2.d0 + death_temp(:,:,4))/6.d0
    Fxp_s = Fxp_s + dt*(Fxp_tmp(:,:,1) + Fxp_tmp(:,:,2)*2.d0 + Fxp_tmp(:,:,3)*2.d0 + Fxp_tmp(:,:,4))/6.d0
    Fxm_s = Fxm_s + dt*(Fxm_tmp(:,:,1) + Fxm_tmp(:,:,2)*2.d0 + Fxm_tmp(:,:,3)*2.d0 + Fxm_tmp(:,:,4))/6.d0
    Fyp_s = Fyp_s + dt*(Fyp_tmp(:,:,1) + Fyp_tmp(:,:,2)*2.d0 + Fyp_tmp(:,:,3)*2.d0 + Fyp_tmp(:,:,4))/6.d0
    Fym_s = Fym_s + dt*(Fym_tmp(:,:,1) + Fym_tmp(:,:,2)*2.d0 + Fym_tmp(:,:,3)*2.d0 + Fym_tmp(:,:,4))/6.d0
    Pxp_s = Pxp_s + dt*(Pxp_tmp(:,:,1) + Pxp_tmp(:,:,2)*2.d0 + Pxp_tmp(:,:,3)*2.d0 + Pxp_tmp(:,:,4))/6.d0
    Pxm_s = Pxm_s + dt*(Pxm_tmp(:,:,1) + Pxm_tmp(:,:,2)*2.d0 + Pxm_tmp(:,:,3)*2.d0 + Pxm_tmp(:,:,4))/6.d0
    Pyp_s = Pyp_s + dt*(Pyp_tmp(:,:,1) + Pyp_tmp(:,:,2)*2.d0 + Pyp_tmp(:,:,3)*2.d0 + Pyp_tmp(:,:,4))/6.d0
    Pym_s = Pym_s + dt*(Pym_tmp(:,:,1) + Pym_tmp(:,:,2)*2.d0 + Pym_tmp(:,:,3)*2.d0 + Pym_tmp(:,:,4))/6.d0
    
    pop_los = pop_los + sum(density_s, MASK=density_s <= 0.d0)
    where (density_s <= 0.d0)
      density_s = 0.d0
    end where

    ! Calculate loss on borders
    
    borderloss = borderloss + sum(density_s(3,3:y-2)) + sum(density_s(x-2,3:y-2)) &
    &+ sum(density_s(3:x-2,3)) + sum(density_s(3:x-2,y-2))
    
    !$omp end workshare
    !print *,borderloss
    density_s(1:3,:) = 0.d0
    density_s(:,1:3) = 0.d0
    density_s(x-2:x,:) = 0.d0
    density_s(:,y-2:y) = 0.d0


    ! Check if potential slice changes with the next step and adjust potential and diffusion accordingly
    if (switch_slices == 1) then
      if (mod(t,time_interval*fraction_time) == 0) then
        !$omp single
        currentslice = currentslice + 1
        !$omp end single
        if (currentslice < (years/time_interval)) then
          !$omp workshare
          !potential = potential_sl(:,:,currentslice)
          potnum = potnum_sl(:,:,currentslice)
          watermask = watermask_sl(:,:,currentslice)
          diffusion = watermask * diffusion
          slope(:,:) = (potnum_sl(:,:,currentslice+1) - potnum_sl(:,:,currentslice)) / (time_interval*fraction_time)
          ! ax - a(x-1) = ax - ax + a = a
          intercept(:,:) = potnum_sl(:,:,currentslice) - slope*time_interval*fraction_time*(currentslice-1)
          !$omp end workshare
        end if
      end if
    end if

    ! Calculate the available potential with the new density
    
    !avail_pot = potential - density_s
    if (switch_slices == 0) then
      !$omp workshare
      !avail_pot = potnum - density_s
      Theta_tmp = Theta_func(density_s,potnum,minnum)
      avail_pot = weibull_func(Theta_tmp,wei_scale)
      !$omp end workshare
    else
      !$omp workshare
      !avail_pot = slope * t + intercept - density_s
      Theta_tmp = Theta_func(density_s,slope*t+intercept,minnum)
      avail_pot = weibull_func(Theta_tmp,wei_scale)
      !$omp end workshare
    end if
    

    !$omp end parallel
    ! Save results every specified save interval
    if (mod(t,save_interval*fraction_time) == 0) then
      currentyear = currentyear + 1
      
      ! Save intermediate results
      call check(nf90_put_var(nc_id_out,numdens_id_out,density_s,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","density_s",t,currentyear)
      density = density_s / (area)
      call check(nf90_put_var(nc_id_out,density_id,density,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","density",t,currentyear)
      call check(nf90_put_var(nc_id_out,availpot_id,avail_pot/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","avail_pot",t,currentyear)
      call check(nf90_put_var(nc_id_out,borderloss_id,borderloss,start=(/currentyear+1/))&
,4,"put_var","borderloss",t,currentyear)
      call check(nf90_put_var(nc_id_out,diffpara_id_out,diffusion,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","diffpara",t,currentyear)
      call check(nf90_put_var(nc_id_out,pop_los_id_out,pop_los,start=(/currentyear+1/))&
,4,"put_var","pop_los",t,currentyear)
      call check(nf90_put_var(nc_id_out,velocity_x_id,vxfield,start=(/1,1,currentyear+1/),count=(/x_1,y_1,1/))&
,4,"put_var","vxfield",t,currentyear)
      call check(nf90_put_var(nc_id_out,velocity_y_id,vyfield,start=(/1,1,currentyear+1/),count=(/x_1,y_1,1/))&
,4,"put_var","vyfield",t,currentyear)
      call check(nf90_put_var(nc_id_out,advection_id,advec_t/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","advec_t",t,currentyear)
      call check(nf90_put_var(nc_id_out,diffusion_id,diffu_t/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","diffu_t",t,currentyear)
      call check(nf90_put_var(nc_id_out,birth_id,birth_t/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","birth_t",t,currentyear)
      call check(nf90_put_var(nc_id_out,death_id,death_t/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","death_t",t,currentyear)
      call check(nf90_put_var(nc_id_out,xi_id_out,xi,start=(/1,1,currentyear+1/),count=(/x_1,y_1,1/))&
,4,"put_var","xi",t,currentyear)
      call check(nf90_put_var(nc_id_out,gradx_id_out,gradxfield,start=(/1,1,currentyear+1/),count=(/x_1,y_1,1/))&
,4,"put_var","gradxfield",t,currentyear)
      call check(nf90_put_var(nc_id_out,grady_id_out,gradyfield,start=(/1,1,currentyear+1/),count=(/x_1,y_1,1/))&
,4,"put_var","gradyfield",t,currentyear)

      call check(nf90_put_var(nc_id_out,Fxp_id_out,Fxp_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,4,"put_var","Fxp_s")
      call check(nf90_put_var(nc_id_out,Fyp_id_out,Fyp_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,3,"put_var","Fyp_s")
      call check(nf90_put_var(nc_id_out,Pxp_id_out,Pxp_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,3,"put_var","Pxp_s")
      call check(nf90_put_var(nc_id_out,Pyp_id_out,Pyp_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,3,"put_var","Pyp_s")
      call check(nf90_put_var(nc_id_out,Fxm_id_out,Fxm_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,3,"put_var","Fxm_s")
      call check(nf90_put_var(nc_id_out,Fym_id_out,Fym_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,3,"put_var","Fym_s")
      call check(nf90_put_var(nc_id_out,Pxm_id_out,Pxm_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,3,"put_var","Pxm_s")
      call check(nf90_put_var(nc_id_out,Pym_id_out,Pym_s/area,start=(/1,1,currentyear+1/),count=(/x,y,1/))&
,3,"put_var","Pym_s")

      if (switch_slices == 1) then
        if (currentslice < 2388) then
          call check(nf90_put_var(nc_id_out,water_id_out,watermask_sl(:,:,currentslice)&
          &,start=(/1,1,currentyear+1/),count=(/x,y,1/)),4,"put_var","watermask_sl",t,currentyear)
          call check(nf90_put_var(nc_id_out,ehep_id_out,potnum_sl(:,:,currentslice)/area&
          &,start=(/1,1,currentyear+1/),count=(/x,y,1/)),4,"put_var","potnum_sl",t,currentyear)
        end if
      end if
      
      ! Print progress to standard output
      print *,"Finished timestep ",t,", saved ",currentyear," times."
      print *,"Current borderloss: ",borderloss
      print *,"Current numerical deathterm loss: ",pop_los
      print *,"Current human total: ",sum(density_s)

      ! Reset the sums that were used before
      borderloss = 0.d0
      pop_los = 0.d0
      advec_t = 0.d0
      diffu_t = 0.d0
      birth_t = 0.d0
      death_t = 0.d0
      Fxp_s(:,:) = 0.d0
      Fxm_s(:,:) = 0.d0
      Fyp_s(:,:) = 0.d0
      Fym_s(:,:) = 0.d0
      Pxp_s(:,:) = 0.d0
      Pxm_s(:,:) = 0.d0
      Pyp_s(:,:) = 0.d0
      Pym_s(:,:) = 0.d0
    end if
    
    
  end do time_do 

! Close the output file and finish the program
call check(nf90_close(nc_id_out),4,"close")
print *,"Finished calculations!"
end program mmf_run

