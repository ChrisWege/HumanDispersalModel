#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 12:58:23 2020

@author: christian

Postprocessing script, calculating different additional values and appends
them for further investigation onto the given .nc file.
Adds:
Dispersal Strength (velocity*population)
Dispersal Strength Density (velocity*density)
Ratio Advection/diffusion (adv/diff term)
Sum of Population, average of non-zero density
Sum of Carrying Capacity, average of non-zero capacity density
Sum growth/death term, -> net gain
Also add:
density EHEP (ehep/area)
"""

from netCDF4 import Dataset
import numpy as np
import sys
from geographiclib.geodesic import Geodesic
from multiprocessing import Pool, Array, RawArray, Lock
import time
from datetime import timedelta
import itertools

var_dict = {}

def init_equi_div(a,b,c,d,e,f,g,h,i):
    var_dict['velx'] = a
    var_dict['vely'] = b
    var_dict['distance_halfx'] = c
    var_dict['distance_halfy'] = d
    var_dict['dlat'] = e
    var_dict['dlon'] = f
    var_dict['equidistant'] = g
    var_dict['distance'] = h
    var_dict['dtime'] = i
    
def equi_div(t,i,j):
    velx = np.frombuffer(var_dict['velx']).reshape((var_dict['dtime'],var_dict['dlat']-1,var_dict['dlon']-1))
    vely = np.frombuffer(var_dict['vely']).reshape((var_dict['dtime'],var_dict['dlat']-1,var_dict['dlon']-1))
    distance_halfx = np.frombuffer(var_dict['distance_halfx']).reshape((var_dict['dlat']-2,var_dict['dlon']-2))
    distance_halfy = np.frombuffer(var_dict['distance_halfy']).reshape((var_dict['dlat']-2,var_dict['dlon']-2))
    
    if not var_dict['equidistant']:
        if (((i == 0) or (j == 0)) or ((i == var_dict['dlat']-1) or (j == var_dict['dlon']-1))):
            diverge_vel = 0.
        else:
            diverge_vel = (velx[t,i,j-1] - velx[t,i,j]) / distance_halfx[i-1,j-1] + (vely[t,i-1,j] - vely[t,i,j]) / distance_halfy[i-1,j-1]
    else:
        if (((i == 0) or (j == 0)) or ((i == var_dict['dlat']-1) or (j == var_dict['dlon']-1))):
            diverge_vel = 0.
        else:
            diverge_vel = (velx[t,i,j-1] - velx[t,i,j]) / var_dict['distance'] + (vely[t,i-1,j] - vely[t,i,j]) / var_dict['distance']
                              
    return diverge_vel

#    
#def init_non_equi_div():
#def non_equi_div(t):
#    
#def init_avg_flux():
#def avg_flux(t):
#    
#def init_d_rho_d_t(t):
#def d_rho_d_t(t):
#    
#def init_add_comp():
#def add_comp(t):

if __name__ == "__main__":
    print("Starting Script, loading data...")
    core_num = 4
    
    #direc = str(sys.argv[1])
    #filen = str(sys.argv[2])
    #print(direc,filen)
    
    #path_in = "/data/sfb806/human_mobility_model/dispersal_model/testcases/ReCalcs/"+direc
    #file_in = filen + ".nc"
    
    #path_in = "/data/sfb806/human_mobility_model/christian/Applic_Auri_Grav/Exp1b/Grav/"
    #file_in = "Dispersal_grav_cwc_2k_open_new.nc"
    
    path_in = "/data/sfb806/Our_Way_Model/HDM/results/new_standard/less_both/"
    file_in = "Dispersal_45k_Rework_run_5hm_less_both.nc"
    
    #path_in = "/data/sfb806/human_mobility_model/dispersal_model/testcases/ReCalcs/Exp9/Results/"
    #file_in = "Equidistant_9a_r.nc"
    
    # change this to True if the grid is equidistant and specify the distance,
    # otherwise the lat/lons are required to calculate these differences
    equidistant = False
    distance = 50. # 50 km
    #if filen == "Equidistant_9i":
    #    distance = 25.
    
    # change this to False if the EHEP is only one time frame long (e.g has no time dimension)
    TimeEHEP = True
    
    data = Dataset(path_in+file_in,mode='r+')
    y = data.dimensions['y']
    dlat = len(y)
    x = data.dimensions['x']
    dlon = len(x)
    lat = data.variables['lat'][:,:]
    lon = data.variables['lon'][:,:]
    time_var = data.variables['time'][:]
    dtime = len(time_var)
    dens = data.variables['Density'][:,:,:]
    borderloss = data.variables['Borderloss'][:]
    velx = data.variables['Velocity_x'][:,:,:]
    vely = data.variables['Velocity_y'][:,:,:]
    adv = data.variables['Advection_term'][:,:,:]
    diff = data.variables['Diffusion_term'][:,:,:]
    birth = data.variables['Birth_term'][:,:,:]
    death = data.variables['Death_term'][:,:,:]
    hnumb = data.variables['hnumb'][:,:,:]
    FluxFxp = data.variables['Adv_Flux_x_pos'][:,:,:]
    FluxFxm = data.variables['Adv_Flux_x_neg'][:,:,:]
    FluxFyp = data.variables['Adv_Flux_y_pos'][:,:,:]
    FluxFym = data.variables['Adv_Flux_y_neg'][:,:,:]
    FluxPxp = data.variables['Diff_Flux_x_pos'][:,:,:]
    FluxPxm = data.variables['Diff_Flux_x_neg'][:,:,:]
    FluxPyp = data.variables['Diff_Flux_y_pos'][:,:,:]
    FluxPym = data.variables['Diff_Flux_y_neg'][:,:,:]
    #FluxPxp = -FluxPxp
    #FluxPxm = -FluxPxm
    #FluxPyp = -FluxPyp
    #FluxPym = -FluxPym
    #mask = data.variables['watermask'][:,:,:]
    if TimeEHEP:
        ehep = data.variables['Ehep'][:,:,:]
    else:
        ehep = data.variables['Ehep'][:,:]
    #ehep = data.variables['Ehep'][:,:]
    area = data.variables['area'][:,:]
    posx_lat = data.variables['posx_lat'][:,:]
    posy_lat = data.variables['posy_lat'][:,:]
    posx_lon = data.variables['posx_lon'][:,:]
    posy_lon = data.variables['posy_lon'][:,:]
    
    print("Initialisation complete!")
    
    # Distance between gridpoints
    if not equidistant:
        print("Calculation distances...")
        distance_fullx = np.zeros([dlat-1,dlon-1])
        distance_fully = np.zeros([dlat-1,dlon-1])
        distance_halfx = np.zeros([dlat-2,dlon-2])
        distance_halfy = np.zeros([dlat-2,dlon-2])
        geod = Geodesic.WGS84
        # calculates the distance between two gridpoints, both for the regular and the 
        # half-step shifted grid
        for i in list(range(dlon-1)):
            for j in list(range(dlat-1)):
                mlon = geod.InverseLine(lat[j,i],lon[j,i],lat[j,i+1],lon[j,i+1])
                distance_fullx[j,i] = mlon.s13 / 1000.
                mlat = geod.InverseLine(lat[j,i],lon[j,i],lat[j+1,i],lon[j+1,i])
                distance_fully[j,i] = mlat.s13 / 1000.
        
        for i in list(range(dlon-2)):
            for j in list(range(dlat-2)):
                mlon = geod.InverseLine(posx_lat[j,i],posx_lon[j,i],posx_lat[j,i+1],posx_lon[j,i+1])
                distance_halfx[j,i] = mlon.s13 / 1000.
                mlat = geod.InverseLine(posy_lat[j,i],posy_lon[j,i],posy_lat[j+1,i],posy_lon[j+1,i])
                distance_halfy[j,i] = mlat.s13 / 1000.
    
    print("Complete!")
    
    
    # Ratio advection/diffusion term
    ratadvdiff = np.zeros([dtime,dlat,dlon])
    ratadvdiff[:,:,:] = -1
    diff_z = np.nonzero(hnumb)
    ratadvdiff[diff_z] = abs(adv[diff_z]/diff[diff_z])
    
    # divergence field
    print("Calculating divergence fields...")
    diverge_vel = np.zeros([dtime,dlat,dlon])
    
    # divergence fields of the velocity. Keep in mind that these
    # are one unit smaller than the parent field! The velocity divergence is
    # surrounded by one row/column of zeros!
    
    start_time = time.monotonic()
    
    raw_velx = RawArray('d',dtime*(dlat-1)*(dlon-1))
    np_velx = np.frombuffer(raw_velx).reshape((dtime,dlat-1,dlon-1))
    np.copyto(np_velx,velx)
    raw_vely = RawArray('d',dtime*(dlat-1)*(dlon-1))
    np_vely = np.frombuffer(raw_vely).reshape((dtime,dlat-1,dlon-1))
    np.copyto(np_vely,vely)
    raw_distance_halfx = RawArray('d',(dlat-2)*(dlon-2))
    np_distance_halfx = np.frombuffer(raw_distance_halfx).reshape((dlat-2,dlon-2))
    np.copyto(np_distance_halfx,distance_halfx)
    raw_distance_halfy = RawArray('d',(dlat-2)*(dlon-2))
    np_distance_halfy = np.frombuffer(raw_distance_halfy).reshape((dlat-2,dlon-2))
    np.copyto(np_distance_halfy,distance_halfy)
    lock=Lock()
    raw_diverge_vel = Array('d',dtime*dlat*dlon,lock=lock)
    np_diverge_vel = np.frombuffer(raw_diverge_vel.get_obj()).reshape((dtime,dlat,dlon))
    np.copyto(np_diverge_vel,diverge_vel)
    
    with Pool(processes=4,initializer=init_equi_div, 
              initargs=(raw_velx,raw_vely,raw_distance_halfx,raw_distance_halfy,
                        dlat,dlon,equidistant,distance,dtime)) as pool:
        result = pool.starmap(equi_div,list(itertools.product(range(dtime),range(dlat),range(dlon))))
    result=np.array(result)
    diverge_vel[:,:,:] = np.reshape(result,(dtime,dlat,dlon))
        
    print("Complete!")
    end_time = time.monotonic()
    print(timedelta(seconds=end_time - start_time))
    
    # average Fluxes and overall Dispersalflux
    print("Calculating average Fluxes...")
    
    
    # averages the Fluxes for each grid point an adds them together for the
    # totalflux, each dimension for itself
    
    avgAdvFlux_x = (FluxFxp + FluxFxm) / 2.
    avgAdvFlux_y = (FluxFyp + FluxFym) / 2.
    avgDiffFlux_x = (FluxPxp + FluxPxm) / 2.
    avgDiffFlux_y = (FluxPyp + FluxPym) / 2.
    avgTotalFlux_x = avgAdvFlux_x + avgDiffFlux_x
    avgTotalFlux_y = avgAdvFlux_y + avgDiffFlux_y
                
    print("Complete!")
    
    
    print("Calcululating drhodt with the individual terms...")
    drhodt = np.zeros([dtime,dlat,dlon])
    
    drhodt[:,:,:] = adv + diff + birth + death
        
    print("Additional computations and saving the results...")
    start_time = time.monotonic()
    # Sums, averages, net gain
    pop_sum = np.zeros([dtime])
    cc_sum = np.zeros([dtime])
    growth_sum = np.zeros([dtime])
    death_sum = np.zeros([dtime])
    net_pop = np.zeros([dtime])
    dens_avg = np.zeros([dtime])
    cc_avg = np.zeros([dtime])
    hnumb_diff = np.zeros([dtime])
    hnumb_diff2 = np.zeros([dtime])
    bd_diff = np.zeros([dtime])
    pop_mask = hnumb > 0.
    death_fix = np.zeros([dtime,dlat,dlon])
    death_fix[:,:,:] = death[:,:,:] * pop_mask[:,:,:]
    drhodtsum = np.zeros([dtime])
    
    
    for t in range(dtime):
        pop_sum[t] = np.sum(hnumb[t,:,:])
        if TimeEHEP:
            cc_sum[t] = np.sum(ehep[t,:,:])#*area[:,:])
        else:
            cc_sum[t] = np.sum(ehep[:,:])#*area[:,:])
        dens_mask = np.ma.masked_equal(dens[t,:,:],0.)
        dens_avg[t] = np.ma.mean(dens_mask)
        if TimeEHEP:
            cc_mask = np.ma.masked_less(ehep[t,:,:],0.0001)
            cc_avg[t] = np.ma.mean(cc_mask)
        growth_sum[t] = np.sum(birth[t,:,:]*area[:,:])
        death_sum[t] = np.sum(death_fix[t,:,:]*area[:,:])
        drhodtsum[t] = np.sum(drhodt[t,:,:]*pop_mask)
        if t == 0:
            hnumb_diff[t] = 0.
            hnumb_diff2[t] = 0.
        else:
            hnumb_diff[t] = pop_sum[t] - pop_sum[t-1]
            hnumb_diff2[t] = hnumb_diff[t] - drhodtsum[t]
        
    net_pop[:] = growth_sum[:] + death_sum[:]
    if not TimeEHEP:
        cc_mask = np.ma.masked_less(ehep[:,:],0.0001)
        cc_avg[:] = np.ma.mean(cc_mask)
    
    end_time = time.monotonic()
    print(timedelta(seconds=end_time - start_time))
    # save/overwrite the new variabels:
    
    try:
        net_drhodt = data.createVariable('d_rho_d_t','f8',('t','y','x'))
    except RuntimeError:
        net_drhodt = data.variables['d_rho_d_t']
    net_drhodt.longname = 'Rate of change per grid cell between each time step'
    net_drhodt.description = 'Calculated by the sum of each dynamic term'
    net_drhodt.unit = 'number of humans/yr/km^2'
    net_drhodt[:,:,:] = drhodt
    
    
    try:
        net_ratadvdiff = data.createVariable('Ratio_adv_Diff','f8',('t','y','x'))
    except RuntimeError:
        net_ratadvdiff = data.variables['Ratio_adv_Diff']
    net_ratadvdiff.longname = 'Ratio of advection/diffusion term'
    net_ratadvdiff.unit = 'humans/km^2 / humans/km^2'
    net_ratadvdiff[:,:,:] = ratadvdiff
    try:
        net_pop_sum = data.createVariable('Pop_Sum','f8',('t'))
    except RuntimeError:
        net_pop_sum = data.variables['Pop_Sum']
    net_pop_sum.longname = 'Sum of Population in the domain'
    net_pop_sum.unit = 'number of humans'
    net_pop_sum[:] = pop_sum
    try:
        net_pop_diff = data.createVariable('Pop_Diff','f8',('t'))
    except RuntimeError:
        net_pop_diff = data.variables['Pop_Diff']
    net_pop_diff.longname = 'Population change for each time step, derived from the total population sum'
    net_pop_diff.unit = 'number of humans'
    net_pop_diff[:] = hnumb_diff
    try:
        net_cc_sum = data.createVariable('CC_Sum','f8',('t'))
    except RuntimeError:
        net_cc_sum = data.variables['CC_Sum']
    net_cc_sum.longname = 'Sum of Carrying Capacity in the domain'
    net_cc_sum.unit = 'humans/km^2'
    net_cc_sum[:] = cc_sum
    try:
        net_growth_sum = data.createVariable('Growth_Sum','f8',('t'))
    except RuntimeError:
        net_growth_sum = data.variables['Growth_Sum']
    net_growth_sum.longname = 'Sum of the Birth term in the domain'
    net_growth_sum.unit = 'number of humans'
    net_growth_sum[:] = growth_sum
    try:
        net_death_sum = data.createVariable('Death_Sum','f8',('t'))
    except RuntimeError:
        net_death_sum = data.variables['Death_Sum']
    net_death_sum.longname = 'Sum of the Death term in the domain'
    net_death_sum.unit = 'number of humans'
    net_death_sum[:] = death_sum
    try:
        net_net_pop = data.createVariable('Net_Pop','f8',('t'))
    except RuntimeError:
        net_net_pop = data.variables['Net_Pop']
    net_net_pop.longname = 'Net growth/decay of Population in the domain, derived from Birth/Death Sum, the latter filtered for populations greater than zero'
    net_net_pop.unit = 'number of humans'
    net_net_pop[:] = net_pop
    try:
        net_dens_avg = data.createVariable('Dens_Avg','f8',('t'))
    except RuntimeError:
        net_dens_avg = data.variables['Dens_Avg']
    net_dens_avg.longname = 'Average of non-zero density in the domain'
    net_dens_avg.unit = 'humans/km^2'
    net_dens_avg[:] = dens_avg
    try: 
        net_cc_avg = data.createVariable('CC_Avg','f8',('t'))
    except RuntimeError:
        net_cc_avg = data.variables['CC_Avg']
    net_cc_avg.longname = 'Average of non-zero carrying capacity density in the domain'
    net_cc_avg.unit = 'humans/km^2'
    net_cc_avg[:] = cc_avg
    try:
        net_overshoot = data.createVariable('Death_Overshoot','f8',('t'))
    except RuntimeError:
        net_overshoot = data.variables['Death_Overshoot']
    net_overshoot.longname = 'Difference between Net_Pop and Pop_Diff, showing the residual overshoot of the filtered death term'
    net_overshoot.unit = 'number of humans'
    net_overshoot[:] = hnumb_diff - net_pop
    
    try:
        net_dispflux_x = data.createVariable('Dispflux_x','f8',('t','y','x'))
    except RuntimeError:
        net_dispflux_x = data.variables['Dispflux_x']
    net_dispflux_x.longname = "Dispersal Flux field, X-direction"
    net_dispflux_x.unit = 'humans/km/10yr'
    net_dispflux_x[:,:,:] = avgTotalFlux_x
    
    try:
        net_dispflux_y = data.createVariable('Dispflux_y','f8',('t','y','x'))
    except RuntimeError:
        net_dispflux_y = data.variables['Dispflux_y']
    net_dispflux_y.longname = "Dispersal Flux field, Y-direction"
    net_dispflux_y.unit = 'humans/km/10yr'
    net_dispflux_y[:,:,:] = avgTotalFlux_y
    
    try: 
        net_diverge_vel = data.createVariable('Diverge_vel','f8',('t','y','x'))
    except RuntimeError:
        net_diverge_vel = data.variables['Diverge_vel']
    net_diverge_vel.longname = "Divergence of the velocity field"
    net_diverge_vel.unit = "1 / yr"
    net_diverge_vel[:,:,:] = diverge_vel
    
    try: 
        net_diverge_vel = data.createVariable('Avg_Adv_Flux_x','f8',('t','y','x'))
    except RuntimeError:
        net_diverge_vel = data.variables['Avg_Adv_Flux_x']
    net_diverge_vel.longname = "Average of the Advection Fluxes of one gridpoint, X-direction"
    net_diverge_vel.unit = 'humans/km/10yr'
    net_diverge_vel[:,:,:] = avgAdvFlux_x
    
    try: 
        net_diverge_vel = data.createVariable('Avg_Adv_Flux_y','f8',('t','y','x'))
    except RuntimeError:
        net_diverge_vel = data.variables['Avg_Adv_Flux_y']
    net_diverge_vel.longname = "Average of the Advection Fluxes of one gridpoint, Y-direction"
    net_diverge_vel.unit = 'humans/km/10yr'
    net_diverge_vel[:,:,:] = avgAdvFlux_y
        
    try: 
        net_diverge_vel = data.createVariable('Avg_Diff_Flux_x','f8',('t','y','x'))
    except RuntimeError:
        net_diverge_vel = data.variables['Avg_Diff_Flux_x']
    net_diverge_vel.longname = "Average of the Diffusion Fluxes of one gridpoint, X-direction"
    net_diverge_vel.unit = 'humans/km/10yr'
    net_diverge_vel[:,:,:] = avgDiffFlux_x
    
    try: 
        net_diverge_vel = data.createVariable('Avg_Diff_Flux_y','f8',('t','y','x'))
    except RuntimeError:
        net_diverge_vel = data.variables['Avg_Diff_Flux_y']
    net_diverge_vel.longname = "Average of the Diffusion Fluxes of one gridpoint, Y-direction"
    net_diverge_vel.unit = 'humans/km/10yr'
    net_diverge_vel[:,:,:] = avgDiffFlux_y
    
    data.close()
    
    print("Script completed!")




