#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:01:54 2020

@author: christian

Script to preprocess the given Input of existance potentials to be useable by the dispersal script
"""

from netCDF4 import Dataset
import numpy as np
import time
import random
from geographiclib.geodesic import Geodesic
from scipy.interpolate import griddata
import itertools

def preprocess(EHEP,lat,lon,x=None,y=None):
    """Preprocesses the results to be useable with the dispersal model."""
    geod = Geodesic.WGS84
    if switch_capa:
        if switch_lineargrowth:
            for t,val in enumerate(run):
                EHEP[t,:,:] = EHEP[t,:,:]*t*growth+EHEP[t,:,:]*startcapa
        else:
            EHEP[:,:,:] = EHEP[:,:,:] * startcapa 
    if gridtype == 'xy':
        density = np.zeros((y,x))
        if starting_dist == 'uniform':
            density[:,:] = uni_dist * water_mask[:,:]
        elif starting_dist == 'gaus':
            st_x0 = 32.
            st_y0 = 13.
            st_sx0 = 1.
            st_sy0 = 1.
            st_ampli = 0.2
            gauß = np.zeros((y,x))
            xarray = np.linspace(0,x,x)
            yarray = np.linspace(0,y,y)
            for i,valx in enumerate(xarray):
                for j,valy in enumerate(yarray):
                    gauß[j,i] = st_ampli* np.exp(-((valx-st_x0)**2/(2.*st_sx0**2)+(valy-st_y0)**2/(2.*st_sy0**2)))
            density[:,:] = gauß * water_mask
        elif starting_dist == 'firstslice':
            density[:,:] = EHEP[0,:,:] * 0.7
    else:
        density = np.zeros((len(lat),len(lon)))
        x = len(lon)
        y = len(lat)
        if starting_dist == 'uniform':
            density[:,:] = uni_dist * water_mask[:,:]
        elif starting_dist == 'gaus':
            st_x0 = 100.
            st_y0 = 8.
            st_sx0 = 2.
            st_sy0 = 2.
            st_ampli = 0.50
            gauß = np.zeros((y,x))
            xarray = np.linspace(0,x,x)
            yarray = np.linspace(0,y,y)
            for i,valx in enumerate(xarray):
                for j,valy in enumerate(yarray):
                    gauß[j,i] = st_ampli* np.exp(-((valx-st_x0)**2/(2.*st_sx0**2)+(valy-st_y0)**2/(2.*st_sy0**2)))
                    
            gauß[gauß < 0.005] = 0.
            density[:,:] = gauß * water_mask
        elif starting_dist == 'firstslice':
            density[:,:] = EHEP[0,:,:] * 0.7
        elif starting_dist == 'random':
            #for i in lat:
            #  for j in lon:
            #    density[i,j] = random.random() * 0.01
            #drx = int(len(lon)/5)
            #dry = int(len(lat)/5)
            #rx = np.zeros(drx)
            #ry = np.zeros(dry)
            rx = lon[::5]
            ry = lat[::5]
            drx = len(rx)
            dry = len(ry)
            rxy = np.array(list(itertools.product(ry,rx)))
            rxyd = np.zeros((dry,drx))
            grid_lat, grid_lon = np.meshgrid(lat,lon)
            for i in range(drx):
              for j in range(dry):
                rxyd[j,i] = random.random() * startcapa/5.
                
            density = griddata(rxy,rxyd.flatten(), (grid_lat,grid_lon),method='linear',rescale=True)
            density = density.T
            density[52:,:] = 0.
            density[:,125:] = 0.
            density = density * water_mask
        
    # create the output file and fill it with structure
    prepro = Dataset(path_output + file_output,mode='w',format="NETCDF4")
    prepro.description = "Accessible EHEP with additional values for "\
"the mobility model."
    prepro.history = "Created: " +time.ctime(time.time())
    
    net_dim_x = prepro.createDimension("x",x)
    net_dim_y = prepro.createDimension("y",y)
    net_dim_x1 = prepro.createDimension("x_1",x-1)
    net_dim_y1 = prepro.createDimension("y_1",y-1)
    net_dim_x2 = prepro.createDimension("x_2",x+1)
    net_dim_y2 = prepro.createDimension("y_2",y+1)
    net_dim_single = prepro.createDimension("single",1)

    net_lat = prepro.createVariable("lat","f8",("y","x"))
    net_lon = prepro.createVariable("lon","f8",("y","x"))
    net_dens = prepro.createVariable("dens","f8",("y","x"))
    net_hnumber = prepro.createVariable("hnumb","f8",("y","x"))
    net_dx = prepro.createVariable("dx","f8",("single"))
    net_dy = prepro.createVariable("dy","f8",("single"))
    net_dx_grid = prepro.createVariable("dx_grid","f8",("y","x_1"))
    net_dy_grid = prepro.createVariable("dy_grid","f8",("y_1","x"))
    net_posx_lat = prepro.createVariable("posx_lat","f8",("y","x_1"))
    net_posx_lon = prepro.createVariable("posx_lon","f8",("y","x_1"))
    net_posy_lat = prepro.createVariable("posy_lat","f8",("y_1","x"))
    net_posy_lon = prepro.createVariable("posy_lon","f8",("y_1","x"))
    net_area = prepro.createVariable("area","f8",("y","x"))
    net_corner_lat = prepro.createVariable("corner_lat","f8",("y_2","x_2"))
    net_corner_lon = prepro.createVariable("corner_lon","f8",("y_2","x_2"))
    
    
    # Adding additional documentation    
    net_lat.description = "Latitude"
    net_lat.units = "Degrees north"
    net_lon.description = "Longitude"
    net_lon.units = "Degrees east"
    net_dx.description = "Mean distance between two gridpoints, x-direction"
    net_dx.unit = "km"
    net_dy.description = "Mean distance between two gridpoints, y-direction"
    net_dy.unit = "km"
    net_dx_grid.description = "Grid specific distance between two points, x-direction"
    net_dx_grid.uni = "km"
    net_dy_grid.description = "Grid specific distance between two points, y-direction"
    net_dy_grid.uni = "km"
    net_posx_lat.description = "Latitude of border between two points, x-direction"
    net_posx_lon.description = "Longitude of border between two points, x-direction"
    net_posy_lat.description = "Latitude of border between two points, y-direction"
    net_posy_lon.description = "Longitude of border between two points, y-direction"
    net_posx_lat.unit = "Degrees north"
    net_posx_lon.unit = "Degrees east"
    net_posy_lat.unit = "Degrees north"
    net_posy_lon.unit = "Degrees east"
    net_dens.description = "Starting distribution of human density"
    net_dens.unit = "Humans/100km^2"
    net_hnumber.description = "Starting distribution of number of humans"
    net_hnumber.unit = "Humans"
    net_area.description = "Area of the grid cell"
    net_area.unit = "km^2"
    
    # save the first values and initialize others
    if gridtype == 'xy':
        net_lat[:,:] = lat
        net_lon[:,:] = lon
    else:
        latm,lonm = np.meshgrid(lat,lon)
        net_lat[:,:] = latm.T    #might need a .T transpose in some cases
        net_lon[:,:] = lonm.T    #might need a .T transpose in some cases
    net_dens[:,:] = density
    latxygrid = np.zeros((x,y))
    lonxygrid = np.zeros((x,y))
    #meterslat = np.zeros((x,y))
    #meterslon = np.zeros((x,y))
    
    # define lat/lon coordinates on curvilinear grid (and transpose them)
    if gridtype == 'xy':
        for i in list(range(x)):
            for j in list(range(y)):
                latxygrid[i,j] = lat[j,i]
                lonxygrid[i,j] = lon[j,i]
    else:
        latxygrid = latm
        lonxygrid = lonm
        
    # transform the distances from lowest left grid point into distances between
    # each gridpoint
    metersx = np.zeros((x-1,y))
    metersy = np.zeros((x,y-1))
    posx_lat = np.zeros((x-1,y))
    posx_lon = np.zeros((x-1,y))
    posy_lat = np.zeros((x,y-1))
    posy_lon = np.zeros((x,y-1))
    for i in list(range(x-1)):
        for j in list(range(y-1)):
            #glon = geod.Inverse(latxygrid[i,j],lonxygrid[i,j],latxygrid[i+1,j],lonxygrid[i+1,j])
            mlon = geod.InverseLine(latxygrid[i,j],lonxygrid[i,j],latxygrid[i+1,j],lonxygrid[i+1,j])
            metersx[i,j] = mlon.s13
            mlon_t = mlon.Position(0.5 * mlon.s13)
            posx_lat[i,j] = mlon_t['lat2']
            posx_lon[i,j] = mlon_t['lon2']
            #glat = geod.Inverse(latxygrid[i,j],lonxygrid[i,j],latxygrid[i,j+1],lonxygrid[i,j+1])
            mlat = geod.InverseLine(latxygrid[i,j],lonxygrid[i,j],latxygrid[i,j+1],lonxygrid[i,j+1]) 
            metersy[i,j] = mlat.s13
            mlat_t = mlat.Position(0.5 * mlat.s13)
            posy_lat[i,j] = mlat_t['lat2']
            posy_lon[i,j] = mlat_t['lon2']
            
    for i in list(range(x-1)):
        mlon = geod.InverseLine(latxygrid[i,-1],lonxygrid[i,-1],latxygrid[i+1,-1],lonxygrid[i+1,-1])
        metersx[i,-1] = mlon.s13
        mlon_t = mlon.Position(0.5 * mlon.s13)
        posx_lat[i,-1] = mlon_t['lat2']
        posx_lon[i,-1] = mlon_t['lon2']
    for j in list(range(y-1)):
        mlat = geod.InverseLine(latxygrid[-1,j],lonxygrid[-1,j],latxygrid[-1,j+1],lonxygrid[-1,j+1])
        metersy[-1,j] = mlat.s13
        mlat_t = mlat.Position(0.5 * mlat.s13)
        posy_lat[-1,j] = mlat_t['lat2']
        posy_lon[-1,j] = mlat_t['lon2']
        
    # take the mean of the distances for dx and dy, introduces small error but 
    # neglectiable for most cases. No error if given grid is equidistant.
    dx = metersx.mean()/1000.
    dy = metersy.mean()/1000.
    
    # save the rest of the initial values
    net_dx[:] = dx
    #net_dx[:] = 50.
    net_dy[:] = dy
    #net_dy[:] = 50.
    net_dx_grid[:,:] = metersx.T/1000.
    net_dy_grid[:,:] = metersy.T/1000.
    net_posx_lat[:,:] = posx_lat.T
    net_posx_lon[:,:] = posx_lon.T
    net_posy_lat[:,:] = posy_lat.T
    net_posy_lon[:,:] = posy_lon.T
    
    # calculate the area that every gridpoint represents
    ## if the grid is equidistant, the corner points are exact. Otherwise, an
    ## approximation is used by estimating the corner through lat/lon direction
    ## independently and averaging over both estimations.
    
    corner_points_lat = np.zeros((x+1,y+1))
    corner_points_lon = np.zeros((x+1,y+1))
    
    ### find the positions of the corners of each gridbox:
    # inner points
    for i in list(range(x-1)):
        for j in list(range(y-1)):
            mlat = geod.InverseLine(posy_lat[i,j],posy_lon[i,j],posy_lat[i+1,j],posy_lon[i+1,j])
            mlat_t = mlat.Position(0.5* mlat.s13)
            
            mlon = geod.InverseLine(posx_lat[i,j],posx_lon[i,j],posx_lat[i,j+1],posx_lon[i,j+1])
            mlon_t = mlon.Position(0.5* mlon.s13)
            #print(mlat_t,mlon_t)
            corner_points_lat[i+1,j+1] = np.mean(np.array([mlat_t['lat2'],mlon_t['lat2']]))
            corner_points_lon[i+1,j+1] = np.mean(np.array([mlat_t['lon2'],mlon_t['lon2']]))
            
    # edges
    for i in list(range(x-1)):
        edge_low = geod.InverseLine(posx_lat[i,0],posx_lon[i,0],posx_lat[i,1],posx_lon[i,1])
        edge_high = geod.InverseLine(posx_lat[i,-2],posx_lon[i,-2],posx_lat[i,-1],posx_lon[i,-1])
        edge_low_t = edge_low.Position(-0.5*edge_low.s13)
        edge_high_t = edge_high.Position(1.5*edge_high.s13)
        corner_points_lat[i+1,0] = edge_low_t['lat2']
        corner_points_lon[i+1,0] = edge_low_t['lon2']
        corner_points_lat[i+1,-1] = edge_high_t['lat2']
        corner_points_lon[i+1,-1] = edge_high_t['lon2']
    
    for j in list(range(y-1)):
        edge_left = geod.InverseLine(posy_lat[0,j],posy_lon[0,j],posy_lat[1,j],posy_lon[1,j])
        edge_right = geod.InverseLine(posy_lat[-2,j],posy_lon[-2,j],posy_lat[-1,j],posy_lon[-1,j])
        edge_left_t = edge_left.Position(-0.5*edge_left.s13)
        edge_right_t = edge_right.Position(1.5*edge_right.s13)
        corner_points_lat[0,j+1] = edge_left_t['lat2']
        corner_points_lon[0,j+1] = edge_left_t['lon2']
        corner_points_lat[-1,j+1] = edge_right_t['lat2']
        corner_points_lon[-1,j+1] = edge_right_t['lon2']
        
    # corners
    ## lower left
    corn_x = geod.InverseLine(corner_points_lat[1,0],corner_points_lon[1,0],
                                  corner_points_lat[2,0],corner_points_lon[2,0])
    corn_y = geod.InverseLine(corner_points_lat[0,1],corner_points_lon[0,1],
                                  corner_points_lat[0,2],corner_points_lon[0,2])
    corn_x_t = corn_x.Position(-1*corn_x.s13)
    corn_y_t = corn_y.Position(-1*corn_y.s13)
    corner_points_lat[0,0] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[0,0] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
    ## lower right
    corn_x = geod.InverseLine(corner_points_lat[-3,0],corner_points_lon[-3,0],
                                     corner_points_lat[-2,0],corner_points_lon[-2,0])
    corn_y = geod.InverseLine(corner_points_lat[-1,1],corner_points_lon[-1,1],
                                     corner_points_lat[-1,2],corner_points_lon[-1,2])
    corn_x_t = corn_x.Position(2*corn_x.s13)
    corn_y_t = corn_y.Position(-1*corn_y.s13)
    corner_points_lat[-1,0] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[-1,0] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
    ## upper left
    corn_x = geod.InverseLine(corner_points_lat[1,-1],corner_points_lon[1,-1],
                                  corner_points_lat[2,-1],corner_points_lon[2,-1])
    corn_y = geod.InverseLine(corner_points_lat[0,-3],corner_points_lon[0,-3],
                                  corner_points_lat[0,-2],corner_points_lon[0,-2])
    corn_x_t = corn_x.Position(-1*corn_x.s13)
    corn_y_t = corn_y.Position(2*corn_y.s13)
    corner_points_lat[0,-1] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[0,-1] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
    ## upper right
    corn_x = geod.InverseLine(corner_points_lat[-3,-1],corner_points_lon[-3,-1],
                                  corner_points_lat[-2,-1],corner_points_lon[-2,-1])
    corn_y = geod.InverseLine(corner_points_lat[-1,-3],corner_points_lon[-1,-3],
                                  corner_points_lat[-1,-2],corner_points_lon[-1,-2])
    corn_x_t = corn_x.Position(2*corn_x.s13)
    corn_y_t = corn_y.Position(2*corn_y.s13)
    corner_points_lat[-1,-1] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[-1,-1] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
     
    net_corner_lat[:,:] = corner_points_lat.T
    net_corner_lon[:,:] = corner_points_lon.T
    ### actual area calculation
    point_area = np.zeros((x,y))
    shape = geod.Polygon()
    for i in list(range(x)):
        for j in list(range(y)):
            shape.AddPoint(corner_points_lat[i+1,j+1],corner_points_lon[i+1,j+1])
            shape.AddPoint(corner_points_lat[i,j+1],corner_points_lon[i,j+1])
            shape.AddPoint(corner_points_lat[i,j],corner_points_lon[i,j])
            shape.AddPoint(corner_points_lat[i+1,j],corner_points_lon[i+1,j])
            number, perimeter, area = shape.Compute()
            point_area[i,j] = area/1000000. #conversion to square km
            shape.Clear()
        
    net_area[:,:] = point_area.T
    
    net_dim_steps = prepro.createDimension("time_s",len(run))
    net_time_steps = prepro.createVariable("time_steps","i4",("time_s"))
    net_watermask = prepro.createVariable("watermask","f8",("time_s","y","x"))
    net_watermask.description = "(Ocean-) Watermask (1=land,0=water)"
    net_watermask.units = "None"
    net_time_steps.description = "Times of the intermediate steps"
    net_time_steps.unit = "years"
    net_Acc_EHEP = prepro.createVariable("EHEP","f8",("time_s","y","x"))
    net_Acc_EHEP.description = "Existence Potential for every major and intermediate step"
    net_Acc_EHEP.unit = "Humans / 100 km^2"
    net_Acc_EHEP_num = prepro.createVariable("EHEPnumb","f8",("time_s","y","x"))
    net_Acc_EHEP_num.description = "Existence Potential for every major and intermediate step, total numbers"
    net_Acc_EHEP_num.unit = "Humans"
    
    net_hnumber[:,:] = density[:,:] * net_area[:,:]
    
    for i,val in enumerate(run):
        net_time_steps[i] = val
        net_Acc_EHEP[i,:,:] = EHEP[i,:,:]
        net_Acc_EHEP_num[i,:,:] = EHEP[i,:,:] * net_area
        net_watermask[i,:,:] = water_mask[:,:]

    prepro.close()
    print("Saved preprocessing for the mobility model")
    
# Main Body of the Program
if __name__ == "__main__":
    
    path_input = "/data/sfb806/human_mobility_model/christian/Applic_Auri_Grav/Exp2/"
    file_input = "HEP_chain_auri_grav_45k_25k_glac_extended.nc"
    path_output = "/data/sfb806/human_mobility_model/christian/for_yaping/plots_paper_cwegener/first_review/simple_cases_calcs/BS07/"
    file_output = "Preprocessed_auri_grav_45k_25k_glac_extended_new_rand.nc"
    switch_capa = True
    switch_lineargrowth = False
    endcapa = 0.05 # capacity estimation, end value for lineargrowth=true
    startcapa = 0.05 # capacity estimation, start value for lineargrowth=true and fixed value for false
    starting_dist = 'random'
    uni_dist = 1.
    gridtype = 'latlon'
    
    data_in = Dataset(path_input+file_input,mode='r')

    
    water_mask = data_in.variables['watermask'][:][:]
    EHEP = data_in.variables['AccHEP'][:][:][:]
    run = data_in.variables['time'][:]
    lat = data_in.variables['lat'][:]
    lon = data_in.variables['lon'][:]
    growth = (endcapa - startcapa)/len(run)
    
    preprocess(EHEP,lat,lon,x=None,y=None)
