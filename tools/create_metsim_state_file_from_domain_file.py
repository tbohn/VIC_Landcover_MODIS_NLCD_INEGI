#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import pandas as pd

def main():
    infile = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:s:o:",["infile=","statedate=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -s <statedate> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -s <statedate> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-s", "--statedate"):
            statedate = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Read input domain file
    ds = xr.open_dataset(infile)

    lat = ds['lat']
    lon = ds['lon']
    mask = ds['mask']

    nLat = len(lat)
    nLon = len(lon)

    prec_value = 3
    tmin_value = 10
    tmax_value = 20
    swe_value = 0

    # Time
    enddate = pd.to_datetime(statedate)
    dt = pd.Timedelta('89 days')
    startdate = enddate - dt
    time = pd.date_range(start=startdate,end=enddate,freq='D')
    nTime = len(time)

    # Variables
    prec = np.empty((nTime,nLat,nLon),dtype=float)
    tmin = np.empty((nTime,nLat,nLon),dtype=float)
    tmax = np.empty((nTime,nLat,nLon),dtype=float)
    swe = np.empty((nTime,nLat,nLon),dtype=float)
    prec[:] = mask * prec_value
    tmin[:] = mask * tmin_value
    tmax[:] = mask * tmax_value
    swe[:] = mask * swe_value

    # Write output file
    ds_out = xr.Dataset(
                        {
                         'prec': (['time','lat','lon'],prec),
                         't_min': (['time','lat','lon'],tmin),
                         't_max': (['time','lat','lon'],tmax),
                         'swe': (['time','lat','lon'],swe),
                        },
                        coords={
                                'lat': (['lat'],lat),
                                'lon': (['lon'],lon),
                                'time': (['time'],time),
                               }
                       )

    for varname in ['lat', 'lon']:
        ds_out[varname].attrs = ds[varname].attrs
        ds_out[varname].encoding = ds[varname].encoding

    ds_out['prec'].attrs['units'] = 'mm'
    ds_out['prec'].attrs['long_name'] = 'Daily Precipitation'
    ds_out['t_min'].attrs['units'] = 'C'
    ds_out['t_min'].attrs['long_name'] = 'Daily Minimum Temperature'
    ds_out['t_max'].attrs['units'] = 'C'
    ds_out['t_max'].attrs['long_name'] = 'Daily Maximum Temperature'
    ds_out['swe'].attrs['units'] = 'mm'
    ds_out['swe'].attrs['long_name'] = 'Snow Water Equivalent'

    print('writing to',outfile)
    ds_out.to_netcdf(outfile)

    ds_out.close()
    ds.close()

if __name__ == "__main__":
    main()    
