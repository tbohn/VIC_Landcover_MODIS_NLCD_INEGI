#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import pandas as pd

def main():
    infile = ''
    domainfile = ''
    varnames = ['Prec','Tmin','Tmax','wind']
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:d:v:o:",["infile=","domainfile=","varnamelist=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -d <domainfile> -v <varnamelist> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -d <domainfile> -v <varnamelist> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-d", "--domainfile"):
            domainfile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Read domain file
    ds = xr.open_dataset(domainfile)

    lat = ds['lat']
    lon = ds['lon']
    mask = ds['mask']

    nLat = len(lat)
    nLon = len(lon)
    cellsize = round(np.asscalar((lat[-1] - lat[0]) / float(nLat - 1)),6)
    minlat = round(np.asscalar(lat[0] - (0.5 * cellsize)),6)
    minlon = round(np.asscalar(lon[0] - (0.5 * cellsize)),6)
    maxlat = round(np.asscalar(lat[-1] + (0.5 * cellsize)),6)
    maxlon = round(np.asscalar(lon[-1] + (0.5 * cellsize)),6)

    # Read infile
    ds_in = xr.open_dataset(infile)

    time = ds_in['time']
    lat_in = ds_in['lat']
    lon_in = ds_in['lon']

    nTime = len(time)
    nLat_in = len(lat_in)
    nLon_in = len(lon_in)
    cellsize_in = round(np.asscalar((lat_in[-1] - lat_in[0]) / float(nLat_in - 1)),6)
    minlat_in = round(np.asscalar(lat_in[0] - (0.5 * cellsize_in)),6)
    minlon_in = round(np.asscalar(lon_in[0] - (0.5 * cellsize_in)),6)
    maxlat_in = round(np.asscalar(lat_in[-1] + (0.5 * cellsize_in)),6)
    maxlon_in = round(np.asscalar(lon_in[-1] + (0.5 * cellsize_in)),6)

    # Open output file
    ds_out = xr.Dataset(
        coords={
            'lat': (['lat'],lat),
            'lon': (['lon'],lon),
            'time': (['time'],time),
        }
    )

    for varname in varnames:
        data_in = ds_in[varname]

        # Map input to output
        data_out = np.empty((nTime,nLat,nLon),dtype=np.single)
        for y in range(nLat):
            ctr_lat = minlat + (y + 0.5) * cellsize
            y_in = int((ctr_lat - minlat_in) / cellsize_in)
            for x in range(nLon):
                ctr_lon = minlon + (x + 0.5) * cellsize
                x_in = int((ctr_lon - minlon_in) / cellsize_in)
                data_out[:,y,x] = data_in[:,y_in,x_in]

        ds_out[varname] = (['time','lat','lon'],data_out)
        ds_out[varname].attrs = ds_in[varname].attrs
        ds_out[varname].encoding = ds_in[varname].encoding

    print('writing to',outfile)
    ds_out.to_netcdf(outfile,engine='scipy')

    ds_in.close()
    ds_out.close()
    ds.close()

if __name__ == "__main__":
    main()    
