#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def read_ascii_gridfile(filename, dtype):

    # Read header
    f = open(filename, 'r')
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    header4 = f.readline()
    header5 = f.readline()
    header6 = f.readline()
    ncols = int(header1.partition(' ')[2])
    nrows = int(header2.partition(' ')[2])
    xllcorner = float(header3.partition(' ')[2])
    yllcorner = float(header4.partition(' ')[2])
    cellsize = float(header5.partition(' ')[2])
    nodata = float(header6.partition(' ')[2])
    f.close()

    # Allocate data
    data = np.empty((nrows,ncols), dtype=dtype)

    # Read data
    tmp = np.loadtxt(filename,dtype=dtype,delimiter =' ', skiprows=6)
    for i in range(nrows):
#        data[nrows-1-i] = tmp[i].copy()
        data[i] = tmp[i].copy()

    llcorner = [yllcorner, xllcorner]

    return (data, nrows, ncols, llcorner, cellsize, nodata)


def main():
    domainfile = ''
    lc_map_file = ''
    lc_table = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hd:m:t:o:",["domainfile=","lc_map_file=","lc_table=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -d <domainfile> -m <lc_map_file> -t <lc_table> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -d <domainfile> -m <lc_map_file> -t <lc_table> -o <outfile>')
            sys.exit()
        elif opt in ("-d", "--domainfile"):
            domainfile = arg
        elif opt in ("-m", "--lc_map_file"):
            lc_map_file = arg
        elif opt in ("-t", "--lc_table"):
            lc_table = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Read domain file
    ds = xr.open_dataset(domainfile)

    nLat = ds['mask'].shape[0]
    nLon = ds['mask'].shape[1]

    lat = np.empty([nLat])
    lon = np.empty([nLon])
    mask = np.empty([nLat,nLon])

    lat[:] = ds['lat'][:]
    lon[:] = ds['lon'][:]
    mask[:] = ds['mask'][:]

    cellsize = round((np.max(lat) - np.min(lat)) / float(nLat - 1),6)
    minlat = round(np.min(lat) - (0.5 * cellsize),6)
    minlon = round(np.min(lon) - (0.5 * cellsize),6)
    maxlat = round(np.max(lat) + (0.5 * cellsize),6)
    maxlon = round(np.max(lon) + (0.5 * cellsize),6)

    # Read LC table
    lc_data = np.loadtxt(lc_table, dtype=bytes, delimiter=',').astype(str)
    nClass = lc_data.shape[0]
    v_of_classID = {}
    classID = np.empty([nClass]).astype(int)
    for v in range(nClass):
        classID[v] = int(lc_data[v,0])
        v_of_classID[classID[v]] = v

    # LC map
    (lc_map, nrows_lc, ncols_lc, llcorner_lc, cellsize_lc, nodata_lc) = \
        read_ascii_gridfile(lc_map_file, np.int32)

    # Aggregate lc classes over domain
    Cv = np.zeros([nClass,nLat,nLon])
    count_tot = np.zeros([nLat,nLon])
    for i in range(nrows_lc):
        lat_lc = llcorner_lc[0] + (i + 0.5) * cellsize_lc
        row = int( (lat_lc - minlat) / cellsize)
        print(i,lat_lc,row,lat[row])
        if row < 0 or row >= nLat:
            continue
        for j in range(ncols_lc):
            lon_lc = llcorner_lc[1] + (j + 0.5) * cellsize_lc
            col = int( (lon_lc - minlon) / cellsize)
            if col < 0 or col >= nLon:
                continue
            if (mask[row,col] == 1 and lc_map[i,j] != nodata_lc and
                lc_map[i,j] in v_of_classID.keys()):
                v = v_of_classID[lc_map[i,j]]
                Cv[v,row,col] += 1
                count_tot[row,col] += 1

    for v in range(nClass):
        Cv[v] /= count_tot

    # Write output file
    ds_out = xr.Dataset(
        {
            'mask': (['lat','lon'], mask),
            'Cv': (['class','lat','lon'], Cv),
        },
        coords = {
            'lat': (['lat'], lat),
            'lon': (['lon'], lon),
            'class': (['class'], classID),
        }
    )

    for varname in ['lat','lon','mask']:
        ds_out[varname].attrs = ds[varname].attrs

    ds_out['Cv'].attrs['long_name'] = 'Area fraction for each land cover class'

    for varname in ['lat','lon','mask','Cv']:
        ds_out[varname].encoding['zlib'] = True

    print('writing to',outfile)
    ds_out.to_netcdf(outfile)

    ds_out.close()
    ds.close()

if __name__ == "__main__":
    main()
