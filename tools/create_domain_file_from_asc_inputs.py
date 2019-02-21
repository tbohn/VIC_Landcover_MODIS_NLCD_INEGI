#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import pandas as pd
import pysal

def fill_gaps(data, gapval, default_val):

    # record where gaps exist
    a = np.where(data == gapval)[0]
    b = np.where(data == gapval)[1]
    ngaps = len(a)
    if ngaps == 0:
        return(data)

    # replace whatever value was used to denote gaps with nan
    data[data == gapval] = np.nan

    # prepare large array to do summing
    shape = [8]
    for i in range(len(data.shape)):
        shape.append(data.shape[i])
    bigdata = np.empty(shape)

    # each layer of the large array is a shift in one of 8 directions
    # (corresponding to each of 8 neighbors)
    # null out the row or column of data that would be rolled to the other end
    data_d = data
    data_d[-1] = np.nan
    data_d = np.roll(data_d, 1, 0)
    data_u = data
    data_u[0] = np.nan
    data_u = np.roll(data_u, -1, 0)
    data_r = data
    data_r[:,-1] = np.nan
    data_r = np.roll(data_r, 1, 1)
    data_l = data
    data_l[:,0] = np.nan
    data_l = np.roll(data_l, -1, 1)
    data_dr = data_r
    data_dr[-1] = np.nan
    data_dr = np.roll(data_dr, 1, 0)
    data_ur = data_r
    data_ur[0] = np.nan
    data_ur = np.roll(data_ur, -1, 0)
    data_dl = data_l
    data_dl[-1] = np.nan
    data_dl = np.roll(data_dl, 1, 0)
    data_ul = data_l
    data_ul[0] = np.nan
    data_ul = np.roll(data_ul, -1, 0)

    bigdata[0] = data_d
    bigdata[1] = data_u
    bigdata[2] = data_r
    bigdata[3] = data_l
    bigdata[4] = data_dr
    bigdata[5] = data_dl
    bigdata[6] = data_ur
    bigdata[7] = data_ul

    # now, compute means, ignoring nans
    meandata = np.nanmean(bigdata, axis=0)

    # replace gaps with the means
    for i,j in zip(a,b):
        data[i,j] = meandata[i,j]

    # replace any unfilled gaps with default value
    for i,j in zip(a,b):
        if np.isnan(data[i,j]):
            data[i,j] = default_val

    return(data)


def map_xy(nrows, ncols, llcorner, cellsize, nrows2, ncols2, llcorner2,
           cellsize2):

    map_y = {}
    for y1 in range(nrows):
        ctr_lat = llcorner[0] + (y1 + 0.5) * cellsize
        if ctr_lat > llcorner2[0]:
            y2 = int((ctr_lat - llcorner2[0]) / cellsize2)
            if y2 < 0 or y2 >= nrows2:
                continue
            map_y[y1] = y2

    map_x = {}
    for x1 in range(ncols):
        ctr_lon = llcorner[1] + (x1 + 0.5) * cellsize
        if ctr_lon > llcorner2[1]:
            x2 = int((ctr_lon - llcorner2[1]) / cellsize2)
            if x2 < 0 or x2 >= ncols2:
                continue
            map_x[x1] = x2

    return(map_y, map_x)


def clip_to_mask(data_in, dtype, shape_in, llcorner_in, cellsize_in,
                 nodata_in, mask, shape, llcorner, cellsize, nodata):

    nrows_in = shape_in[-2]
    ncols_in = shape_in[-1]
    nrows = shape[-2]
    ncols = shape[-1]

    [map_y, map_x] = map_xy(nrows, ncols, llcorner, cellsize, nrows_in,
                            ncols_in, llcorner_in, cellsize_in)

    y_min = min(map_y.keys())
    y_max = max(map_y.keys())
    y_in_min = map_y[y_min]
    y_in_max = map_y[y_max]
    x_min = min(map_x.keys())
    x_max = max(map_x.keys())
    x_in_min = map_x[x_min]
    x_in_max = map_x[x_max]

    data = np.full(shape, nodata, dtype=dtype)
    data[..., y_min:y_max+1, x_min:x_max+1] = data_in[..., y_in_min:y_in_max+1, x_in_min:x_in_max+1]
    data[..., mask != 1] = nodata

    return data


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
        data[nrows-1-i] = tmp[i].copy()

    llcorner = [yllcorner, xllcorner]

    return (data, nrows, ncols, llcorner, cellsize, nodata)


def main():
    maskfile = ''
    fracmaskfile = ''
    demfile = ''
    prcp_param_dir = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hm:f:d:p:o:",["maskfile=","fracmaskfile=","demfile=","prcp_param_dir=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -m <maskfile> -f <fracmaskfile> -d <demfile> [-p <prcp_param_dir>] -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -m <maskfile> -f <fracmaskfile> -d <demfile> [-p <prcp_param_dir>] -o <outfile>')
            sys.exit()
        elif opt in ("-m", "--maskfile"):
            maskfile = arg
        elif opt in ("-f", "--fracmaskfile"):
            fracmaskfile = arg
        elif opt in ("-d", "--demfile"):
            demfile = arg
        elif opt in ("-p", "--prcp_param_dir"):
            prcp_param_dir = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Constants
    const_MINPHOUR = 60
    const_HOURPDAY = 24
    const_MINPDAY = const_MINPHOUR * const_HOURPDAY
    fill_value_mask = np.int32(0)
    fill_value_int = np.int32(-9999)
    fill_value_float = 9.96920996838687e+36
    default_value_DEM = 0
    default_value_frac = 0
    default_value_dur = 240 # minutes
    default_value_t_pk = 960 # hours

    # Read mask file
    (mask, nrows, ncols, llcorner, cellsize, nodata) = \
        read_ascii_gridfile(maskfile, np.int32)
    shape_mask = [nrows, ncols]
    a = np.where(mask != 1)[0]
    b = np.where(mask != 1)[1]

    # Define coord vars
    lat = np.asarray([round(llcorner[0] + (i+0.5)*cellsize, 6) for i in
                      range(nrows)], dtype=np.float)
    lon = np.asarray([round(llcorner[1] + (i+0.5)*cellsize, 6) for i in
                      range(ncols)], dtype=np.float)

    # Fractional mask - if not supplied, set equal to mask
    if (fracmaskfile == 'null' or fracmaskfile == ''):
        frac = mask.astype(np.single)
        shape_frac = shape_mask
    else:
        (frac_tmp, nrows_frac, ncols_frac, llcorner_frac, cellsize_frac,
            nodata_frac) = read_ascii_gridfile(fracmaskfile, np.single)
        shape_frac = [nrows_frac, ncols_frac]

        # mask the frac data
        frac = clip_to_mask(frac_tmp, np.single, shape_frac, llcorner_frac,
                            cellsize_frac, nodata_frac, mask, shape_mask,
                            llcorner, cellsize, np.nan)
        frac = fill_gaps(frac, nodata_frac, default_value_frac)

    # Read DEM
    (DEM_in, nrows_DEM, ncols_DEM, llcorner_DEM, cellsize_DEM, nodata_DEM) = \
        read_ascii_gridfile(demfile, np.single)
    shape_DEM = [nrows_DEM, ncols_DEM]

    # mask the elev data
    elev = clip_to_mask(DEM_in, np.single, shape_DEM, llcorner_DEM,
                        cellsize_DEM, nodata_DEM, mask, shape_mask, llcorner,
                        cellsize, np.nan)
    elev = fill_gaps(elev, nodata_DEM, default_value_DEM)

    if prcp_param_dir != '':
        # Read prcp param files
        month = np.arange(12,dtype=int);
        dur = np.empty((12,nrows,ncols), dtype=np.single)
        t_pk= np.empty((12,nrows,ncols), dtype=np.single)
        for m in month:
            monthstr = '{:02d}'.format(m+1)
            file = prcp_param_dir + '/MeanRainDur_Month' + monthstr + '.asc'
            (dur_tmp, nrows_tmp, ncols_tmp, llcorner_tmp, cellsize_tmp,
                nodata_tmp) = read_ascii_gridfile(file, np.single)
            shape_tmp = [nrows_tmp, ncols_tmp]
            dur_tmp[dur_tmp != nodata_tmp] *= const_MINPHOUR
            dur_tmp[np.where((dur_tmp != nodata_tmp) & (dur_tmp < 1))] = 1
            dur_tmp[np.where((dur_tmp != nodata_tmp) &
                             (dur_tmp > const_MINPDAY - 1))] = const_MINPDAY
            dur[m] = clip_to_mask(dur_tmp, np.single, shape_tmp, llcorner_tmp,
                                  cellsize_tmp, nodata_tmp, mask, shape_mask,
                                  llcorner, cellsize, np.nan)
            dur[m] = fill_gaps(dur[m], nodata_tmp, default_value_dur)
            file = prcp_param_dir + '/TimePeakDiurnalCycle_Month' + monthstr + '.asc'
            (t_pk_tmp, nrows_tmp, ncols_tmp, llcorner_tmp, cellsize_tmp,
                nodata_tmp) = read_ascii_gridfile(file, np.single)
            shape_tmp = [nrows_tmp, ncols_tmp]
            t_pk_tmp[t_pk_tmp != nodata_tmp] *= const_MINPHOUR
            t_pk_tmp = np.where(t_pk_tmp != nodata_tmp, t_pk_tmp %
                                const_MINPDAY, t_pk_tmp)
            t_pk[m] = clip_to_mask(t_pk_tmp, np.single, shape_tmp, llcorner_tmp,
                                  cellsize_tmp, nodata_tmp, mask, shape_mask,
                                  llcorner, cellsize, np.nan)
            t_pk[m] = fill_gaps(t_pk[m], nodata_tmp, default_value_t_pk)
#        for i,j in zip(a,b):
#            dur[:,i,j] = np.nan
#            t_pk[:,i,j] = np.nan

    # Compute area
    area = np.empty([nrows,ncols],dtype=np.single)
    for i in range(nrows):
        dx = pysal.cg.sphere.arcdist((lon[0] - (cellsize / 2),lat[i]),(lon[0] +
                                    (cellsize / 2),lat[i])) * 1000
        dy = pysal.cg.sphere.arcdist((lon[0],lat[i] - (cellsize / 2)),(lon[0],
                                    lat[i] + (cellsize / 2))) * 1000
        area[i] = dx * dy
    area = area.copy() * mask
    area = np.where(area == 0, np.nan, area).astype(np.single)

    # Write output file
    ds = xr.Dataset(
        {
            'mask': (['lat','lon'], mask),
            'frac': (['lat','lon'], frac),
            'area': (['lat','lon'], area),
            'elev': (['lat','lon'], elev),
        },
        coords = {
            'lat': (['lat'], lat),
            'lon': (['lon'], lon),
        }
    )

    ds['lat'].attrs['long_name'] = 'Grid cell center latitude'
    ds['lat'].attrs['standard_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'degrees North'
    ds['lat'].attrs['axis'] = 'Y'
    ds['lon'].attrs['long_name'] = 'Grid cell center longitude'
    ds['lon'].attrs['standard_name'] = 'longitude'
    ds['lon'].attrs['units'] = 'degrees East'
    ds['lon'].attrs['axis'] = 'X'
    ds['mask'].attrs['long_name'] = 'Domain mask (1=inside domain, 0=outside)'
    ds['mask'].attrs['units'] = 'unitless'
    ds['mask'].attrs['missing_value'] = fill_value_mask
    ds['mask'].encoding['_FillValue'] = fill_value_mask
    ds['frac'].attrs['long_name'] = 'Fraction of grid cell that is within domain'
    ds['frac'].attrs['units'] = 'unitless'
    ds['elev'].attrs['long_name'] = 'Grid cell elevation above sea level'
    ds['elev'].attrs['units'] = 'm'
    ds['area'].attrs['long_name'] = 'Grid cell area'
    ds['area'].attrs['units'] = 'm2'

    for varname in ['lat','lon','mask','frac','elev','area']:
        ds[varname].encoding['zlib'] = True

    if prcp_param_dir != '':
        ds['month'] = (['month'], month)
        ds = ds.set_coords('month')
        ds['dur'] = (['month','lat','lon'], dur)
        ds['dur'].attrs['long_name'] = 'Precipitation duration (minutes)'
        ds['dur'].attrs['units'] = 'min'
        ds['t_pk'] =(['month','lat','lon'], t_pk)
        ds['t_pk'].attrs['long_name'] = 'Time of peak precipitation intensity (minutes since beginning of day)'
        ds['t_pk'].attrs['units'] = 'min'

        for varname in ['dur','t_pk']:
            ds[varname].encoding['zlib'] = True

    print('writing to',outfile)
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()
