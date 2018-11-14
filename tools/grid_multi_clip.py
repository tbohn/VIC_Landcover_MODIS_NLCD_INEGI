#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt


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
    infile = ''
    maskfile = ''
    prefix = ''
    outdir = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:m:p:o:",["infile=","maskfile=","prefix=","outdir="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -m <maskfile> -p <prefix> -o <outdir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -m <maskfile> -p <prefix> -o <outdir>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-m", "--maskfile"):
            maskfile = arg
        elif opt in ("-p", "--prefix"):
            prefix = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg


    # Define fill values
    fill_value_mask = np.int32(0)
    fill_value_int = np.int32(-1)
    fill_value_str = ''
    fill_value_float = 9.96920996838687e+36

    # Compute suffix
    basename = os.path.basename(infile)
    suffix = basename[len(prefix):]
    if suffix[0] == '_' or suffix[0] == '.':
      suffix = suffix[1:]

    # Read input dataset
    ds = xr.open_dataset(infile)

    latvar_tmp = ds['lat']
    lonvar_tmp = ds['lon']

    nLat = len(latvar_tmp)
    nLon = len(lonvar_tmp)
    minLat = round(np.asscalar(latvar_tmp[0]),6)
    maxLat = round(np.asscalar(latvar_tmp[-1]),6)
    minLon = round(np.asscalar(lonvar_tmp[0]),6)
    maxLon = round(np.asscalar(lonvar_tmp[-1]),6)
    resolution = (maxLat-minLat)/(nLat-1)
    minLat -= 0.5*resolution
    maxLat += 0.5*resolution
    minLon -= 0.5*resolution
    maxLon += 0.5*resolution
    llcorner = [minLat, minLon]

    # Hack for masks that don't have fillvalues
    if 'mask' in ds.variables:
        ds['mask'].encoding['_FillValue'] = 0

    # Read mask file
    (mask, nLat_mask, nLon_mask, llcorner_mask, resolution_mask, nodata) = \
        read_ascii_gridfile(maskfile, np.int)

    # Create new lat/lon dimensions large enough to hold the union of
    # the input dataset and the mask
    minLat_new = min(minLat, llcorner_mask[0])
    maxLat_new = max(maxLat, llcorner_mask[0]+nLat_mask*resolution_mask)
    minLon_new = min(minLon, llcorner_mask[1])
    maxLon_new = max(maxLon, llcorner_mask[1]+nLon_mask*resolution_mask)
    nLat_new = int((maxLat_new-minLat_new)/resolution)
    nLon_new = int((maxLon_new-minLon_new)/resolution)
    lat_new = np.asarray([round(minLat_new + (i+0.5)*resolution, 6) for i in
                          range(nLat_new)], dtype=np.float)
    lon_new = np.asarray([round(minLon_new + (i+0.5)*resolution, 6) for i in
                          range(nLon_new)], dtype=np.float)
    mask_new = np.ones([nLat_new, nLon_new], dtype=np.int)
    llcorner_new = [minLat_new, minLon_new]

    # Assemble new coord vars
    coords_new = {
                  'lat': (['lat'], lat_new),
                  'lon': (['lon'], lon_new),
                 }
    for varname in ds.coords:
        if varname not in ['lat','lon']:
            coords_new[varname] = ds[varname]

    # Create new dataset
    ds_new = xr.Dataset(coords=coords_new)
    for varname in ['lat','lon']:
        ds_new[varname].attrs = ds[varname].attrs
        ds_new[varname].encoding = ds[varname].encoding

    # Create new variables to store old variables (with padding)
    for varname in ds.variables:
        if varname not in ds_new.coords:
            shape = ds[varname].shape
            shape_new = []
            for dim in ds[varname].dims:
                if dim in ['lat', 'lon']:
                    shape_new.append(ds_new.dims[dim])
                else:
                    shape_new.append(ds.dims[dim])
            dtype = ds[varname].dtype
            if '_FillValue' in ds[varname].encoding.keys():
                fill_value = ds[varname].encoding['_FillValue']
            elif dtype == int or dtype == np.int16 or dtype == np.int32:
                fill_value = fill_value_int
            else:
                fill_value = fill_value_float
            if 'lat' in ds[varname].dims and 'lon' in ds[varname].dims:
                tmp = clip_to_mask(ds[varname], dtype, shape, llcorner,
                                   resolution, fill_value, mask_new,
                                   shape_new, llcorner_new, resolution,
                                   fill_value)
            else:
                tmp = ds[varname]
            ds_new[varname] = (ds[varname].dims, tmp)
            ds_new[varname].attrs = ds[varname].attrs
            ds_new[varname].encoding = ds[varname].encoding

    # loop over mask and clip ds_new for each tile
    for y in range(nLat_mask):
        south = llcorner_mask[0] + (y * resolution_mask)
        north = south + resolution_mask
        for x in range(nLon_mask):
            if mask[y,x]:
                west = llcorner_mask[1] + (x * resolution_mask)
                east = west + resolution_mask
                y0 = int((south-minLat_new)/resolution)
                y1 = int((north-minLat_new)/resolution)
                x0 = int((west-minLon_new)/resolution)
                x1 = int((east-minLon_new)/resolution)

                ds_out = ds_new.isel(lat=slice(y0,y1),lon=slice(x0,x1))

                # Write to output file
                outfile = outdir + '/' + prefix + '.' + '{:.0f}'.format(south) + '_' + '{:.0f}'.format(north) + 'n.' + '{:.0f}'.format(west) + '_' + '{:.0f}'.format(east) + 'e.' + suffix
                print('writing to',outfile)
                ds_out.to_netcdf(outfile)

    ds.close()
    ds_new.close()

if __name__ == "__main__":
    main()
