#!/usr/bin/env python

import os
import re
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


def copy_data(data_in, dtype, shape_in, llcorner_in, cellsize_in,
              nodata_in, data, shape, llcorner, cellsize, nodata):

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

    data[..., y_min:y_max+1, x_min:x_max+1] = data_in[..., y_in_min:y_in_max+1, x_in_min:x_in_max+1]

    return data


def main():
    domainfile = ''
    indir = ''
    prefix = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hd:i:p:o:",["domainfile=","indir=","prefix=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -d <domainfile> -i <indir> -p <prefix> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -d <domainfile> -i <indir> -p <prefix> -o <outfile>')
            sys.exit()
        elif opt in ("-d", "--domainfile"):
            domainfile = arg
        elif opt in ("-i", "--indir"):
            indir = arg
        elif opt in ("-p", "--prefix"):
            prefix = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg


    # Define default fill values
    fill_value_mask = np.int32(0)
    fill_value_int = np.int32(-1)
    fill_value_str = ''
    fill_value_float = 9.96920996838687e+36

    # Read domain file
    ds_dom = xr.open_dataset(domainfile)

    lat_mask = ds_dom['lat']
    lon_mask = ds_dom['lon']
    mask = ds_dom['mask']

    nLat_mask = len(lat_mask)
    nLon_mask = len(lon_mask)
    minLat_mask = round(np.asscalar(lat_mask[0]),6)
    maxLat_mask = round(np.asscalar(lat_mask[-1]),6)
    minLon_mask = round(np.asscalar(lon_mask[0]),6)
    maxLon_mask = round(np.asscalar(lon_mask[-1]),6)
    resolution_mask = (maxLat_mask-minLat_mask)/(nLat_mask-1)
    minLat_mask -= 0.5*resolution_mask
    maxLat_mask += 0.5*resolution_mask
    minLon_mask -= 0.5*resolution_mask
    maxLon_mask += 0.5*resolution_mask
    llcorner_mask = [minLat_mask, minLon_mask]

    # Loop over infiles
    out_dict = {}
    out_attrs = {}
    out_encoding = {}
    first = True
    filelist = os.listdir(indir)
    for file in filelist:

        if not re.search(prefix, file):
            continue

        # build filename
        file_with_path = indir + '/' + file

        # open file
        print('reading',file_with_path)
        ds = xr.open_dataset(file_with_path)
        lat = ds['lat']
        lon = ds['lon']
        nLat = len(lat)
        nLon = len(lon)
        minLat = round(np.asscalar(lat[0]),6)
        maxLat = round(np.asscalar(lat[-1]),6)
        minLon = round(np.asscalar(lon[0]),6)
        maxLon = round(np.asscalar(lon[-1]),6)
        resolution = (maxLat-minLat)/(nLat-1)
        minLat -= 0.5*resolution
        maxLat += 0.5*resolution
        minLon -= 0.5*resolution
        maxLon += 0.5*resolution
        llcorner = [minLat, minLon]

        if first:

            first = False

            # Assemble new coord vars
            coords_out = {
                          'lat': (['lat'], lat_mask),
                          'lon': (['lon'], lon_mask),
                         }
            for varname in ds.coords:
                if varname not in ['lat','lon','lats','lons','LandMask','mask']:
                    coords_out[varname] = ds[varname]

            # Create output dataset
            print('creating ds_out')
            ds_out = xr.Dataset(coords=coords_out)
            for varname in ds_out.coords:
                if varname in ['lat','lon']:
                    ds_out[varname].attrs = ds_dom[varname].attrs
                    ds_out[varname].encoding = ds_dom[varname].encoding
                else:
                    ds_out[varname].attrs = ds[varname].attrs
                    ds_out[varname].encoding = ds[varname].encoding
            ds_out['mask'] = (['lat','lon'], mask)
            ds_out['mask'].attrs = ds_dom['mask'].attrs
            ds_out['mask'].encoding = ds_dom['mask'].encoding
            # this statement seems to be necesary to avoid an error in writing
            print('ds_out',ds_out)

            # Allocate space for output variables
            for varname in ds.variables:
                if varname not in ds.coords:
                    shape = ds[varname].shape
                    shape_tmp = []
                    for dim in ds[varname].dims:
                        if dim in ['lat', 'lon']:
                            shape_tmp.append(ds_dom.dims[dim])
                        else:
                            shape_tmp.append(ds.dims[dim])
                    dtype = ds[varname].dtype
                    if '_FillValue' in ds[varname].encoding.keys():
                        fill_value = ds[varname].encoding['_FillValue']
                    elif dtype == int or dtype == np.int16 or dtype == np.int32:
                        fill_value = fill_value_int
                    else:
                        fill_value = fill_value_float
                    out_dict[varname] = np.full(shape_tmp, fill_value, dtype=dtype)

        # Create tmp variables to store old variables (with padding)
        for varname in ds.variables:
            if varname not in ds.coords and varname != 'mask':
                print('processing',varname)
                shape = ds[varname].shape
                shape_out = []
                for dim in ds[varname].dims:
                    if dim in ['lat', 'lon']:
                        shape_out.append(ds_dom.dims[dim])
                    else:
                        shape_out.append(ds.dims[dim])
                dtype = ds[varname].dtype
                if '_FillValue' in ds[varname].encoding.keys():
                    fill_value = ds[varname].encoding['_FillValue']
                elif dtype == int or dtype == np.int16 or dtype == np.int32:
                    fill_value = fill_value_int
                else:
                    fill_value = fill_value_float
                if 'lat' in ds[varname].dims and 'lon' in ds[varname].dims:
                    out_dict[varname] = copy_data(ds[varname], dtype, shape,
                                                  llcorner, resolution,
                                                  fill_value, out_dict[varname],
                                                  shape_out, llcorner_mask,
                                                  resolution_mask, fill_value)
                else:
                    out_dict[varname] = ds[varname]

                out_attrs[varname] = ds[varname].attrs
                out_encoding[varname] = ds[varname].encoding

        ds.close()

    # Assign to output
    for varname in ds.variables:
        if varname not in ds.coords and varname != 'mask':
            print('writing',varname)
#            out_dict[varname][...,mask == 0] = out_encoding[varname]['_FillValue']
            ds_out[varname] = (ds[varname].dims, out_dict[varname])
            ds_out[varname].attrs = out_attrs[varname]
            ds_out[varname].encoding = out_encoding[varname]

    # Write to output file
    print('writing to',outfile)
    ds_out.to_netcdf(outfile)

    ds_dom.close()
    ds_out.close()

if __name__ == "__main__":
    main()
