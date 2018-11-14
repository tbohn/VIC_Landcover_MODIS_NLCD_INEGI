#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import math
import time

# Function for simple linear interpolation of temporal gaps
def interp_gaps(y, maxlen, cycle):

    N = len(y)
    maxlen = min(maxlen, N)
    maxlen = max(maxlen, 1)

    gap = {}
    istart = np.full(N,-1,dtype=int)

    # Identify gaps
    for i in range(N):
        if np.isnan(y[i]):
            if i > 0 and istart[i-1] >= 0:
                istart[i] = istart[i-1]
                gap[istart[i]] += 1
            else:
                istart[i] = i
                gap[istart[i]] = 1

    if cycle:
        if ( 0 in gap.keys()
             and gap[0] > 0
             and istart[N-1] > 0 ):
            gap[istart[N-1]] += gap[0]
            gap[0] = 0

    # Interpolate
    for i in (gap.keys()):
        if gap[i] > 0 and gap[i] < maxlen:
            x0 = i-1
            x1 = i + gap[i]
            if x0 >= 0:
                x0_adj = x0
            else:
                x0_adj = x0 + N
            if x1 < N:
                x1_adj = x1
            else:
                x1_adj = x1 - N
            y0 = y[x0_adj]
            y1 = y[x1_adj]
            slope = (y1 - y0) / (x1 - x0)
            for j in range(i,i+gap[i]):
                x = j % N
                dist = x - x0
                if dist < 0:
                    dist += N
                y[x] = y0 + dist * slope

    return y


# function to make a gaussian kernel
def make_kernel(rad_time, rad_space, sigma_time, sigma_space):
    ktlen = 2*rad_time+1
    kxlen = 2*rad_space+1
    kylen = 2*rad_space+1
    kernel = np.zeros([ktlen,kxlen,kylen])
    sum = 0
    for t in range(ktlen):
        for i in range(kxlen):
            for j in range(kylen):
              if (sigma_space > 0):
                  space_term = ((i-rad_space)**2+(j-rad_space)**2) / (2*sigma_space**2)
              else:
                  space_term = 0
              if (sigma_time > 0):
                  time_term = (t-rad_time)**2/(2*sigma_time**2)
              else:
                  time_term = 0
              kernel[t,i,j] = np.exp( -( space_term + time_term ) )
              sum += kernel[t,i,j]
    for t in range(ktlen):
        for i in range(kxlen):
            for j in range(kylen):
                kernel[t,i,j] /= sum
    return kernel

# gapfilling with a gaussian kernel
def gapfill (mydata, opmask, data_with_buffer, data_count_with_buffer, rt, rs, st, ss):

    # Note: ideally we could simplify this function

    tlen = mydata.shape[0]
    xlen = mydata.shape[1]
    ylen = mydata.shape[2]
    tlen_wbuf = data_with_buffer.shape[0]
    xlen_wbuf = data_with_buffer.shape[1]
    ylen_wbuf = data_with_buffer.shape[2]
    buft = int(0.5 * (tlen_wbuf - tlen))
    bufs = int(0.5 * (ylen_wbuf - ylen))

    # make kernel
    kernel = make_kernel(rt, rs, st, ss)
    [ktlen,kylen,kxlen] = kernel.shape

    # loop for filling isolated gaps with kernel-weighted spatio/temporal average
    for t in np.arange(tlen):

        t0 = t-rt
        t1 = t+rt+1
        t0_wbuf = t0+buft
        t1_wbuf = t1+buft

        # Loop over isolated gaps within this time slice
        a = np.where(opmask[t,:,:] == True)[0]
        b = np.where(opmask[t,:,:] == True)[1]
        for i,j in zip(a,b):

            # Fill tmpdata and tmpcount arrays with data in neighborhood surrounding this gap
            i0 = i-rs
            i1 = i+rs+1
            j0 = j-rs
            j1 = j+rs+1
            i0_wbuf = i0+bufs
            if i0_wbuf < 0:
                i0_offset = -i0_wbuf
                i0_wbuf = 0
            else:
                i0_offset = 0
            i1_wbuf = i1+bufs
            if i1_wbuf > ylen_wbuf:
                i1_wbuf = ylen_wbuf
            j0_wbuf = j0+bufs
            if j0_wbuf < 0:
                j0_offset = -j0_wbuf
                j0_wbuf = 0
            else:
                j0_offset = 0
            j1_wbuf = j1+bufs
            if j1_wbuf > xlen_wbuf:
                j1_wbuf = xlen_wbuf
            ilen = i1_wbuf-i0_wbuf
            jlen = j1_wbuf-j0_wbuf

            tmpdata = np.empty(kernel.shape,dtype=np.single)
            tmpdata[:,i0_offset:i0_offset+ilen,j0_offset:j0_offset+jlen] = data_with_buffer[t0_wbuf:t1_wbuf,i0_wbuf:i1_wbuf,j0_wbuf:j1_wbuf]
            tmpcount = np.empty(kernel.shape,dtype=np.single)
            tmpcount[:,i0_offset:i0_offset+ilen,j0_offset:j0_offset+jlen] = data_count_with_buffer[t0_wbuf:t1_wbuf,i0_wbuf:i1_wbuf,j0_wbuf:j1_wbuf]

            # Mask out any invalid data points
            tmpmask = 1 - np.isnan(tmpdata).astype(int)
            tmpdata = tmpdata * tmpmask.astype(float)
            tmpcount = tmpcount * tmpmask

            # If we have valid points to interpolate between, do interpolation
            weights = kernel * tmpcount
            wtdata = weights * tmpdata
            sumwts = np.nansum(weights)
            sumwtdata = np.nansum(wtdata)
            if (sumwts > 0):
                mydata[t,i,j] = sumwtdata / sumwts

    return mydata


def main():
    infile = ''
    outfile = ''
    maxlen = -1
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:v:o:t:c:m:",["infile=","varnamelist=","outfile=","type=","count_thresh=","maxlen="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outfile> -t <type> -c <count_thresh> -m <maxlen>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outfile> -t <type> -c <count_thresh> -m <maxlen>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outfile"):
            outfile = arg
        elif opt in ("-t", "--type"):
            stat_type = arg
        elif opt in ("-c", "--count_thresh"):
            count_thresh = int(arg)
        elif opt in ("-m", "--maxlen"):
            maxlen = int(arg)

    # Default values for gaps that can't otherwise be filled
    defaults = {}
    for tmpvar in ['LAI','NDVI','fcanopy','albedo']:
        defaults[tmpvar] = {}
    defaults['LAI']['mean'] = 0
    defaults['LAI']['std'] = 0
    defaults['LAI']['anom'] = 0
    defaults['NDVI']['mean'] = 0.1
    defaults['NDVI']['std'] = 0
    defaults['NDVI']['anom'] = 0
    defaults['fcanopy']['mean'] = 0.01
    defaults['fcanopy']['std'] = 0
    defaults['fcanopy']['anom'] = 0
    defaults['albedo']['mean'] = 0.15
    defaults['albedo']['std'] = 0
    defaults['albedo']['anom'] = 0

    ds = xr.open_dataset(infile)

    if (stat_type == 'clim'):
        timevarname = 'day_of_year'
    else:
        timevarname = 'time'
    timevar = ds[timevarname]
    latvar = ds['lat']
    lonvar = ds['lon']
    classvar = ds['veg_class']
    classname = ds['class_name']
    mask = ds['LandMask']
    Cv = ds['Cv']

    nTime = len(timevar)
    nClasses = len(classvar)
    nLat = len(latvar)
    nLon = len(lonvar)

    # stat_type-dependent parameters
    cycle = False
    set_boundaries = False
    if (stat_type == 'clim'):
        types = ['mean','std']
        cycle = True
    elif (stat_type == 'anom'):
        types = ['anom']
        count_thresh = 1
        set_boundaries = True
        boundary_val = 0
    else:
        types = ['']

    if maxlen < 1:
        maxlen = nTime - 2

    # dims of data with buffer
    buft = 1
    if cycle:
        bufs = nLat
    else:
        bufs = 1
    nTime_wbuf = nTime+2*buft
    nLat_wbuf = nLat+2*bufs
    nLon_wbuf = nLon+2*bufs

    # set offsets, ordered starting at upper left, moving left-right and top-down
    slices_buf = np.empty(8,dtype=slice)
    slices_buf[0] = [slice(None,None,1),slice(None,None,1),slice(-bufs,None,1),slice(-bufs,None,1)]
    slices_buf[1] = [slice(None,None,1),slice(None,None,1),slice(None,None,1),slice(-bufs,None,1)]
    slices_buf[2] = [slice(None,None,1),slice(None,None,1),slice(None,bufs,1),slice(-bufs,None,1)]
    slices_buf[3] = [slice(None,None,1),slice(None,None,1),slice(-bufs,None,1),slice(None,None,1)]
    slices_buf[4] = [slice(None,None,1),slice(None,None,1),slice(None,bufs,1),slice(None,None,1)]
    slices_buf[5] = [slice(None,None,1),slice(None,None,1),slice(-bufs,None,1),slice(None,bufs,1)]
    slices_buf[6] = [slice(None,None,1),slice(None,None,1),slice(None,None,1),slice(None,bufs,1)]
    slices_buf[7] = [slice(None,None,1),slice(None,None,1),slice(None,bufs,1),slice(None,bufs,1)]
    slices_data_wbuf = np.empty(8,dtype=slice)
    slices_data_wbuf[0] = [slice(None,None,1),slice(None,None,1),slice(None,bufs,1),slice(None,bufs,1)]
    slices_data_wbuf[1] = [slice(None,None,1),slice(None,None,1),slice(bufs,-bufs,1),slice(None,bufs,1)]
    slices_data_wbuf[2] = [slice(None,None,1),slice(None,None,1),slice(-bufs,None,1),slice(None,bufs,1)]
    slices_data_wbuf[3] = [slice(None,None,1),slice(None,None,1),slice(None,bufs,1),slice(bufs,-bufs,1)]
    slices_data_wbuf[4] = [slice(None,None,1),slice(None,None,1),slice(-bufs,None,1),slice(bufs,-bufs,1)]
    slices_data_wbuf[5] = [slice(None,None,1),slice(None,None,1),slice(None,bufs,1),slice(-bufs,None,1)]
    slices_data_wbuf[6] = [slice(None,None,1),slice(None,None,1),slice(bufs,-bufs,1),slice(-bufs,None,1)]
    slices_data_wbuf[7] = [slice(None,None,1),slice(None,None,1),slice(-bufs,None,1),slice(-bufs,None,1)]

    # Loop over variables
    for varname in varnames:

        print('processing',varname)

        # data count
        tmpvar_count = varname + '_count'
        data_count = ds[tmpvar_count]
        data_count = np.where(np.isnan(data_count),0,data_count)

        for type in (types):

            if (type != ''):
                tmpvar = varname + '_' + type
            else:
                tmpvar = varname

            print('processing',tmpvar)

            # read data from infile
            data = ds[tmpvar]

            # allocate gapfill_flag
            gapfill_flag = np.full([nTime,nClasses,nLat,nLon],0,dtype=int)

            # Null out data for which count is below threshold
            data[data_count < count_thresh] = np.nan
            should_exist = np.where(Cv > 0, 1, 0)
            exist_before = np.where(~np.isnan(data), 1, 0)

            # Temporal Gaps

            # Set up tmp arrays with optional buffer for boundaries
            if set_boundaries:
                data_wbuf_time = np.full([nTime_wbuf,nClasses,nLat,nLon],boundary_val,dtype=np.single)
                data_wbuf_time[buft:-buft] = data[:]
            else:
                data_wbuf_time = np.full([nTime,nClasses,nLat,nLon],0,dtype=np.single)
                data_wbuf_time[:] = data[:]

            # Do gapfilling
            data_wbuf_time = np.apply_along_axis(interp_gaps, 0, data_wbuf_time, maxlen, cycle)

            # Save back in data array
            if set_boundaries:
                data[:] = data_wbuf_time[buft:-buft]
            else:
                data[:] = data_wbuf_time[:]

            # Spatial Gaps (those unfilled by temporal gapfill)

            # Adjust data_count for any gaps that have been filled
            data_count_tmp = np.zeros(data_count.shape, dtype=int)
            data_count_tmp[:] = data_count[:].astype(int)
            exist_after = np.where(~np.isnan(data), 1, 0)
            data_count_tmp = np.maximum(data_count_tmp, exist_after * count_thresh)

            # Spatial buffer
            data_wbuf_space = np.full([nTime,nClasses,nLon_wbuf,nLat_wbuf],np.nan,dtype=np.single)
            data_count_wbuf_space = np.full([nTime,nClasses,nLon_wbuf,nLat_wbuf],0,dtype=np.int)
            for i in range(8):
                if (cycle):
                    data_wbuf_space[slices_data_wbuf[i][0],
                                    slices_data_wbuf[i][1],
                                    slices_data_wbuf[i][2],
                                    slices_data_wbuf[i][3]] \
                        = data[slices_buf[i][0],
                               slices_buf[i][1],
                               slices_buf[i][2],
                               slices_buf[i][3]].copy()
                    data_count_wbuf_space[slices_data_wbuf[i][0],
                                          slices_data_wbuf[i][1],
                                          slices_data_wbuf[i][2],
                                          slices_data_wbuf[i][3]] \
                        = data_count_tmp[slices_buf[i][0],
                                         slices_buf[i][1],
                                         slices_buf[i][2],
                                         slices_buf[i][3]].copy()
                elif (set_boundaries):
                    data_wbuf_space[slices_data_wbuf[i][0],
                                    slices_data_wbuf[i][1],
                                    slices_data_wbuf[i][2],
                                    slices_data_wbuf[i][3]] = boundary_val
                    data_count_wbuf_space[slices_data_wbuf[i][0],
                                          slices_data_wbuf[i][1],
                                          slices_data_wbuf[i][2],
                                          slices_data_wbuf[i][3]] = count_thresh

            # Do gapfilling
            opmask = should_exist - exist_after
            for c in range(nClasses):

                rt = 0
                st = 0
                rs = bufs
                ss = 1
                data[:,c] = gapfill(data[:,c],opmask[:,c],data_wbuf_space[:,c],data_count_wbuf_space[:,c],rt,rs,st,ss)


            # Assign defaults to gaps that weren't filled
            print('assigning defaults')
            exist_after = np.where(~np.isnan(data), 1, 0)
            opmask = should_exist - exist_after
            for c in range(nClasses):
                a = np.where(opmask[0,c] == 1)[0]
                b = np.where(opmask[0,c] == 1)[1]
                for i,j in zip(a,b):
                    data[:,c,i,j] = defaults[varname][type] 
            ds[tmpvar].encoding['zlib'] = True

            # Update gapfill_flag
            exist_after = np.where(~np.isnan(data), 1, 0)
            gapfill_flag += exist_after - exist_before

            # Save gapfill_flag
            tmpvar_gap = tmpvar + '_gapfill_flag'
            ds[tmpvar_gap] = ((timevarname,'veg_class','lat','lon'),gapfill_flag)
            ds[tmpvar_gap].encoding['zlib'] = True


    # Write to output file
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

