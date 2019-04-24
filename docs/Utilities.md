# Utility Scripts for the VIC_Landcover_MODIS_NLCD_INEGI Project

Some useful utility scripts:
 - To create a domain file for a new domain:

   `create_domain_file_from_asc_inputs.py -m $MASK -f $FRACMASK -d $DEM -p $PRCP_PARAMDIR -o $OUTFILE`

   where

   `$MASK` = ESRI ascii-format mask delineating the domain (1 = in domain; 0 = outside)  
   `$FRACMASK` = ESRI ascii-format grid showing the fraction of each cell's area that lies within the domain  
   `$DEM` = ESRI ascii-format grid of elevation values (Digital Elevation Model)  
   `$PRCP_PARAMDIR` = directory containing ESRI ascii-format grid files of precipitation parameters (2 parameters, "duration" and "peak_time", with 12 monthly mean values each, one file per month)  
   `$OUTFILE` = output path/filename

 - To create initial state files for MetSim:

   `create_metsim_state_file_from_domain_file.py -i $INFILE -s $STATEDATE -o $OUTFILE`

   where

   `$INFILE` = domain file  
   `$STATEDATE` = date of 1 day prior to the desired start; format `$YYYY-$MM-$SDD`  
   `$OUTFILE` = output MetSim state file  

 - To clip out a spatial subset (lat/lon box) of a netcdf file:

   `gridclip.py -i $INFILE -s $SOUTH -n $NORTH -w $WEST -e $EAST -o $OUTFILE`

   where

   `$INFILE` = input file  
   `$SOUTH`,`$NORTH`,`$WEST`,`$EAST` = geographic boundaries of subset  
   `$OUTFILE` = output file  

 - To subsample VIC parameters to a finer resolution AND/OR clip to a new irregularly-shaped domain (via the mask of a domain file):

   `subsample_vic_params_over_domain.py -i $INFILE -d $DOMAINFILE [-n] [-s] [-f] -o $OUTFILE`

   where

   `$INFILE` = input file  
   `$DOMAINFILE` = domain file to clip to (at either the same or finer resolution)  
   `-n` = optional; if specified, use default fill ("nodata") values  
   `-s` = optional; if specified, `$INFILE` contains snow bands  
   `-f` = optional; if specified, `$INFILE` contains Fcanopy  
   `$OUTFILE` = output file  

 - To subsample VIC gridded daily forcings to a finer resolution AND/OR clip to a new irregularly-shaped domain (via the mask of a domain file):

   `subsample_forcings_over_domain.py -i $INFILE -d $DOMAINFILE -v $VARNAMELIST -o $OUTFILE`

   where

   `$INFILE` = input file  
   `$DOMAINFILE` = domain file to clip to (at either the same or finer resolution)  
   `$VARNAMELIST` = (optional) comma-separated list of variables to process; default = "Prec,Tmin,Tmax,wind"  
   `$OUTFILE` = output file  

 - To set the "run_cell" variable equal to the mask of a forcing file (to ensure that VIC doesn't try to run where no forcings are available):

   `set_run_cell.py -p $PARAMFILE -r $RUNCELLFILE -v $VARNAME -o $OUTFILE`

   where

   `$PARAMFILE` = VIC parameter file  
   `$RUNCELLFILE` = file from which to get the values for setting "run_cell", e.g., a domain file, or maybe another VIC parameter file  
   `$VARNAME` = name of variable in `$RUNCELLFILE` the values of which will be assigned to "run_cell"  
   `$OUTFILE` = output file
