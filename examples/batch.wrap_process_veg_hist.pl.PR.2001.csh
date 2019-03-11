# To run this script, replace $PROJECT with the path to your clone of this GitHub repo; and $AGGROOT with the location of the top-level directory for processing aggregated MODIS observations
wrap_process_veg_hist.pl $AGGROOT/NLCD_INEGI_MODIS $AGGROOT/NLCD_INEGI_MODIS 2011 veg_hist.30 0,1,2,3,4,5,6,7 null 2011 3600 20 1 >& log.wrap_process_veg_hist.pl.USMX.NLCD_INEGI.2011.30_40.0-7.txt
!!! make note in procedure.md that we include the .30 in prefix to make sure we only process the 30_40n files
!!! point people to L2015 params - they'll need that before stage 7 and 9

!!! Do we do this here? Can we include as another step in the process* script?

set_run_cell.py -p $AGGROOT/NLCD_INEGI_MODIS/2011/vic_params.allyears/params.USMX.NLCD_INEGI.2011.2000_2016.10_20n.-70_-60e.nc -r /home/tjbohn/data/PR/domain/domain.PR.10_20n.-70_-60e.nc -v mask -o /media/tjbohn/BigData/data/VegHist/NLCD_INEGI_MODIS/2001.PR/vic_params.allyears/params.PR.S2006.2001.PR.2000_2016.10_20n.-70_-60e.nc.tmp
mv /media/tjbohn/BigData/data/VegHist/NLCD_INEGI_MODIS/2001.PR/vic_params.allyears/params.PR.S2006.2001.PR.2000_2016.10_20n.-70_-60e.nc.tmp /media/tjbohn/BigData/data/VegHist/NLCD_INEGI_MODIS/2001.PR/vic_params.allyears/params.PR.S2006.2001.PR.2000_2016.10_20n.-70_-60e.nc

wrap_process_veg_hist.pl $AGGROOT/NLCD_INEGI_MODIS $AGGROOT/NLCD_INEGI_MODIS 2011 veg_hist.30 8,9 null 2011 3600 10 1 >& log.wrap_process_veg_hist.pl.USMX.NLCD_INEGI.2011.30_40.8-9.txt

set_run_cell.py -p /media/tjbohn/BigData/data/VegHist/NLCD_INEGI_MODIS/2001.PR/vic_params.allyears/params.PR.S2006.2001.PR.2001_2001.10_20n.-70_-60e.nc -r /home/tjbohn/data/PR/domain/domain.PR.10_20n.-70_-60e.nc -v mask -o /media/tjbohn/BigData/data/VegHist/NLCD_INEGI_MODIS/2001.PR/vic_params.allyears/params.PR.S2006.2001.PR.2001_2001.10_20n.-70_-60e.nc.tmp
mv /media/tjbohn/BigData/data/VegHist/NLCD_INEGI_MODIS/2001.PR/vic_params.allyears/params.PR.S2006.2001.PR.2001_2001.10_20n.-70_-60e.nc.tmp /media/tjbohn/BigData/data/VegHist/NLCD_INEGI_MODIS/2001.PR/vic_params.allyears/params.PR.S2006.2001.PR.2001_2001.10_20n.-70_-60e.nc

