#!/usr/bin/perl

$rootdir = shift;
$archivedir = shift;
$lcid = shift;
$file = shift;
$stagelist = shift;
$lcid_cv = shift;
$lcid_out = shift;
$clean = shift; # 1 = delete files from previous steps as we go
$data_root = shift;
$domain_pfx = shift;
$param_pfx = shift;
$lcscheme = shift;

@stages = split /,/, $stagelist;
foreach $stage (@stages) {
  $do_stage{$stage} = 1;
}

# Fill class = open shrubland
# have to give index of class, not classID (for NLCD, they're not equal)
if ($lcid_out =~ /(mode|IGBP|MODIS)/) {
  $fill_class = 7;
}
else {
  $fill_class = 10;
}

@vars = ("LAI","NDVI","fcanopy","albedo");
@vars_for_vic = ("LAI","fcanopy","albedo");
$varnamelist = join ",", @vars;
$varnamelist_for_vic = join ",", @vars_for_vic;

if ($lcid_out eq "" || $lcid_out eq "null") {
  die "$0: ERROR: lcid_out must be defined\n";
}


if ($file =~ /\.([\d\.-]+)_([\d\.-]+)n\./) {
  ($minlat,$maxlat) = ($1,$2);
}

# Extra filtering of bad albedo values that are difficult to do in aggregation
if ($do_stage{0}) {

  $subdir_in = "aggregated";
  $indir = "$rootdir/$lcid/$subdir_in";
  $subdir_out = "aggregated.filter_albedo";
  $outdir = "$archivedir/$lcid/$subdir_out";
  `mkdir -p $outdir`;
  $infile = "$indir/$file";
  $outfile = "$outdir/$file";
  $cmd = "cleanup_albedo.py -i $infile -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


# Replace the Cv with that of the specified full file
if ($do_stage{101}) {

  if ($lcid_cv eq "null" || $lcid_cv eq "") { exit; }

  $subdir_in = "aggregated.filter_albedo";
  $indir = "$archivedir/$lcid/$subdir_in";
  $cvdir = "$archivedir/$lcid_cv/$subdir_in";
  $subdir_out = "$subdir_in";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;
  $infile = "$indir/$file";
  $cvfile = "$cvdir/$file";
  $cvfile =~ s/\.$lcid\./.$lcid_cv./;
  $outfile = "$outdir/$file";
  $outfile =~ s/\.$lcid\./.$lcid_out./;
  $cmd = "replace_cv.py -i $infile -c $cvfile -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


# Compute climatology of 8-day observations
if ($do_stage{1}) {

  $subdir_in = "aggregated.filter_albedo";
  $indir = "$archivedir/$lcid_out/$subdir_in";
  $subdir_out = "climatology";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;
  $infile = "$indir/$file";
  $infile =~ s/\.$lcid\./.$lcid_out./;
  $outfile = "$outdir/$file";
  $outfile =~ s/\.$lcid\./.$lcid_out./;
  $outfile =~ s/nc$/clim.nc/;
  $cmd = "compute_clim_veg_hist.py -i $infile -v $varnamelist -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


# Fill gaps in climatology via gaussian kernel
if ($do_stage{2}) {

  $subdir_in = "climatology";
  $indir = "$archivedir/$lcid_out/$subdir_in";
  $subdir_out = "clim_gapfill";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;
  $infile = "$indir/$file";
  $infile =~ s/\.$lcid\./.$lcid_out./;
  $infile =~ s/nc$/clim.nc/;
  $outfile = "$outdir/$file";
  $outfile =~ s/\.$lcid\./.$lcid_out./;
  $outfile =~ s/nc$/clim_gapfill.nc/;
  $cmd = "gapfill_veg_hist.py -i $infile -t clim -v $varnamelist -c 5 -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

  if ($clean) {
    $cmd = "rm $infile";
    print "$cmd\n";
    (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
  }

}


# Compute anomalies (z-transform) of 8-day observations relative to climatology
if ($do_stage{3}) {

  $subdir_in = "aggregated.filter_albedo";
  $indir = "$archivedir/$lcid_out/$subdir_in";
  $subdir_out = "anomaly";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;
  $infile = "$indir/$file";
  $infile =~ s/\.$lcid\./.$lcid_out./;
  $outfile = "$outdir/$file";
  $outfile =~ s/\.$lcid\./.$lcid_out./;
  $outfile =~ s/nc$/anom.nc/;
  $outfile1 = "$outdir/tmp.nc";
  $subdir_clim = "clim_gapfill";
  $climdir = "$archivedir/$lcid_out/$subdir_clim";
  $climfile = "$climdir/$file";
  $climfile =~ s/\.$lcid\./.$lcid_out./;
  $climfile =~ s/nc$/clim_gapfill.nc/;
  $cmd = "compute_anom_veg_hist.py -i $infile -c $climfile -v $varnamelist -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


# Gapfill anomalies via gaussian kernel
if ($do_stage{4}) {

  $subdir_in = "anomaly";
  $indir = "$archivedir/$lcid_out/$subdir_in";
  $subdir_out = "anom_gapfill";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;
  $infile = "$indir/$file";
  $infile =~ s/\.$lcid\./.$lcid_out./;
  $infile =~ s/nc$/anom.nc/;
  $outfile = "$outdir/$file";
  $outfile =~ s/\.$lcid\./.$lcid_out./;
  $outfile =~ s/nc$/anom_gapfill.nc/;
  $cmd = "gapfill_veg_hist.py -i $infile -t anom -v $varnamelist -c 1 -m 40 -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

  if ($clean) {
    $cmd = "rm $infile";
    print "$cmd\n";
    (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
  }

}


# Recombine climatology and anomalies
if ($do_stage{5}) {

  $subdir_in = "anom_gapfill";
  $indir = "$archivedir/$lcid_out/$subdir_in";
  $subdir_clim = "clim_gapfill";
  $climdir = "$archivedir/$lcid_out/$subdir_clim";
  $subdir_out = "recombine";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;
  if ($varname =~ /ndvi/i) { next; }
  $infile = "$indir/$file";
  $infile =~ s/\.$lcid\./.$lcid_out./;
  $infile =~ s/nc$/anom_gapfill.nc/;
  $climfile = "$climdir/$file";
  $climfile =~ s/nc$/clim_gapfill.nc/;
  $climfile =~ s/\.$lcid\./.$lcid_out./;
  $outfile = "$outdir/$file";
  $outfile =~ s/\.$lcid\./.$lcid_out./;
  $outfile =~ s/nc$/recombine.nc/;
  $cmd = "recombine_clim_anom_veg_hist.py -a $infile -c $climfile -v $varnamelist -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


# Interpolate to daily, agg to monthly, and compute monthly climatology 
if ($do_stage{6}) {

  $subdir_in = "recombine";
  $indir = "$archivedir/$lcid_out/$subdir_in";
  $subdir_out_daily = "daily";
  $outdir_daily = "$archivedir/$lcid_out/$subdir_out_daily";
  `mkdir -p $outdir_daily`;
  $subdir_out_monthly = "monthly";
  $outdir_monthly = "$archivedir/$lcid_out/$subdir_out_monthly";
  `mkdir -p $outdir_monthly`;
  $subdir_out_monthly_clim = "monthly_clim";
  $outdir_monthly_clim = "$archivedir/$lcid_out/$subdir_out_monthly_clim";
  `mkdir -p $outdir_monthly_clim`;
  $infile = "$indir/$file";
  $infile =~ s/\.$lcid\./.$lcid_out./;
  $infile =~ s/nc$/recombine.nc/;
  $outfile_monthly = "$outdir_monthly/$file";
  $outfile_monthly =~ s/\.$lcid\./.$lcid_out./;
  $outfile_monthly =~ s/nc$/monthly.nc/;
  $outfile_monthly_clim = "$outdir_monthly_clim/$file";
  $outfile_monthly_clim =~ s/\.$lcid\./.$lcid_out./;
  $outfile_monthly_clim =~ s/nc$/monthly_clim.nc/;
  $prefix = $file;
  $prefix =~ s/.nc$//;
  $prefix = "$prefix.daily";
  $cmd = "interp_and_agg2monthly_veg_hist.py -i $infile -v $varnamelist_for_vic -o $outdir_daily -p $prefix -m $outfile_monthly -c $outfile_monthly_clim";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

  if ($clean) {
    $cmd = "rm -rf $infile";
    print "$cmd\n";
    (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
  }

}

# Combine these veg params with existing soil params 
if ($do_stage{7}) {

  $domain_dir = "$data_root/$domain_pfx";
  $param_dir = "$data_root/$param_pfx";
  $subdir_clim = "monthly_clim";
  $climdir = "$archivedir/$lcid_out/$subdir_clim";

  $subdir_out = "vic_params.allyears";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;

  if ($file =~ /veg_hist\.([-\d\.]+_[-\d\.]+n\.[-\d\.]+_[-\d\.]+e)\.(\d+_\d+)\./) {
    ($locstr,$year_range) = ($1,$2);
  }
  $domainfile = "$domain_dir/$domain_pfx.$locstr.nc";
  $paramfile = "$param_dir/$param_pfx.$locstr.nc";
  $libfile = "$data_root/veg_lib.$domain.$lcscheme.txt";
  $rootfile = "$data_root/root_zones.$domain.$lcscheme.txt";
  $climfile = "$climdir/$file";
  $climfile =~ s/\.$lcid\./.$lcid_out./;
  $climfile =~ s/nc$/monthly_clim.nc/;
  $outfile = "$outdir/$param_pfx.$lcid_out.$year_range.$locstr.nc";
  $outfile =~ s/\.orig//;

  $albedo_flag = "";
  $snowband_flag = "-s";

  $cmd = "replace_vegparams_with_veghist_clim.py -d $domainfile -p $paramfile -b $libfile -r $rootfile -c $climfile $albedoflag $snowband_flag -w 0 -f $fill_class -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}

# Compute climatology based on subset of years, from monthly
if ($do_stage{8}) {

  if ($lcid_out =~ /^s?(\d\d\d\d)/) {
    $year = $1;
    if ($year < 2001) {
      $year = 2000;
    }
    $startyear = $year;
    $endyear = $year;
  }
  else {
    die "$0: ERROR: start/end years not clear for lcid $lcid_out\n";
  }
  $startyear_idx = $startyear - 2000;
  $endyear_idx = $endyear - 2000;

  $subdir_in = "monthly";
  $indir = "$archivedir/$lcid_out/$subdir_in";
  $subdir_out = "monthly_clim.$startyear\_$endyear";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;
  $infile = "$indir/$file";
  $infile =~ s/\.$lcid\./.$lcid_out./;
  $infile =~ s/nc$/monthly.nc/;
  $outfile = "$outdir/$file";
  $outfile =~ s/\.$lcid\./.$lcid_out./;
  $outfile =~ s/nc$/monthly_clim.nc/;
  $outfile =~ s/2000_2016/$startyear\_$endyear/;
  $cmd = "compute_clim_from_monthly.py -i $infile -v $varnamelist_for_vic -s $startyear_idx -e $endyear_idx -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}

# Combine these veg params with existing soil params 
if ($do_stage{9}) {

  if ($lcid_out =~ /^s?(\d\d\d\d)/) {
    $year = $1;
    if ($year < 2001) {
      $year = 2000;
    }
    $startyear = $year;
    $endyear = $year;
  }
  else {
    die "$0: ERROR: start/end years not clear for lcid $lcid_out\n";
  }

  $domain_dir = "$data_root/domain";
  $param_dir = "$data_root/$param_pfx";
  $subdir_clim = "monthly_clim.$startyear\_$endyear";
  $climdir = "$archivedir/$lcid_out/$subdir_clim";

  $subdir_out = "vic_params.$startyear\_$endyear";
  $outdir = "$archivedir/$lcid_out/$subdir_out";
  `mkdir -p $outdir`;

  if ($file =~ /veg_hist\.([-\d\.]+_[-\d\.]+n\.[-\d\.]+_[-\d\.]+e)\.(\d+_\d+)\./) {
    ($locstr,$year_range) = ($1,$2);
  }
  $year_range = "$startyear\_$endyear";
  $domainfile = "$domain_dir/$domain_pfx.$locstr.nc";
  $paramfile = "$param_dir/$param_pfx.$locstr.nc";
  $libfile = "$data_root/veg_lib_$lcscheme";
  $rootfile = "$data_root/root_zones_$lcscheme";
  $climfile = "$climdir/$file";
  $climfile =~ s/\.$lcid\./.$lcid_out./;
  $climfile =~ s/nc$/monthly_clim.nc/;
  $climfile =~ s/2000_2016/$startyear\_$endyear/;
  $outfile = "$outdir/$param_pfx.$lcid_out.$year_range.$locstr.nc";
  $outfile =~ s/\.orig//;

  $albedo_flag = "";
  $snowband_flag = "-s";

  $cmd = "replace_vegparams_with_veghist_clim.py -d $domainfile -p $paramfile -b $libfile -r $rootfile -c $climfile $albedoflag $snowband_flag -l 0 -o $outfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}

