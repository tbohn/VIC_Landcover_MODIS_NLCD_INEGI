#!/usr/bin/perl

$config = shift; # config file defining some user-specific terms for data download; can be 'null' if skipping the download stage
$lctype = shift; # 'modis' or 'asc'
$lcdir = shift;
$lcpfx = shift;
$lcid = shift;
$lc_table = shift;
$startyear = shift;
$endyear = shift;
$vmin = shift;
$vmax = shift;
$hmin = shift;
$hmax = shift;
$minlat = shift;
$maxlat = shift;
$minlon = shift;
$maxlon = shift;
$resolution = shift;
$landmask = shift;
$domain = shift;
$outdir = shift;
$outpfx = shift;
$force = shift;
$stagelist = shift; # comma-separated

@stages = split /,/, $stagelist;
foreach $stage (@stages) {
  $do_stage{$stage} = 1;
}

$lcmaplist_file = "lcmaplist.$lctype.$lcid.$minlat.$maxlat.$minlon.$maxlon.$startyear-$endyear.txt";

$filelist_file = "filelist.$lctype.$lcid.$minlat.$maxlat.$minlon.$maxlon.$startyear-$endyear.txt";

`mkdir -p $outdir`;

if ($lctype eq "modis") {

  # MODIS tiles
  opendir(DIR,$lcdir) or die "$0: ERROR: cannot open dir $lcdir for reading\n";
  @lcfiles = grep /^$lcpfx/, readdir(DIR);
  closedir(DIR);
  $j=0;
  for ($v=$vmin; $v<=$vmax; $v++) {
    for ($h=$hmin; $h<=$hmax; $h++) {
      $hvstr = sprintf "h%02dv%02d", $h, $v;
      $found = 0;
      foreach $file (@lcfiles) {
        if ($file =~ /$hvstr/) {
          push @lcmapfiles, "$lcdir/$file";
          $found = 1;
        }
      }
      if (!$found) {
        push @lcmapfiles, "null"
      }
      $j_of_v_and_h{$v}{$h} = $j;
      $j++;
    }
  }
  $nTiles = $j;
  $hv_boundlist = "$hmin,$hmax,$vmin,$vmax";
  $lcmaplist = join ",", @lcmapfiles;

}
elsif ($lctype eq "asc") {

  # Ascii grid files
  opendir(DIR,$lcdir) or die "$0: ERROR: cannot open dir $lcdir for reading\n";
  @lcfiles = grep /^$lcpfx/, readdir(DIR);
  closedir(DIR);
  foreach $file (@lcfiles) {
    if ($file =~ /([\d-]+\.?\d*)_([\d\.-]+)n.([\d\.-]+)_([\d\.-]+)(e|w)/) {
      ($lat0,$lat1,$lon0,$lon1,$ew) = ($1,$2,$3,$4,$5);
      if ($ew eq "w") {
        $lon0 *= -1;
        $lon1 *= -1;
      }
      if ($lat0 >= $minlat && $lat1 <= $maxlat && $lat0 < $lat1 && $lon0 >= $minlon && $lon1 <= $maxlon && $lon0 < $lon1) {
        push @lcmapfiles, "$lcdir/$file";
      }
    }
  }
  $lcmaplist = join ",", @lcmapfiles;

  # MODIS tiles
  $j=0;
  for ($v=$vmin; $v<=$vmax; $v++) {
    for ($h=$hmin; $h<=$hmax; $h++) {
      $j_of_v_and_h{$v}{$h} = $j;
      $j++;
    }
  }
  $nTiles = $j;
  $hv_boundlist = "$hmin,$hmax,$vmin,$vmax";

}
$cmd = "echo $lcmaplist > $lcmaplist_file";
(system($cmd)==0) or die "$0: ERROR: $cmd failed\n";


if ($do_stage{1}) {

  # read user-specific USGS access info from config file
  open(FILE, $config) or die "$0: ERROR: cannot open file $config for reading\n";
  foreach (<FILE>) {
    chomp;
    @fields = split /\s+/;
    if ($fields[0] =~ /USERNAME/i) {
      $username = $fields[1];
    }
    elsif ($fields[0] =~ /PASSWORD/i) {
      $password = $fields[1];
    }
    elsif ($fields[0] =~ /MODISROOT/i) {
      $modisroot = $fields[1];
    }
  }
  close(FILE);

  # Define MODIS directories etc
  $rootdir_MODIS = "$modisroot/MODIS";
  $url_MODIS = "e4ftl01.cr.usgs.gov";
  $collection = "006";
  $mission_LAI1 = "MOLT";
  $prefix_LAI1 = "MOD15A2H";
  $subdir_LAI1 = "$mission_LAI1/$prefix_LAI1.$collection";
  $mission_LAI2 = "MOTA";
  $prefix_LAI2 = "MCD15A2H";
  $subdir_LAI2 = "$mission_LAI2/$prefix_LAI2.$collection";
  $mission_NDVI = "MOLT";
  $prefix_NDVI = "MOD13A1";
  $subdir_NDVI = "$mission_NDVI/$prefix_NDVI.$collection";
  $mission_albedo = "MOTA";
  $prefix_albedo = "MCD43A3";
  $subdir_albedo = "$mission_albedo/$prefix_albedo.$collection";

  # Get listing of subdirs in both LAI dirs
  $dir_LAI1 = "$rootdir_MODIS/LAI/$url_MODIS/$subdir_LAI1";
  opendir(DIR,$dir_LAI1) or die "$0: ERROR: cannot open dir $dir_LAI1 for reading\n";
  @subdirs1 = grep /^\d\d\d\d\.\d\d\.\d\d$/, readdir(DIR);
  closedir(DIR);
  $dir_LAI2 = "$rootdir_MODIS/LAI/$url_MODIS/$subdir_LAI2";
  opendir(DIR,$dir_LAI2) or die "$0: ERROR: cannot open dir $dir_LAI2 for reading\n";
  @subdirs2 = grep /^\d\d\d\d\.\d\d\.\d\d$/, readdir(DIR);
  closedir(DIR);
  @subdirs = (@subdirs1,@subdirs2);
#  $first_dir_LAI2 = $subdirs2[0];

  # Fill gaps in date list
  foreach $subdir (sort(@subdirs)) {
    ($year,$month,$day) = split /\./, $subdir;
    $monthday = "$month.$day";
    if ($year % 4 == 0) {
      $valid_monthday_combos_leap{$monthday} = 1;
    }
    else {
      $valid_monthday_combos{$monthday} = 1;
    }
  }
  @subdirs_new = ();
  for ($year=$startyear; $year<=$endyear; $year++) {
    if ($year % 4 == 0) {
      foreach $monthday (sort(keys(%valid_monthday_combos_leap))) {
        push @subdirs_new, "$year.$monthday";
      }
    }
    else {
      foreach $monthday (sort(keys(%valid_monthday_combos))) {
        push @subdirs_new, "$year.$monthday";
      }
    }
  }
  @subdirs = @subdirs_new;

  # Loop over LAI subdirs
  $use_subdir_LAI2 = 0;
  $count = 0;
  $first_date = 1;
  $year_save = -1;
  foreach $subdir (sort(@subdirs)) {
    ($year,$month,$day) = split /\./, $subdir;
    if ($year >= $startyear && $year <= $endyear) {
      $jday = &compute_jday($year,$month,$day);
      $jday = sprintf "%03d", $jday;
      # Set LAI dir to first or second based on date
#      if ($subdir eq $subdir_LAI2) {
      if ($year == 2002 && $month eq "07" && $day eq "04") {
        $use_subdir_LAI2 = 1;
      }
      if ($use_subdir_LAI2) {
        $subdir_LAI = $subdir_LAI2;
        $dir_LAI = $dir_LAI2;
        $prefix_LAI = $prefix_LAI2;
        $mission_LAI = $mission_LAI2;
      }
      else {
        $subdir_LAI = $subdir_LAI1;
        $dir_LAI = $dir_LAI1;
        $prefix_LAI = $prefix_LAI1;
        $mission_LAI = $mission_LAI1;
      }
      $dir_NDVI = "$rootdir_MODIS/NDVI/$url_MODIS/$subdir_NDVI";
      $dir_albedo = "$rootdir_MODIS/albedo/$url_MODIS/$subdir_albedo";
      @dirs = ($dir_LAI,$dir_NDVI,$dir_albedo);
      @pfxs = ($prefix_LAI,$prefix_NDVI,$prefix_albedo);
      @missions = ($mission_LAI,$mission_NDVI,$mission_albedo);

      # Download files
      for ($i=0; $i<3; $i++) {
        $index = "$dirs[$i]/$subdir/index.html";
        for ($j=0; $j<$nTiles; $j++) {
          $files[$i][$j] = "null";
        }
        if (-e $index) {
          open (FILE,$index) or die "$0: ERROR: cannot open file $index for reading\n";
          foreach (<FILE>) {
            if (/\"($pfxs[$i].A$year\d\d\d\.h\d\dv\d\d\.$collection\.\d+\.hdf)\"/) {
              $file = $1;
              if ($file =~ /$pfxs[$i].A$year(\d\d\d)\.h(\d\d)v(\d\d)\.$collection\.\d+\.hdf/) {
                ($jday,$h,$v) = ($1,$2,$3);
                $v_int = $v*1;
                $h_int = $h*1;
                if ($v_int >= $vmin && $v_int <= $vmax && $h_int >= $hmin && $h_int <= $hmax) {
                  $j = $j_of_v_and_h{$v_int}{$h_int};
                  print "downloading $file\n";
                  $outfile = "$dirs[$i]/$subdir/$file"; # this will be created by the wget command
                  if ($force && -e $outfile) {
                    $cmd = "rm $outfile";
                    print "$cmd\n";
                    system($cmd);
                  }
                  if (!-e $outfile) {
                    $cmd = "wget -nv -r -P $dirs[$i]/$subdir --http-user=$username --http-password=$password http://$url_MODIS/$missions[$i]/$pfxs[$i].$collection/$subdir/$file";
                    print "$cmd\n";
                    system($cmd);
                  }
                  if (-e $outfile) {
                    $files[$i][$j] = $outfile;
                  }
                  else {
                    print "$0: WARNING: attempted to download $file but $outfile does not exist\n";
                  }
                }
              }
            }
          }
          close(FILE);
        }
      }

      $filelist_LAI = join ",", @{$files[0]};
      $filelist_NDVI = join ",", @{$files[1]};
      $filelist_albedo = join ",", @{$files[2]};

      # Append date and filelists to input file for python processing script
      if ($first_date) {
        $redirect = ">";
        $first_date = 0;
      }
      else {
        $redirect = ">>";
      }
      $info = "$year$jday $filelist_LAI $filelist_NDVI $filelist_albedo";
      $cmd = "echo \"$info\" $redirect $filelist_file";
      print "$cmd\n";
      (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

      $year_save = $year;

#$count++;
#if($count == 1) {
#last;
#}
    }
  }

}


if ($do_stage{2}) {

  # Call python script to read all appropriate subdirs for this date and agg
#  $lon0 = -($minlon);
#  $lon1 = -($maxlon);
#  $outfile = "$outdir/$outpfx.$minlat" . "_" . "$maxlat". "n.$lon0" ."_$lon1" . "w.$startyear\_$endyear";
  $outfile = "$outdir/$outpfx.$minlat" . "_" . "$maxlat". "n.$minlon" ."_$maxlon" . "e.$startyear\_$endyear";
  $logfile = "log.join_and_agg_MODIS_over_landcover.py.$lcid.$minlat.$maxlat.$minlon.$maxlon.$startyear-$endyear.txt";
  $cmd = "join_and_agg_MODIS_over_landcover.py -c $lctype -i $lcid -l $lcmaplist_file -t $lc_table -m $landmask -d $domain -f $filelist_file -b $hv_boundlist -s $minlat -n $maxlat -w $minlon -e $maxlon -r $resolution -o $outfile > $logfile 2> $logfile";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


sub compute_jday {
  my $year = shift @_;
  my $month = shift @_;
  my $day = shift @_;
  my @month_days = (31,28,31,30,31,30,31,31,30,31,30,31);
  my $days_in_month;
  my $jday;
  $jday = $day;
  $month--;
  while ($month >= 1) {
    $days_in_month = $month_days[$month-1];
    if ($year % 4 == 0 && $month == 2) {
      $days_in_month++;
    }
    $jday += $days_in_month;
    $month--;
  }
  return $jday;
}
