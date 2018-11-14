#!/usr/bin/perl

$rootdir = shift;
$archivedir = shift;
$lcid = shift;
$prefix = shift;
$stagelist = shift;
$lcid_cv = shift;
$lcid_out = shift;
$time_offset = shift; # seconds
$nParallel = shift; # number of tiles to process in parallel
$clean = shift; # 1 = delete files from previous steps as we go
$data_root = shift;
$domain_pfx = shift;
$param_pfx = shift;
$lcscheme = shift;

@vars = ("LAI","NDVI","fcanopy","albedo");
@vars_for_vic = ("LAI","fcanopy","albedo");
$varnamelist = join ",", @vars;
$varnamelist_for_vic = join ",", @vars_for_vic;

if ($lcid_out eq "" || $lcid_out eq "null") {
  $lcid_out = "$lcid.$lcid_cv";
}

# Get list of tiles to process; assumes the "aggregated" directory exists
$subdir_in = "aggregated";
$indir = "$rootdir/$lcid/$subdir_in";
opendir(DIR,"$indir") or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
closedir(DIR);
@files = sort(@files);

$cmd = "ps -ef | grep process_veg_hist.single_file.pl | grep -v wrap | grep perl";
@info = `$cmd`;
$nRunning = @info;
foreach $file (@files) {

  while ($nRunning >= $nParallel) {
    sleep 60;
    $cmd = "ps -ef | grep process_veg_hist.single_file.pl | grep -v wrap | grep perl";
    @info = `$cmd`;
    $nRunning = @info;
  }

  $logfile = "log.process_veg_hist.single_file.pl.$file.txt";
  $cmd = "process_veg_hist.single_file.pl $rootdir $archivedir $lcid $file $stagelist $lcid_cv $lcid_out $clean $data_root $domain_pfx $param_pfx $lcscheme > $logfile 2> $logfile &";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
  sleep $time_offset;

  $cmd = "ps -ef | grep process_veg_hist.single_file.pl | grep -v wrap | grep perl";
  @info = `$cmd`;
  $nRunning = @info;

}

