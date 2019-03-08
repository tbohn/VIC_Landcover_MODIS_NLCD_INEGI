#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$startdate = shift; # yyyy.mm.dd
$enddate = shift; # yyyy.mm.dd
$hmin = shift;
$hmax = shift;
$vmin = shift;
$vmax = shift;
$outdir = shift;
$outpfx = shift;

`mkdir -p $outdir`;
$tmpdir = "$outdir.tmp";
`mkdir -p $outdir.tmp`;

# Make list of hv strings identifying the desired tiles
@hvlist = ();
#%hash_of_filelists = {};
for ($h=$hmin; $h<=$hmax; $h++) {
  for ($v=$vmin; $v<=$vmax; $v++) {
    $hvstr = sprintf "h%02dv%02d", $h, $v;
    push @hvlist, $hvstr;
#    @{$hash_of_filelists{$hvstr}} = ();
  }
}

# Get list of all date subdirectories
opendir(DIR,$indir) or die "$0: ERROR: cannot open directory $indir for reading\n";
@subdirs = grep /^\d/, readdir(DIR);
closedir(DIR);

# Copy the desired tiles to a single directory
foreach $subdir (sort(@subdirs)) {
  if ($subdir < $startdate || $subdir > $enddate) {
    next;
  }
  opendir(DIR,"$indir/$subdir") or die "$0: ERROR: cannot open directory $indir/$subdir for reading\n";
  @files = grep /^$prefix/, readdir(DIR);
  closedir(DIR);

  foreach $file (sort(@files)) {
    for $hvstr (@hvlist) {
      if ($file =~ /$hvstr/) {
        $cmd = "cp $indir/$subdir/$file $tmpdir/";
        (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
      }
    }
  }
}

# For each tile, assemble the list of files from all acqusition dates
# for that tile and run find_mode_MODIS_PFT.py
for $hvstr (@hvlist) {
  # Get list of all files corresponding to this tile
  opendir(DIR,$tmpdir) or die "$0: ERROR: cannot open directory $tmpdir for reading\n";
  @files = grep /$hvstr/, grep /^$prefix/, readdir(DIR);
  closedir(DIR);

  # Run the script
  $filelist = join ",", sort(@files);
  $outfile = "$outdir/$outpfx.$hvstr.hdf";
  $cmd = "find_mode_MODIS_PFT.py -i $tmpdir -f $filelist -o $outfile";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
