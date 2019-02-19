#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$maskfile = shift;
$outdir = shift;
$outpfx = shift;
$verbose = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
closedir(DIR);

foreach $file (sort(@files)) {

  $cmd = "grid_multi_clip.py -i $indir/$file -m $maskfile -p $outpfx -o $outdir";
  if ($verbose) { print "$cmd\n"; }
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


