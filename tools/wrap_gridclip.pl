#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$south = shift;
$north = shift;
$west = shift;
$east = shift;
$outdir = shift;
$outpfx = shift;
$verbose = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
closedir(DIR);

foreach $file (sort(@files)) {

  $outfile = $file;
  $outfile =~ s/$prefix/$outpfx/g;

  $cmd = "gridclip.py -i $indir/$file -s $south -n $north -w $west -e $east -o $outdir/$outfile.tmp";
  if ($verbose) { print "$cmd\n"; }
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

  $cmd = "make_time_unlimited_and_conv_nc4.py -i $outdir/$outfile.tmp -o $outdir/$outfile";
  if ($verbose) { print "$cmd\n"; }
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

  $cmd = "rm $outdir/$outfile.tmp";
  if ($verbose) { print "$cmd\n"; }
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


