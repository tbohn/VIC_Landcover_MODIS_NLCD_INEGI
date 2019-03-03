#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$outdir = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
`mkdir -p $outdir`;
foreach $file (sort(@files)) {
  $outfile = $file;
  $outfile =~ s/monthly/hourly/g;
  if ($file =~ /(\d\d\d\d).nc/) {
    $startyear = $1;
    $endyear = $1;
  }
  $cmd = "python disagg_veghist_monthly2hourly_nc.py -i $indir/$file -s $startyear -e $endyear -o $outdir/$outfile";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
closedir(DIR);
