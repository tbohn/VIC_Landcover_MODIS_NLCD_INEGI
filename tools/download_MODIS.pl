#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$startdate = shift; # yyyy.mm.dd
$enddate = shift; # yyyy.mm.dd
$hmin = shift;
$hmax = shift;
$vmin = shift;
$vmax = shift;
$username = shift;
$password = shift;

@hvlist = ();
for ($h=$hmin; $h<=$hmax; $h++) {
  for ($v=$vmin; $v<=$vmax; $v++) {
    $hvstr = sprintf "h%02dv%02d", $h, $v;
    push @hvlist, $hvstr;
  }
}

opendir(DIR,$indir) or die "$0: ERROR: cannot open directory $indir for reading\n";
@subdirs = grep /^\d/, readdir(DIR);
closedir(DIR);

foreach $subdir (sort(@subdirs)) {
  if ($subdir < $startdate || $subdir > $enddate) {
    next;
  }
  $index_file = "$indir/$subdir/index.html";
  open(FILE,$index_file) or die "$0: ERROR: cannot open file $index_file for reading\n";
  foreach (<FILE>) {
    if (/href=\"($prefix\.A.+\.hdf)\"/) {
      $filename = $1;
      for $hvstr (@hvlist) {
        if ($filename =~ /$hvstr/) {
          $cmd = "wget -P $indir/$subdir --http-user=$username --http-password=$password https://$indir/$subdir/$filename";
          (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
        }
      }
    }
  }
  close(FILE);
}
