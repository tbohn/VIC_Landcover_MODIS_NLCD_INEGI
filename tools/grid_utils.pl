#!/usr/bin/perl

sub read_header {

# This function reads the header of an ascii arc/info grid file

# Input arguments
my $filename = shift;
my $ncols;
my $nrows;
my $xllcorner;
my $yllcorner;
my $cellsize;
my $nodata;

# Open file, read header, close file
open (FILE, $filename) or die "$0: ERROR: cannot open file $filename\n";
foreach (<FILE>) {
  if (/^ncols\s+(\S+)/i) {
    $ncols = $1;
  }
  elsif (/^nrows\s+(\S+)/i) {
    $nrows = $1;
  }
  elsif (/^xllcorner\s+(\S+)/i) {
    $xllcorner = $1;
  }
  elsif (/^yllcorner\s+(\S+)/i) {
    $yllcorner = $1;
  }
  elsif (/^cellsize\s+(\S+)/i) {
    $cellsize = $1;
  }
  elsif (/^nodata(_value)?\s+(\S+)/i) {
    $nodata = $2;
  }
}
close(FILE);

# Return header info
return ($ncols,$nrows,$xllcorner,$yllcorner,$cellsize,$nodata);

}

sub latlon2ease {

# This function converts latitude and longitude to the corresponding
# row and column in the northern hemisphere EASE-grid mapping system

# Command-line arguments
my $lat = shift;  # Latitude
my $lon = shift;  # Longitude
my $C   = shift;  # Cell edge length (km) (common values are 25.067525 and 100.2701)
my $row0 = shift; # Row index of North Pole
my $col0 = shift; # Column index of North Pole
               # North pole is typically located at (360.0,360.0) (25km cell length)
               #                                 or (89.625,89.625) (100km cell length)
my $exact = shift;# 1 or 0
               # 1 = row and column numbers will be decimals
               #     (i.e. if the given location is on the border between two rows,
               #      its row will end in .5)
               # 0 = round the computed row and column to integers
               #     (i.e. what people normally think of as row and column)

# Constants
my $R = 6371.228; # Radius of Earth (km)
my $PI = 3.1415927;

# Convert deg to rad
my $phi = $lat*$PI/180;
my $lambda = $lon*$PI/180;

# Convert lat and lon to EASE-grid row and col
my $col = 2*$R/$C * sin($lambda) * sin($PI/4 - $phi/2) + $col0;
my $row = 2*$R/$C * cos($lambda) * sin($PI/4 - $phi/2) + $row0;

# Aspect ratio
my $h = cos($PI/4 - $phi/2);
my $k = 1/cos($PI/4 - $phi/2);
my $asp_ratio = $h/$k;

# Round if desired
if (!$exact) {
  $col = int($col + 0.5);
  $row = int($row + 0.5);
}

return ($row,$col);

}

sub ease2latlon {

# This function converts northern hemisphere EASE-grid row and column
# to the corresponding latitude and longitude

# Command-line arguments
my $row = shift;  # EASE-grid row
my $col = shift;  # EASE-grid column
my $C   = shift;  # Cell edge length (km) (common values are 25.067525 and 100.2701)
my $row0 = shift; # Row index of North Pole
my $col0 = shift; # Column index of North Pole
               # North pole is typically located at (360.0,360.0) (25km cell length)
               #                                 or (89.625,89.625) (100km cell length)

# Constants
my $R = 6371.228; # Radius of Earth (km)
my $PI = 3.1415927;

# Compute lambda and phi
my $tmp = 0.5*$C/$R*sqrt(($col-$col0)*($col-$col0) + ($row-$row0)*($row-$row0));
if ($tmp > 1) {
  $tmp = 1;
}
my $phi = 2*($PI/4 - atan2($tmp, sqrt(1-$tmp*$tmp)));
my $lambda = atan2(($col-$col0),($row-$row0));

# Convert rad to deg
my $lat = $phi*180/$PI;
my $lon = $lambda*180/$PI;

return ($lat,$lon);

}

sub utm2latlon {

# This function converts UTM Easting and Northing (with zone) to
# latitude and longitude (assuming WGS84 datum)

# Command-line arguments
my $x = shift;  # Easting
my $y = shift;  # Northing
my $zone   = shift;  # Longitudinal UTM zone

# Constants
my $PI = 3.1415927;

# Datum-specific constants (assuming WGS84 datum)
my $x_offset = 500000; # 500 km false Easting
my $k0 = 0.9996;
my $a = 6378137; # Equatorial radius of Earth (m)
my $b = 6356752.3142; # Polar radius of Earth (m)

# Derived constants
my $e = sqrt(1-($b*$b)/($a*$a)); # eccentricity
my $eSq = $e*$e; # eccentricity squared
my $ePmSq = ($eSq)/(1-$eSq); # e prime squared
my $e1 = (1 - sqrt(1-$eSq))/(1 + sqrt(1-$eSq));
my $J1 = (3*$e1/2 - 27*$e1*$e1*$e1/32);
my $J2 = (21*$e1*$e1/16 - 55*$e1*$e1*$e1*$e1/256);
my $J3 = (151*$e1*$e1*$e1/96);
my $J4 = (1097*$e1*$e1*$e1*$e1/512);

# Intermediate terms
my $lon0 = $zone*6 - 183;  # Central meridian of zone (deg)
my $M = $y/$k0;  # Meridional arc
my $mu = $M/($a*(1-$eSq/4 - 3*$eSq*$eSq/64 - 5*$eSq*$eSq*$eSq/256));
my $fp = $mu + $J1*sin(2*$mu) + $J2*sin(4*$mu) + $J3*sin(6*$mu) + $J4*sin(8*$mu); # footprint latitude
my $C1 = $ePmSq*cos($fp)*cos($fp);
my $tanfp = sin($fp)/cos($fp);
my $T1 = $tanfp*$tanfp;
my $tmp = sqrt(1-$eSq*sin($fp)*sin($fp));
my $R1 = $a*(1-$eSq)/($tmp*$tmp*$tmp); # Radius of curvature of Earth in meridional plane
my $N1 = $a/$tmp; # Radius of curvature of Earth perpendicular to meridional plane
my $D = ($x-$x_offset)/($N1*$k0);
my $Q1 = $N1*$tanfp/$R1;
my $Q2 = $D*$D/2;
my $Q3 = (5 + 3*$T1 + 10*$C1 - 4*$C1*$C1 - 9*$ePmSq)*$D*$D*$D*$D/24;
my $Q4 = (61 + 90*$T1 + 298*$C1 + 45*$T1*$T1 - 3*$C1*$C1 - 252*$ePmSq)*$D*$D*$D*$D*$D*$D/720;
my $Q5 = $D;
my $Q6 = (1 + 2*$T1 + $C1)*$D*$D*$D/6;
my $Q7 = (5 - 2*$C1 + 28*$T1 - 3*$C1*$C1 + 8*$ePmSq + 24*$T1*$T1)*$D*$D*$D*$D*$D/120;

# Finally, lat and lon
my $phi = $fp - $Q1*($Q2 - $Q3 + $Q4);
my $lat = $phi*180/$PI;
my $dlon = ($Q5 - $Q6 + $Q7)/cos($fp);
$dlon *=180/$PI;
my $lon = $lon0 + $dlon;

return ($lat,$lon);

}

sub latlon2utm {

# This function converts latitude and longitude (with UTM zone) to
# UTM Easting and Northing (assuming WGS84 datum)

# Command-line arguments
my $lat = shift;  # Latitude
my $lon = shift;  # Longitude
my $zone   = shift;  # Longitudinal UTM zone

# Constants
my $PI = 3.1415927;

# Datum-specific constants (assuming WGS84 datum)
my $x_offset = 500000; # 500 km false Easting
my $k0 = 0.9996;
my $a = 6378137; # Equatorial radius of Earth (m)
my $b = 6356752.3142; # Polar radius of Earth (m)

# Derived constants
my $lon0 = $zone*6 - 183;  # Central meridian of zone (deg)
#my $dlon = ($lon - $lon0)*$PI/180;  # Offset between input longitude and centeral meridian (radians)
my $dlon = ($lon - $lon0)*3600/10000;  # Offset between input longitude and centeral meridian (1e4*sec)
my $phi = $lat*$PI/180;  # input latitude in radians
my $e = sqrt(1-($b*$b)/($a*$a)); # eccentricity
my $eSq = $e*$e; # eccentricity squared
my $ePmSq = ($eSq)/(1-$eSq); # e prime squared
my $n = ($a-$b)/($a+$b);
my $n2 = $n*$n;
my $n3 = $n2*$n;
my $n4 = $n3*$n;
my $n5 = $n4*$n;

# Compute meridional arc
my $A0 = $a*(1 - $n + (5/4)*($n2-$n3) + (81/64)*($n4-$n5));
my $B0 = (3*$a*$n/2)*(1 - $n + (7/8)*($n2-$n3) + (55/64)*($n4-$n5));
my $C0 = (15*$a*$n2/16)*(1 - $n + (3/4)*($n2-$n3));
my $D0 = (35*$a*$n3/48)*(1 - $n + (11/16)*($n2-$n3));
my $E0 = (315*$a*$n4/51)*(1 - $n);
my $M = $A0*$phi - $B0*sin(2*$phi) + $C0*sin(4*$phi) - $D0*sin(6*$phi) + $E0*sin(8*$phi);

# More intermediate terms
my $tmp = sqrt(1-$eSq*sin($phi)*sin($phi));
my $R1 = $a*(1-$eSq)/($tmp*$tmp*$tmp);
my $N1 = $a/$tmp;
my $tanphi = sin($phi)/cos($phi);
my $sin1sec = sin(1/3600*$PI/180);
my $K1 = $M*$k0;
my $K2 = $k0*$sin1sec*$sin1sec*$N1*sin($phi)*cos($phi)*1e8/2;
my $K3 = ($k0*$sin1sec*$sin1sec*$sin1sec*$sin1sec*$N1*sin($phi)*cos($phi)*cos($phi)*cos($phi)/24)*(5 - $tanphi*$tanphi + 9*$ePmSq*cos($phi)*cos($phi) + 4*$ePmSq*$ePmSq*cos($phi)*cos($phi)*cos($phi)*cos($phi))*1e16;
my $K4 = $k0*$sin1sec*$N1*cos($phi)*1e4;
my $K5 = ($k0*$sin1sec*$sin1sec*$sin1sec*$N1*cos($phi)*cos($phi)*cos($phi)/6)*(1 - $tanphi*$tanphi + $ePmSq*cos($phi)*cos($phi))*1e12;

# Finally, x and y
my $y = $K1 + $K2*$dlon*$dlon + $K3*$dlon*$dlon*$dlon*$dlon;
my $x = $x_offset + $K4*$dlon + $K5*$dlon*$dlon*$dlon;

return ($x,$y);

}

sub great_circle_distance {

# This function computes the great circle distance (km) between two points, given their lat/lon coordinates

# Command-line arguments
my $lat1 = shift;  # Latitude of 1st point
my $lon1 = shift;  # Longitude of 1st point
my $lat2 = shift;  # Latitude of 2nd point
my $lon2 = shift;  # Longitude of 2nd point

# Constants
my $R = 6371.228; # Radius of Earth (km)
my $PI = 3.1415927;
my $deg2rad = $PI/180.0;

# Compute theta and phi
my $theta1 = $deg2rad*$lon1;
my $phi1 = $deg2rad*$lat1;
my $theta2 = $deg2rad*$lon2;
my $phi2 = $deg2rad*$lat2;

# Compute distance
my $term1 = cos($phi1)*cos($theta1)*cos($phi2)*cos($theta2);
my $term2 = cos($phi1)*sin($theta1)*cos($phi2)*sin($theta2);
my $term3 = sin($phi1)*sin($phi2);
my $temp = $term1+$term2+$term3;
if ($temp > 1.0) {
  $temp = 1.0;
}
my $dist = $R*atan2(sqrt(1-$temp**2), $temp);

return ($dist);

}

sub modis_sinusoidal_to_latlon {

  my $h = shift @_;
  my $v = shift @_;
  my $npix_per_deg_y = shift @_;
  my $row = shift @_;
  my $col = shift @_;

  my $PI = 3.14159265358979328462;
  my $minlat = 80 - 10*$v;
  my $maxlat = 90 - 10*$v;
  my $cellsize_y = 1/$npix_per_deg_y;
  my $npix_across = 10*$npix_per_deg_y;
  my $latcenter = $maxlat - ($row/$npix_per_deg_y + 0.5*$cellsize_y);
  my $lat_for_x;
  if ($v < 9) {
    $lat_for_x = $latcenter - 0.5*$cellsize_y;
  }
  else {
    $lat_for_x = $latcenter + 0.5*$cellsize_y;
  }
  my $ncols_map_half_width = 180*$npix_per_deg_y*cos($PI/180*$lat_for_x); # number of columns between 0 longitude and 180 longitude for this row
  my $npix_per_deg_x = $npix_per_deg_y*cos($PI/180*$lat_for_x);
  my $cellsize_x = 1/($npix_per_deg_x);

  my $minlon;
  my $loncenter;

  if ($h < 18) {
    $minlon = 10*($h-18)*$npix_per_deg_y*$cellsize_x;
    $loncenter = $minlon + ($col+0.5)*$cellsize_x;
  }
  else {
    $minlon = 10*($h-18)*$npix_per_deg_y*$cellsize_x;
    $loncenter = $minlon + ($col+0.5)*$cellsize_x;
  }

  return ($latcenter,$loncenter,$cellsize_y,$cellsize_x);

}

sub latlon_to_modis_sinusoidal {
  my $lat = shift @_;
  my $lon = shift @_;
  my $npix_per_deg_y = shift @_;

  my $PI = 3.14159265358979328462;
  my $v = int((90-$lat)/10);
  my $row = int((90-($v*10)-$lat)*$npix_per_deg_y);
  my $npix_per_deg_x = $npix_per_deg_y*cos($PI/180*$lat);
  my $cellsize_x = 1/($npix_per_deg_x);
  my $cellsize_y = 1/($npix_per_deg_y);
  my $aspect_ratio = $cellsize_x/$cellsize_y;
  my $ncols_map_half_width = int(180*$npix_per_deg_y*cos($PI/180*$lat)); # number of columns between 0 longitude and 180 longitude for this row
  my $col;
  my $h;
  if ($lon >= 0) {
    $col = int($ncols_map_half_width*$lon/180 + 0.5);
  }
  else {
    $col = -1 - int($ncols_map_half_width*abs($lon)/180 + 0.5);
  }
  $col += 180*$npix_per_deg_y;
  $h = int($col/($npix_per_deg_y*10));
  $col = $col-$h*10*$npix_per_deg_y;

  return($h,$v,$row,$col);

}

1; # This line is necessary for this file to be included in perl scripts via "require"
