#!/usr/local/bin/perl

if ( $#ARGV < 1 ) {
    print( "runalard.pl: template image -st saturation Temp (60000) -si saturation Image (60000) -FWHMt (FWHM Temp) -FWHMi (FWHM Image) -nsx nstampsx(4) -nsy nstampsy(4) -sx imagepiecesx(1) -sy imagepiecesy(1) -min mingoodpix (5)
-ms  minimumforstamp(130) -order spatialkernal order(2) -reverse -hms meshsize (9) -stampsbyxy (FILE) -scaletotemp (scale sub image to template scale) -removeconv \n-iterkernalsig [2] \n\n");
    die;
}

$snpath = $ENV{ 'SNPATH' };
$snsearch = $ENV{ 'SNSEARCH' };
($template,$img)=@ARGV;

$_=$template;
if (/(\S+).fits/) {
    $template=$1;
}
$_=$img;
if (/(\S+).fits/) {
    $img=$1;
}

$iterkernalsig=2.;
$removeconv=0;
$FWHMi=0.;
$FWHMt=0.;
$saturation1=60000;
$saturation2=60000;
$nsx=9;
$nsy=9;
$sx=1;
$sy=1;
$min=5;
$ms=130;
$order=2;
$reverseflag=0;
$hms=9;
$stampsbyxy=0;
$starlist="";
$scaletotemp=0;

for ($i = 2; $i <= $#ARGV; $i++) {
    if ($ARGV[$i] eq '-st') {
	$i++;
	$saturation2=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-si') {
	$i++;
	$saturation1=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-FWHMt') {
	$i++;
	$FWHMt=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-FWHMi') {
	$i++;
	$FWHMi=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-nsx') {
	$i++;
	$nsx=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-nsy') {
	$i++;
	$nsy=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-sx') {
	$i++;
	$sx=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-sy') {
	$i++;
	$sy=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-min') {
	$i++;
	$min=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-ms') {
	$i++;
	$ms=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-order') {
	$i++;
	$order=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-starlist') {
	$i++;
	$starlist="-instarfile $ARGV[$i]";
    }
    elsif ($ARGV[$i] eq '-reverse') {
	$reverseflag=1;
    }
    elsif ($ARGV[$i] eq '-scaletotemp') {
	$scaletotemp=1;
    }
    elsif ($ARGV[$i] eq '-removeconv') {
	$removeconv=1;
    }
    elsif ($ARGV[$i] eq '-stampsbyxy') {
	$i++;
	$stampsbyxy=$ARGV[$i];
	system("cp $stampsbyxy STAMPS");
	$stampsbyxy=1;
    }
    elsif ($ARGV[$i] eq '-hms') {
       $i++;
       $hms=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-iterkernalsig') {
	$i++;
	$iterkernalsig=$ARGV[$i];
    }
    else {
	die "Argument $i [$ARGV[$i] is not understood\n";
    }
    
}

if ($FWHMi ==0 && $FWHMt ==0) {
    if (! -e "$img.fits.stars") {
	system("$snsearch/runsex.pl $img.fits 10");
    }
    if (! -e "$template.fits.stars") {
	system("$snsearch/runsex.pl $template.fits 10");
    }
    $_=`$snsearch/getseeingratio.pl $img.fits.stars $template.fits.stars`;
    ($ratio)=/(\S+)/;
    printf("Image to template seeing ratio is %7.3f\n",$ratio);
    if ($ratio <1) {
	$reverse=0; #template worse than image
    }
    else {
	$reverse=1; #template better than image
    }
}
else {
    if ($FWHMi ==0) {
	$_=`$snsearch/getstar.pl $img.fits $starlist`;
	($FWHMi)=/FWHM:\s*(\S+)/;
    }
    if ($FWHMt==0) {
	$_=`$snsearch/getstar.pl $template.fits $starlist`;
	($FWHMt)=/FWHM:\s*(\S+)/;
    }
    
    print "FWHM image:         $FWHMi\nFWHM template:      $FWHMt\n";        
    
    if ($FWHMi < $FWHMt) {$reverse=0;}
    else {$reverse=1;}
}

if ($reverseflag ==1) {
    if ($reverse==1) {
	$reverse=0;
    }
    else {
	$reverse=1;
    }
}

if ($reverse==1) {
    $tmp=$saturation1;
    $saturation1=$saturation2;
    $saturation2=$tmp;
}

if ($scaletotemp) { #see main.c in ISIS for some help on this...
    $reverse +=2;
}


$hss=$hms+3;
open (F,">default_config");
print  F <<EOT;
nstamps_x         $nsx       
nstamps_y         $nsy      
sub_x             $sx       
sub_y             $sy 
half_mesh_size    $hms      
half_stamp_size   $hss     
deg_bg            1      
saturation1       $saturation1   
saturation2       $saturation2   
pix_min           $min 
min_stamp_center  $ms  
ngauss            3       
deg_gauss1        6       
deg_gauss2        4       
deg_gauss3        2       
sigma_gauss1      0.7     
sigma_gauss2      1.5     
sigma_gauss3      2.0     
deg_spatial       $order  
reverse           $reverse
stampsbyxy        $stampsbyxy
iter_kernal_sig   $iterkernalsig
EOT

close(F);

if ($reverse ==0 || $reverse == 2 ) {
    print "$snpath/alardsub $img.fits $template.fits (reverse =$reverse)\n";
    system("$snpath/alardsub $img.fits $template.fits");
    system("mv conv.fits $img.sub.fits");
    if ($removeconv) {
	system("rm conv0.fits");
    }
    else { 
	system("mv conv0.fits $img.conv.fits");
    }
}
else {
    print("$snpath/alardsub $template.fits $img.fits (reverse =$reverse)\n");
    system("$snpath/alardsub $template.fits $img.fits");
    system("mv conv.fits $img.sub.fits"); 
    if ($removeconv) {
	system("rm conv0.fits");
    }
    else {
	system("cp $img.fits $img.conv.fits");
    }
}
system("mv sum_kernel $img.sum_kernel");
unlink "kernel_coeff0" if -e "kernel_coeff0";
unlink "kernel_table" if -e "kernel_table";
unlink "toto.bmp" if -e "toto.bmp";




