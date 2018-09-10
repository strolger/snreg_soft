#!/usr/local/bin/perl

@img=@ARGV;
for ($i=0; $i<=$#img; $i++){
    $_=`gethead $img[$i] naxis1 naxis2 bitpix imagetyp`;
    @c=split(/\s+/);
    $bi="unk";
    if ($c[2]=="-32"){$bit="real";}
    if ($c[2]=="16"){$bit="ushort";}
    $c[3]= lc ($c[3]);
    $filter=`gethead $img[$i] filter`;
    chomp($filter);
    $object=`gethead $img[$i] object`;
    chomp($object);
    $object=`gethead $img[$i] object`;
    chomp($object);
    $T=`gethead $img[$i] trim`;
    if ($T ne ""){$T=1;}else{$T=0;}
    $O=`gethead $img[$i] overscan`;
    if ($O ne ""){$O=1;}else{$O=0;}
    $Z=`gethead $img[$i] zerocor`;
    if ($Z ne ""){$Z=1;}else{$Z=0;}
    $F=`gethead $img[$i] flatcor`;
    if ($F ne ""){$F=1;}else{$F=0;}

    print "$img[$i]\[$c[0],$c[1]\]\[$bit\]\[$c[3]\]\[$filter\]\:$object\n";
}
