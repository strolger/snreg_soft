#! /bin/csh -f

 set dir        = `grep IM_DIR process_config    | awk '{print $2}'`
 set dir_mrj    = `grep MRJ_DIR process_config   | awk '{print $2}'`
 set ref_file   = `grep REFERENCE process_config | awk '{print $2}'`
 set date_file  = `grep INFILE  process_config   | awk '{print $2}'`
 set deg_reg    = `grep DEGREE  process_config   | awk '{print $2}'`
 set dir_config = `grep CONFIG_DIR process_config| awk '{print $2}'`
 set cthresh    = `grep COSMIC_THRESH process_config| awk '{print $2}'`



 echo $dir $dir_mrj $ref_file $date_file

  cd $dir

 $dir_mrj"/bin/extract"  -i $ref_file -s $cthresh  -c $dir_config"/phot_config"
 mv bright.data  good.data

 set liste = `awk '{print $1}' $date_file`

 #set liste = `echo file1.fits`
 set nb    = $#liste
 set i     = 1





 while($i <= $nb)


  $dir_mrj"/bin/extract" -i $liste[$i]  -s $cthresh -c $dir_config"/phot_config"
  mv bright.data  bad.data
  
  $dir_mrj"/bin/hist2" -o good.data -i bad.data -q

  set rad = `cat bin.data|awk '{print $1}'`

  echo rad: $rad

 set j=1
 set sigma = 1
 while($j < 30 && $sigma)
  $dir_mrj"/bin/cross" -i bad.data -j good.data -r $rad -f -q -s $rad
 echo pouet $j

$dir_mrj"/bin/fitn" -i dayfile -d 1 -q|tail -1| awk '{if(($2+$4)/2.0<0.4) print 0; else print 1}'


  set sigma = `$dir_mrj"/bin/fitn" -i dayfile -d 1 -q|tail -1| awk '{if(($2+$4)/2.0<0.4) print 0; else print 1}'`
 echo pouet1 
 echo sigma $sigma
 echo pouet2
  @ j += 1
 end




  set j=1
  while($j < 4)
   $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
   $dir_mrj"/bin/cross" -i bad.data -j good.data -r 1.0 -f -q
   @ j += 1
  end

 

 set deg = 1

  while($deg < $deg_reg) 
   @ deg += 1
   $dir_mrj"/bin/fitn" -i dayfile -d $deg -q
   $dir_mrj"/bin/cross" -i bad.data -j good.data -r 1.0 -f -q
   echo deg: $deg
  end

  set f = `$dir_mrj"/bin/fitn" -i dayfile -d $deg -q| tail -1`
  echo $f $liste[$i] >> log_interp2
  
 $dir_mrj"/bin/fitn" -i dayfile -d $deg_reg -q
 $dir_mrj"/bin/interp" -i  temp.fits
 echo  interp_$liste[$i]
 mv interp.fits interp_$liste[$i]
 @ i += 1

 end
