#! /bin/csh -f

 set dir       = `grep IM_DIR process_config|awk '{print $2}'`
 set dir_mrj   = `grep MRJ_DIR process_config|awk '{print $2}'`
 set ref_file  = `grep REFERENCE process_config|awk '{print $2}'`
 set date_file = `grep INFILE  process_config|awk '{print $2}'`


 echo $dir $dir_mrj $ref_file $date_file


  cd $dir

 $dir_mrj"/bin/extract"  -i $ref_file   -c $dir_mrj"/register/phot_config"
 sort +2n bright.data |tail -500 > good.data

 set liste = `awk '{print $1}' $date_file |grep -v $ref_file`


 set nb    = $#liste
 set i     = 199




 


 while($i <= $nb)

  $dir_mrj"/bin/extract" -i $liste[$i]   -c $dir_mrj"/register/phot_config"
   sort +2n bright.data |tail -500 > bad.data



  $dir_mrj"/bin/hist2" -o good.data -i bad.data -q



 $dir_mrj"/bin/cross" -i bad.data -j good.data -r 7.0 -f -q
 $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
 $dir_mrj"/bin/cross" -i bad.data -j good.data -r 7.0 -f -q
 $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
 $dir_mrj"/bin/cross" -i bad.data -j good.data -r 7.0 -f -q
 $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
 $dir_mrj"/bin/cross" -i bad.data -j good.data -r 2.5 -f -q
 $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
 $dir_mrj"/bin/cross" -i bad.data -j good.data -r 2.5 -f -q
 $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
 $dir_mrj"/bin/cross" -i bad.data -j good.data -r 1.0 -f -q
 $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
 $dir_mrj"/bin/cross" -i bad.data -j good.data -r 1.0 -f -q
 $dir_mrj"/bin/fitn" -i dayfile -d 1 -q
 
 $dir_mrj"/bin/interp" -i $liste[$i]   
 mv interp.fits interp_$liste[$i]
 @ i += 1

 end
