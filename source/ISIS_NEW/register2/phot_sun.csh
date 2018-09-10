#! /bin/csh -f

 set dir       = `grep IM_DIR process_config|awk '{print $2}'`
 set dir_mrj   = `grep MRJ_DIR process_config|awk '{print $2}'`
 set ref_file  = `grep REFERENCE process_config|awk '{print $2}'`
 set date_file = `grep INFILE  process_config|awk '{print $2}'`
 set phot_file = `grep VARIABLES  process_config|awk '{print $2}'`

  cd $dir


  $dir_mrj"/bin/Bphot" -i $ref_file   -c  $dir_mrj"/register/phot_config"

  set list = `awk '{print "interp_"$1}' $date_file`
  set dates = `awk '{print $2}' $date_file`



 echo > $dir_mrj"/register/lc.data"
 set i  = 1
 set nb = $#list


 while($i <= $nb)

   $dir_mrj"/bin/mrj_phot" $ref_file $list[$i]   -c $dir_mrj"/register/default_config"
   $dir_mrj"/bin/Aphot" -i conv0.fits -j ref.data -k conv.fits   -c $dir_mrj"/register/phot_config" -o $phot_file -l $list[$i] 


  set phot = `head -1 aper.data `
  echo $list[$i] $dates[$i] $phot>> $dir_mrj"/register/lc.data"
  @ i += 1

 end
