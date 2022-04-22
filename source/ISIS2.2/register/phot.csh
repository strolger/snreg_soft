#! /bin/csh -f

 set dir           = `grep IM_DIR process_config|awk '{print $2}'`
 set dir_mrj       = `grep MRJ_DIR process_config|awk '{print $2}'`
 set ref_file      = `grep REF_SUB process_config|awk '{print $2}'`
 set date_file     = `grep INFILE  process_config|awk '{print $2}'`
 set phot_file     = `grep VARIABLES  process_config|awk '{print $2}'`
 set dir_config    = `grep CONFIG_DIR process_config|awk '{print $2}'`

  cd $dir
 
  

  $dir_mrj"/bin/Bphot" -i $ref_file   -c  $dir_config"/phot_config"

  set list  = `awk '{print $1}' $date_file`
  set dates = `awk '{print $2}' $date_file`



 echo > $dir_mrj"/register/lc.data"
 set i  = 1
 set nb = $#list


 while($i <= $nb)

   $dir_mrj"/bin/Cphot" -i "conv_"$list[$i] -j "interp_"$list[$i]   -o  $dir_config"/phot.data" -e $dates[$i] -c $dir_config"/default_config" -d  $dir_config"/phot_config"

  @ i += 1

 end
