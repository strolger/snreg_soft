#! /bin/csh -f

 set dir           = `grep IM_DIR process_config|awk '{print $2}'`
 set dir_mrj       = `grep MRJ_DIR process_config|awk '{print $2}'`
 set ref_file      = `grep REF_SUB process_config|awk '{print $2}'`
 set date_file     = `grep INFILE  process_config|awk '{print $2}'`
 set phot_file     = `grep VARIABLES  process_config|awk '{print $2}'`
 set dir_config    = `grep CONFIG_DIR process_config|awk '{print $2}'`

 set list = `awk '{print $1}' $date_file`

 cd $dir
 

 set i  = 1
 set nb = $#list

 while($i <= $nb)
 
  set f = `$dir_mrj"/bin/mrj_phot" $ref_file "interp_"$list[$i]   -c $dir_config"/default_config" | grep scatter`
  mv conv.fits "conv_"$list[$i]
  echo $f $list[$i]>>log_subtract

 @ i += 1

 end
