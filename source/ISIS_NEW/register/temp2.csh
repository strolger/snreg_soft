#! /bin/csh -f

 
 set dir           = `grep IM_DIR process_config     | awk '{print $2}'`
 set dir_mrj       = `grep MRJ_DIR process_config    | awk '{print $2}'`
 set ref_file      = `grep REFERENCE process_config  | awk '{print $2}'`
 set date_file     = `grep INFILE  process_config    | awk '{print $2}'`
 set phot_file     = `grep VARIABLES  process_config | awk '{print $2}'`
 set dir_config    = `grep CONFIG_DIR process_config | awk '{print $2}'`


 set list = `awk '{print "interp_"$1}'  $dir_config"/ref_list"`

 cd $dir


 set i  = 1
 set nb = $#list

  



 
 


 while($i <= $nb)
 
  $dir_mrj"/bin/mrj_phot" $list[$i] t2.fits  -c $dir_config"/default_config"
  
  echo "conv_"$list[$i] "conv0_"$list[$i]

  mv conv.fits "conv2_"$list[$i]
  mv conv0.fits "conv20_"$list[$i]

  @ i += 1

 end


 $dir_mrj"/bin/stack" `awk '{print "conv20_interp_"$1}'  $dir_config"/ref_list"` -o ref.fits -c $dir_config"/default_config"
