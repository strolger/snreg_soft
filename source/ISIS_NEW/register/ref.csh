#! /bin/csh -f

 
 set dir           = `grep IM_DIR process_config     | awk '{print $2}'`
 set dir_mrj       = `grep MRJ_DIR process_config    | awk '{print $2}'`
 set ref_file      = `grep REFERENCE process_config  | awk '{print $2}'`
 set ref_file2     = `grep REF_STACK process_config  | awk '{print $2}'`
 set date_file     = `grep INFILE  process_config    | awk '{print $2}'`
 set phot_file     = `grep VARIABLES  process_config | awk '{print $2}'`
 set dir_config    = `grep CONFIG_DIR process_config | awk '{print $2}'`


 set list = `awk '{print "interp_"$1}'  $dir_config"/ref_list"`

 cd $dir


 set i  = 1
 set nb = $#list

 
  

 while($i <= $nb)
 
 #echo  $list[$i] $ref_file 
 #exit

  $dir_mrj"/bin/mrj_phot" $list[$i] $ref_file2  -c $dir_config"/default_config"
  
  echo "conv_"$list[$i] "conv0r_"$list[$i]

  mv conv.fits  "convr_"$list[$i]
  mv conv0.fits "conv0r_"$list[$i]

  @ i += 1

 end

 $dir_mrj"/bin/stack" `awk '{print "conv0r_interp_"$1}'  $dir_config"/ref_list"` -o ref.fits -c $dir_config"/default_config"



