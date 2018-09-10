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

  


 $dir_mrj"/bin/stack" `awk '{print "interp_"$1}'  $dir_config"/ref_list"` -o ref.fits -c $dir_config"/default_config" -m



