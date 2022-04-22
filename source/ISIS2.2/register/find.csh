#! /bin/csh -f

 set dir           = `grep IM_DIR process_config|awk '{print $2}'`
 set dir_mrj       = `grep MRJ_DIR process_config|awk '{print $2}'`
 set ref_file      = `grep REF_SUB process_config|awk '{print $2}'`
 set date_file     = `grep INFILE  process_config|awk '{print $2}'`
 set phot_file     = `grep VARIABLES  process_config|awk '{print $2}'`
 set dir_config    = `grep CONFIG_DIR process_config|awk '{print $2}'`
 set thresh        = `grep SIG_THRESH process_config|awk '{print $2}'`

  cd $dir
 
 set list = `awk '{print "conv_"$1}' $date_file`
 
 
 #$dir_mrj"/bin/abs" $list -o var.fits -c ../register/phot_config 

 $dir_mrj"/bin/detect" -i var.fits -o $dir_config"/phot.data" -r abs.fits -t $thresh -c $dir_config"/default_config"
