#! /bin/csh -f

 set dir           = `grep IM_DIR process_config|awk '{print $2}'`
 set dir_mrj       = `grep MRJ_DIR process_config|awk '{print $2}'`
 set ref_file      = `grep REF_SUB process_config|awk '{print $2}'`
 set date_file     = `grep INFILE  process_config|awk '{print $2}'`
 set phot_file     = `grep VARIABLES  process_config|awk '{print $2}'`
 set dir_config    = `grep CONFIG_DIR process_config|awk '{print $2}'`

 
  cd $dir

  set liste = `awk '{print $1}' $date_file`

 set n=$#liste
 set i=1

 echo > temp

 while($i <= $n)
 
 set s = ` $dir_mrj"/bin/phot_ref" -i $liste[$i]  -c   $dir_config"/phot_config" | awk '{print $2,$4,$6,$8}' `

 echo $liste[$i] $s >> temp
 echo $liste[$i] $s

  @ i += 1
 end

  sort +1n temp | head -20 | awk '{print $1}' > $dir_config"/ref_list"




 
