#! /usr/local/bin/tcsh -f

 set liste = `ls ../images/*brr*.fits`
 set nb = $#liste
 set i=1

 while($i <= $nb)

  
  set day=`fold -80 $liste[$i] | grep "DATE-OBS" | awk '{print substr($2,2,2)}'`
  set month=`fold -80 $liste[$i] | grep "DATE-OBS" | awk '{print substr($2,5,2)}'`
  set year=`fold -80 $liste[$i] | grep "DATE-OBS" | awk '{print substr($2,8,2)}'`


  set tu_h = `fold -80 $liste[$i] | grep "TU-START"| awk '{print substr($2,2,2)}'`
  set tu_m = `fold -80 $liste[$i] | grep "TU-START"| awk '{print substr($2,5,2)}'`

  set tm_h = `fold -80 $liste[$i] | grep "TM-START"| awk '{print substr($2,2,2)}'`

  if($tm_h<18) then
   set day=`echo $day|awk '{print $1+1.0}'`
  endif


  set year = `echo $year | awk '{print $1+1900.0}'`

  echo $year>temp2
  echo $month>>temp2
  echo $day>>temp2
  echo $tu_h>>temp2
  echo $tu_m>>temp2

  set jd = `jul<temp2| tail -1`
  echo $liste[$i] $jd
 
   @ i += 1
 end
