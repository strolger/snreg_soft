#! /bin/csh

set liste1 = `awk '{print $1}' dates`
set liste2 = `awk '{print $2}' dates`

set nb = $#liste1

set i=1
while($i <= $nb)
 set t = `echo $i | awk '{print $1/6.0-int($1/6.0)+1}'`
 if($t == 1) then 
  echo "conv_"$liste1[$i]
 endif
 @ i += 1
end
