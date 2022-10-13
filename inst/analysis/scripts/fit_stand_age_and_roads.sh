sp=bomvos
#for reg in south midnorth north; do
#  for res in 30m 60m 90m; do
#    ./fit_stand_age_and_roads.R --species $sp --region $reg --resolution $res &>_${reg}.${sp}.${res}.log
#  done
#done

#example
reg=midnorth
res=90m
./fit_stand_age_and_roads.R --species $sp --region $reg --resolution $res --grid_size 5 &>_${reg}.${sp}.${res}.log
