sp=bomvos
#for reg in south midnorth north; do
#  for res in 30m 60m 90m; do
#    ../scripts/fit_stand_age.R --species $sp --region $reg --resolution $res &>_${reg}.${sp}.${res}.log
#  done
#done

reg=south
res=30m
../scripts/fit_stand_age.R --species $sp --region $reg --resolution $res --block_size 1000 &>_${reg}.${sp}.${res}.log2

reg=midnorth
res=90m
../scripts/fit_stand_age.R --species $sp --region $reg --resolution $res --block_size 1000 &>_${reg}.${sp}.${res}.log2
