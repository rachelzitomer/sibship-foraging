sp=bomvos
for reg in south midnorth north; do
  for res in 30m 60m 90m; do
    ./fit_stand_age.R --species $sp --region $reg --resolution $res &>_${reg}.${sp}.${res}.log
  done
done
