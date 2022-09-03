sp=bomvos
for reg in north south; do
  for res in 30m; do
    ../scripts/fit_stand_age.R --species $sp --region $reg --resolution $res &>_${reg}.${sp}.${res}.log
  done
done
