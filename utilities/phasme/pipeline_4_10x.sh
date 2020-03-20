

#!/bin/sh


tr=20
ty='10x'
#for tr in 5 10 20 50 100 ; do 
#echo 'threshold' ${tr}
#(



#for ty in 'ont' 'ccs' 'clr'; do
echo -n '' > ${ty}_${tr}_swerS.txt
echo -n '' > ${ty}_${tr}_swerG.txt
echo -n '' > ${ty}_${tr}_len.txt

for chr in $(seq 1 22)  ; do
cat results_${tr}/${ty}_${chr}_improved_90.txt| grep "smaller  "| sed 's/smaller  //'| tr '\n' , >>${ty}_${tr}_swerS.txt
cat results_${tr}/${ty}_${chr}_improved_90.txt| grep "greater  "| sed 's/greater  //'| tr '\n' , >>${ty}_${tr}_swerG.txt
cat results_${tr}/${ty}_${chr}_improved_90.txt| grep "length "| sed 's/length //'| tr '\n' , >>${ty}_${tr}_len.txt
done

echo '' >> ${ty}_${tr}_swerS.txt
echo '' >> ${ty}_${tr}_swerG.txt
echo '' >> ${ty}_${tr}_len.txt

for chr in $(seq 1 22)  ; do
cat results_${tr}/${ty}_${chr}.txt | grep "smaller  " | sed 's/smaller  //' | tr '\n' , >>${ty}_${tr}_swerS.txt
cat results_${tr}/${ty}_${chr}.txt | grep "greater  " | sed 's/greater  //' | tr '\n' , >>${ty}_${tr}_swerG.txt
cat results_${tr}/${ty}_${chr}.txt | grep "length " | sed 's/length //' | tr '\n' ,  >>${ty}_${tr}_len.txt
done
#done
#) &
#done
