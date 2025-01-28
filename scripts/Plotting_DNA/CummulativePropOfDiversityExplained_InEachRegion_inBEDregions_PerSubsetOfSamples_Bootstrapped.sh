#!/bin/bash


GenoFiles=$1
BED=$2
SamplesOrderInVCF=$3
output=$4
samples=$5

##############
# To debug / develop
##############


rm ${output} 2> ~/null
mkdir -p $(dirname ${output})
toprint=()
samplelines=()
echo ${GenoFiles}
SimOrNot=$(echo ${GenoFiles} | grep 'SLiM'  | wc -l)
echo ${SimOrNot}
if [ ${SimOrNot} -gt 0 ]
then
	samplelines=$(echo ${samples} | sed 's/ /\n/g' | awk '{print $1+4}' | shuf | awk '{str=str","$1}END{print str}' | sed 's/^,//g')
else
	samplelines=$(awk '{if(NR==FNR){a[$1]=1;next}if(a[$1]){print FNR+4}}' <(echo ${samples} | sed 's/ /\n/g') <(cat ${SamplesOrderInVCF} | sed 's/ /\n/g') | shuf | awk '{str=str","$1}END{print str}' | sed 's/^,//g')
fi
echo "Num lines and order: "${samplelines}
numsamples=$(echo ${samples} | sed 's/ /\n/g' | wc -l)
totalvar=$(bedtools intersect -a <(zcat ${GenoFiles}) -b <(zcat ${BED}) | awk -v sls=${samplelines} '{split(sls, a, ","); str=""; for(i=1; i<=length(a); i++){str=str"\t"$a[i]} print str}' | sed 's/^\t//g' | sed 's/:/\t/g' | awk '{for(i=1;i<=NF;i++){a[$i]=1}if(length(a)>1){print $0} delete(a)}' | wc -l)
toprint+=(1)
for numchr in $(seq $((${numsamples}-1)) -1 1)
do
	var=$(bedtools intersect -a <(zcat ${GenoFiles}) -b <(zcat ${BED}) | awk -v sls=${samplelines} '{split(sls, a, ","); str=""; for(i=1; i<=length(a); i++){str=str"\t"$a[i]} print str}' | sed 's/^\t//g' | cut -f-${numchr} | sed 's/:/\t/g' | awk '{for(i=1;i<=NF;i++){a[$i]=1}if(length(a)>1){print $0} delete(a)}' | wc -l)
	toprint+=($(echo ${var}}/${totalvar}  | awk '{split($1, a, "/"); print a[1]/a[2]}'))
done

echo "${toprint[*]}" | sed 's/ /\n/g' | tac | awk '{str=str" "$1}END{print str}' | sed 's/^ //g' > ${output}
echo "${toprint[*]}" | sed 's/ /\n/g' | tac | awk '{str=str" "$1}END{print str}' | sed 's/^ //g'










