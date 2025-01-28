#!/bin/bash

inputTreeFiles=$1
outputEvsIRatio=$2

if [[ -s ${outputEvsIRatio} ]]
then
	rm ${outputEvsIRatio} 2> ~/.null
fi
mkdir -p $(dirname ${outputEvsIRatio})
echo "chrsize mu Ne Rep TypeOfAlleles externalmean externalsd externalmedian externaltotal internalmean internalsd internalmedian internaltotal" | sed 's/\s\+/\t/g' > ${outputEvsIRatio}
for i in ${inputTreeFiles}; do
	SimParams=$(echo ${i} | sed 's/.*SLiM_Sample_//g' | sed 's/.treefile//g' | sed 's/_/ /g');
	externalmean=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -P '^[0-9]' | cut -f2 -d':' | awk '{sum+=$1}END{sum=sum/NR; print sum}');
	externalsd=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -P '^[0-9]' | cut -f2 -d':' | awk -v mean=${externalmean} '{sum += ($1 - mean)^2} END {print sum/(NR-1)}');
	externalmedian=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -P '^[0-9]' | cut -f2 -d':' | awk '{a[i++]=$1}END{if(i%2==0){print (a[i/2]+a[i/2-1])/2}else{print a[(i-1)/2]}}');
	externaltotal=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -P '^[0-9]' | cut -f2 -d':' | awk '{sum+=$1}END{print sum}');
	internalmean=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -Pv '^[0-9]' | cut -f2 -d':' | awk '{sum+=$1}END{sum=sum/NR; print sum}');
	internalsd=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -Pv '^[0-9]' | cut -f2 -d':' | awk -v mean=${internalmean} '{sum += ($1 - mean)^2} END {print sum/(NR-1)}');
	internalmedian=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -Pv '^[0-9]' | cut -f2 -d':' | awk '{a[i++]=$1}END{if(i%2==0){print (a[i/2]+a[i/2-1])/2}else{print a[(i-1)/2]}}');
	internaltotal=$(cat ${i} | sed 's/\([0-9)]\+:[0-9].[0-9]\+\)/\n\1/g' | grep ':' | sed 's/\([0-9)]\+:[0-9].[0-9]\+\).\+/\1/g' | grep -Pv '^[0-9]' | cut -f2 -d':' | awk '{sum+=$1}END{print sum}');
	echo ${SimParams}" "${externalmean}" "${externalsd}" "${externalmedian}" "${externaltotal}" "${internalmean}" "${internalsd}" "${internalmedian}" "${internaltotal} | sed 's/\s\+/\t/g' >> ${outputEvsIRatio}
done
