#!/bin/bash

chrsExtraCallableRegions=$1
chrsNonExtraCallableRegions=$2
ExtraCallableRegions=${3%.gz}
NonExtraCallableRegions=${4%.gz}
chrs=$5

chrs=( $(echo ${chrs} | sed 's/ /\n/g') );
chrsExtraCallableRegions=( $(echo ${chrsExtraCallableRegions} | sed 's/ /\n/g') );
chrsNonExtraCallableRegions=( $(echo ${chrsNonExtraCallableRegions} | sed 's/ /\n/g') );
rm ${ExtraCallableRegions}* ${NonExtraCallableRegions}* 2> ~/null
mkdir -p $(dirname ${ExtraCallableRegions})

for c in ${!chrs[@]}
do
	echo ${chrs[${c}]}
	zcat ${chrsExtraCallableRegions[${c}]} >> ${ExtraCallableRegions}
	zcat ${chrsNonExtraCallableRegions[${c}]} >> ${NonExtraCallableRegions}
done
gzip ${ExtraCallableRegions} ${NonExtraCallableRegions}

