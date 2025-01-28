#!/bin/bash

chrmincovs=$1
chrcovgenomes=$2
mincov=$3
covgenome=$4
chrs=$5




chrmincovs=( $(echo ${chrmincovs} | sed 's/ /\n/g') )
chrcovgenomes=( $(echo ${chrcovgenomes} | sed 's/ /\n/g') )
chrs=( $(echo ${chrs} | sed 's/ /\n/g') )
mkdir -p $(dirname ${mincov})
rm ${mincov} ${covgenome} 2> ~/null


echo "### mincov"; >&2 echo "### mincov"
for c in ${!chrs[@]}
do
	echo ${chrs[$c]}
	zcat ${chrmincovs[$c]} >> ${mincov%.gz}
	#rm ${chrmincovs[$c]} 2> ~/null
done
gzip ${mincov%.gz}
echo "### covgenome"; >&2 echo "### covgenome"
for c in ${!chrs[@]}
do
	echo ${chrs[$c]}
	zcat ${chrcovgenomes[$c]} >> ${covgenome%.gz}
	#rm ${chrcovgenomes[$c]} 2> ~/null
done
gzip ${covgenome%.gz}













