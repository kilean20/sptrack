#!/bin/sh
#usage: ./runQ.sh
while read label seed nux nuz QERR scale
do  
	if [ ${label:0:1} != '#' ]
	then
		./Qerrgen $label $seed $nux $nuz $QERR $scale
	fi
done << HERE
#label seed nux nuz QERR scale
0  2 2.6 2.9 0.001 0
1  2 2.6 2.9 0.001 1
5  2 2.6 2.9 0.001 5
HERE
