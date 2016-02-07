#!/bin/sh
while read casenum epsx epsz N0 N_particle N_TURN N_INJTURN APx APz Q0 Q1 TBT_OUT RAW_OUT SBL_OUT EBE_OUT DNU_OUT ENV_OUT 
do  
	if [ ${casenum:0:1} != '#' ]
	then
		todolist=$todolist"$casenum $epsx $epsz $N0 $N_particle $N_TURN $N_INJTURN $APx $APz $Q0 $Q1 $TBT_OUT $RAW_OUT $SBL_OUT $EBE_OUT $DNU_OUT $ENV_OUT "
#the last space is important, can't be ignored
	fi
done < $1
parallel -j 16 -n 17 ./loadQ_track -- $todolist 

