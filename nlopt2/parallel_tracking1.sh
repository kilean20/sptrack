#!/bin/bash
SERVER1=apiucf1
SERVER2=apiucf2
SERVER3=apiucf3
SERVER8=apiucf8
SERVER9=apiucf9

mkdir -p data_track fig
while read label epsx epsz N0 N_particle N_TURN N_INJTURN APx APz Q0 Q1 TBT_OUT RAW_OUT SBL_OUT EBE_OUT DNU_OUT ENV_OUT 
do  
	if [ ${label:0:1} != '#' ]
	then
		todolist=$todolist"$label $epsx $epsz $N0 $N_particle $N_TURN $N_INJTURN $APx $APz $Q0 $Q1 $TBT_OUT $RAW_OUT $SBL_OUT $EBE_OUT $DNU_OUT $ENV_OUT "
#the last space is important, can't be ignored
		labellist=$labellist"$label "
	fi
done < table.in


run_core1() {
	#echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17}
	./loadQ_track $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17}
	mv tbt_$1.dat dnu_$1.dat env_$1.dat data_track/
}


export -f run_core1
parallel --joblog log --progress -n 17 --env run_core1 --workdir . -S 16/$SERVER1,16/$SERVER2,16/$SERVER3,8/$SERVER8,8/$SERVER9 --sshdelay 0.2 run_core1 ::: $todolist 
#parallel --joblog log --progress -n 17 --env run_core1 --workdir . -S 16/$SERVER1,16/$SERVER2,16/$SERVER3 --sshdelay 0.2 run_core1 ::: $todolist 
#parallel --joblog log --progress --env run_core1 -n 17 run_core1 ::: $todolist 

post_process() {
	./plot.sh emit $1
	./plot.sh histu $1
	convert -density 100 -flatten fig/emit_$1.eps fig/emit_$1.png
	convert -density 100 -flatten fig/histu_$1.eps fig/histu_$1.png
	#./movie.sh $1
}
export -f post_process
#parallel --joblog log --progress --env post_process -n 1 post_process ::: $labellist 

