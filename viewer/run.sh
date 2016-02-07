#!/bin/sh
while read A B C D E F
do  
	if [ ${A:0:1} != '#' ] && [ ${A:0:1} != '%' ]
	then
		./svd_case1 $A $B $C $D $E $F
		./plot.sh $A
	fi
done << HERE
%label bareQ errQ trimQ objective_parameter_file number_of_SV_truncated
#a Q0.dat Q1.dat Trim9.dat obj_param1.dat 2
#a2 Q0.dat Q1.dat Trim12.dat obj_param1.dat 2
#aa Q0.dat Q1.dat Trim9.dat obj_param1.dat 1
#b Q0.dat Q1.dat Trim9.dat obj_param2.dat 0
%test of changing tunes by 2 quads
#c Q0.dat Q1.dat Trim2.dat obj_param3.dat 0
%iterative test -> OK!
#test1 Q0.dat Q1.dat Trim2.dat obj_test1.dat 0
#test2 Q0.dat Q2_test1.dat Trim2.dat obj_test1.dat 0
%five times bigger errors
#d Q0.dat Q_1_5.dat Trim12.dat obj_param2.dat 0
#d2 Q0.dat Q2_d.dat Trim12.dat obj_param2.dat 0
HERE
