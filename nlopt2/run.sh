#!/bin/sh
while read A B C D E F
do  
	if [ ${A:0:1} != '#' ] && [ ${A:0:1} != '%' ]
	then
		./nlopt2 $A $B $C $D $E $F
		mv dqba_$A.dat harm_$A.dat dnxz_$A.dat FODO0_$A.twiss FODO1_$A.twiss FODO2_$A.twiss Q2_$A.dat data/
		./plot.sh $A
		#mv bbeat_$A.dat quad_$A.dat data/
	fi
done << HERE
%label bareQ errQ trimQ objective_parameter_file method
%method should be one of the following LN_COBYLA LN_BOBYQA LN_NEWUOA LN_PRAXIS LN_NELDERMEAD LN_SBPLX
% 56
#a1 Q0.in Q1.in Trm_10.in Obj_a.in LN_SBPLX
#a2 Q0.in Q1.in Trm_10.in Obj_a.in LN_NELDERMEAD
#a3 Q0.in Q1.in Trm_10.in Obj_a.in LN_PRAXIS
#a4 Q0.in Q1.in Trm_10.in Obj_a.in LN_NEWUOA
#a5 Q0.in Q1.in Trm_10.in Obj_a.in LN_BOBYQA
#a6 Q0.in Q1.in Trm_10.in Obj_a.in LN_COBYLA
% 456
#b1 Q0.in Q1.in Trm_14.in Obj_b.in LN_SBPLX
#b2 Q0.in Q1.in Trm_14.in Obj_b.in LN_NELDERMEAD
#b3 Q0.in Q1.in Trm_14.in Obj_b.in LN_PRAXIS
#b4 Q0.in Q1.in Trm_14.in Obj_b.in LN_NEWUOA
#b5 Q0.in Q1.in Trm_14.in Obj_b.in LN_BOBYQA
#b6 Q0.in Q1.in Trm_14.in Obj_b.in LN_COBYLA
% 45
#c1 Q0.in Q1.in Trm_10.in Obj_c.in LN_SBPLX
#c2 Q0.in Q1.in Trm_10.in Obj_c.in LN_NELDERMEAD
#c3 Q0.in Q1.in Trm_10.in Obj_c.in LN_PRAXIS
#c4 Q0.in Q1.in Trm_10.in Obj_c.in LN_NEWUOA
#c5 Q0.in Q1.in Trm_10.in Obj_c.in LN_BOBYQA
#c6 Q0.in Q1.in Trm_10.in Obj_c.in LN_COBYLA
#c11 Q0.in Q1.in Trm_10.in Obj_c1.in LN_SBPLX
#c12 Q0.in Q1.in Trm_10.in Obj_c1.in LN_NELDERMEAD
#c13 Q0.in Q1.in Trm_10.in Obj_c1.in LN_PRAXIS
#c14 Q0.in Q1.in Trm_10.in Obj_c1.in LN_NEWUOA
#c15 Q0.in Q1.in Trm_10.in Obj_c1.in LN_BOBYQA
#c16 Q0.in Q1.in Trm_10.in Obj_c1.in LN_COBYLA
% 345
#d1 Q0.in Q1.in Trm_14.in Obj_d.in LN_SBPLX
#d2 Q0.in Q1.in Trm_14.in Obj_d.in LN_NELDERMEAD
#d3 Q0.in Q1.in Trm_14.in Obj_d.in LN_PRAXIS
#d4 Q0.in Q1.in Trm_14.in Obj_d.in LN_NEWUOA
#d5 Q0.in Q1.in Trm_14.in Obj_d.in LN_BOBYQA
#d6 Q0.in Q1.in Trm_14.in Obj_d.in LN_COBYLA
% 56, Q=5
#e1 Q0.in Q5.in Trm_10.in Obj_a.in LN_SBPLX
#e2 Q0.in Q5.in Trm_10.in Obj_a.in LN_NELDERMEAD
#e3 Q0.in Q5.in Trm_10.in Obj_a.in LN_PRAXIS
#e4 Q0.in Q5.in Trm_10.in Obj_a.in LN_NEWUOA
#e5 Q0.in Q5.in Trm_10.in Obj_a.in LN_BOBYQA
#e6 Q0.in Q5.in Trm_10.in Obj_a.in LN_COBYLA
% 456, Q=5
#f1 Q0.in Q5.in Trm_14.in Obj_b.in LN_SBPLX
#f2 Q0.in Q5.in Trm_14.in Obj_b.in LN_NELDERMEAD
#f3 Q0.in Q5.in Trm_14.in Obj_b.in LN_PRAXIS
#f4 Q0.in Q5.in Trm_14.in Obj_b.in LN_NEWUOA
#f5 Q0.in Q5.in Trm_14.in Obj_b.in LN_BOBYQA
#f6 Q0.in Q5.in Trm_14.in Obj_b.in LN_COBYLA
% 45, Q=5
#g1 Q0.in Q5.in Trm_10.in Obj_c.in LN_SBPLX
#g2 Q0.in Q5.in Trm_10.in Obj_c.in LN_NELDERMEAD
#g3 Q0.in Q5.in Trm_10.in Obj_c.in LN_PRAXIS
#g4 Q0.in Q5.in Trm_10.in Obj_c.in LN_NEWUOA
#g5 Q0.in Q5.in Trm_10.in Obj_c.in LN_BOBYQA
#g6 Q0.in Q5.in Trm_10.in Obj_c.in LN_COBYLA
% 345, Q=5
#h1 Q0.in Q5.in Trm_14.in Obj_d.in LN_SBPLX
#h2 Q0.in Q5.in Trm_14.in Obj_d.in LN_NELDERMEAD
#h3 Q0.in Q5.in Trm_14.in Obj_d.in LN_PRAXIS
#h4 Q0.in Q5.in Trm_14.in Obj_d.in LN_NEWUOA
#h5 Q0.in Q5.in Trm_14.in Obj_d.in LN_BOBYQA
#h6 Q0.in Q5.in Trm_14.in Obj_d.in LN_COBYLA
%yy Q0.in Q1.in Trm_12s.in Obj_d.in GN_DIRECT
%ba LN_SBPLX Q0b.dat Q1b.dat Trim_9.dat obj_param3.dat 
HERE
