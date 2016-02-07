#!/bin/sh
#usage: ./get_table_out.sh > table.out 
#rm -f table_out
awk '{printf "%s %s ",$1,$2 }' harm_a1.dat
tail -n1 tbt_k1_1.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k1_1.dat | avg.awk 7
tail -n1 tbt_k2_1.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k2_1.dat | avg.awk 7
tail -n1 tbt_k3_1.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k3_1.dat | avg.awk 7
tail -n1 tbt_k4_1.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k4_1.dat | avg.awk 7
echo

awk '{printf "%s %s ",$1,$2 }' harm_e1.dat
tail -n1 tbt_k1_5.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k1_5.dat | avg.awk 7
tail -n1 tbt_k2_5.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k2_5.dat | avg.awk 7
tail -n1 tbt_k3_5.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k3_5.dat | avg.awk 7
tail -n1 tbt_k4_5.dat | awk '{printf "%s ",$2 }' 
sed -n '992, 1001p' tbt_k4_5.dat | avg.awk 7
echo

#for i in a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 b6 c1 c2 c3 c4 c5 c6 d1 d2 d3 d4 d5 d6 e1 e2 e3 e4 e5 e6 f1 f2 f3 f4 f5 f6 g1 g2 g3 g4 g5 g6 h1 h2 h3 h4 h5 h6
for i in a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 b6 c1 c2 c3 c4 c5 c6 d1 d2 d3 d4 d5 d6 
do
	#echo -n "$i "
	awk '{printf "%s %s ",$3,$4 }' harm_$i.dat
	tail -n1 tbt_k1_$i.dat | awk '{printf "%s ",$2 }' 
	sed -n '992, 1001p' tbt_k1_$i.dat | avg.awk 7
	tail -n1 tbt_k2_$i.dat | awk '{printf "%s ",$2 }' 
	sed -n '992, 1001p' tbt_k2_$i.dat | avg.awk 7
	tail -n1 tbt_k3_$i.dat | awk '{printf "%s ",$2 }' 
	sed -n '992, 1001p' tbt_k3_$i.dat | avg.awk 7
	tail -n1 tbt_k4_$i.dat | awk '{printf "%s ",$2 }' 
	sed -n '992, 1001p' tbt_k4_$i.dat | avg.awk 7
	echo
done
