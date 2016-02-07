set term x11 1
#set xlabel "turn number"
#set ylabel "{/Symbol e}_x ({/Symbol m}m)" offset 2,0
set yrange [1:6]
#set xtics nomirror 1000
p 'tbt_k1_1.dat' u 1:($7*1e6) every 20 lt 3 lc 1 lw 3 w l notitle ,\
'tbt_k1_a1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_a2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_a3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_a4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_a5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_a6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_b1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_b2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_b3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_b4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_b5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_b6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_c1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_c2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_c3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_c4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_c5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_c6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_d1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_d2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_d3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_d4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_d5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_d6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k2_1.dat' u 1:($7*1e6) every 20 lt 3 lc 2 lw 3 w l notitle ,\
'tbt_k2_a1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_a2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_a3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_a4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_a5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_a6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_b1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_b2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_b3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_b4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_b5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_b6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_c1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_c2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_c3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_c4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_c5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_c6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_d1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_d2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_d3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_d4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_d5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_d6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k3_1.dat' u 1:($7*1e6) every 20 lt 3 lc 3 lw 3 w l notitle ,\
'tbt_k3_a1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_a2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_a3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_a4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_a5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_a6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_b1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_b2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_b3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_b4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_b5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_b6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_c1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_c2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_c3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_c4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_c5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_c6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_d1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_d2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_d3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_d4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_d5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_d6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k4_1.dat' u 1:($7*1e6) every 20 lt 3 lc 4 lw 3 w l notitle ,\
'tbt_k4_a1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_a2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_a3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_a4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_a5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_a6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_b1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_b2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_b3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_b4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_b5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_b6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_c1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_c2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_c3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_c4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_c5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_c6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_d1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_d2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_d3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_d4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_d5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_d6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k1_5.dat' u 1:($7*1e6) every 20 lt 3 lc 1 lw 3 w l notitle ,\
'tbt_k1_e1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_e2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_e3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_e4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_e5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_e6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_f1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_f2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_f3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_f4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_f5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_f6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_g1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_g2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_g3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_g4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_g5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_g6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_h1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_h2.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_h3.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_h4.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_h5.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k1_h6.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k2_5.dat' u 1:($7*1e6) every 20 lt 3 lc 2 lw 3 w l notitle ,\
'tbt_k2_e1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_e2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_e3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_e4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_e5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_e6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_f1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_f2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_f3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_f4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_f5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_f6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_g1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_g2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_g3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_g4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_g5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_g6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_h1.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_h2.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_h3.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_h4.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_h5.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k2_h6.dat' u 1:($7*1e6) every 20 lt 1 lc 2 w l notitle,\
'tbt_k3_5.dat' u 1:($7*1e6) every 20 lt 3 lc 3 lw 3 w l notitle ,\
'tbt_k3_e1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_e2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_e3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_e4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_e5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_e6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_f1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_f2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_f3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_f4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_f5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_f6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_g1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_g2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_g3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_g4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_g5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_g6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_h1.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_h2.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_h3.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_h4.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_h5.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k3_h6.dat' u 1:($7*1e6) every 20 lt 1 lc 3 w l notitle,\
'tbt_k4_5.dat' u 1:($7*1e6) every 20 lt 3 lc 4 lw 3 w l notitle ,\
'tbt_k4_e1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_e2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_e3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_e4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_e5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_e6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_f1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_f2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_f3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_f4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_f5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_f6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_g1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_g2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_g3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_g4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_g5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_g6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_h1.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_h2.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_h3.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_h4.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_h5.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\
'tbt_k4_h6.dat' u 1:($7*1e6) every 20 lt 1 lc 4 w l notitle,\

reset
set term x11 2
set xrange [0:1000]
set yrange [1:2.2]
set xlabel "Turn number"
set ylabel "{/Symbol e}_x ({/Symbol m}m)" offset 2,0
set label "N=2e9" at 60, 1.83
set label "N=1e9" at 65, 1.23
set label "N=6e8" at 70, 1.08
p '../case18/track3/tbt_a0031.dat' u 1:($7*1e6) every 10 lt 1 lc rgb"#009900" lw 2 w l t "no err",\
'tbt_k2_1.dat' u 1:($7*1e6) every 10 lt 1 lc 1 lw 2 w l t "w/ err" ,\
'tbt_k2_d6.dat' u 1:($7*1e6) every 10 lt 1 lc 3 lw 2 w l t "w correction",\
'tbt_k3_1.dat' u 1:($7*1e6) every 10 lt 1 lc 1 lw 2 w l notitle,\
'tbt_k3_d6.dat' u 1:($7*1e6) every 10 lt 1 lc 3 lw 2 w l notitle,\
'../case18/track3/tbt_a0051.dat' u 1:($7*1e6) every 10 lt 1 lc rgb"#009900" lw 2 w l notitle,\
'tbt_k4_1.dat' u 1:($7*1e6) every 10 lt 1 lc 1 lw 2 w l notitle,\
'tbt_k4_d6.dat' u 1:($7*1e6) every 10 lt 1 lc 3 lw 2 w l notitle,\
'../case18/track3/tbt_a0101.dat' u 1:($7*1e6) every 10 lt 1 lc rgb"#009900" lw 2 w l notitle

set term post eps enhanced color font 'Helvetica,28'
set key samplen 2 spacing 0.90 font ',28' at 900,1.6 
set output 'emit_all.eps'
replot

set term x11 3
#set xlabel "turn number"
#set ylabel "{/Symbol e}_x ({/Symbol m}m)" offset 2,0
set yrange [1:2.2]
set xrange [0:1000]
p 'tbt_k4_a1.dat' u 1:($7*1e6) every 20 lt 1 lc 1 w l notitle,\
'tbt_k4_a2.dat' u 1:($7*1e6) every 10 lt 1 lc 1 w l notitle,\
'tbt_k4_a3.dat' u 1:($7*1e6) every 10 lt 1 lc 1 w l notitle,\
'tbt_k4_a4.dat' u 1:($7*1e6) every 10 lt 1 lc 1 w l notitle,\
'tbt_k4_a5.dat' u 1:($7*1e6) every 10 lt 1 lc 1 w l notitle,\
'tbt_k4_a6.dat' u 1:($7*1e6) every 10 lt 1 lc 1 w l notitle,\
'tbt_k4_b1.dat' u 1:($7*1e6) every 10 lt 1 lc 2 w l notitle,\
'tbt_k4_b2.dat' u 1:($7*1e6) every 10 lt 1 lc 2 w l notitle,\
'tbt_k4_b3.dat' u 1:($7*1e6) every 10 lt 1 lc 2 w l notitle,\
'tbt_k4_b4.dat' u 1:($7*1e6) every 10 lt 1 lc 2 w l notitle,\
'tbt_k4_b5.dat' u 1:($7*1e6) every 10 lt 1 lc 2 w l notitle,\
'tbt_k4_b6.dat' u 1:($7*1e6) every 10 lt 1 lc 2 w l notitle,\
'tbt_k4_c1.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
'tbt_k4_c2.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
'tbt_k4_c3.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
'tbt_k4_c4.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
'tbt_k4_c5.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
'tbt_k4_c6.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
'tbt_k4_d1.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
'tbt_k4_d2.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
'tbt_k4_d3.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
'tbt_k4_d4.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
'tbt_k4_d5.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
'tbt_k4_d6.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\


reset
set term x11 4
#set xlabel "turn number"
#set ylabel "{/Symbol e}_x ({/Symbol m}m)" offset 2,0
set yrange [1:2.2]
set xrange [0:1000]
p 'tbt_k3_c11.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
  'tbt_k3_c12.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
  'tbt_k3_c13.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
  'tbt_k3_c14.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
  'tbt_k3_c15.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
  'tbt_k3_c16.dat' u 1:($7*1e6) every 10 lt 1 lc 3 w l notitle,\
  'tbt_k4_c11.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
  'tbt_k4_c12.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
  'tbt_k4_c13.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
  'tbt_k4_c14.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
  'tbt_k4_c15.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\
  'tbt_k4_c16.dat' u 1:($7*1e6) every 10 lt 1 lc 4 w l notitle,\

