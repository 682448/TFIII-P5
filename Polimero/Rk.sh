#!/bin/bash
for N in 4 8 16 32 64;do
for T in 0.01 0.1 1 10;do
noum="DataMeanWithForce_F=${T}_N=${N}.csv"
for cte in `LANG=en_US seq -5 0.05 5`;do
sed -i -e '9,9d;9d' Rk.c
sed -i "8a #define N $N"  Rk.c
sed -i -e '11,11d;11d' Rk.c
sed -i "10a #define cte (float)($cte)"  Rk.c
gcc Rk.c -o Rk.out -lm;./Rk.out <<EOF
$T
$noum
EOF
done
sed -i "1a N,F,T,B,EtaOnM,KOnM,Tiempo,dt,R_G,Ree,Cinetica,Potencial"  DataMeanWithForce_F=${T}_N=${N}.csv
done
done

rm Rk.out


