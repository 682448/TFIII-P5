#!/bin/bash
for N in 9;do
for T in 1;do
for K in 1000;do
noum="Data/Ree_T=${T}_KOnM=${K}_N=${N}.csv"
for cte in `LANG=en_US seq -5 0.1 5`;do
#for cte in -3 -4 -5;do
sed -i -e '9,9d;9d' Rk.c
sed -i "8a #define N $N"  Rk.c
sed -i -e '11,11d;11d' Rk.c
sed -i "10a #define cte (float)($cte)"  Rk.c
gcc Rk.c -o Rk.out -lm;./Rk.out <<EOF
$T
$K
$noum
EOF
done
sed -i "1a N,F,T,B,EtaOnM,KOnM,Tiempo,dt,Ree" Data/Ree_T=${T}_KOnM=${K}_N=${N}.csv
done
done
done

rm Rk.out


