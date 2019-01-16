#!/bin/bash
for N in 5;do
for T in 1;do
for K in 100;do
#noum="Data/Ree_T=${T}_KOnM=${K}_N=${N}.csv"
noum="Basura.csv"
noum2="Data/RadioGiro.csv"
#for cte in `LANG=en_US seq 5.6 0.2 10`;do
for cte in 0;do
sed -i -e '9,9d;9d' Rk.c
sed -i "8a #define N $N"  Rk.c
sed -i -e '11,11d;11d' Rk.c
sed -i "10a #define cte (float)($cte)"  Rk.c
gcc Rk.c -o Rk.out -lm;./Rk.out <<EOF
$T
$K
$noum
$noum2
EOF
done
sed -i "1a N,F,T,B,EtaOnM,KOnM,Tiempo,dt,Ree,Ree_best,err,err1" Data/Ree_T=${T}_KOnM=${K}_N=${N}.csv
done
done
done

rm Rk.out

# 9 17 33 65
