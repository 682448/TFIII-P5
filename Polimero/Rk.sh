#!/bin/bash

for T in 10 0.01;do
for N in 4 8 16 32 64 128 256;do
sed -i -e '9,9d;9d' Rk.c
sed -i "8a #define N $N"  Rk.c
gcc Rk.c -o Rk.out -lm;./Rk.out <<EOF
$T
EOF
done
done

rm Rk.out


