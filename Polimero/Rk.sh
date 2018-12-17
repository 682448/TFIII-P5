#!/bin/bash

for N in 4 8 16 32 64 128 256;do
sed -i -e '9,9d;9d' Rk.c
sed -i "8a #define N $N"  Rk.c
gcc Rk.c -o Rk.out -lm;./Rk.out
done

rm Rk.out


