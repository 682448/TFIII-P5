#!/bin/bash
T=0.2
chmod +x ./Velocidades.py
chmod +x ./Posiciones.py
chmod +x ./FasEvol.py

for B in 1 2 5;do
for Eta in 0.1 1 5;do
Posiciones="Posicion_T=${T}_B=${B}_Eta=${Eta}.png"
noum2="Velocidad_T=${T}_B=${B}_Eta=${Eta}.png"
noum3="Ecolucion_T=${T}_B=${B}_Eta=${Eta}.png"
noum4="Diagrama_Fases_T=${T}_B=${B}_Eta=${Eta}.png"
gcc Rk.c -o Rk.out -lm;./Rk.out <<EOF
$T
$B
$Eta
EOF

#./Posiciones.py <<EOF
#$Posiciones
#EOF

#./Velocidades.py <<EOF
#$noum2
#EOF

#./FasEvol.py <<EOF
#$noum3
#$noum4
#EOF
done
done

rm Rk.out
rm Data_Basico.csv

