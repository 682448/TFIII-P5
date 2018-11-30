#!/bin/bash
T=1
chmod +x ./Velocidades.py
chmod +x ./Posiciones.py
chmod +x ./FasEvol.py

for dt in 1E-1 1E-2 1E-3 1E-4 ;do
for B in 0.5;do
for Eta in 0;do
Posiciones="Oscilador_IMG/Posicion_dt=${dt}_T=${T}_B=${B}_Eta=${Eta}.png"
noum2="Oscilador_IMG/Velocidad_dt=${dt}_T=${T}_B=${B}_Eta=${Eta}.png"
#noum3="Oscilador_IMG/Evolucion_dt=${dt}_T=${T}_B=${B}_Eta=${Eta}.png"
#noum4="Oscilador_IMG/Diagrama_Fases_dt=${dt}_T=${T}_B=${B}_Eta=${Eta}.png"
gcc Rk.c -o Rk.out -lm;./Rk.out <<EOF
$T
$B
$Eta
$dt
EOF

./Posiciones.py <<EOF
$Posiciones
EOF

./Velocidades.py <<EOF
$noum2
EOF

#./FasEvol.py <<EOF
#$noum3
#$noum4
#EOF
done
done
done

rm Rk.out
rm Data_Basico.csv

