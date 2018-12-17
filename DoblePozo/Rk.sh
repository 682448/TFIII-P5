#!/bin/bash
chmod +x ./Velocidades.py
chmod +x ./Posiciones.py
chmod +x ./FasEvol.py

for T in 100;do
for B in 1 2 5;do
for Eta in 0.1 1 5;do
Posiciones="DoblePozo_IMG/Posicion_T=${T}_B=${B}_Eta=${Eta}_Termalizado.png"
noum2="DoblePozo_IMG/Velocidad_T=${T}_B=${B}_Eta=${Eta}.png"
noum3="DoblePozo_IMG/Evolucion_T=${T}_B=${B}_Eta=${Eta}_Termalizado.png"
noum4="DoblePozo_IMG/Diagrama_Fases_T=${T}_B=${B}_Eta=${Eta}.png"
echo T=${T} B=${B} Eta=${Eta}
gcc Rk.c -o Rk.out -lm; ./Rk.out <<EOF
$T
$B
$Eta
EOF

./Posiciones.py <<EOF
$Posiciones
$T
EOF

./Velocidades.py <<EOF
$noum2
EOF

./FasEvol.py <<EOF
$noum3
$noum4
EOF
done
done
done

rm Rk.out
rm Data_Basico.csv

