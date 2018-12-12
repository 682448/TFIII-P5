#!/bin/bash
chmod +x ./DistribucionTiempos.py
cte=0

for T in 1;do
for cte in `LANG=en_US seq -5 0.05 5` ;do
for B in 2 ;do
for Eta in 1;do
noum1="DoblePozo_IMG/DisTiempo/DisTiempo_pos_T=${T}_F=${cte}_B=${B}_Eta=${Eta}.png"
noum2="DoblePozo_IMG/DisTiempo/DisTiempo_neg_T=${T}_F=${cte}_B=${B}_Eta=${Eta}.png"
echo F=$cte T=${T} B=${B} Eta=${Eta}
gcc Rk2.c -o Rk2.out -lm; ./Rk2.out <<EOF
$T
$B
$Eta
$cte
EOF

./DistribucionTiempos.py <<EOF
$noum1
$noum2
EOF

done
done
done
done
rm Rk2.out
rm Distribucion_t_neg.csv
rm Distribucion_t_pos.csv

