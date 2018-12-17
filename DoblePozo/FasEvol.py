#!/usr/bin/env python
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
noum1=input()# Para el diagrama de evolución
noum2=input()# Para el diagrama de fases
'''
Data=pd.read_csv('Data_Basico.csv')

# Esta parte se encarga de hacer el diagrama de evolución
plt.figure(figsize=(25,10))
plt.xlim(0,10000)
plt.plot(Data.loc[:,"Tiempo"],Data.loc[:,"Posicion"])
plt.xlabel('Tiempo',fontsize=40)
plt.ylabel('Posición',fontsize=40)
plt.savefig(noum1)
'''
'''
# Esta parte se encarga de hacer el diagrama de fases
plt.figure(figsize=(10,10))
plt.scatter(Data.loc[:,"Posicion"],Data.loc[:,"Velocidad"],s=0.001)
plt.xlim(-2,2)
plt.xlabel('Posición',fontsize=20)
plt.ylabel('Velocidad',fontsize=20)
plt.savefig(noum2)
'''


