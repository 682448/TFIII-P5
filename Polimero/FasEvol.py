#!/usr/bin/env python
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
noum1=input()
noum2=input()
Data=pd.read_csv('Data_Basico.csv')
plt.figure(figsize=(25,10))
plt.xlim(0,10000)
plt.plot(Data.loc[:,"Tiempo"],Data.loc[:,"Posicion"])
plt.xlabel('Tiempo',fontsize=40)
plt.ylabel('Posición',fontsize=40)
plt.savefig(noum1)
plt.figure(figsize=(10,10))
plt.scatter(Data.loc[:,"Posicion"],Data.loc[:,"Velocidad"],s=0.000003)
plt.xlabel('Posición',fontsize=20)
plt.ylabel('Velocidad',fontsize=20)
plt.savefig(noum2)



