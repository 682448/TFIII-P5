#!/usr/bin/env python
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
noum=input()
'''
# Cargo los datos
Data=pd.read_csv('Data_Basico.csv')
# Calculo varianza y media para el ajuste a gausiana
(mu, sigma) = norm.fit(Data.loc[:,"Velocidad"])
# Calculo la función gaussiana
dist = norm(mu, sigma)
# Hago el histograma
n, bins, patches = plt.hist(Data.loc[:,"Velocidad"], 60, density=True, facecolor='b', alpha=0.2,label='data')
# Genero los valores de x que tomará la gaussiana
x = np.linspace(min(bins), max(bins), 1000)
# Represento la gaussiana
plt.plot(x, dist.pdf(x), 'b', linewidth = 3,label='model')
plt.xlabel('Velocidad')
plt.ylabel('Densidad de probabilidad')
plt.grid(True)
plt.legend()
# Guardo la gráfica
plt.savefig(noum)
'''

