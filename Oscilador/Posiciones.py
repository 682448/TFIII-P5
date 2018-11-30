#!/usr/bin/env python
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd

Data=pd.read_csv('Data_Basico.csv')
(mu, sigma) = norm.fit(Data.loc[:,"Posicion"])
dist = norm(mu, sigma)
noum=input()
n, bins, patches = plt.hist(Data.loc[:,"Posicion"], 60, density=True, facecolor='g', alpha=0.2,label='data')
#x = np.linspace(min(bins), max(bins), 1000)
#plt.plot(x, dist.pdf(x), 'g', linewidth = 3,label='model')
plt.xlabel('Posición')
plt.ylabel('Densidad de probabilidad')
plt.grid(True)
plt.legend()
plt.savefig(noum)
