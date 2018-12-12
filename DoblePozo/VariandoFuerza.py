#!/usr/bin/env python
# coding: utf-8

# In[71]:


import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pylab import *
# Cargo los datos
Data=pd.read_csv('Data/VariandoFuerza/Probabilidades.csv')
# Defino la funcion
def f(x,A):
    return 1/(1+np.exp(A*x))
# Hago el ajuste
popt,pcov=curve_fit(f,Data.loc[:,"F"],Data.loc[:,"+x"])
# Defino el tamaño de la gráfica
plt.figure(figsize=(10,10))
# Defino los puntos que tomara la funcion en el modelo
xfine=np.linspace(-5,5,500)
# Dibujo la funcion y hago scatteringn de los datos
plot(xfine,f(xfine,*popt),color='b',lw=1,label='model')
plt.scatter(Data.loc[:,"F"],Data.loc[:,"+x"],alpha=1,s=10,label='data',color='red')
# Coloco el nombre a los ejes, el grid y la leyenda
plt.xlabel(r"Fuerza",fontsize='20')
plt.ylabel(r"Probabilidad de ocupación $x>0$",fontsize='20')
plt.grid()
plt.legend(prop={'size':'20'})
# Guardo los datos
plt.savefig("DoblePozo_IMG/VariandoFuerza/ProbabilidadFrenteFuerza+x.png")
plt.show()
# Hago otra vez el ajuste para x negativo
popt,pcov=curve_fit(f,Data.loc[:,"F"],Data.loc[:,"-x"])
# Dibujo la gráfica
plt.figure(figsize=(10,10))
plot(xfine,f(xfine,*popt),label='model')
plt.scatter(Data.loc[:,"F"],Data.loc[:,"-x"],label='data',color='red',s=10)
plt.grid()
plt.ylabel(r"Probabilidad de ocupación x<0",fontsize='20')
plt.xlabel(r"Fuerza",fontsize='20')
plt.legend(prop={'size':'20'})
plt.savefig("DoblePozo_IMG/VariandoFuerza/ProbabilidadFrenteFuerza-x.png")
plt.show()

