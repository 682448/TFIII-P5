#!/usr/bin/env python
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from pylab import *
# El nombre con el que se guarda
noum1=input()
# El nombre con el que se guarda
noum2=input()

# Cargo el fichero con la distribucion positiva
Data=pd.read_csv('Distribucion_t_pos.csv')

plt.figure(figsize=(10,10))
numbins=102
# Saco la gráfica con el histograma
n, bins, patches = plt.hist(Data.loc[:,"t"], numbins, density=True, facecolor='g', alpha=0.3,label='data') 
# Utiles para normalizar y hacer el curve_fit
hist,bins = np.histogram(Data.loc[:,"t"],bins=numbins)
Paso=bins[1]-bins[0]
xData=np.arange(min(bins)+Paso/2,max(bins),Paso)
# Calculo el Area del histograma
Area=0
for i in range(0,len(bins)-1,1):
    Area=Area+Paso*hist[i]
# Defino la funcion para el curve_fit
def f(x,A,B):
    return A*np.exp(-B*x)
# Hago el ajuste
popt,pcov = curve_fit(f,xData,hist/Area)
# Los valores de x que tomará la funcion
xfine = np.linspace(0,120,500) 
# Grafica de la funcion
plot(xfine,f(xfine,*popt),color='g',lw=3,label='model')
plt.legend()
plt.xlim(0,30)
plt.xlabel('Tiempo')
plt.ylabel('Distribución de tiempos')
plt.grid(True)
plt.savefig(noum1)



# Cargo el fichero con la distribucion negativa
Data2=pd.read_csv('Distribucion_t_neg.csv')

plt.figure(figsize=(10,10))
numbins=102
# Saco la gráfica con el histograma
n, bins, patches = plt.hist(Data2.loc[:,"t"], numbins, density=True, facecolor='g', alpha=0.3,label='data')
# Utiles para hacer el curve_fit
#  hist la frecuencia, bins guarda la posicion del bin a la que se dan las frecuencias
hist,bins=np.histogram(Data2.loc[:,"t"],bins=numbins)
#  Calculo el ancho de cada bin
Paso=bins[1]-bins[0]
xData=np.arange(min(bins)+Paso/2,max(bins),Paso)
# Calculo el Area para poder normalizar 
Area=0
#  Recorro cada bin sumando su área
for i in range(0,len(bins)-1,1):
    Area=Area+Paso*hist[i]
# La función f ya esta definida de antes
# Hago el ajuste
popt,pcov=curve_fit(f,xData,hist/Area)
# Los valores que tomará f en x
xfine=np.linspace(0,120,500)
# Hago un plot de la funcion
plot(xfine,f(xfine,*popt),color='g',lw=3,label='model')
plt.legend()
plt.xlim(0,30)
plt.xlabel('Tiempo')
plt.ylabel('Distribución de tiempos')
plt.grid(True)
plt.savefig(noum2)

