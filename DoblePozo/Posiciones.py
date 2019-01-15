#!/usr/bin/env python
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import scipy.integrate as integrate
from pylab import *

# El nombre con el que se guarda el histograma de posción
noum=input()
# Pongo k_b*T=0.2
T=float(input())
B=float(input())
#print(T)
#print(B)
# Cargo los datos
Data=pd.read_csv('Data_Basico.csv')
numbins=120
# Saco la gráfica con el histograma
n, bins, patches = plt.hist(Data.loc[:,"Posicion"], numbins, density=True, facecolor='g', alpha=0.2,label='data') 
# Utiles para normalizar y hacer el curve_fit
hist,bins = np.histogram(Data.loc[:,"Posicion"],bins = numbins) 
#Paso=bins[1]-bins[0]
#xData=np.arange(min(bins)+Paso/2,max(bins),Paso)
#xData=xData.astype(np.float)
#hist=hist.astype(np.float)
#x=xData
# Calculo el area de histograma para poder normalizar
#Area=0
#for i in range(0,len(bins)-1,1):
#    Area=Area+Paso*hist[i]
# Defino la funcion para el curve_fit
#def gauss(x,A):
#    return A*np.exp(-(x**2-1)**2/T)

def f(x):
    return np.exp(-0.5*B*(x**2-1)**2/T)

AA=integrate.quad(lambda x: f(x),-6,6)
print("Constante de normalizacion:%f"%(AA[0]))
# Hago el ajuste
#popt, pcov = curve_fit(gauss,xData,hist/Area)
# Los valores de x que tomará la funcion
xfine = np.linspace(-2,2,500)  
# Grafica de la funcion
plot(xfine,f(xfine)/AA[0],color='g',lw=3,label='model')
#plt.ylim(0,2.5)
plt.xlabel('Posición')
plt.ylabel('Densidad de probabilidad')
plt.grid(True)
plt.legend()
# Guardo la imagen
plt.savefig(noum)

