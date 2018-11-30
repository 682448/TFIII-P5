#!/usr/bin/env python
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from pylab import *

noum=input()
Data=pd.read_csv('Data_Basico.csv')
numbins=120
n, bins, patches = plt.hist(Data.loc[:,"Posicion"], numbins, density=True, facecolor='g', alpha=0.2,label='data') 
hist,bins = np.histogram(Data.loc[:,"Posicion"],bins = numbins) 
Paso=bins[1]-bins[0]
xData=np.arange(min(bins)+Paso/2,max(bins),Paso)
xData=xData.astype(np.float)
hist=hist.astype(np.float)
x=xData
Area=0
for i in range(0,len(bins)-1,1):
    Area=Area+Paso*hist[i]

def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bigauss(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)
expected=(-1,.2,.7,1,.2,.7) 
popt, pcov = curve_fit(bigauss,xData,hist/Area,expected)

xfine = np.linspace(min(bins)+Paso/2,max(bins),1000)  # define values to plot the function for
plot(xfine,bigauss(xfine,*popt),color='g',lw=3,label='model')
plt.xlabel('Posición')
plt.ylabel('Densidad de probabilidad')
plt.grid(True)
plt.legend()
plt.savefig(noum)

