

```python
import matplotlib.pyplot as plt
import pandas as pd


PosVel=pd.read_csv('Rk_Hist_dt00001.csv')
n, bins, patches = plt.hist(PosVel.loc[:,"z"], 95, density=True, facecolor='g', alpha=0.75)
plt.xlabel('P(z)')
plt.ylabel('z')
plt.title('Distribución Gaussiana')
plt.axis([-.25, .25, 0, 10])
plt.grid(True)
plt.show()
n, bins, patches = plt.hist(PosVel.loc[:,"Posicion"], 60, density=True, facecolor='g', alpha=0.75)
plt.xlabel('Posicion')
plt.ylabel('P(x)')
plt.title('Distribución de posiciones')
plt.axis([-4, 4, 0, 0.7])
plt.grid(True)
plt.show()
n, bins, patches = plt.hist(PosVel.loc[:,"Velocidad"], 100, density=True, facecolor='g', alpha=0.75)
plt.xlabel('Velocidad')
plt.ylabel('P(v)')
plt.title('Distribución de velocidades')
plt.axis([-4, 4, 0, 0.5])
plt.grid(True)
plt.show()
```


![png](README_files/README_0_0.png)



![png](README_files/README_0_1.png)



![png](README_files/README_0_2.png)



```python
Pos=pd.read_csv('Rk_Eta0.csv')
plt.figure(figsize=(25,3))
plt.plot(Pos.loc[:,"Tiempo"],Pos.loc[:,"Posicion"])
plt.figure(figsize=(5,5))
plt.scatter(Pos.loc[:,"Posicion"],Pos.loc[:,"Velocidad"],s=0.001)
```




    <matplotlib.collections.PathCollection at 0x7fa08b9f1080>




![png](README_files/README_1_1.png)



![png](README_files/README_1_2.png)



```python
Pos=pd.read_csv('Rk_Eta001.csv')
plt.figure(figsize=(25,3))
plt.plot(Pos.loc[:,"Tiempo"],Pos.loc[:,"Posicion"])
plt.figure(figsize=(5,5))
plt.scatter(Pos.loc[:,"Posicion"],Pos.loc[:,"Velocidad"],s=0.001)
```




    <matplotlib.collections.PathCollection at 0x7fa08b7bd240>




![png](README_files/README_2_1.png)



![png](README_files/README_2_2.png)



```python
Pos=pd.read_csv('Rk_Eta01.csv')
plt.figure(figsize=(25,3))
plt.plot(Pos.loc[:,"Tiempo"],Pos.loc[:,"Posicion"])
plt.figure(figsize=(5,5))
plt.scatter(Pos.loc[:,"Posicion"],Pos.loc[:,"Velocidad"],s=0.001)
```




    <matplotlib.collections.PathCollection at 0x7fa08bb12748>




![png](README_files/README_3_1.png)



![png](README_files/README_3_2.png)



```python
Pos=pd.read_csv('Rk_Eta1.csv')
plt.figure(figsize=(25,3))
plt.plot(Pos.loc[:,"Tiempo"],Pos.loc[:,"Posicion"])
plt.figure(figsize=(5,5))
plt.scatter(Pos.loc[:,"Posicion"],Pos.loc[:,"Velocidad"],s=0.001)
```




    <matplotlib.collections.PathCollection at 0x7fa08bfa9828>




![png](README_files/README_4_1.png)



![png](README_files/README_4_2.png)



```python
Pos=pd.read_csv('Rk_Eta10.csv')
plt.figure(figsize=(25,3))
plt.plot(Pos.loc[:,"Tiempo"],Pos.loc[:,"Posicion"])
plt.figure(figsize=(5,5))
plt.scatter(Pos.loc[:,"Posicion"],Pos.loc[:,"Velocidad"],s=0.001)
```




    <matplotlib.collections.PathCollection at 0x7fa08b9d3128>




![png](README_files/README_5_1.png)



![png](README_files/README_5_2.png)



```python

```
