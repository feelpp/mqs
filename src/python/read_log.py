import numpy as np
from math import sqrt
from scipy.stats import linregress
import scipy.sparse as sp
import scipy.sparse.linalg as spsp
import matplotlib.pyplot as plt

fichier=open("log","r")

t=[]
IV0=[]
IV1=[]
Bz0=[]
ligne = fichier.readline().rstrip("\n\r")
while (ligne != "J_induct:coil"):
    ligne = fichier.readline().rstrip("\n\r")

ligne=fichier.readline().rstrip("\n\r")
while( ligne.split(":")[0] != "[env] Time "):
    ligne=fichier.readline().rstrip("\n\r")
    data=ligne.split(",")
    t.append(float(data[0].split("=")[1]))
    Bz0.append(float(data[5].split("}")[0]))
    ligne=fichier.readline().rstrip("\n\r")
    data=ligne.split(",")
    IV0.append(float(data[1].split("=")[1]))
    IV1.append(float(data[2].split("=")[1]))
    ligne=fichier.readline().rstrip("\n\r")
fichier.close()

minB = min(Bz0)
max0 = max(IV0)
#max1 = min(IV1)

Bz0 = [i/minB for i in Bz0]
IV0 = [i/max0 for i in IV0]
#IV1 = [i/max1 for i in IV1]

plt.plot(t,Bz0,'.',label="Bz(t)")
plt.plot(t,IV0,'+',label='I(t)')
#plt.plot(t,IV1,'x')
plt.legend()
plt.show()
