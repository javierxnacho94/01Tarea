#Tarea 1

import numpy as np
import matplotlib.pyplot as plt
import math as m
from astropy import constants as co
from astropy import units as un
from scipy import integrate as inte
import time as t

#parte 1

datos=np.loadtxt('sun_AM0.dat')
long_onda=datos[:,0]
flujo=datos[:,1]
long_onda=long_onda*un.nm #Nanometros
flujo=flujo*un.J*(un.m**-2)*(un.nm**-1)*(un.s**-1) # J s-1 m-2 nm-1
long_onda=long_onda.to('um') #pasamos a micron
flujo=flujo.to('erg/(s cm2 um)') #pasamos a erg s-1 cm-2 um-1

fig1 = plt.figure()
ax=fig1.add_subplot(111)
ax.fill_between(long_onda,flujo,edgecolor='k',facecolor='w')
ax.set_xlim(0, 5)
ax.set_ylim(0, 2.5e6)
plt.xlabel('Longitud de onda [$um$]')
plt.ylabel('Flujo [$erg$ $s^{-1}cm^{-2}um^{-1}$]')
plt.title('Espectro del Sol')
plt.savefig('espectro_sol.png')

#parte 2
tiempo1=t.time()
Integral=0
for i in range(0,len(flujo)-2):
    Integral=Integral+(flujo[i+1]+flujo[i])*(long_onda[i+1]-long_onda[i])/2

tiempo1=t.time()-tiempo1


print Integral #Esta en unidades de erg s-1 cm-2
tiempo2=t.time()
Integral2=inte.trapz(flujo, long_onda)
tiempo2=t.time()-tiempo2
print Integral2

#parte 3

def fint(x):   #funcion de la integral
    return (np.tan(x)**3) / ((np.cos(x)**2)*((np.exp(np.tan(x)))-1))


a=0.05
b=m.pi/2-0.05
n=2000
dx=(b-a)/n
x=np.linspace(a,b,n)
x=np.array(x)
y=fint(x)
integral=0
tiempo3=t.time()
for i in range(0,len(y)-1):
    integral=integral+(y[i+1]+y[i])*dx/2
tiempo3=t.time()-tiempo3
T=5778*un.K

factor=((2*np.pi*co.h)/((co.c)**2)) * (co.k_B*T/(co.h))**4

integral=integral*factor
integral=integral.to('erg/(s cm2)')
print integral
tiempo4=t.time()
integral2=inte.trapz(y,x)
tiempo4=t.time()-tiempo4
print integral2*factor.to('erg/(s cm2)')

R2=(Integral/integral)/(4*np.pi)
R=np.sqrt(R2)
print R
print tiempo1, tiempo2, tiempo3, tiempo4
