#!/usr/bin/python3

from math import *
import numpy as np

gamma=2*pi*29.1646943 # rad/s/uT

b1=pi/2/gamma/2 # uT for a 2 s long b1 pulse
print(b1)
b0=1.0 # uT

t=180.0 # s, free precession time
tau=pi/gamma/b1/2 # induce perfect spin reversal on resonance


def ramsey(w,b1):

    # implementation of Ramsey eq. (12)
    # N.F. Ramsey, Phys. Rev. 78, 695 (1950).
    # This is the probability to find the particle in the flipped state.

    # intermediate equal signs are to relate to notations of Ramsey
    

    w0=gamma*b0 # eq (13)
    w1=twob=gamma*b1 # eq (13)

    deltaw=lam=w0-w # eq (5)
    wprime=a=np.sqrt(deltaw**2+w1**2) # eq (5)

    costheta=deltaw/wprime # eq (5)
    sintheta=w1/wprime # eq (5)
    
    return 4*w1**2/wprime**2*np.sin(0.5*wprime*tau)**2*(np.cos(0.5*deltaw*t)*np.cos(0.5*wprime*tau)-deltaw/wprime*np.sin(0.5*deltaw*t)*np.sin(0.5*wprime*tau))**2 # eq (12)


f0=gamma*b0/2/pi # central frequency to graph
deltaf=10*gamma*b1/2/pi # range of frequencies to graph
n=100000 # number of points to graph

f=np.linspace(f0-deltaf,f0+deltaf,n)

# Study to see what happens if b1 is set improperly

p1=ramsey(2*pi*f,b1)
print(ramsey(2*pi*f0,b1))

p2=ramsey(2*pi*f,b1*1.1)
print(ramsey(2*pi*f0,b1*1.1))

p3=ramsey(2*pi*f,b1*0.9)
print(ramsey(2*pi*f0,b1*0.9))

import matplotlib.pyplot as plt

plt.plot(f,p1,label=r'$B_1=\pi/(2\gamma\tau)$')
plt.plot(f,p2,label=r'$B_1=\pi/(2\gamma\tau)\times 1.1$')
plt.plot(f,p3,label=r'$B_1=\pi/(2\gamma\tau)\times 0.9$')
plt.axvline(f0,color='r')
plt.axvline(f0+gamma*b1/2/pi,color='r',linestyle='--')
plt.axvline(f0-gamma*b1/2/pi,color='r',linestyle='--')
plt.legend()
plt.show()

