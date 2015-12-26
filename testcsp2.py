#!/usr/bin/env python
# -*- coding:utf-8 -*-
#スライドする

import scipy as sp
from scipy.io.wavfile import read
from scipy.fftpack import fft, ifft
import math
import numpy as np
import matplotlib.pyplot as plt

wavfL = "left.wav"
wavfR = "right.wav"

fsL, dL = read(wavfL)
fsR, dR = read(wavfR)

print "L sampling rate:", fsL
print "L shape", dL.shape 

print "R sampling rate:", fsR
print "R shape", dR.shape 

print dL
print dR


ds, = dR.shape #data size
ws = 2 ** 16 #window size

plt.plot(dL[:ws],"r")
plt.plot(dR[:ws],"b")
plt.show()

csps = []
win = sp.hamming(ws)
for i in range(ds-ws):

    #dl = dL[i:ws+i]
    #dr = dR[:ws]
    dl = dL[:ws]
    dr = dR[i:ws+i]

    dl *= win
    dr *= win

    spL = fft(dl)
    spR = fft(dr)
    spL = spL[:ws/2+1]
    spR = spR[:ws/2+1]

    #print "corr:",np.corrcoef(spL,spR)
    #振幅スペクトル
    ampL = np.sqrt( np.dot(spL,np.conj(spL).T) )
    ampR = np.sqrt( np.dot(spR,np.conj(spR).T) )
    #ampL = np.array([np.sqrt(c.real ** 2 + c.imag ** 2) for c in spL])
    #ampR = np.array([np.sqrt(c.real ** 2 + c.imag ** 2) for c in spR])
    #ampL = np.sqrt(sum(c ** 2 for c in spL))
    #ampR = np.sqrt(sum(c ** 2 for c in spR))
    #print "ampL:",ampL
    #print "ampR:",ampR
    #nL = spL / ampL
    #nR = spR / ampR
    #print "sL*sR.T:",np.dot(spL,spR.T)
    #C = np.dot(spL,spR.T) / (ampL*ampR)
    C = np.corrcoef(spL,spR)
    csp = float(C[0:1,1:2])
    print i,"---C:",csp
    #csp = ifft(C)

    #csp_coef=max(csp)
    #print "csp:",csp_coef
    csps.append(csp)

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(311)
ax1.plot(csps,"r-")
ax1.grid(True)

maxf = csps.index(max(csps))

ax2 = fig1.add_subplot(312)
dl = dL[:ws]
dr = dR[:ws]

ax2.plot(dl,"b-")
ax2.plot(dr,"r-")
ax2.grid(True)

ax3 = fig1.add_subplot(313)
#dl = dL[:ws]
#dr = dR[:ws]

ax3.plot(dL,"b-")
ax3.plot(dR,"r-")
ax3.grid(True)


plt.show()

"""
r, = dR.shape


fig1 = plt.figure(1)
ax1 = fig1.add_subplot(211)

for i in range(10):
    print "i:",i

    #dL = np.delete(dL,i,1)
    #dR = np.delete(dL,r-i-1,1)
    dl = dL[i:]
    dr = dR[:r-i]
    #dr = dR[i:]
    #dl = dL[:r-i]

    
     np.corrcoef(dl,dr)

plt.plot(dL,"ro-")
plt.plot(dR,"bo-")
plt.grid(True)
plt.title("Left")
plt.xlim((4050,4100))
plt.ylim((-5000,5000))

"""

"""
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(211)
ax1.plot(dL,"ro-")
ax1.plot(dR,"bo-")
ax1.grid(True)
ax1.set_title("Left")
ax1.set_xlim((4050,4100))
ax1.set_ylim((-5000,5000))
"""

"""
ax2 = fig1.add_subplot(212)
ax2.plot(dR,"bo-")
ax2.grid(True)
ax2.set_title("Right")
ax2.set_xlim((4050,4100))
ax2.set_ylim((-5000,5000))
"""

#plt.show()
