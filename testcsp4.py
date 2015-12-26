#!/usr/bin/env python
# -*- coding:utf-8 -*-
#スライドする

import scipy as sp
import scipy.signal as sig 
from scipy.io.wavfile import read
from scipy.fftpack import fft, ifft
import math
import numpy as np
import matplotlib.pyplot as plt

wavfL = "left.wav"
wavfR = "right.wav"

fsL, dL = read(wavfL)
fsR, dR = read(wavfR)

ds, = dR.shape #data size
ws = ds - 10 #2 ** 16 #window size
win = sp.hamming(ws)

#plt.plot(dL[:ws],"r")
#plt.plot(dR[:ws],"b")
#plt.show()

csps = []
rang = ds-ws
print rang


for i in range(-rang, rang):
    
    if i < 0:
        dl = dL[:ws].copy()
        dr = dR[abs(i):ws+abs(i)].copy()
    else:
        dl = dL[i:ws+i].copy()
        dr = dR[:ws].copy()

    #dl *= win
    #dr *= win

    spL = fft(dl)
    spR = fft(dr)

    spL = spL[:ws/2+1]
    spR = spR[:ws/2+1]

    ampL = np.array([np.sqrt(c.real ** 2 + c.imag ** 2) for c in spL])
    ampR = np.array([np.sqrt(c.real ** 2 + c.imag ** 2) for c in spR])
    #C = np.corrcoef(spL,spR)
    #C = np.corrcoef(dl,dr)
    #C = sig.correlate(dl,dr)
    C = spL * np.conj(spR) / ampL*ampR
    csp = ifft(C)
    plt.plot(C)
    plt.show()
    #print "C:",C
    print "c max index:",np.argmax(C),", max:",max(C)
    #csp = float(C[0:1,1:2])
    #print i,"---C:",csp

    #csps.append(csp)


"""
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(311)
ax1.plot(csps,"r-")
ax1.grid(True)

maxf = csps.index(max(csps))
print "maxf:",maxf

ax2 = fig1.add_subplot(312)
dr2 = dR[:ws].copy()
dl2 = dL[:ws].copy()
ax2.plot(dl2,"b-")
ax2.plot(dr2,"r-")
ax2.grid(True)
ax2.set_xlim((0,20))
ax2.set_ylim((-2000,2000))

#maxf = 3
ax3 = fig1.add_subplot(313)
dr3 = dR[:ws].copy()
dl3 = dL[maxf:ws+maxf].copy()
ax3.plot(dr3,"b-")
ax3.plot(dl3,"r-")
ax3.grid(True)
ax3.set_xlim((0,20))
ax3.set_ylim((-2000,2000))

plt.show()
"""
