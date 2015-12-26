#!/usr/bin/env python
# -*- coding:utf-8 -*-

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
r, = dR.shape
csps = []
for i in range(1):

    dl = dL[i:]
    dr = dR[:r-i]
    rL, = dl.shape
    rR, = dr.shape

    win = sp.hamming(rL)
    #dl *= win
    #dr *= win
    #sp = fft(dL*win)
    
    #FFT(left)
    spL = fft(dl)
    #plt.plot(spL)
    #plt.show()
    spL = spL[:rL/2+1]
    #FFT(right)
    spR = fft(dr)
    spR = spR[:rR/2+1]

    #print "corr:",np.corrcoef(spL,spR)
    #振幅スペクトル
    ampL = np.sqrt(spL * np.conj(spL))
    ampR = np.sqrt(spR * np.conj(spR))
    
    nL = spL / ampL
    #nR = spR / ampR
    
    C = (spL*np.conj(spR)) / (ampL*ampR)
    
    csp = ifft(C)

    #csp_coef=max(csp)
    #print "csp:",csp_coef
    #csps.append(csp_coef)
    #plt.plot(C)
    plt.plot(dl,"r--")
    plt.plot(dr,"b-")

    #plt.xlim((0,500))
    #plt.ylim((-5000,5000))
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
