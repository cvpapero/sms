#!/usr/bin/env python
# -*- coding:utf-8 -*-

import scipy as sp
from scipy.io.wavfile import read
from scipy.fftpack import fft, ifft

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


rL, = dL.shape
rR, = dR.shape

win = sp.hamming(rL)
#sp = fft(dL*win)

#FFT(left)
spL = fft(dL)
spL = spL[:rL/2+1]
#FFT(right)
spR = fft(dR)
spR = spR[:rR/2+1]

#振幅スペクトル
#ampL = np.array([np.sqrt(c.real ** 2 + c.imag **2) for c in spL])
#ampR = np.array([np.sqrt(c.real ** 2 + c.imag **2) for c in spR])
#ampL = np.array([np.abs(c) for c in spL])
#ampR = np.array([np.abs(c) for c in spR])
ampL = np.sqrt(spL * np.conj(spL))
ampR = np.sqrt(spR * np.conj(spR))

nL = spL / ampL
nR = spR / ampR

#C = nL * np.conj(nR)

C = (spL*np.conj(spR)) / (ampL*ampR)

print C
csp = ifft(C)
#print "sp_h:",sp_h

#plt.plot(C)
plt.plot(csp)
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
