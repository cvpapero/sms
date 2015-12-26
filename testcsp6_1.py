#!/usr/bin/env python
# -*- coding:utf-8 -*-
#スライドする

#

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
ws = ds - 100 #2 ** 16 #window size
#win = sp.hamming(ws)

#plt.plot(dL[:ws],"r")
#plt.plot(dR[:ws],"b")
#plt.show()

#csps = []
#rang = ds-ws
#print rang
def CCF(array1,array2):
    #print "array1",array1
    s1, = array1.shape
    s2, = array2.shape

    if s1 < s2:
        array1,array2 = array2,array1
        s1, = array1.shape
        s2, = array2.shape    
    
    ccf = []
    
    for i in range(s1 - s2):
        temp = 0
        for j in range(s2):
            temp += array2[j] * array1[j - i]
            #print j,",",temp
        ccf.append(temp / s2)

    #print ccf

    return ccf

spL = fft(dL)
spR = fft(dR)

spL = spL[:ds/2+1]
spR = spR[:ds/2+1]

#振幅スペクトルらしい
ampL = np.array([np.sqrt(c.real ** 2 + c.imag ** 2) for c in spL])
ampR = np.array([np.sqrt(c.real ** 2 + c.imag ** 2) for c in spR])
nspL = spL/ampL
nspR = spR/ampR

sds, = spL.shape
ws = sds - 100

#C1 = nspL[:ws]np.conj(nspR[:sds])
#C2 = CCF(nspL[:sds],np.conj(nspR[:ws]))

#C = CCF(spL[:ws],spR[:sds])
C1 = sig.correlate(nspR[:sds],nspL[:ws], mode="valid")
C2 = sig.correlate(nspL[:sds],nspR[:ws], mode="valid")
csp1 = ifft(C1)
csp2 = ifft(C2)

#print "ds:",sds,", ws:",ws 
#print "c1 max index:",C1.index(max(C1)),", max:",max(C1)
#print "c2 max index:",C2.index(max(C2)),", max:",max(C2)
print "csp1 max index:",np.argmax(csp1),", max:",max(csp1)
print "csp2 max index:",np.argmax(csp2),", max:",max(csp2)
#
#print "c.shape:",C.shape
#print "dL.shape:",dL.shape
#print "spL.shape:",spL.shape
plt.plot(csp1,"r")
plt.plot(csp2,"b")
#plt.plot(C2,"b")


plt.show()



"""
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

    
    #C = np.corrcoef(spL,spR)
    #C = np.corrcoef(dl,dr)
    C = sig.correlate(dl,dr)
    print "C:",C
    #csp = float(C[0:1,1:2])
    #print i,"---C:",csp

    #csps.append(csp)


"""
"""
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(311)
ax1.plot(csp,"r-")
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
"""
#plt.show()

