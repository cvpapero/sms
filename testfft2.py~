#!/usr/bin/env python
# -*- coding:utf-8 -*-

from scipy import arange, hamming, sin, pi 
from scipy.fftpack import fft, ifft
#import numpy as np
from matplotlib import pylab as plt


fs = 16000 #サンプリングレート 
L = 1024 #長さ

sin200 = sin(2.*pi*arange(L)*200./fs)
sin600 = sin(2.*pi*arange(L)*600./fs)
sin900 = sin(2.*pi*arange(L)*900./fs)

sig = sin200 + sin600 + sin900

#win function
win = hamming(L)

#フーリエ変換
spec = fft(sig)
spec_w = fft(sig*win)
h_spec = abs(spec[:L/2+1])
h_spec_w = abs(spec_w[:L/2+1])

#フーリエ逆変換
r_sig = ifft(spec_w)
r_sig /= win

# 図を表示
fig = plt.figure()
fig.add_subplot(411)
plt.plot(sig)
plt.xlim([0, L])
plt.title("1. Input signal", fontsize = 20)
fig.add_subplot(412)
plt.plot(h_spec_w)
plt.xlim([0, len(h_spec_w)])
plt.title("2. Spectrum (no window)", fontsize = 20)
fig.add_subplot(413)
plt.plot(h_spec)
plt.xlim([0, len(h_spec)])
plt.title("3. Spectrum (with window)", fontsize = 20)
fig.add_subplot(414)
plt.plot(r_sig)
plt.xlim([0, L])
plt.title("4. Resynthesized signal", fontsize = 20)

plt.show()
