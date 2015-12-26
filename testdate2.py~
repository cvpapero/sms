#!/usr/bin/env python
# -*- coding:utf-8 -*-
import wave
import numpy as np
import scipy as sp
import scipy.fftpack
from scipy.io.wavfile import read
from pylab import *

def main():
    #wavファイルをインプットしてみる

    #wf = wave.open("left.wav" , "r" )
    #fs = wf.getframerate()  # サンプリング周波数
    #x = wf.readframes(wf.getnframes()) #初期化
    #x = frombuffer(x, dtype= "int16") / max(x)  # -1 - +1に正規化
    #wf.close()

    wavfL = "left.wav"
    wavfR = "right.wav"
    
    fsL, dL = read(wavfL)
    fsR, dR = read(wavfR)

    rL, = dL.shape
    rR, = dR.shape
    #fs = 8000 # Sampling rate
    #L = 256 # Signal length

    # 440[Hz]のサイン波を作る。
    #sine_440 = 5* sp.sin(2. * sp.pi * sp.arange(L) * 500. / fs)
    # 600[Hz]のサイン波を作る。
    #sine_600 = 2 * sp.sin(2. * sp.pi * sp.arange(L) * 1000. / fs)
    # 800[Hz]のサイン波を作る。
    #sine_800 = 3 * sp.sin(2. * sp.pi * sp.arange(L) * 1500. / fs)

    # 全部足す
    #x = sine_440 + sine_600 + sine_800
    win = sp.hamming(fsL)
    dL *= win
    dL /= max(x)

    start = 0  # サンプリングする開始位置
    N = 256    # FFTのサンプル数

    X = np.fft.fft(x[start:start+N])  # FFT
    #    X = scipy.fftpack.fft(x[start:start+N])         # scipy版
    
    freqList = np.fft.fftfreq(N, d=1.0/fs)  # 周波数軸の値を計算
    #    freqList = scipy.fftpack.fftfreq(N, d=1.0/ fs)  # scipy版
    
    amplitudeSpectrum = [np.sqrt(c.real ** 2 + c.imag ** 2) for c in X]  # 振幅スペクトル
    phaseSpectrum = [np.arctan2(int(c.imag), int(c.real)) for c in X]    # 位相スペクトル
    
    # 波形を描画
    subplot(311)  # 3行1列のグラフの1番目の位置にプロット
    plot(range(start, start+N), x[start:start+N])
    axis([start, start+N, -1.0, 1.0])
    xlabel("time [sample]")
    ylabel("amplitude")

    # 振幅スペクトルを描画
    subplot(312)
    #plot(freqList, amplitudeSpectrum, marker= 'o', linestyle='-')
    plot(freqList, abs(X), marker= 'o', linestyle='-')
    axis([0, fs/2, 0, 50])
    xlabel("frequency [Hz]")
    ylabel("amplitude spectrum")

    # 位相スペクトルを描画
    subplot(313)
    plot(freqList, phaseSpectrum, marker= 'o', linestyle='-')
    axis([0, fs/2, -np.pi, np.pi])
    xlabel("frequency [Hz]")
    ylabel("phase spectrum")

    show()



if __name__ == "__main__":
    main()
