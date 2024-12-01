
# Assignment 1 - Analog Communication Systems

## Overview of the Assignment
This assignment involves simulating and analyzing an analog communication system using MATLAB. Key aspects include:
1. Data generation and signal processing.
2. Amplitude Modulation (AM) and Frequency Modulation (FM).
3. Simulating the effect of noise on transmitted signals.
4. Demodulation and signal reconstruction.

The assignment consists of analytical questions and MATLAB simulations.

---

## Questions and Solutions

### Question 1: AM Modulation - MATLAB (50%)

#### **1.1 Data Generation**
1. **Steps:**
   - Load the provided `.wav` file, extract the left channel, and create a one-dimensional signal \( v_m(t) \).
   - Calculate the Fourier Transform \( V_m(f) \) using MATLAB's `fft` function.
2. **Output:**
   - Plots of \( v_m(t) \) and its Fourier transform \( V_m(f) \).
3. **Bandwidth Evaluation:**
   - Use MATLAB's `powerbw` function to calculate the signal's bandwidth.

---

#### **1.2 Modulator**
1. **Steps:**
   - Simulate AM modulation using `ammod` with a carrier frequency \( f_c = 15 \text{kHz} \) and modulation index \( k_{AM} = 0.02 \).
   - Calculate the Fourier Transform \( V_{AM}(f) \) of the modulated signal.
2. **Output:**
   - Plot of \( V_{AM}(f) \).
   - Listen to the modulated signal using MATLAB's `sound` function.

---

#### **1.3 Channel**
1. **Steps:**
   - Add Gaussian noise \( z(t) \) with \( N_0 = 2 \times 0.02^2 \) to the AM signal.
   - Compute and plot the Fourier Transform \( X_R(f) \) of the received noisy signal.
2. **Output:**
   - Plot of \( X_R(f) \).

---

#### **1.4 Demodulator**
1. **Steps:**
   - Filter the noisy signal \( x_r(t) \) with a bandpass filter (`bandpass`) to obtain \( x_L(t) \).
   - Demodulate \( x_L(t) \) using `amdemod` to reconstruct the original signal \( x_d(t) \).
   - Compare \( x_d(t) \) and \( v_m(t) \) in both time and frequency domains.
2. **Output:**
   - Time domain and frequency domain plots.
   - Listen to the demodulated signal using `sound`.
   - Calculate correlation between \( v_m(t) \) and \( x_d(t) \) using `xcorr`.

---

#### **1.5 Different Noise and Modulations**
1. **Steps:**
   - Repeat steps with \( N_0 = 2 \times 0.1^2 \).
   - Simulate FM modulation and demodulation using `fmmod` and `fmdemod`.
   - Compare correlations for AM and FM demodulated signals.
2. **Output:**
   - Plots of signals and their transforms under different noise levels.
   - Correlation coefficients for AM and FM demodulated signals.

---

### Key MATLAB Implementation Highlights

#### **Data Signal and Fourier Transform**
```matlab
[v_m, fs] = audioread("in-the-air.wav");
v_m = v_m(:,1);
Ts = 1/fs;
N = length(v_m);
t = 0:Ts:(N-1)*Ts;
f = linspace(-fs/2, fs/2, N);
V_m = fftshift(fft(v_m)) / sqrt(N);
```

#### **AM Modulation**
```matlab
fc = 15e3;
kAM = 0.02;
v_AM = ammod(v_m, fc, fs, 0, kAM);
V_AM = fftshift(fft(v_AM)) / sqrt(N);
```

#### **Gaussian Noise Addition**
```matlab
N0 = 2 * (0.02)^2;
z = sqrt(N0/2) * randn(size(v_AM));
x_r = v_AM + z;
X_R = fftshift(fft(x_r)) / sqrt(N);
```

#### **Demodulation**
```matlab
x_L = bandpass(x_r, [fc - 5e3, fc + 5e3], fs);
x_d = amdemod(x_L, fc, fs, 0, kAM);
x_d = lowpass(x_d, 5e3, fs);
corr = xcorr(x_d, v_m, 0, "coeff");
```

---

## Summary
This assignment demonstrates the practical implementation of analog communication principles. By simulating an AM and FM communication system in MATLAB, the effects of noise and modulation techniques on signal quality and recovery are analyzed. Correlation analysis confirms the efficiency of both modulation schemes under different noise conditions.
