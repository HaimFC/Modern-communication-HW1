clc;
close all;

%Q.1 - Data Generation
%%

[v_m,fs] = audioread("in-the-air.wav");
%sound(v_m, fs);
v_m = v_m(:,1);
Ts = 1/fs;
N=length(v_m);
t = 0:Ts:(N-1)*Ts;
f = linspace(-fs/2,fs/2,N);
V_m = fftshift(fft(v_m)) / sqrt(N);


% Plot the data signal and its Fourier transform
figure;

% Plot the data signal vm(t)
subplot(2, 1, 1);
plot(t, v_m);
title('Data Signal v_m(t)');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform Vm(f)
subplot(2, 1, 2);
plot(f, abs(V_m));
title('Fourier Transform V_m(f)');
xlabel('Frequency');
ylabel('Amplitude');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);


% Evaluate the bandwidth of the information signal
BW = powerbw(V_m, fs);
fprintf('BW of the information signal: %.4f [kHz]\n', BW(1));

%Q.2 - Modulator
%%

% Choose the carrier frequency fc and modulation index kAM
fc = 15 * 10^3; % Carrier frequency in Hz
kAM = 0.02; % Modulation index

% Create the AM modulated signal v_AM(t)
v_AM = ammod(v_m, fc, fs, 0, kAM);

% Compute the Fourier transform of v_AM(t)
V_AM = fftshift(fft(v_AM)) / sqrt(N);

% Plot the Fourier transform V_AM(f)
figure;
plot(f, abs(V_AM));
title('Fourier Transform V_{AM}(f)');
xlabel('Frequency');
ylabel('Amplitude');

% Listen to the modulated signal v_AM(t)
sound(v_AM, fs);

%Q.3 - Channel
%%
N0 = 2*(0.02)^2;
z = sqrt(N0/2) * randn(1, N).';
x_r = v_AM + z;

% Fourier transform of x_r
X_R = fftshift(fft(x_r)) / sqrt(N);


figure;
% Plot the data signal xr(t)
subplot(2, 1, 1);
plot(t, x_r);
title('Data Signal x_r(t)');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform Xr(f)
subplot(2, 1, 2);
plot(f, abs(X_R));
title('Fourier Transform X_R(f)');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);

%Q.4 - Demodulator
%%
fm = 5*10^3;
x_L = bandpass(x_r,[fc-fm fc+fm],fs);

% Compute the Fourier transform of Lxt
X_L = fftshift(fft(x_L))/sqrt(N);

figure(1);
% Plot the data signal xl(t)
subplot(2, 1, 1);
plot(t, x_L);
title('Data Signal x_L(t)');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform XL(f)
subplot(2, 1, 2);
plot(f, abs(X_L));
title('Fourier Transform X_L(f)');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);


% Demodulate the received signal x_L(t)
x_d = amdemod(x_L, fc, fs, 0, kAM); 
x_d = lowpass(x_d,fm,fs);
X_d = fftshift(fft(x_d))/sqrt(N); 

figure(2);
% Plot the data signal xd and original vm(t)
subplot(2,1,1);
plot(t,x_d)
hold on
plot(t,v_m)
title('v_m(t) and x_d(t) - High Noise')
legend('x_d(t)','v_m(t)')
xlabel('Time'); ylabel('Amplitude'); grid on

% Plot the Fourier transform Xd and Vm(f)
subplot(2,1,2);
plot(f,abs(X_d))
hold on
plot(f,abs(V_m))
title("|V_m(f)| and |X_d(f)|- Noise")
legend('|X_d(f)|', '|V_m(f)|')
xlabel('Frequency'); ylabel('|Signal(f)|'); grid on

% Listen to the decoded signal x_d(t)
sound(x_d, fs);

%Correlation calculation
corr = xcorr(x_d, v_m, 0, "coeff");
fprintf('The correlation between vm(t) and xd(t) is: %.4f\n', corr);

%Q.5 (1.5.1)
%%

N0 = 2*(0.1)^2;
z = sqrt(N0/2) * randn(1, N).';
x_r = v_AM + z;

% Fourier transform of x_r
X_R = fftshift(fft(x_r)) / sqrt(N);


figure;
% Plot the data signal xr(t)
subplot(2, 1, 1);
plot(t, x_r);
title('Data Signal x_r(t)');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform Xr(f)
subplot(2, 1, 2);
plot(f, abs(X_R));
title('Fourier Transform X_R(f)');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);

fm = 5*10^3;
x_L = bandpass(x_r,[fc-fm fc+fm],fs);

% Compute the Fourier transform of Lxt
X_L = fftshift(fft(x_L))/sqrt(N);

figure(1);
% Plot the data signal xl(t)
subplot(2, 1, 1);
plot(t, x_L);
title('Data Signal x_L(t)');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform XL(f)
subplot(2, 1, 2);
plot(f, abs(X_L));
title('Fourier Transform X_L(f)');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);


% Demodulate the received signal x_L(t)
x_d = amdemod(x_L, fc, fs, 0, kAM); 
x_d = lowpass(x_d,fm,fs);
X_d = fftshift(fft(x_d))/sqrt(N); 

figure(2);
% Plot the data signal xd and original vm(t)
subplot(2,1,1);
plot(t,x_d)
hold on
plot(t,v_m)
title('v_m(t) and x_d(t) - High Noise')
legend('x_d(t)','v_m(t)')
xlabel('Time'); ylabel('Amplitude'); grid on

% Plot the Fourier transform Xd and Vm(f)
subplot(2,1,2);
plot(f,abs(X_d))
hold on
plot(f,abs(V_m))
title("|V_m(f)| and |X_d(f)|- Noise")
legend('|X_d(f)|', '|V_m(f)|')
xlabel('Frequency'); ylabel('|Signal(f)|'); grid on

% Listen to the decoded signal x_d(t)
sound(x_d, fs);

%Correlation calculation
corr = xcorr(x_d, v_m, 0, "coeff");
fprintf('The correlation between vm(t) and xd(t) is: %.4f\n', corr);

%Q.5 (1.5.2)
%%  ---------------------Noise of N0 = 0.02--------------------------------

fd = 10*10^3;
v_FM = fmmod(v_m,fc,fs,fd); 

%-------Channel--------

N0 = 2*(0.02)^2;
z = sqrt(N0/2) * randn(1, N).';
x_r = v_FM + z;

% Fourier transform of x_r
X_R = fftshift(fft(x_r)) / sqrt(N);

figure(1);
% Plot the data signal xr(t)
subplot(2, 1, 1);
plot(t, x_r);
title('Data Signal x_r(t) - Noise of 0.02');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform Xr(f)
subplot(2, 1, 2);
plot(f, abs(X_R));
title('Fourier Transform X_R(f) - Noise of 0.02');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);

fm = 5*10^3;
x_L = bandpass(x_r,[fc-fm fc+fm],fs);

% Compute the Fourier transform of Lxt
X_L = fftshift(fft(x_L))/sqrt(N);

figure(2);
% Plot the data signal xl(t)
subplot(2, 1, 1);
plot(t, x_L);
title('Data Signal x_L(t) - Noise of 0.02');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform XL(f)
subplot(2, 1, 2);
plot(f, abs(X_L));
title('Fourier Transform X_L(f) - Noise of 0.02');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);

%-------Demodulator--------

% Demodulate the received signal x_L(t)
x_d = fmdemod(x_L, fc, fs, fd); 
x_d = lowpass(x_d,fm,fs);
X_d = fftshift(fft(x_d))/sqrt(N); 

figure(3);
% Plot the data signal xd and original vm(t)
subplot(2,1,1);
plot(t,x_d)
hold on
plot(t,v_m)
title('v_m(t) and x_d(t) - High Noise(0.02)')
legend('x_d(t)','v_m(t)')
xlabel('Time'); ylabel('Amplitude'); grid on

% Plot the Fourier transform Xd and Vm(f)
subplot(2,1,2);
plot(f,abs(X_d))
hold on
plot(f,abs(V_m))
title("|V_m(f)| and |X_d(f)|- Noise (0.02)")
legend('|X_d(f)|', '|V_m(f)|')
xlabel('Frequency'); ylabel('|Signal(f)|'); grid on

% Listen to the decoded signal x_d(t)
sound(x_d, fs);

%Correlation calculation
corr = xcorr(x_d, v_m, 0, "coeff");
fprintf('The correlation between vm(t) and xd(t) is: %.4f\n', corr);

%------------------------Noise of N0 = 0.1---------------------------------

fd = 10*10^3;
v_FM = fmmod(v_m,fc,fs,fd); 

%-------Channel--------

N0 = 2*(0.1)^2;
z = sqrt(N0/2) * randn(1, N).';
x_r = v_FM + z;

% Fourier transform of x_r
X_R = fftshift(fft(x_r)) / sqrt(N);

figure(4);
% Plot the data signal xr(t)
subplot(2, 1, 1);
plot(t, x_r);
title('Data Signal x_r(t) - Noise of 0.1 ');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform Xr(f)
subplot(2, 1, 2);
plot(f, abs(X_R));
title('Fourier Transform X_R(f) - Noise of 0.1');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);

fm = 5*10^3;
x_L = bandpass(x_r,[fc-fm fc+fm],fs);

% Compute the Fourier transform of Lxt
X_L = fftshift(fft(x_L))/sqrt(N);

figure(5);
% Plot the data signal xl(t)
subplot(2, 1, 1);
plot(t, x_L);
title('Data Signal x_L(t) - Noise of 0.1');
xlabel('Time');
ylabel('Amplitude');

% Plot the Fourier transform XL(f)
subplot(2, 1, 2);
plot(f, abs(X_L));
title('Fourier Transform X_L(f) - Noise of 0.1');
xlabel('Frequency');
ylabel('|Signal(f)|');

% Adjust the subplot positions for better visibility
subplot(2, 1, 1);
pos = get(gca, 'Position');
pos(2) = pos(2) + 0.08;
pos(4) = pos(4) - 0.08;
set(gca, 'Position', pos);

%-------Demodulator--------

% Demodulate the received signal x_L(t)
x_d = fmdemod(x_L, fc, fs, fd); 
x_d = lowpass(x_d,fm,fs);
X_d = fftshift(fft(x_d))/sqrt(N); 

figure(6);
% Plot the data signal xd and original vm(t)
subplot(2,1,1);
plot(t,x_d)
hold on
plot(t,v_m)
title('v_m(t) and x_d(t) - High Noise (0.1)')
legend('x_d(t)','v_m(t)')
xlabel('Time'); ylabel('Amplitude'); grid on

% Plot the Fourier transform Xd and Vm(f)
subplot(2,1,2);
plot(f,abs(X_d))
hold on
plot(f,abs(V_m))
title("|V_m(f)| and |X_d(f)|- Noise (0.1)")
legend('|X_d(f)|', '|V_m(f)|')
xlabel('Frequency'); ylabel('|Signal(f)|'); grid on

% Listen to the decoded signal x_d(t)
sound(x_d, fs);

%Correlation calculation
corr = xcorr(x_d, v_m, 0, "coeff");
fprintf('The correlation between vm(t) and xd(t) is: %.4f\n', corr);

