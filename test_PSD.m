close all; clear; clc;

% --- Parameters ---
Fs = 50e6;              
dt = 1/Fs;
Nt = 2000;              
NT = 500;               
N_total = Nt * NT;
t = (0:N_total-1)' * dt;

% --- Signal Generation (Filtered Noise) ---
f_s1 = 10e6;
BW = 3e6;
noise_raw = (randn(N_total, 1) + 1i*randn(N_total, 1));

f_axis = (-N_total/2 : N_total/2-1)' * (Fs/N_total);
filter_gauss = exp(-((f_axis - f_s1).^2) / (2 * (BW/2.355)^2)); 

Es_fft_full = fftshift(fft(noise_raw)) .* filter_gauss;
Es = ifft(ifftshift(Es_fft_full)); 

% --- Analysis ---
% 1. Single Batch FFT (Normalized for comparison)
Es_batch = Es(1:Nt); 
fft_batch = fftshift(fft(Es_batch)) / Nt; 

% 2. PSD (Welch Method) - Using 'psd' to get Power/Hz
[pxx, f_psd] = pwelch(Es, Nt, Nt/2, Nt, Fs, 'centered');

% 3. Autocorrelation (The "hidden" step in PSD)
[autocorr, lags] = xcorr(Es(1:Nt), 'biased');
tau = lags * dt;

% --- Plotting ---
figure('Position', [50, 50, 1200, 400]);

% Plot 1: Raw FFT (Very Jagged)
subplot(1,3,1);
plot(f_axis(1:100:end)/1e6, 20*log10(abs(Es_fft_full(1:100:end)/N_total)), 'b');
title('Full Signal FFT (Raw)');
xlabel('MHz'); ylabel('dB'); grid on;
axis([0 20 -120 -40]); % Forced axis to see the signal

% Plot 2: PSD (Smooth Envelope)
subplot(1,3,2);
plot(f_psd/1e6, 10*log10(pxx), 'r', 'LineWidth', 1.5);
title('PSD (Averaged)');
xlabel('MHz'); ylabel('dB/Hz'); grid on;
axis([0 20 -120 -40]); % Match the FFT axis

% Plot 3: Autocorrelation (The "Bridge")
subplot(1,3,3);
plot(tau*1e6, abs(autocorr), 'k');
title('Autocorrelation Magnitude');
xlabel('Lag (\mu s)'); ylabel('|R(\tau)|'); grid on;