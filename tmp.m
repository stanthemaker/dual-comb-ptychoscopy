close all;
clear;
clc;

c = 299792458; 
% --- Comb Parameters ---
frep_1 = 27e3;
frep_2 = 27.2e3;
d_frep     = abs(frep_2 - frep_1);
d_tau = abs(1/frep_2 - 1/frep_1);

n_indices = 1:20;

% --- sampling setup ---
Fs = 1400e3;
dt = 1 / Fs;
N_avg = 100 + 2;             %% number of averages ; +2 for buffer
T_igm = 1 / d_frep;        %% time duration for an entire inteferogram
T = N_avg * T_igm;
Nt = T / dt;
df = 1 / T;
t = (0:Nt-1) * dt;
N_tau = round(1/frep_1 / d_tau);

% --- Create Comb ---
E_c1 = zeros(size(t));
E_c2 = zeros(size(t));

for n = n_indices
    
    fn1 = n * frep_1;
    fn2 = n * frep_2;
    
    E_c1 = E_c1 + exp(1i * 2*pi * fn1 * t);
    E_c2 = E_c2 + exp(1i * 2*pi * fn2 * t);
end
E_c1 = E_c1 / max(abs(E_c1));
E_c2 = E_c2 / max(abs(E_c2));

% --- params incoherent signal ---
f_s = 250e3;
BW = 60e3;
f_axis = (-Nt/2 : Nt/2-1)' * (Fs/Nt);
filter_gauss = exp(-((f_axis - f_s).^2) / (2 * (BW/2.355)^2)); 

noise_raw = (randn(Nt, 1) + 1i*randn(Nt, 1));
Es_fft = fftshift(fft(noise_raw)) .* filter_gauss;
E_s = ifft(ifftshift(Es_fft));
E_s = E_s / max(abs(E_s));
E_s = E_s(:).';
plot_setup(E_c1, E_c2, E_s,Fs);

%% Ground truth
N_igm = round(T_igm * Fs); % Number of elements in one T_ig window
E_s_mat = reshape(E_s, N_igm, N_avg);
tau_c = 12.5e-6;                             % Estimated coherence time
alpha_max = round(4 * tau_c * Fs);

R_ss_sum = zeros(2 * alpha_max + 1, 1);
for k = 1:N_avg

    chunk = E_s_mat(:, k);
    [R_chunk, lags] = xcorr(chunk, alpha_max);
    R_ss_sum = R_ss_sum + R_chunk;
end
R_ss_expected = R_ss_sum / N_avg;
alpha_raw = lags(:) * (1/Fs);


% --- 1. The pwelch Ground Truth ---
% Use N_igm (the length of one chunk) instead of the full length Nt
N_fft = 2^nextpow2(N_igm);
[pxx, f_psd] = pwelch(E_s, rectwin(N_igm), N_igm/2, N_fft, Fs, 'centered');

%--- essentially what pwelch is doing ---%
% N_igm = round(T_igm * Fs); % Number of elements in one T_ig window
% E_s_mat = reshape(E_s, N_igm, N_avg);
% tau_c = 12.5e-6;                             % Estimated coherence time
% alpha_max = round(4 * tau_c * Fs);
% 
% R_ss_sum = zeros(2 * alpha_max + 1, 1);
% for k = 1:N_avg
% 
%     chunk = E_s_mat(:, k);
%     [R_chunk, lags] = xcorr(chunk, alpha_max);
%     R_ss_sum = R_ss_sum + R_chunk;
% end
% R_ss_expected = R_ss_sum / N_avg;
% alpha_raw = lags(:) * (1/Fs);

% --- 1. Extract exactly one interferogram period ---
% This perfectly captures Comb 2 sliding past Comb 1 by one full pulse period
% --- 1. The CPSD Built-in Method ---
% VERY IMPORTANT: The argument order matters!
% cpsd(X, Y) puts the complex conjugate on Y. 
% We want conj(E_c1) * E_c2(shifted), so X = E_c2 and Y = E_c1
N_fft = N_igm; 
[P12_cpsd, f_comb] = cpsd(E_c2, E_c1, rectwin(N_igm), N_igm/2, N_fft, Fs, 'centered');

% Normalize the magnitude
P12_cpsd_norm = abs(P12_cpsd) / max(abs(P12_cpsd));

% --- 2. The Manual xcorr + FFT Method ---
% Extract exactly one interferogram period (one full relative slip of the combs)
E_c1_chunk = E_c1(1 : N_igm);
E_c2_chunk = E_c2(1 : N_igm);

% xcorr(A, B) computes sum(A(shifted) * conj(B)). 
% To get conj(E_c1), A = E_c2_chunk, B = E_c1_chunk
[R12_raw, lags] = xcorr(E_c2_chunk, E_c1_chunk);

% To make the xcorr FFT match the cpsd frequency bins exactly, we need to extract 
% the center N_fft points of the cross-correlation (or zero-pad appropriately).
% Since xcorr is 2*N_igm-1 long, we take the FFT of the shifted sequence:
R12_padded = zeros(N_fft, 1);
center_idx = N_igm; % Center lag of xcorr output
% Wrap the xcorr output into the N_fft window to match the periodic assumption of cpsd
R12_padded(1:N_igm/2) = R12_raw(center_idx : center_idx + N_igm/2 - 1);
R12_padded(end-N_igm/2+1:end) = R12_raw(center_idx - N_igm/2 : center_idx - 1);

% Take the FFT and shift
P12_xcorr_raw = fftshift(fft(R12_padded));
P12_xcorr_norm = abs(P12_xcorr_raw) / max(abs(P12_xcorr_raw));

% --- 3. Plot and Compare ---
figure;
plot(f_comb, 10*log10(P12_cpsd_norm), 'b', 'LineWidth', 2, 'DisplayName', 'cpsd(E_{c2}, E_{c1})');
hold on;
plot(f_comb, 10*log10(P12_xcorr_norm), 'r--', 'LineWidth', 2, 'DisplayName', 'FFT of xcorr');
title('Cross-Spectral Density Comparison');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude (dB)');
legend;
grid on;