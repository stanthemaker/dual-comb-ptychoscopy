close all;
clear;
clc;

c = 299792458; 
lambda_center = 1550e-9;      
f_center      = c/lambda_center;

% --- Comb Parameters ---
frep_1 = 27.2e9;
frep_2 = 27.3e9;
d_frep     = frep_2 - frep_1;   % 100 MHz difference
comb2_shift = 80e6;
n_indices  = -3:3;              % Indices centered at 0

% --- sampling batch setup ---
% to prevent aliasing, 1 / dt > 2 f_max = 2 ( 0.5 frep) 
% NT * Nt *dt > 1 / d_frep = 1ms

Fs = 200e9;
dt = 1 / Fs;
Nt = 800; % number of time step in a batch
NT = 250; % number of time batch 
T_batch = Nt * dt; % sampling time for a batch
T = NT * T_batch;
t = 0:dt:T-dt;

df_batch = 1 / T_batch; %spectral resolution
df = 1 / T;

% --- Create Combs ---
E_c1 = zeros(size(t));
E_c2 = zeros(size(t));
% E_s = zeros(size(t));
for n = n_indices
    fn1 = n * frep_1;
    fn2 = n * frep_2 + comb2_shift;
    E_c1 = E_c1 + exp(1i * 2*pi * fn1 * t);
    E_c2 = E_c2 + exp(1i * 2*pi * fn2 * t);
end
% Normalize
E_c1 = E_c1 / max(abs(E_c1));
E_c2 = E_c2 / max(abs(E_c2));

% --- Create Incoherent Signal ---
%% noisy single wavelength with linewdith
% width_sig_Hz  = 100e6;           % 10 MHz FWHM
% noise = (randn(size(t)) + 1i*randn(size(t)));
% envelope = noise; 
% [b, a] = butter(1, width_sig_Hz/(Fs/2));
% envelope_filtered = filter(b, a, envelope);
% signal = envelope_filtered .* exp(1i * 2*pi * f_sig_rel * t);
% signal = exp(1i * 2*pi * f_sig_rel * t);
% E_s = 2 * signal / mean(abs(signal)); 
%% Hann
% f_sig_rel = 40.6e9 ;
% % signal = exp(1i * 2*pi * f_sig_rel * t);
% f_sig = f_center + f_sig_rel;
% BW = 30e9;
% tau = 1 / BW; 
% t_center = t(end)/2; 
% envelope = exp(-( (t - t_center).^2 ) / (2 * tau^2));
% signal = envelope .* exp(1i * 2*pi * f_sig_rel * (t - t_center));
% E_s = 1e5 * signal / max(abs(signal));
data = load('BroadbandSignal.mat');
E_s = data.E_s;
f_sig_rel = data.f_sig_rel;
f_sig = f_center + f_sig_rel;

plot_setup(E_c1,E_c2,E_s,Fs);
%%
I1 = abs(E_s + E_c1).^2 - abs(E_c1).^2 - abs(E_s).^2;
I2 = abs(E_s + E_c2).^2 - abs(E_c2).^2 - abs(E_s).^2;

cutoff_freq = frep_1 /2;
I1 = lowpass(I1, cutoff_freq, Fs);
I2 = lowpass(I2, cutoff_freq, Fs);

I1 = reshape(I1, Nt, NT);
I2 = reshape(I2, Nt, NT);

F1 = fftshift(fft(I1,[],1),1);
F2 = fftshift(fft(I2,[],1),1);

freq_axis = (-Nt/2 : Nt/2-1) * (Fs / Nt);
mask = abs(freq_axis) <= cutoff_freq;
F1 = F1(mask, :);
F2 = F2(mask, :);
%% paper implementation ( need to know the complex field
% S1 = conj(E_c1) .* E_s + E_c1 .* conj(E_s);
% S2 = conj(E_c2) .* E_s + E_c2 .* conj(E_s);
% 
% f_cutoff = frep_1 / 2; 
% [b_lp, a_lp] = butter(6, f_cutoff / (Fs/2));
% S1 = filtfilt(b_lp, a_lp, S1);
% S2 = filtfilt(b_lp, a_lp, S2);
% 
% plot_complex(S1,Fs);

% 
% S1 = reshape(S1, Nt, NT);
% S2 = reshape(S2, Nt, NT);
% 
% f_RF = (-Nt/2: Nt/2 - 1) * d_fRF;
% omega_k = 2 * pi * f_RF;
% 
% ti = (0:Nt-1).' * dt;
% F1 = (1/Nt) * (exp(-1i * omega_k.' * ti.') * S1);
% F2 = (1/Nt) * (exp(-1i * omega_k.' * ti.') * S2);

%%
max_freq_span = (length(n_indices)+2) * frep_1; % buffer
total_bins    = ceil(max_freq_span / df_batch);


% Initialize global spectrum
spec_recon = zeros(total_bins, 1);
center_idx = floor(total_bins / 2);
global_freq_axis = ((1:total_bins).' - center_idx) * df_batch + f_center;

f = (-size(F1,1)/2:size(F1,1)/2-1) * (Fs/Nt);
for n = 1
    Cn = zeros(size(F1,1), 1);

    delta_n = (n * d_frep) + comb2_shift;
    shift_bins = round(delta_n / df_batch);

    for j = 1:NT
        F1_n = F1(:,j);
        F2_n_shifted = circshift(F2(:,j), shift_bins, 1);
        Tj = (j-1) * T_batch;
        phasor = exp(1i * 2 * pi * delta_n * Tj);
        Cn  = Cn + F2_n_shifted .* conj(F1_n) .* phasor;
    
        % plot(f/1e9 , abs(F1_n));
        % hold on;
        % plot(f/1e9 , abs(F2(:,j)));
        % plot(f/1e9 , abs(F2_n_shifted));
    end
    Cn = Cn / NT;

   
    freq_offset = n * frep_1;
    bin_offset_global = round(freq_offset / df_batch);
    target_center_idx = center_idx + bin_offset_global;

    % Recalculate global indices
    Cn_L = length(Cn);
    global_start_idx = target_center_idx - floor(Cn_L/2);
    global_end_idx   = global_start_idx + Cn_L - 1;
    
    if global_start_idx > 0 && global_end_idx <= total_bins
        spec_recon(global_start_idx:global_end_idx) = spec_recon(global_start_idx:global_end_idx) + Cn;
    else
        warning('Comb line n=%d is outside the spec_recon array bounds.', n);
    end
end

getError(spec_recon , E_s, Fs, global_freq_axis, f_center);

spec_recon = 20*log10(abs(spec_recon));
figure;
plot(global_freq_axis / 1e12, spec_recon, "LineWidth",1.5);
title('Stitched Dual-Comb Spectrum');
xlabel('Optical Frequency (THz)');
ylabel('Signal Power (dB)')
ylim([max(spec_recon)-60 max(spec_recon)+10]);
xline(f_center / 1e12,'--b', "Center freq");
xline(f_sig / 1e12,'--r', "Signal freq");
set(gca,"Fontsize",18)
