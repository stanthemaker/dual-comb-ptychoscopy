close all;
clear;
clc;

c = 299792458; 
lambda_center = 1550e-9;      
f_center      = c/lambda_center;
f_s = 40e9;
f_s_global = f_center + f_s;

% --- Comb Parameters ---
frep_1 = 27.2e9;
frep_2 = 27.21e9;
d_frep     = frep_2 - frep_1;   % 10 MHz difference
aom_shift = 0;

n_indices  = 1:3;
delta_n_max = abs(d_frep) * max(n_indices) + aom_shift;
delta_n_min = abs(d_frep) * min(n_indices) + aom_shift;

Fs = 200e9;
dt = 1 / Fs;
T_batch = 2e-9; 
Nt = round(T_batch / dt);
d_fRF = 1 / T_batch;

T_ifgm = 1 / abs(d_frep); 
k = 50;                      
T_total = k * T_ifgm;
df = 1/ T_total;
NT = round(T_total / T_batch);

t = 0:dt:T_total-dt;
N_total = length(t);

% --- Create Combs ---
E_c1 = zeros(size(t));
E_c2 = zeros(size(t));
for n = n_indices
    fn1 = n * frep_1;
    fn2 = n * frep_2 + aom_shift;
    E_c1 = E_c1 + exp(1i * 2*pi * fn1 * t);
    E_c2 = E_c2 + exp(1i * 2*pi * fn2 * t);
end
% Normalize
E_c1 = E_c1 / max(abs(E_c1));
E_c2 = E_c2 / max(abs(E_c2));

% --- Create coherent single-wavelength signal ---
% E_s = exp(1i * 2*pi * f_s * t);
% E_s = E_s / max(abs(E_s));

% rect_width = 2e9;
% N_total = length(t); 
% df_exact = Fs / N_total;
% idx_start = round(f_sig_rel / df_exact) + 1;
% idx_end   = round((f_sig_rel + rect_width) / df_exact) + 1;
% 
% S_f = zeros(1, N_total);
% num_bins = idx_end - idx_start + 1;
% S_f(idx_start:idx_end) = exp(1i * 2 * pi * rand(1, num_bins));
% E_s = ifft(S_f);
% E_s = E_s / max(abs(E_s));

% --- Create Incoherent Signal ---

BW = 10e9;
f_axis = (-N_total/2 : N_total/2-1)' * (Fs/N_total);
filter_gauss = exp(-((f_axis - f_s).^2) / (2 * (BW/2.355)^2)); 
noise_raw = (randn(N_total, 1) + 1i*randn(N_total, 1));
Es_fft = fftshift(fft(noise_raw)) .* filter_gauss;
E_s = ifft(ifftshift(Es_fft));
E_s = E_s / max(abs(E_s));
E_s = E_s(:).';

plot_setup(E_c1, E_c2, E_s,Fs,f_center);
%%

I1 = abs(E_s + E_c1).^2 - abs(E_c1).^2 - abs(E_s).^2;
I2 = abs(E_s + E_c2).^2 - abs(E_c2).^2 - abs(E_s).^2;

% I1 = lowpass(I1, cutoff_freq, Fs);
% I2 = lowpass(I2, cutoff_freq, Fs);

I1 = reshape(I1, Nt, NT);
I2 = reshape(I2, Nt, NT);

%% getting expected Cn
% 1. Define resolution based on the FFT padding FIRST
N_pad = 2^nextpow2(abs(1/d_frep/dt));
% N_pad = 4e2;
df_pad = (Fs/N_pad);
f_batch = (-N_pad/2 : N_pad/2-1) * df_pad;

% 2. Initialize global spectrum using the df_pad resolution
max_n = max(abs(n_indices));
max_freq_span = 2 * (max_n + 1) * frep_1; 
total_bins    = ceil(max_freq_span / df_pad);

spec_recon = zeros(total_bins, 1);
center_idx = floor(total_bins / 2) + 1; % Global baseband DC
global_freq_axis = ((1:total_bins).' - center_idx) * df_pad + f_center;

% Bins needed for exactly +/- 0.5 * frep_1
span_bins = round(0.5 * frep_1 / df_pad); 
Cn_center = N_pad/2 + 1; % 0 Hz bin in the padded baseband FFT
Cn_start_idx = Cn_center - span_bins;
Cn_end_idx   = Cn_center + span_bins - 1;

cutoff_freq = 0.6 * frep_1;
cutoff_bins = round(cutoff_freq / df_pad);
freq_mask = zeros(N_pad, 1);
freq_mask(Cn_center - cutoff_bins : Cn_center + cutoff_bins) = 1;

for n = n_indices
    Cn = zeros(N_pad, 1);
    
    % delta_n is defined as f2 - f1
    delta_n = 2 * pi * ((n * d_frep) + aom_shift);
    tau = (0:Nt-1).' * dt;

    for j = 1:NT
        F1 = fftshift(fft(I1(:,j), N_pad)) .* freq_mask;
        
        Tj = (j-1) * T_batch;
        I2_shift = I2(:,j) .* exp(1i* delta_n * tau);
        F2 = fftshift(fft(I2_shift, N_pad)) .* freq_mask;
        Cn  = Cn + F2 .* conj(F1) .* exp(1i* delta_n * Tj);
        % if j == 16
        %     figure
        %     plot(f_batch/1e9 , abs(F1));hold on;
        %     plot(f_batch/1e9 , abs(F2));
        %     pause;
        % end
        % F2 = fftshift(fft(I2(:,j), N_pad)) .* freq_mask;
        % F2_shift = circshift(F2, -shift_bins, 1);
        % Tj = (j-1) * T_batch;
        % phasor = exp(-1i * 2 * pi * delta_n * Tj);
        % Cn  = Cn + F2 .* conj(F1) * phasor;
        
    end    
    Cn = Cn / NT;
    % figure
    % plot(f_batch/1e9 , Cn);
    % pause;
    % close all;

    freq_offset_opt = n * frep_1;
    bin_offset_global = round(freq_offset_opt / df_pad);
    target_center_idx = center_idx + bin_offset_global;
    
    global_start_idx = target_center_idx - span_bins;
    global_end_idx   = target_center_idx + span_bins - 1;
    
    if global_start_idx > 0 && global_end_idx <= total_bins
        spec_recon(global_start_idx:global_end_idx) = Cn(Cn_start_idx:Cn_end_idx);
    else
        warning('Comb line n=%d is outside the spec_recon array bounds.', n);
    end
end

% Plotting
spec_recon_dB = 20*log10(abs(spec_recon));

figure;
plot(global_freq_axis / 1e12, spec_recon_dB, "LineWidth", 1.5);
title('Stitched Dual-Comb Spectrum');
xlabel('Optical Frequency (THz)');
ylabel('Signal Power (dB)')
ylim([max(spec_recon_dB)-60 max(spec_recon_dB)+10]);
xline(f_center / 1e12, '--b', "Baseband Center", 'LabelVerticalAlignment', 'bottom');
xline(f_s_global / 1e12, '--r', "Signal freq", 'LabelVerticalAlignment', 'bottom');
set(gca, "Fontsize", 14)