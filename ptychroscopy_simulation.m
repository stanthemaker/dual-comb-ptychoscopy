close all;
clear;
clc;

opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [30, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "Var2"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable("./0403_AOM/15528-2715.CSV", opts);
data = table2array(data);

c = 299792458;
center_wl = 1552.8e-9;
center_freq_THz = c / center_wl / 1e12;
wl = data(:,1);
freq = c ./ wl * 1e-3; % Evaluates to THz
p = data(:,2);         % Power in dBm


[freq_sorted, sort_idx] = sort(freq);
p_sorted = p(sort_idx);

[p_peaks_dBm, locs] = findpeaks(p_sorted, 'MinPeakProminence', 3, 'MinPeakHeight', -50);
comb_peaks_THz = freq_sorted(locs);
[~, center_idx_comb] = min(abs(comb_peaks_THz - center_freq_THz));
f0_THz = comb_peaks_THz(center_idx_comb);

p_peaks_linear = 10.^(p_peaks_dBm / 10);
E_amp_comb = sqrt(p_peaks_linear);
E_amp_comb = E_amp_comb / max(E_amp_comb);
N_comblines = length(E_amp_comb);
frep_1 = 27.3e9;
frep_2 = 27.298e9;
aom_shift = 80e6;

% Simulation parameters
Fs = 2e12; % frep_1 * N_comblines/2 
dt = 1/Fs;

d_frep  = frep_2 - frep_1;
T_ifgm  = 1 / abs(d_frep);
T_batch = 1e-9; % T_batch <= 1 / 20 / n / ∆frep
Nt = round(T_batch / dt);
rotation_factor = 2 * N_comblines /2 * T_batch * abs(d_frep); %

df_batch = 1 / T_batch;
k = 5;
T_total = k * T_ifgm;
t = 0:dt: T_total - dt;
N = length(t);
NT = round(T_total / T_batch);

%% Generate Combs
n_indices = (1:N_comblines) - center_idx_comb;
E_c1 = zeros(1, N);
E_c2 = zeros(1, N);

for n = n_indices
    freq1 = n * frep_1; 
    freq2 = aom_shift + n * frep_2;
    
    E_c1 = E_c1 + E_amp_comb(n+center_idx_comb) * exp(1j * 2 * pi * freq1 * t);
    E_c2 = E_c2 + E_amp_comb(n+center_idx_comb) * exp(1j * 2 * pi * freq2 * t);
end

%% Generate Incoherent Source

BW = 100e9;
gaussian_center = -80e9; % from center_freq
f_axis = (-N/2 : N/2-1)' * (Fs/N);
filter_gauss = exp(-((f_axis - gaussian_center).^2) / (2 * (BW/2.355)^2)); 
noise_raw = (randn(N, 1) + 1i*randn(N, 1));
Es_fft = fftshift(fft(noise_raw)) .* filter_gauss;
E_s = ifft(ifftshift(Es_fft));
E_s = E_s / max(abs(E_s));
E_s = E_s(:).';
[S_ss_dB, f_Es_PSD] = getPSD(E_s, Fs, 10e9);
f_Es_PSD = f_Es_PSD + center_freq_THz*1e12;
S_ss_dB = S_ss_dB - max(S_ss_dB);

plot_setup(E_c1, E_c2, E_s,Fs,center_freq_THz * 1e12);
% data = readtable("./0320_dualbeatnote/EDFA_spectrum", opts);
% data = table2array(data);
% 
% wl = data(:,1);
% freq = c ./ wl * 1e-3;
% p = data(:,2);   
% [freq, sort_idx] = sort(freq);
% p = p(sort_idx);
% p_linear = 10.^(p / 10);
% E_amp = sqrt(p_linear);
% E_amp = E_amp/ max(E_amp);
% 
% f_axis_Hz = (-N/2 : N/2 - 1) * (Fs / N);              
% f_axis_THz = (f_axis_Hz / 1e12) + center_freq_THz; 
% H_f = interp1(freq, E_amp, f_axis_THz, 'linear', 0);
% 
% noise_raw = (randn(1, N) + 1i * randn(1, N)) / sqrt(2);
% Es_fft = fftshift(fft(noise_raw)) .* (H_f(:).');
% E_s = ifft(ifftshift(Es_fft));
% E_s = E_s / max(abs(E_s));
% 
% [S_xx_dB, f_axis_Hz] = getPSD(E_s, Fs, 10e9);
% S_xx_dB = S_xx_dB - max(S_xx_dB);
% f_THz = (f_axis_Hz / 1e12) + center_freq_THz;
% plot(f_THz, S_xx_dB, 'y', 'LineWidth', 1.5, 'DisplayName', 'Simulated PSD (Autocorrelation)');


% plot(f_axis_THz , fftshift(fft(E_s))); hold on
% plot(f_axis_THz , H_f);
% 
% plot_setup(E_c1, E_c2, E_s,Fs);
%%

I1 = abs(E_s + E_c1).^2 - abs(E_c1).^2 - abs(E_s).^2;
I2 = abs(E_s + E_c2).^2 - abs(E_c2).^2 - abs(E_s).^2;
I1 = I1(1:Nt * NT);
I2 = I2(1:Nt * NT);
I1 = reshape(I1, Nt, NT);
I2 = reshape(I2, Nt, NT);

%%
N_pad = 2^nextpow2(Nt);
df_pad = (Fs/N_pad);
f_batch = (-N_pad/2 : N_pad/2-1) * df_pad;

% 2. Initialize global spectrum using the df_pad resolution

max_n = max(abs(n_indices));
freq_range_global = 2 * (max_n + 1) * frep_1; 
total_bins    = ceil(freq_range_global / df_pad);
spec_recon = zeros(total_bins, 1);
center_idx_global = floor(total_bins / 2) + 1;

% Global freq axis in hz
global_freq_axis = ((1:total_bins).' - center_idx_global) * df_pad  + f0_THz*1e12;
S_ss_global = interp1(f_Es_PSD, S_ss_dB, global_freq_axis, 'linear', 0);

% For reconstructing info from n-comb 
span_bins = round(0.5 * frep_1 / df_pad);
Cn_center = N_pad/2 + 1;
Cn_start_idx = Cn_center - span_bins;
Cn_end_idx   = Cn_center + span_bins - 1;

cutoff_freq = 0.6 * frep_1;
cutoff_bins = round(cutoff_freq / df_pad);
freq_mask = zeros(N_pad, 1);
freq_mask(Cn_center - cutoff_bins : Cn_center + cutoff_bins) = 1;


for n = -10:10
    Cn = zeros(N_pad, 1);
    
    % delta_n is defined as f2 - f1
    delta_n = 2 * pi * ((n * d_frep) + aom_shift);
    tau = (0:Nt-1).' * dt;

    for j = 1:NT
        % I1 = ;
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
        % figure
        % plot(f_batch, F2);
        % pause;
        
    end
    freq_offset_opt = n * frep_1;

    %% regularization
    comb_power = E_amp_comb(n+center_idx_comb)^2;

    % % Tikhonov regularization
    % Thv_reg = 0.05 * max(E_amp_comb.^2); 
    % Cn = Cn / NT / (comb_power + Thv_reg) ;

    % Wiener regularization

    alpha = 0.05;         % Default: 1.0. Lower it (e.g., 0.1) to make it LESS aggressive.
    SNR_peak_dB = 10;    % Assumed noise floor. Raise it (e.g., 15) to make it LESS aggressive.
    min_signal = 0.02;   % Prevents the penalty from blowing up at the far edges.
    
    local_signal_dB = interp1(global_freq_axis, S_ss_global, f0_THz*1e12 + freq_offset_opt, 'linear', -100);
    P_signal = 10^(local_signal_dB / 10);
    P_noise = 10^(-SNR_peak_dB / 10);
    wiener_reg = alpha * P_noise / (P_signal + 1e-12);

    Cn = Cn / NT / (comb_power + wiener_reg) ;
    % Cn = Cn / NT ;    
    
    % figure
    % plot(f_batch/1e9 , Cn); hold on
    
    % pause;
    % close all;

    
    bin_offset_global = round(freq_offset_opt / df_pad);
    target_center_idx = center_idx_global + bin_offset_global;
    
    global_start_idx = target_center_idx - span_bins;
    global_end_idx   = target_center_idx + span_bins - 1;
    
    if global_start_idx > 0 && global_end_idx <= total_bins
        spec_recon(global_start_idx:global_end_idx) = Cn(Cn_start_idx:Cn_end_idx);
    else
        warning('Comb line n=%d is outside the spec_recon array bounds.', n);
    end

end
spec_recon_dB = 10*log10(abs(spec_recon));
spec_recon_dB = spec_recon_dB - max(spec_recon_dB);
disp("Comb number: " + n)
figure;
plot(global_freq_axis / 1e12, spec_recon_dB, "LineWidth", 1.5); hold on
plot(global_freq_axis /1e12, S_ss_global, "LineWidth", 1.5)
title('Stitched Dual-Comb Spectrum');
xlabel('Optical Frequency (THz)');
ylabel('Signal Power (dB)')
ylim([max(spec_recon_dB)-60 max(spec_recon_dB)+10]);
set(gca, "Fontsize", 14)

% Plotting
% spec_recon_dB = 20*log10(abs(spec_recon));
% spec_recon_dB = spec_recon_dB - max(spec_recon_dB);
% 
% figure;
% plot(global_freq_axis / 1e12, spec_recon_dB, "LineWidth", 1.5); hold on
% plot(global_freq_axis /1e12, S_ss_global, "LineWidth", 1.5)
% title('Stitched Dual-Comb Spectrum');
% xlabel('Optical Frequency (THz)');
% ylabel('Signal Power (dB)')
% ylim([max(spec_recon_dB)-60 max(spec_recon_dB)+10]);
% set(gca, "Fontsize", 14)