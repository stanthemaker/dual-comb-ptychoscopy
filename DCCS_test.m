close all;
clear;
clc;

c = 299792458; 
lambda_center = 1550e-9;      % Reference Center 
lambda_sig    = 1550.6e-9;      % Incoherent Source Center
width_sig_Hz  = 10e6;         % 10 MHz FWHM
f_center      = c/lambda_center;
f_sig     = c/lambda_sig;
f_sig_rel = f_sig - f_center;

% --- Comb Parameters ---
frep_1 = 20.2e9;
frep_2 = 20.3e9;
d_frep     = frep_2 - frep_1; % 1 MHz difference
N_combline    = 11;               % Number of comb lines
n_indices  = -5:5;              % Indices centered at 0

% --- sampling batch setup ---
% to prevent aliasing, 1 / dt > 2 f_max = 2 ( 0.5 frep) 
% NT * Nt *dt = 1 / d_frep = 1ms
dt = 4e-12;
Fs = 1 / dt;
Nt = 10000; % number of time step in a batch
NT = 250; % number of time batch 
T_batch = Nt *dt; % sampling time for a batch
d_fRF = 1 / T_batch; %spectral resolution
T = NT * Nt *dt; % total sampling time
t = 0:dt:T-dt;

% --- Create Combs ---
comb1 = zeros(size(t));
comb2 = zeros(size(t));
for n = n_indices
    fn1 = n * frep_1;
    fn2 = n * frep_2;
    comb1 = comb1 + exp(1i * 2*pi * fn1 * t);
    comb2 = comb2 + exp(1i * 2*pi * fn2 * t);
end
% Normalize
E_c1 = comb1 / max(abs(comb1));
E_c2 = comb2 / max(abs(comb2));

% --- Create Incoherent Signal ---
noise = (randn(size(t)) + 1i*randn(size(t)));
envelope = noise; 
[b, a] = butter(1, width_sig_Hz/(Fs/2));
envelope_filtered = filter(b, a, envelope);
signal = envelope_filtered .* exp(1i * 2*pi * f_sig_rel * t);
E_s = 2 * signal / mean(abs(signal)); 

S1 = conj(E_c1) .* E_s;
S2 = conj(E_c2) .* E_s;

S1 = reshape(S1, Nt, NT);
S2 = reshape(S2, Nt, NT);

f_RF = (-Nt/2 : Nt/2 - 1) * d_fRF;
omega_k = 2 * pi * f_RF;

ti = (0:Nt-1).' * dt;
F1 = (1/Nt) * (exp(-1i * omega_k.' * ti.') * S1);
F2 = (1/Nt) * (exp(-1i * omega_k.' * ti.') * S2);
%%
max_freq_span = (max(n_indices) - min(n_indices)) * frep_1 + 2*Fs; % buffer
total_bins    = ceil(max_freq_span / d_fRF);
N_bins_linespan = ceil(frep_1 / d_fRF);
% Initialize global spectrum
spec_recon = zeros(total_bins, 1);
center_idx = floor(total_bins / 2);
global_freq_axis = ((1:total_bins).' - center_idx) * d_fRF + f_center;
%%

for n = n_indices
    Cn = zeros(Nt, 1);

    delta_n = n * d_frep;
    shift_bins = round(delta_n / d_fRF);

    for j = 1:NT
        F1_n = F1(:,j);
        F2_n = circshift(F2(:,j), shift_bins, 1);
        Tj = (j-1) * T_batch;
        phasor = exp(1i * 2*pi * delta_n * Tj);
        Cn  = Cn + F2_n .* conj(F1_n) .* phasor;
    end
    Cn = Cn / NT;

   
    freq_offset = n * frep_1;
    bin_offset_global = round(freq_offset / d_fRF);
    target_center_idx = center_idx + bin_offset_global;
    global_start_idx = target_center_idx - floor(N_bins_linespan/2);
    global_end_idx   = global_start_idx + N_bins_linespan - 1;

    Cn_start_idx = Nt/2 - floor(N_bins_linespan/2);
    Cn_end_idx   = Cn_start_idx + N_bins_linespan - 1;
    
    if global_start_idx > 0 && global_end_idx <= total_bins
        spec_recon(global_start_idx:global_end_idx) = spec_recon(global_start_idx:global_end_idx) + Cn(Cn_start_idx: Cn_end_idx);
    else
        warning('Comb line n=%d is outside the spec_recon array bounds.', n);
    end
end
%%
figure;
plot(global_freq_axis / 1e12, 20*log10(abs(spec_recon)));
title('Stitched Dual-Comb Spectrum');
xlabel('Optical Frequency (THz)');
% ylim([-120 -40])
xline(f_center / 1e12,'--b', "Center freq");
xline(f_sig / 1e12,'--r', "signal freq");
set(gca,"Fontsize",18)