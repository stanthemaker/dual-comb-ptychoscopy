close all;
clear;
clc;

c = 299792458; 
lambda_center = 1550e-9;      % Reference Center 
lambda_sig    = 1550.6e-9;      % Incoherent Source Center
width_sig_Hz  = 10e6;         % 10 MHz FWHM
f_center      = c/lambda_center;
f_sig_abs     = c/lambda_sig;
f_sig_rel     = f_sig_abs - f_center; % Relative freq

% --- Comb Parameters ---
frep_1 = 27.2e9;
frep_2 = 27.201e9;
d_frep     = frep_2 - frep_1; % 1 MHz difference
n_lines    = 21;                  % Total lines
n_indices  = -10:10;              % Indices centered at 0

% --- sampling batch setup ---
% to prevent aliasing, 1 / dt > 2 f_max = 2 ( 0.5 frep) 
% NT * Nt *dt = 1 / d_frep = 1ms
dt = 4e-12;
Fs = 1 / dt;
Nt = 1000; % number of time step in a batch
NT = 250; % number of time batch 
T_batch = Nt *dt; % sampling time for a batch
d_fRF = 1 / T_batch;
T = NT * Nt *dt; % total sampling time
t = 0:dt:T-dt;

%% Generate Combs

comb1 = zeros(size(t));
comb2 = zeros(size(t));

for n = n_indices
    fn1 = n * frep_1;
    fn2 = n * frep_2;
    comb1 = comb1 + exp(1i * 2*pi * fn1 * t);
    comb2 = comb2 + exp(1i * 2*pi * fn2 * t);
end

comb1 = comb1 / max(abs(comb1));
comb2 = comb2 / max(abs(comb2));

%% Generate Incoherent Source
noise = (randn(size(t)) + 1i*randn(size(t)));
signal = noise .* exp(1i * 2*pi * f_sig_rel * t);
envelope = signal .* exp(-1i * 2*pi * f_sig_rel * t);
fc_lpf = width_sig_Hz; 
[b, a] = butter(1, fc_lpf/(Fs/2));
envelope_filtered = filter(b, a, envelope);
signal = envelope_filtered .* exp(1i * 2*pi * f_sig_rel * t);
signal = 0.5 * signal / mean(abs(signal)); % scaling

%% Interference

E1 = comb1 + signal;
E2 = comb2 + signal;
plot_beatnotes(comb1,comb2, signal,dt);

I1 = abs(E1).^2;
I2 = abs(E2).^2;

% Remove DC (High-pass filter or subtract mean) to keep only beats
I1 = I1 - mean(I1);
I2 = I2 - mean(I2);

%% create frequency x time spectrograms, each time step is a sampling batch
I1_mat = reshape(I1, Nt, NT);
I2_mat = reshape(I2, Nt, NT);

% FFT each time batch to get spectrograms
F1 = fftshift(fft(I1_mat, [], 1),1);
F2 = fftshift(fft(I2_mat, [], 1),1);


%% Ptychoscopy inversion
% Define RF Frequency Axis
f_RF = (-Nt/2 : Nt/2 - 1) * d_fRF;
batch_index = 1 : NT;


% Optical Axis: Wide enough to see multiple ghosts (-200 to +100 GHz)
f_opt_axis = linspace(-200e9, 200e9, 100000);
spec_recon = zeros(size(f_opt_axis));


for n = n_indices
    Cn = zeros(size(F1(:,1)));
    delta_n = n * d_frep;
    shift_bins = round(delta_n / d_fRF); 

    for j = batch_index
        F1_n = F1(:,j);
        F2_n = circshift(F2(:,j), shift_bins, 1);
        Tj = (j-1) * T_batch;
        phasor = exp(1i * 2*pi * delta_n * Tj);
        Cn  = Cn + F2_n .* conj(F1_n) .* phasor;
    end
    Cn = Cn / NT;
    
    f_comb_line = n * frep_1;
    f_current_grid = f_comb_line + f_RF;
    
   
    spec_recon = spec_recon + interp1(f_current_grid, abs(Cn), f_opt_axis, 'linear', 0);
end

%% Plotting
close all;
figure('Color', 'w', 'Position', [100 100 1000 400]);

% Plot Reconstruction
plot(f_opt_axis/1e9, spec_recon, 'k', 'LineWidth', 1.5); hold on;

% Plot Ground Truth (The Signal)
xline(f_sig_rel/1e9, '-r', 'True Signal', 'LineWidth', 2);

title('Reconstruction result');
xlabel('Optical Frequency (GHz)');
ylabel('Accumulated Magnitude');
xlim([-130, 130]); 
set(gca, "Fontsize",18);
grid on;