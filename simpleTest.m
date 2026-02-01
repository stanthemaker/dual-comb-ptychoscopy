close all; clear; clc;

% --- 1. Setup Parameters ---
Fs = 4000;          % High sample rate to simulate analog world
dt = 1 / Fs;
T = 2;              % Duration (longer = sharper peaks)
t = 0:dt:T-dt;

% Source Signal (Target)
fs_target = 220;    % Signal at 220 Hz

% Comb Parameters
frep1 = 100;        % Comb 1 spacing
frep2 = 110;        % Comb 2 spacing

% Hardware Limit (The Filter you requested)
% We only look at beats < Frep/2. Anything higher is filtered out.
cutoff_freq = 60;   

% --- 2. Generate Signals ---
Es = exp(1i * 2*pi * fs_target * t);

% Create Combs (Sums of discrete lines)
lines1 = 0:10;      % 0, 100, 200... 1000
lines2 = 0:10;      % 0, 130, 260... 1300

comb1 = zeros(size(t));
comb2 = zeros(size(t));

for n = lines1
    comb1 = comb1 + exp(1i * 2*pi * (n*frep1) * t);
end
for n = lines2
    comb2 = comb2 + exp(1i * 2*pi * (n*frep2) * t);
end

% --- 3. Measurement (Intensity) ---
% Assuming pass through BPD
I1_raw = abs(Es + comb1).^2 - abs(comb1).^2 - abs(Es).^2;
I2_raw = abs(Es + comb2).^2 - abs(comb2).^2 - abs(Es).^2;

I1 = lowpass(I1_raw, cutoff_freq, Fs);
I2 = lowpass(I2_raw, cutoff_freq, Fs);

% --- 5. Get RF Spectrum (Positive Only) ---
L = length(t);
f_RF = (0:L/2-1)*(Fs/L); % Positive frequency axis

% FFT and take first half
Spec1 = abs(fft(I1)/L);
Spec1 = Spec1(1:L/2);

Spec2 = abs(fft(I2)/L);
Spec2 = Spec2(1:L/2);

% --- 6. The Reconstruction (Back-Projection) ---
% Scan Optical Frequencies
f_opt_grid = linspace(0, 500, 5000); 

% Calculate EXPECTED beat frequencies for every optical point
% Formula: Beat = | f_opt - Nearest_Comb_Tooth |
beat_pred_1 = min(abs(f_opt_grid - (lines1'*frep1)), [], 1);
beat_pred_2 = min(abs(f_opt_grid - (lines2'*frep2)), [], 1);

% --- 7. Sample the Spectrum ---
% We interpolate the measured spectrum at the predicted beat locations
Recon1 = interp1(f_RF, Spec1, beat_pred_1, 'linear', 0);
Recon2 = interp1(f_RF, Spec2, beat_pred_2, 'linear', 0);

% IMPORTANT: Enforce the Filter logic in Reconstruction
% If the model predicts a beat > cutoff, we shouldn't have seen it.
% Set those points to zero to avoid reading noise.
Recon1(beat_pred_1 > cutoff_freq) = 0;
Recon2(beat_pred_2 > cutoff_freq) = 0;

Product = Recon1 .* Recon2;

% --- Plotting ---
figure('Position', [100, 100, 800, 800]);

subplot(4,1,1);
plot(f_RF, Spec1); xlim([0 150]);
title('Measured RF Spectrum 1 (Filtered)'); xlabel('RF Frequency (Hz)');

subplot(4,1,2);
plot(f_opt_grid, Recon1); 
title(['Recon 1 (Ambiguous) - True Signal and Ghost at ', num2str(200-20), 'Hz']); 
grid on;

subplot(4,1,3);
plot(f_opt_grid, Recon2); 
title('Recon 2 (Ambiguous)'); grid on;

subplot(4,1,4);
plot(f_opt_grid, Product, 'r', 'LineWidth', 2);
xline(fs_target, 'g--');
title('Product (Final Result)'); 
xlabel('Optical Frequency (Hz)'); grid on;