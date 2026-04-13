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
center_freq = c / center_wl / 1e12;
wl = data(:,1);
freq = c ./ wl * 1e-3; % Evaluates to THz
p = data(:,2);         % Power in dBm


[freq_sorted, sort_idx] = sort(freq);
p_sorted = p(sort_idx);

[p_peaks_dBm, locs] = findpeaks(p_sorted, 'MinPeakProminence', 3, 'MinPeakHeight', -50);
comb_peaks_THz = freq_sorted(locs);
[~, center_idx] = min(abs(comb_peaks_THz - center_freq));
f0_THz = comb_peaks_THz(center_idx);

p_peaks_linear = 10.^(p_peaks_dBm / 10);
E_amp_comb = sqrt(p_peaks_linear);
E_amp_comb = E_amp_comb / max(E_amp_comb);
N_comblines = length(E_amp_comb);



%%
% Simulation parameters
Fs = 3e12;
dt = 1/Fs;

frep_1 = 27.3e9;
frep_2 = 27.298e9;
aom_shift = 80e6;
d_frep  = frep_1 - frep_2;
T_ifgm  = 1 / abs(d_frep);
k = 10;
T_total = k * T_ifgm;
t = 0:dt: T_total - dt;
df = 1 / T_total;

N = length(t);

E_c1 = zeros(1, N);
E_c2 = zeros(1, N);

for k = 1:N_comblines
    m = k - center_idx;
    freq1 = m * frep_1; 
    freq2 = aom_shift + m * frep_2;
    
    E_c1 = E_c1 + E_amp_comb(k) * exp(1j * 2 * pi * freq1 * t);
    E_c2 = E_c2 + E_amp_comb(k) * exp(1j * 2 * pi * freq2 * t);
end

% ==========================================
% --- Dual comb beatnotes simulation ---
% ==========================================

% I = abs(E_c1 + E_c2).^2 - abs(E_c1).^2 - abs(E_c2).^2;
I = 2 * real(E_c1 .* conj(E_c2));

% 1. Define the frequency axis FIRST
f_axis = (-N/2 : N/2 - 1) * (Fs / N);

% 2. Calculate the FFT of Comb 1 (with a tiny offset to prevent log10(0))
P_c1_linear = abs(fftshift(fft(E_c1))/N);
P_c1_dB = 20*log10(P_c1_linear + 1e-12); 

% 3. Normalize both arrays so their peaks align perfectly at 0 dB
p_norm = p - max(p);
P_c1_norm = P_c1_dB - max(P_c1_dB);

%%
figure('Position', [100, 100, 900, 700]);

% ==========================================
% --- TOP PANEL: Optical Spectrum Overlay ---
% ==========================================
subplot(2,1,1)
% Plot Experimental Data
plot(freq, p_norm, 'b', 'LineWidth', 1.5, 'DisplayName', 'Experimental OSA'); 
hold on;

% Plot Simulated Data (Shifted by f0_THz for perfect alignment)
plot((f_axis/1e12) + f0_THz, P_c1_norm, 'r', 'LineWidth', 1, 'DisplayName', 'Simulated Comb 1');

title('Optical Spectrum: Experimental vs. Simulated');
xlabel('Absolute Frequency (THz)');
ylabel('Normalized Power (dB)');
legend('Location', 'best');
grid on;

% Lock the Y-axis so the -300 dB simulated nulls don't ruin the view
ylim([-60, 5]); 
xlim([min(freq), max(freq)]);


% ==========================================
% --- BOTTOM PANEL: RF Beat Notes ---
% ==========================================
subplot(2,1,2)

f_max = 200e6; % 200 MHz limit
crop_idx = (f_axis >= 0) & (f_axis <= f_max);

% Calculate RF FFT (with offset)
I_fft = abs(fftshift(fft(I)));
I_fft = I_fft / max(I_fft);
I_fft_log = 20*log10(I_fft + 1e-12);

% Plot the cropped RF spectrum
plot(f_axis(crop_idx)/1e6, I_fft_log(crop_idx), 'b', 'LineWidth', 1.5);

title('Dual-Comb RF Beat Note Spectrum');
xlabel('RF Frequency (MHz)');
ylabel('Normalized Amplitude (dB)');
grid on;
ylim([-80, 5]); % Keep the Y-axis clean