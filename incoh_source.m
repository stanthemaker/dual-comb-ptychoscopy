% --- Your Original Initialization ---
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [30, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "Var2"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable("./0320_dualbeatnote/EDFA_spectrum", opts);
data = table2array(data);
clear opts

c = 299792458;
N = 2^22;                 % Number of sampling points (use power of 2 for fast FFT)
dt = 0.1e-12;             % Time step (e.g., 0.1 picoseconds)
Fs = 1 / dt;              % Sampling frequency
t = (0:N-1) * dt;         % Time vector
df = Fs / N;    
center_wl = 1552.8e-9;
center_freq_THz = c / center_wl / 1e12;

wl = data(:,1);
freq = c ./ wl * 1e-3;    % Evaluates to THz (assuming wl is in nm)
p = data(:,2);


%%
% ==========================================
% --- COMPLETION OF THE SCRIPT STARTS HERE ---
% ==========================================
% 

[freq, sort_idx] = sort(freq);
p = p(sort_idx);
p_linear = 10.^(p / 10);
E_amp = sqrt(p_linear);
E_amp = E_amp/ max(E_amp);

f_axis_Hz = (-N/2 : N/2 - 1) * (Fs / N);              
f_axis_THz = (f_axis_Hz / 1e12) + center_freq_THz; 
H_f = interp1(freq, E_amp, f_axis_THz, 'linear', 0);

noise_raw = (randn(1, N) + 1i * randn(1, N)) / sqrt(2);
Es_fft = fftshift(fft(noise_raw)) .* (H_f(:).');
E_s = ifft(ifftshift(Es_fft));
E_s = E_s / max(abs(E_s));

% crop_idx = (f_axis_THz >= 190) & (f_axis_THz <= 195);
% Y = abs(fftshift(fft(E_s)));
% plot(f_axis_THz(crop_idx), Y(crop_idx))

% ==========================================
% --- COMPARISON PLOT ---
% ==========================================
figure('Position', [100, 100, 900, 500]);

% Plot 1: Original Experimental Data
p_normalized = p - max(p);
plot(freq, p_normalized, 'b', 'LineWidth', 2, 'DisplayName', 'Original OSA Data');
hold on;

% Plot 2: Simulated PSD from Autocorrelation
RBW = 10e9;
[S_xx_dB, f_axis_Hz] = getPSD(E_s, Fs, RBW);
f_auto_THz = (f_axis_Hz / 1e12) + center_freq_THz;
S_xx_normalized = S_xx_dB - max(S_xx_dB);
plot(f_auto_THz, S_xx_normalized, 'r', 'LineWidth', 1.5, 'DisplayName', 'Simulated PSD (pwelch method)');

RBW = 10e9;
[S_xx_dB, f_axis_Hz] = getPSD_lag(E_s, Fs, RBW);
f_auto_THz = (f_axis_Hz / 1e12) + center_freq_THz;
S_xx_normalized = S_xx_dB - max(S_xx_dB);
plot(f_auto_THz, S_xx_normalized, 'y', 'LineWidth', 1.5, 'DisplayName', 'Simulated PSD (Autocorrelation)');


title('Wiener-Khinchin Theorem: Experimental vs. Simulated PSD');
xlabel('Absolute Frequency (THz)');
ylabel('Normalized Power (dB)');
legend('Location', 'best');
grid on;

% Zoom in strictly on the signal region to avoid the empty noise floor
xlim([min(freq) - 0.2, max(freq) + 0.2]);
ylim([-60, 5]);