close all;
clear;
clc;

c = 299792458;  
f_center      = 193.151e12;
f_sig         = 193.181e12;
f_sig_rel     =  f_center - f_sig;
% --- Comb Parameters ---
frep_1 = 27.320e9;
frep_2 = 27.220e9;
d_frep     = frep_2 - frep_1;   % 100 MHz difference
comb2_shift = 80e6;

% --- read inteferogram ---
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "Var2"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

data = readtable("./interfer1.Wfm.csv", opts);
data = table2array(data);
t = data(:,1);
I1 = data(:,2);

data = readtable("./interfer2.Wfm.csv", opts);
data = table2array(data);
I2 = data(:,2);

dt = t(2) - t(1);           % Time step
Fs = 1/dt;                  % Sampling frequency
L = length(I1);              % Length of signal
T = L * dt;




Nt = 508;
T_batch = Nt * dt;
df = 1 / T_batch;
NT = 3;
f = (-Nt/2:Nt/2-1) * (Fs/Nt);
%%
cutoff_freq = frep_1 /2;
% I1 = lowpass(I1, cutoff_freq, Fs);
% I2 = lowpass(I2, cutoff_freq, Fs);

I1 = reshape(I1, Nt, NT);
I2 = reshape(I2, Nt, NT);

F1 = fftshift(fft(I1,[],1),1);
F2 = fftshift(fft(I2,[],1),1);

bins = round(39.8e6 / df);
center_idx = Nt/2;
F2_neg = F2(1:center_idx, :);
F2_pos = F2(center_idx+1:end, :);
F2_neg_shifted = circshift(F2_neg, -bins, 1);
F2_pos_shifted = circshift(F2_pos, bins, 1);
F2 = [F2_neg_shifted; F2_pos_shifted];
%%
n_indices = -3:3;
max_freq_span = (length(n_indices)+2) * frep_1;
total_bins    = ceil(max_freq_span / df);

% Initialize global spectrum
spec_recon = zeros(total_bins, 1);
center_idx = floor(total_bins / 2);
global_freq_axis = ((1:total_bins).' - center_idx) * df + f_center;

for n = n_indices
    Cn = zeros(Nt, 1);

    delta_n = (n * d_frep) + comb2_shift;
    shift_bins = round(delta_n / df);

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
    % plot(f/1e9 , Cn);
    
    freq_offset = n * frep_1;
    bin_offset_global = round(freq_offset / df);
    target_center_idx = center_idx + bin_offset_global;

    Cn_length = length(Cn);
    global_start_idx = target_center_idx - floor(Cn_length/2);
    global_end_idx   = global_start_idx + Cn_length - 1;
    
    if global_start_idx > 0 && global_end_idx <= total_bins
        spec_recon(global_start_idx:global_end_idx) = spec_recon(global_start_idx:global_end_idx) + Cn(1: end);
    else
        warning('Comb line n=%d is outside the spec_recon array bounds.', n);
    end
end
spec_recon = 20*log10(abs(spec_recon));
figure;
plot(global_freq_axis / 1e12, spec_recon, "LineWidth",1.5);

xlabel('Optical Frequency (THz)');
ylabel('Signal Power (dB)')
ylim([max(spec_recon)-60 max(spec_recon)+10]);
xline(f_center / 1e12,'--b', "Center freq");
xline(f_sig / 1e12,'--r');
legend("reconstructed signal","original signal");
set(gca,"Fontsize",18)