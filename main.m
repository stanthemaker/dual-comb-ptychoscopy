close all;
clear;
clc;

c = 299792458; 
% --- Comb Parameters ---
frep_1 = 10e3;
frep_2 = 10.08e3;
d_frep     = abs(frep_2 - frep_1);
d_tau = abs(1/frep_2 - 1/frep_1);

n_indices = 8:42;

% --- sampling setup ---
Fs = 1400e3;
dt = 1 / Fs;
N_avg = 100 + 2;             %% number of averages ; +2 for buffer
T_igm = 1 / d_frep;        %% time duration for an entire inteferogram
T = N_avg * T_igm;
Nt = ceil(T / dt);
df = 1 / T;
t = (0:Nt-1) * dt;
N_tau = round(1/frep_1 / d_tau) -1;

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
N_fft = 2^nextpow2(N_igm);
[Pss, f_psd] = pwelch(E_s, rectwin(N_igm), N_igm/2, N_fft, Fs, 'centered');
Pss_norm = Pss / max(abs(Pss));

[P12, ] = cpsd(E_c2, E_c1, rectwin(N_igm), N_igm/2, N_fft, Fs, 'centered');
P12_norm = P12 / max(abs(P12));

Gt = P12_norm .* Pss_norm;
Gt = Gt / max(abs(Gt));

figure

plot(f_psd/1e3, 10*log10(Pss_norm)); hold on
plot(f_psd/1e3, 10*log10(P12_norm));
% plot(f_psd/1e3, 10*log10(Gt));
ylim([-150 , 0])
legend("Pss norm", "P12 norm")
xlabel("Frequency (kHz)")
set(gca, "Fontsize", 18,"Linewidth",1.5)

%% Measurement
E_c1s = E_c1 .* conj(E_s);
E_c2s = E_c2 .* conj(E_s);

M1M2_avg = zeros(N_tau, 1);

% Define a fixed margin to capture the pulse widths (e.g., 15% of a period).
pulse_margin = round(0.1 * Fs / frep_1); 

for n_avg = 1:N_avg-1
    T_sweep = T_igm * (n_avg);
    
    for n_tau = 1:N_tau
        % 1. Position of the stationary E_c1 pulse
        T_pulse1 = (n_tau-1) * (1/frep_1) + T_sweep;
        
        % 2. Find the NEAREST E_c2 pulse. 
        m_nearest = round((T_pulse1 - T_sweep) * frep_2);
        T_pulse2 = T_sweep + m_nearest * (1/frep_2);
        
        % 3. Convert both exact times to indices
        idx_1 = round(T_pulse1 * Fs) + 1;
        idx_2 = round(T_pulse2 * Fs) + 1;
        
        % 4. INDEPENDENT WINDOWS: Tight bounding box around EACH pulse
        idx1_start = max(1, idx_1 - pulse_margin);
        idx1_end   = min(length(E_c1s), idx_1 + pulse_margin);
        
        idx2_start = max(1, idx_2 - pulse_margin);
        idx2_end   = min(length(E_c2s), idx_2 + pulse_margin);
        
        % 5. Sum and multiply independently
        M1 = sum(E_c1s(idx1_start:idx1_end));
        M2 = sum(E_c2s(idx2_start:idx2_end));
        M1M2_avg(n_tau) = M1M2_avg(n_tau) + (conj(M1) * M2);
        
        % --- plotting to verify ---
        % Create a static viewing window around E_c1 to watch E_c2 walk through it
        % view_margin = round(20e-6 * Fs); % 20 us view margin
        % idx_view_start = max(1, idx_1 - view_margin);
        % idx_view_end   = min(length(t), idx_1 + view_margin);
        % 
        % t_window = t(idx_view_start:idx_view_end) * 1e6; 
        % E1_window = abs(E_c1(idx_view_start:idx_view_end)); 
        % E2_window = abs(E_c2(idx_view_start:idx_view_end));
        % Es_window = abs(E_s(idx_view_start:idx_view_end));
        % 
        % clf; 
        % plot(t_window, E1_window, 'b', 'LineWidth', 1.5); hold on;
        % plot(t_window, E2_window, 'r', 'LineWidth', 1.5);
        % plot(t_window, Es_window, 'Color',[.7 .7 .7], 'LineWidth', 1.5);
        % 
        % % Draw lines for BOTH calculated centers
        % xline(T_pulse1 * 1e6, 'b--', 'E_{c1} Center', 'LineWidth', 1.5);
        % xline(T_pulse2 * 1e6, 'r--', 'E_{c2} Nearest', 'LineWidth', 1.5);
        % 
        % % Highlight the independent integration windows with dotted lines
        % xline(t(idx1_start)*1e6, 'b--'); xline(t(idx1_end)*1e6, 'b--');
        % xline(t(idx2_start)*1e6, 'r--'); xline(t(idx2_end)*1e6, 'r--');
        % 
        % xlabel('Laboratory Time (\mu s)');
        % ylabel('Electric Field Envelope |E|');
        % legend('E_{c1}', 'E_{c2}', 'E_s', 'Location', 'northeast');
        % grid on;
        % 
        % % Lock the X-axis width 
        % xlim([T_pulse1*1e6 - 20, T_pulse1*1e6 + 20]);
        % 
        % pause;
    end
end

%time domain -> delay time domain
%%
tau_wrapped = zeros(N_tau, 1);
for n_tau = 1:N_tau
    if n_tau <= ceil(N_tau/2)
        tau_wrapped(n_tau) = -(n_tau - 1) * d_tau; % 0 to -T/2
    else
        tau_wrapped(n_tau) = (N_tau - n_tau + 1) * d_tau; % +T/2 to 0
    end
end
[tau_sorted, sort_idx] = sort(tau_wrapped);

M1M2_avg = M1M2_avg / (N_avg-1);
% 1. sort_idx unscrambles the data
% 2. flip() reverses the time axis to match MATLAB's cpsd convention
M1M2_avg = flip(M1M2_avg(sort_idx));
N_pad = length(M1M2_avg) * 10;

M1M2_fft = fftshift(fft(M1M2_avg, N_pad));
M1M2_fft_norm = M1M2_fft / max(abs(M1M2_fft));

Fs_M1M2 = 1 / d_tau;
df_M1M2 = 1 / (N_pad * d_tau);
f_M1M2 = ((0 : N_pad-1) - floor(N_pad/2)) * df_M1M2;
% df_M1M2 = 1 / (N_tau * d_tau);
% f_M1M2 = ((0 : N_tau-1) - floor(N_tau/2)) * df_M1M2;
%%
figure
plot(f_psd/1e3, 10*log10(Pss_norm)); hold on 0
% plot(f_psd/1e3, 10*log10(P12_norm));
plot(f_M1M2/1e3, 10*log10(M1M2_fft_norm));
ylim([-70 ,0])
xlim([-700,700])
% --- 3. Stacked Subplot Comparison ---
% figure('Name', 'Dual-Comb Spectrum Validation', 'Position', [100, 100, 800, 600]);
% 
% % Top Plot: Simulated Extraction
% subplot(2, 1, 1);
% plot(f_M1M2, 10*log10(M1M2_fft_norm), 'b', 'LineWidth', 1.5);
% title('Simulated Extraction: FFT of \langle M_1 M_2^* \rangle');
% ylabel('Magnitude (dB)');
% grid on;
% % Limit the X-axis to the relevant RF band to match visual scales if necessary
% % xlim([-frep_1/2, frep_1/2]); 
% ylim([-200 5]);
% 
% % Bottom Plot: Theoretical Ground Truth
% subplot(2, 1, 2);
% plot(f_psd, 10*log10(Gt), 'r', 'LineWidth', 1.5);
% title('Theoretical Ground Truth: P_{ss} \times P_{12}');
% xlabel('RF Frequency (Hz)');
% ylabel('Magnitude (dB)');
% grid on;
% % Match the X-axis limits of the top plot for a 1-to-1 visual comparison
% xlim([min(f_M1M2) max(f_M1M2)]); 
% ylim([-200 5]);

