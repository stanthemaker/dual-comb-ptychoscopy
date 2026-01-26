
%% Visualization: True Optical Fields (Prior to Interference)
% --- Visualization Parameters ---
dt_vis = 1e-12;        % 1 ps time step (1 THz Bandwidth, No Aliasing)
Fs_vis = 1/dt_vis;     
T_vis  = 20/frep_1;    % Simulate just 20 pulse periods for the plot
t_vis  = 0:dt_vis:T_vis-dt_vis;
L_vis  = length(t_vis);

fprintf('--- Generating High-Res Visualization ---\n');
fprintf('Vis Sampling Rate: %.1f THz (Captures all comb lines)\n', Fs_vis/1e12);

% --- 1. Generate Comb 1 (High Res) ---
E_c1_vis = zeros(size(t_vis));
for n = n_indices
    f_n = n * frep_1;
    E_c1_vis = E_c1_vis + exp(1i * 2*pi * f_n * t_vis);
end
E_c1_vis = E_c1_vis / max(abs(E_c1_vis));

% --- 2. Generate Signal (High Res) ---
% Signal is at f_sig_rel (~ -125 GHz). 
% We generate it directly at this frequency.
noise_vis = (randn(size(t_vis)) + 1i*randn(size(t_vis)));

% Apply Spectral Width (Filter)
% We work in baseband 0Hz then shift, to ensure shape is correct
[b_vis, a_vis] = butter(1, (width_sig_Hz)/(Fs_vis/2));
env_vis = filter(b_vis, a_vis, noise_vis);

% Shift to -125 GHz
E_sig_vis = env_vis .* exp(1i * 2*pi * f_sig_rel * t_vis);
E_sig_vis = 0.05 * E_sig_vis / mean(abs(E_sig_vis));

% --- 3. Compute Spectra for Plotting ---
f_axis_vis = (-L_vis/2 : L_vis/2 - 1) * (Fs_vis/L_vis);
Spec_C1 = fftshift(fft(E_c1_vis)) / L_vis;
Spec_Sig = fftshift(fft(E_sig_vis)) / L_vis;

% --- 4. Plot ---
figure('Name', 'True Optical Input Fields', 'Color', 'w', 'Position', [100 100 1000 700]);

% Subplot 1: Frequency Domain (The "Map")
subplot(2,1,1);
% Plot Comb
stem(f_axis_vis/1e9, abs(Spec_C1).^2, 'b', 'LineWidth', 1.5, 'Marker', 'none', 'DisplayName', 'Comb Lines'); 
hold on;
% Plot Signal
plot(f_axis_vis/1e9, abs(Spec_Sig).^2, 'r', 'LineWidth', 1, 'DisplayName', 'Incoherent Signal');

% Formatting
xline(f_sig_rel/1e9, '--k', 'Target Freq');
title('Optical Spectrum (Relative to 1550 nm)');
xlabel('Frequency Offset (GHz)');
ylabel('Power (a.u.)');
legend;
grid on;
xlim([-300, 300]); % Show full range
ylim([0, max(abs(Spec_C1).^2)*1.1]);

% Subplot 2: Time Domain (The "Pulses")
subplot(2,1,2);
% Plot Comb Intensity (Short snippet)
plot(t_vis*1e12, abs(E_c1_vis).^2, 'b', 'DisplayName', 'Comb Intensity'); hold on;
% Plot Signal Intensity (Magnified)
plot(t_vis*1e12, 10 * abs(E_sig_vis).^2 + 0.5, 'r', 'DisplayName', 'Signal (Offset +0.5 and Amp x10)');

title('Time Domain Snapshot');
xlabel('Time (ps)');
ylabel('Intensity');
legend;
grid on;
xlim([0, 200]); % Zoom in on first 200ps (a few pulses)

fprintf('Visualization Complete. Proceeding to Ptychoscopy Simulation...\n');