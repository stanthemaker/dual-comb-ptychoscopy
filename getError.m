function getError(spec_recon_linear, E_s, Fs, global_freq_axis, f_center)
    % spec_recon_linear: The magnitude of spec_recon (NOT in dB)
    % E_s: The original time-domain signal
    % Fs: Sampling frequency
    % global_freq_axis: The optical frequency axis from your stitching
    % f_center: The reference optical center frequency
    
    % --- 1. Calculate Original Spectrum ---
    N = length(E_s);
    % Get the FFT of the original signal
    F_orig = fftshift(fft(E_s));
    % Create the frequency axis for the original signal (relative to f_center)
    freq_orig = (-N/2 : N/2-1) * (Fs/N) + f_center;
    
    % --- 2. Align Grids ---
    % Interpolate the original magnitude onto the global_freq_axis
    % We use the magnitude (abs) because recon is often power/magnitude based
    orig_mag_interp = interp1(freq_orig, abs(F_orig), global_freq_axis, 'linear', 0);
    
    % --- 3. Normalization ---
    % Normalize both to a peak of 1 for a fair shape comparison
    recon_norm = spec_recon_linear / max(spec_recon_linear(:));
    orig_norm  = orig_mag_interp / max(orig_mag_interp(:));
    
    % --- 4. Error Metrics ---
    % Mean Squared Error (MSE)
    mse = mean((recon_norm - orig_norm).^2);
    
    % Root Mean Squared Error (RMSE)
    rmse = sqrt(mse);
    
    % Spectral Overlap (Fidelity)
    % 1 means perfect match, 0 means no overlap
    fidelity = sum(recon_norm .* orig_norm) / (norm(recon_norm) * norm(orig_norm));

    % --- 5. Plot Comparison ---
    figure;
    set(gca, "Fontsize", 18);
    subplot(2,1,1);
    plot(global_freq_axis/1e12, orig_norm, 'k', 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(global_freq_axis/1e12, recon_norm, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed');
    title('Normalized Spectral Comparison');
    ylabel('Normalized Amplitude');
    legend; grid on;
    
    subplot(2,1,2);
    plot(global_freq_axis/1e12, abs(recon_norm - orig_norm));
    title(sprintf('Absolute Error (RMSE: %.4f, Fidelity: %.4f)', rmse, fidelity));
    xlabel('Optical Frequency (THz)');
    ylabel('Error');
    grid on;
    
    fprintf('--- Reconstruction Accuracy ---\n');
    fprintf('RMSE: %.4f\n', rmse);
    fprintf('Fidelity (Overlap): %.4f\n', fidelity);
    
end