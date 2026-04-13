function getError(f_recon, spec_recon, E_s, Fs, f_center)
    % f_recon: The frequency axis for the reconstructed spectrum
    % spec_recon: The reconstructed power spectrum
    % E_s: The target time-domain incoherent signal
    % Fs: Sampling frequency
    % N_igm: (Optional) Window size for pwelch. Defaults to length(E_s)/8.
    % N_fft: (Optional) Number of FFT points for pwelch. Defaults to next power of 2.
    
    N = length(E_s);
    N_window = round (N / 1000);
    N_pad = 2^nextpow2(N_window);

    % --- 1. Calculate Original Spectrum (PSD) ---
    % Use pwelch to get the Power Spectral Density of the incoherent signal
    [Pss, f_psd] = pwelch(E_s, rectwin(N_window), floor(N_window/2), N_pad, Fs, 'centered');
    f_psd = f_psd + f_center;

    % --- 2. Align Grids ---
    % Interpolate the ground truth PSD onto the recon frequency axis.
    % This correctly "downsamples" and aligns the high-res pwelch output to your f_recon grid.
    % We use '0' as the extrapolation value for any out-of-bounds frequencies.
    orig_power_interp = interp1(f_psd, Pss, f_recon, 'linear', 0);
    
    % --- 3. Normalization ---
    % Normalize both to a peak of 1 for a fair shape comparison
    recon_norm = spec_recon / max(spec_recon(:));
    orig_norm  = orig_power_interp / max(orig_power_interp(:));
    
    % --- 4. Error Metrics ---
    % % Mean Squared Error (MSE)
    % mse = mean((recon_norm - orig_norm).^2);
    % 
    % % Root Mean Squared Error (RMSE)
    % rmse = sqrt(mse);
    % 
    % % Spectral Overlap (Fidelity)
    % % 1 means perfect match, 0 means no overlap
    % fidelity = sum(recon_norm .* orig_norm) / (norm(recon_norm) * norm(orig_norm));
    % 
    % --- 5. Plot Comparison ---
    figure;
    
    % Top plot: Overlaid Spectra
    % subplot(2,1,1);
    % Note: Assumes f_recon is in Hz. Divides by 1e12 for THz plotting.
    plot(f_recon/1e12, orig_norm, 'k', 'LineWidth', 2, 'DisplayName', 'Original (pwelch)');
    hold on;
    plot(f_recon/1e12, recon_norm, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed');
    title('Normalized Spectral Power Comparison');
    ylabel('Normalized Power');
    % xline(f_s_global / 1e12, '--r', "Signal freq", 'LabelVerticalAlignment', 'bottom');
    set(gca, 'FontSize', 16);
    legend; grid on;
    
    % Bottom plot: Absolute Error
    % subplot(2,1,2);
    % plot(f_recon/1e12, abs(recon_norm - orig_norm), 'b', 'LineWidth', 1.5);
    % title(sprintf('Absolute Error (RMSE: %.4f, Fidelity: %.4f)', rmse, fidelity));
    % xlabel('Frequency (THz)');
    % ylabel('Absolute Error');
    % set(gca, 'FontSize', 14);
    % grid on;
    % 
    % Print metrics to command window
    % fprintf('--- Reconstruction Accuracy ---\n');
    % fprintf('RMSE: %.4f\n', rmse);
    % fprintf('Fidelity (Overlap): %.4f\n', fidelity);
    
end