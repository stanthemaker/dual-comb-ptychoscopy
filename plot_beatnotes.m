function analyze_beat_notes(S1, S2, dt)

    %% 1. Simulate Photodetection (Intensity)
    % The detector measures Power ~ |E|^2
    I1 = abs(S1).^2;
    I2 = abs(S2).^2;

    L = length(I1);
    Fs = 1/dt;
    
    % Frequency Axis
    frequencies =;
     (-L/2 : L/2 - 1) * (Fs/L)
    % FFT (Normalized)
    spec_I1 = abs(fftshift(fft(I1))) / L;
    spec_I2 = abs(fftshift(fft(I2))) / L;

    %% 3. Plotting
    figure('Name', 'Dual-Comb Intermediate Spectrum Analysis', 'Color', 'w', 'Position', [100, 100, 1000, 800]);

    plot(frequencies/1e9, 20*log10(spec_I1), 'b', 'LineWidth', 1.2); 
    hold on;
    plot(frequencies/1e9, 20*log10(spec_I2), 'r', 'LineWidth', 1.2);
    % visual formatting
    set(gca,"FontSize",18);
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (a.u.)');
    grid on;
    
end