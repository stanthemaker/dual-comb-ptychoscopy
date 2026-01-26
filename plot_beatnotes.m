function analyze_beat_notes(comb1, comb2, signal,dt)

    %% 1. Simulate Photodetection (Intensity)
    % The detector measures Power ~ |E|^2
    I1 = abs(comb1+signal).^2;
    I11 = abs(comb1).^2;
    I2 = abs(comb2+signal).^2;
    I22 = abs(comb2).^2;

    % Remove DC component (Mean) to visualize AC beat notes clearly
    I1 = I1 - mean(I1);
    I2 = I2 - mean(I2);
    I11 = I11 - mean(I11);
    I22 = I22 - mean(I22);

    %% 2. Compute FFT
    L = length(I1);
    Fs = 1/dt;
    
    % Frequency Axis
    frequencies = (-L/2 : L/2 - 1) * (Fs/L);
    
    % FFT (Normalized)
    spec_I1 = abs(fftshift(fft(I1)) - fftshift(fft(I11))) / L;
    spec_I2 = abs(fftshift(fft(I2))- fftshift(fft(I22))) / L;

    %% 3. Plotting
    figure('Name', 'Dual-Comb Intermediate Spectrum Analysis', 'Color', 'w', 'Position', [100, 100, 1000, 800]);

    plot(frequencies/1e9, spec_I1, 'b', 'LineWidth', 1.2); 
    hold on;
    plot(frequencies/1e9, spec_I2, 'r', 'LineWidth', 1.2);
    % visual formatting
    set(gca,"FontSize",18);
    title('E1');
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (a.u.)');
    grid on;
    xlim([-12 12])
    
end