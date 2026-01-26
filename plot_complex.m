function plot_complex(S1, Fs)
    figure('Color','w'); 
    S1_flat = S1(:);
    L = length(S1_flat);
    f_axis_shifted = (-L/2 : L/2-1) * (Fs/L); % Frequency axis from -Fs/2 to +Fs/2
    Y_shifted = fftshift(fft(S1_flat));       % Shift zero freq to center
    
    plot(f_axis_shifted/1e9, 20*log10(abs(Y_shifted/L)));
    xlabel('RF Frequency (GHz)'); ylabel('Magnitude (dB)');
    title('Two-Sided Spectrum of S1');
    xlim([-40 40]); % Zoom in on the central region
    grid on;
    xline(0, 'k');
end