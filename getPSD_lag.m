function [S_xx_dB, f_axis] = getPSD_lag(signal, Fs, RBW)
    % getPSD: Calculates the Power Spectral Density of a complex signal
    % mimicking the Resolution Bandwidth (RBW) of a hardware spectrum analyzer.
    %
    % Inputs:
    %   signal - Time domain electric field or voltage
    %   Fs     - Sampling frequency in Hz
    %   RBW    - Desired Resolution Bandwidth in Hz (e.g., 10e9 for 10 GHz)
    %
    % Outputs:
    %   pxx_dB - Power Spectral Density in dB
    %   f_axis - frequency axis in Hz

    % 1. Determine the window size to achieve the desired RBW.
    % A standard Hamming window has a noise bandwidth factor of ~1.36.
    % N_window = (Fs / RBW) * 1.36
    [Rxx, lags] = xcorr(signal, 'biased');
    M = length(Rxx); % This evaluates to 2N - 1

    window_width = round(Fs / RBW);
    lag_window = exp(-0.5 * (lags / (window_width/2)).^2); % Gaussian lag window

    Rxx_windowed = Rxx .* lag_window;
    S_xx_linear = abs(fftshift(fft(Rxx_windowed)));
    S_xx_dB = 10 * log10(S_xx_linear + eps);
    f_axis = (-(M-1)/2 : (M-1)/2) * (Fs / M);
end