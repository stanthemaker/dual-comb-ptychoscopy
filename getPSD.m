function [pxx_dB, f_axis] = getPSD(signal, Fs, RBW)
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
    window_length = round((Fs / RBW) * 1.36);
    
    % Ensure the window isn't longer than the signal itself
    if window_length > length(signal)
        warning('RBW is too narrow for the given signal length. Using maximum possible length.');
        window_length = length(signal);
    end
    
    % 2. Standard 50% overlap for Welch's method
    noverlap = round(window_length / 2);
    
    % 3. Calculate PSD using pwelch
    % 'centered' automatically shifts the zero-frequency to the middle
    [pxx_linear, f_axis] = pwelch(signal, window_length, noverlap, window_length, Fs, 'centered');
    
    % 4. Convert to Decibels, adding eps to prevent log10(0)
    pxx_dB = 10 * log10(pxx_linear + eps);
end