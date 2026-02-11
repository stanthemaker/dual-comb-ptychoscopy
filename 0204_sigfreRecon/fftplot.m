function fftplot(t, I)
    dt = t(2) - t(1);           % Time step
    Fs = 1/dt;                  % Sampling frequency
    L = length(I);              % Length of signal
    f = (0:(L/2)) * (Fs/L);     % Frequency vector
    
    %% 2. Perform FFT
    % Compute the two-sided spectrum P2. 
    % Then compute the single-sided spectrum P1 based on P2 and the signal length L.
    Y = fft(I);
    Y = Y(1:L/2+1);
    
    P1 = abs(Y/L);
    P1(2:end-1) = 2*P1(2:end-1);
    % eps = 1e-9;
    P1 = 20*log10(P1); % Magnitude in dB
    
    %% 2. Plot with Linear Frequency Axis (Standard for RF/Comb)
    figure;
    % plot(f, Y, 'LineWidth', 1.5) 
    plot(f, P1, 'LineWidth', 1.5)
    
    xlabel('Frequency (THz)')
    ylabel('Magnitude (dB)')
   
end