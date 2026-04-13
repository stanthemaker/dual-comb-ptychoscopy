%% Input: time domain signal
function plotsetup(comb1 , comb2, signal, Fs, f_center)
    L = length(comb1);
    f = (0:(L/2)) * (Fs/L);

    Y = fft(comb1);
    Y = Y(1:L/2+1);
    P = abs(Y/L);
    P(2:end-1) = 2*P(2:end-1);
    P = 20*log10(P); 
    figure;
    plot((f+f_center)/1e12, P, 'LineWidth', 1.5);
    hold on;
    
    Y = fft(comb2);
    Y = Y(1:L/2+1);
    P = abs(Y/L);
    P(2:end-1) = 2*P(2:end-1);
    P = 20*log10(P);
    plot((f+f_center)/1e12, P, 'LineWidth', 1.5);

    RBW = 10e9;
    [S_xx_dB, f] = getPSD_lag(signal, Fs, RBW);
    plot((f+f_center)/1e12, S_xx_dB ,'LineWidth', 1.5);
    

    set(gca, "Fontsize", 18,"Linewidth",1.5)
    xlabel("Frequency (THz)")
    ylabel("PWD (dB)")
    legend("comb1", "comb2", "signal");
end