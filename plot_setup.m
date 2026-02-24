%% Input: time domain signal
function plotsetup(comb1 , comb2, signal, Fs)
    L = length(comb1);
    f = (0:(L/2)) * (Fs/L);

    Y = fft(comb1);
    Y = Y(1:L/2+1);
    P = abs(Y/L);
    P(2:end-1) = 2*P(2:end-1);
    P = 20*log10(P); 
    figure;
    plot(f, P, 'LineWidth', 1.5);
    hold on;
    
    Y = fft(comb2);
    Y = Y(1:L/2+1);
    P = abs(Y/L);
    P(2:end-1) = 2*P(2:end-1);
    P = 20*log10(P);
    plot(f, P, 'LineWidth', 1.5);

    Y = fft(signal);
    Y = Y(1:L/2+1);
    P = abs(Y/L);
    P(2:end-1) = 2*P(2:end-1);
    P = 20*log10(P);
    plot(f, P, 'LineWidth', 1.5);

    set(gca, "Fontsize", 18,"Linewidth",1.5)
    xlabel("Frequency")
    ylabel("PWD (dB)")
    legend("comb1", "comb2", "signal");
end