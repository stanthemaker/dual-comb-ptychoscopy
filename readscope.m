opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [22, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable("./0403_AOM/AOM_beatnote_longest.csv", opts);
data = table2array(data);

clear opts
c = 299792458;
t = data(:,1);
v = data(:,2);

plot(t * 1e9, v*1e3);
xlabel("Time (ns)")
ylabel("voltage (mV)")
set(gca, "Fontsize", 18,"Linewidth",1.5)
%%

dt = t(2) - t(1);
pulse_freq = 100e6;
pulse_T = 1 / pulse_freq;
pulse_N = round(pulse_T / dt); % Added round() to prevent indexing errors
v_N = v(1:pulse_N);

% 1. Calculate comb power
V_square_sum = sum(v_N.^2);
pulse_energy = V_square_sum * dt / 50; 
pulse_power_watts = pulse_energy / pulse_T; 
pulse_power_dBm = 10 * log10(pulse_power_watts * 1000);



%%
dt = t(2) - t(1);
Fs = 1/dt;
L = length(v);
f = (0:(L/2)) * (Fs/L);

Y = fft(v);
Y = Y(1:L/2+1);
P = abs(Y/L);
P(2:end-1) = 2*P(2:end-1);
P = 20*log10(P);
plot(f/1e6,P)

xlabel("frequency (MHz)")
ylabel("Power (dB)")
set(gca, "Fontsize", 18,"Linewidth",1.5)