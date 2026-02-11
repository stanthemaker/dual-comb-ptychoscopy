
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [30, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["CSV", "VarName2"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable("./combspec_thz.CSV", opts);
data = table2array(data);

clear opts

x = data(:,1);
y = data(:,2);
p = 10.^(y/10);
c = 299792458; 
%%
%all freq are in THz, time in ps
f_offset = 190; % for simulation purpose
f = c ./ lambda / 1e3 - f_offset; 
f_center = c / 1552.1 / 1e3 - f_offset;
f_sig = 193.212 - f_offset;
frep = 0.027;
fmax = 195 - f_offset ; % max freq in comb

Fs = 2 * fmax;
dt = 1/ Fs; % in ps
Nt = 2^16;
t = 0:dt:(Nt-1)*dt;
T = Nt * dt; % in ps
df = 1 / T; % in THz

%%
E_c1 = zeros(size(t));
range = round(length(f)/2)-400: round(length(f)/2);
for i = range
    E_c1 = E_c1 + sqrt(p(i)) * exp(1i * 2*pi * f(i) * t);
end
% fftplot(t, E_c1);
%%

% E_c1 = exp(1i * 2*pi * 3.206 * t);
E_s = exp(1i * 2*pi * f_sig * t);

I1 = abs(E_s + E_c1).^2 - abs(E_c1).^2 - abs(E_s).^2;

cutoff_freq = frep / 2;
% I1 = lowpass(I1, cutoff_freq, Fs);
I1_low = lowpass(I1, cutoff_freq, Fs, 'ImpulseResponse', 'iir', 'Steepness', 0.95);

%%
% fftplot(t, I1);
fftplot(t, I1_low);
% xlim([0, cutoff_freq]);