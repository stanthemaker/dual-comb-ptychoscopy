% --- Parameters ---
Fs = 200e9;
dt = 1 / Fs;
Nt = 800; % number of time step in a batch
NT = 250; % number of time batch 
T_batch = Nt * dt; % sampling time for a batch
T = NT * T_batch;
t = 0:dt:T-dt;
df_batch = 1 / T_batch; %spectral resolution
df = 1 / T;

% --- Gaussian Signal Parameters ---
f_sig_rel = 27e9;    % Center frequency
BW_sig    = 10e9;    % "Span" (we will treat this as 4*sigma)
f_spacing = df;

% Define frequency vector
sigma = BW_sig / 4;  % Standard deviation for the Gaussian
f_vec = (f_sig_rel - 2*BW_sig) : f_spacing : (f_sig_rel + 2*BW_sig);

% --- Generate the Signal ---
E_s = zeros(size(t));

% Seed the random generator for reproducibility if desired
rng(42); 

for f_tone = f_vec
    weight = exp(-(f_tone - f_sig_rel)^2 / (2 * sigma^2));
    phi = 2 * pi * rand(); 
    E_s = E_s + weight * exp(1i * (2 * pi * f_tone * t + phi));
end

% Normalize
E_s = E_s / max(abs(E_s));
plot(t, E_s);

%%
save('BroadbandSignal.mat', 'E_s', 't','f_sig_rel');