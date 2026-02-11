
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
data = readtable("./comb1.CSV", opts);
data = table2array(data);


c = 299792458; 
x = data(:,1);
f = c ./ x / 1e3; % in THz
comb1 = data(:,2);

data = readtable("./comb2.CSV", opts);
data = table2array(data);
comb2 = data(:,2);

data = readtable("./signal.CSV", opts);
data = table2array(data);
x = data(:,1);
f_sig = c ./ x / 1e3;
signal = data(:,2);

plot(f,comb1);
hold on;
plot(f,comb2);
plot(f_sig,signal);
set(gca, "Fontsize", 18, "Linewidth",1.5);
xline(c/1552.1/1e3, '--r')
legend("comb1", "comb2", "signal","center freq");
clear opts

