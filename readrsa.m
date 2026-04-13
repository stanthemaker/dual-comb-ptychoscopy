
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["FrequencyHz", "DatadBm"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable("./0403_AOM/AOM_beatnotes_rsa_zoomin.csv", opts);
data = table2array(data);
clear opts

f = data(:,1);
p = data(:,2);

%%
plot(f/1e6,p)
xlabel("Frequency (MHz)")
ylabel("Power (dB)")
set(gca, "Fontsize", 18,"Linewidth",1.5)