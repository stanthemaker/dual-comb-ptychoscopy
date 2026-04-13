
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [30, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable("./0403_AOM/15528-27226-AOM-EDFA.CSV", opts);
data = table2array(data);

clear opts
c = 299792458;
wl = data(:,1);
freq = c ./ wl *1e-3;
p = data(:,2);

