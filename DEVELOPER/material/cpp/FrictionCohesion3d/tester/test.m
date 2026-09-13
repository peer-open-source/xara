clear all
clc
close all








filename = '..\..\..\..\..\Win64\bin\node_displ.out';
delimiter = ' ';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

tau   = dataArray{:, 1};
gammaY    = dataArray{:, 2};
gammaZ = dataArray{:, 3};
eps = dataArray{:, 4};

clearvars filename delimiter formatSpec fileID dataArray ans;

if length(tau)>length(eps)
    tau = tau(1:length(eps))
end

plot(gammaZ, tau)