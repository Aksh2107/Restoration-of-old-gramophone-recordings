function [residual] = getResidual(data, coeffs)
% Generate the residual over the entire block
% put your code here
coeffs=[1,coeffs];
residual = filter(coeffs, 1, data);
end

% Load data and plot it for sanity
load referenceARSignal.mat
figure(1); plot(data); title('AR3 process');
% Model the data
model_order = 3;
coeffs = [-2.4 2.3 -0.9];
avg = mean(data);
% Calculate the residual here (note the normalisation)
res = getResidual(data - avg, coeffs);
