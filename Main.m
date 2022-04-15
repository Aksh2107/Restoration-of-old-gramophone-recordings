clc;
clear all;
close all;
% Initialize model order, block size
model_order = 16;
block_size = 2400;
threshold  = 0.02;
% Audio Read of wav
[samples,fs] = audioread("Media1.wav");
% Taking the first 5 seconds of the audio
if mod(fs,block_size) == 0
    t = fs * 20;
else
% To make the no.of samples divisible by block size    
    temp = fs;
    temp = temp + block_size - mod(fs,block_size);
    t = temp * 20;
end    
raw_data = samples(1:t,1);        
restored1 = [];
restored = [];
detected_missing = zeros(1,block_size);
% Spliting the data into blocks
data = buffer(raw_data,block_size,block_size/2,'nodelay');
% Estimate AR coeffiecients, prediction error and interpolation for all
for j = 1:size(data,2)
    block = data(:,j);
    % Detect the clicks and remove from the inner block 
    inner_block   = block( block_size * 0.25 + 1 : block_size * 0.75);
    % from samples
    missing_samples = find(inner_block == 0);
    % from clicks
    missing_data = find(abs(inner_block) > threshold);
    % concatenate all missing values
    if isempty(missing_samples) == 0
        missing_data = vertcat(missing_samples,missing_data);
        missing_data = sort(missing_data);
    end    
    % Function call to estimate AR coefficients
    [coeffs, avg] = estimateARcoeffs(block,model_order);
    % Function call to get residual
    residual_data = getResidual(block - avg, coeffs);
    if isempty(missing_data) == 1
        restored2 = block( block_size * 0.25 + 1 : block_size * 0.75);
    else
        %detected_missing(1:block_size * 0.25) = zeros(1,block_size * 0.25);
        detected_missing((block_size * 0.25) + missing_data) = 1;
        block((block_size * 0.25) + missing_data) = 0;
        %detected_missing = zeros(1,block_size * 0.25);
        [restored2, Ak2, Au2, ik2] = interpolateAR(block, detected_missing, coeffs');
        restored2 = restored2( block_size * 0.25 + 1 : block_size * 0.75);
    end    
    restored1 = [restored1 restored2'];
    missing_data = [];
    detected_missing = zeros(1,block_size);
end
restored = [zeros(1,block_size * 0.25),restored1];
rem = audioplayer(restored,fs);
play(rem);
% *********************************************************************** %
% Funtion to estimate the AR coefficients 
function [coeffs, avg] = estimateARcoeffs(data, model_order)
    r = zeros(model_order,1); 
    R = zeros(model_order,model_order); 
    B = size(data,1);
    for j = 1:model_order   
        for k = 1:B-1
            if k - j <= 0             
                continue          
            end 
            r(j,1) = r(j,1) + ( data(k) * data(k - j));
        end  
    end
    clear j;
    for p = 1:model_order
        for j = 1:model_order
            for k = 1:B-1 
                if k - p <= 0
                    continue
                end    
                if k - j <= 0
                    continue
                end    
                R(p,j) = R(p,j) + (data(k - p) * data(k - j));    
            end    
        end
    end
    coeffs = - (inv(R) * r);
    avg = abs(mean(coeffs));
end
% *********************************************************************** %
% Funtion to get the prediction error %
function [residual] = getResidual(data, coeffs)
% Generate the residual over the entire block
    coeffs = [1 coeffs'];
    [residual] = filter(coeffs,1,data);
end
% *********************************************************************** %
% Funtion to estimate the interpolated values %
function [restored2, Ak2, Au2, ik2] = interpolateAR(block, detected_missing, coeffs)
    K = flip(coeffs);
    L = not(detected_missing);
    N = length(block);
    M = N - length(coeffs);
    R = [K, 1, zeros(1,N-(length(K)+1))];
    C = [K(1), zeros(1,M-1)];
    % Building main matrix A
    A = toeplitz (C,R);
    % Finding Au
    Au2 = A .* detected_missing;
    Au2(:, (all(Au2 == 0,1))) = []; 
    % Finding Ak
    Ak2 = A .* L;
    Ak2(:, (all(Ak2 == 0,1))) = [];
    % Finding yk
    yk = nonzeros(block(length(block) * 0.25 + 1 : length(block) * 0.75));
    yk = vertcat(block(1:length(block)*0.25),yk,block(length(block)*0.75 + 1:length(block)));
    ik2 = yk - mean(yk);
    % Finding yu
    yu = (-1 * inv(Au2.' * Au2) * (Au2.' * Ak2 * ik2));
    % Denormalising yu to get restored2
    yu = yu +mean(yk);
    % Provides the index value of the missing
    F = find(detected_missing); 
    % Values in the detected_missing matrix
    for i = (1: length(yu))
        block(F(i)) = yu(i);
    end
    restored2 = block;
end
