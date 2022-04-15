
load referenceARSignal.mat
figure(1); plot(data); title('AR3 process');
model_order = 3;
[coeffs, avg] = estimateARcoeffs(data, model_order);

%Estimate AR coefficient 
function [coeffs, avg] = estimateARcoeffs(data, model_order)
%length of data and its mean value
N=length(data);
mean_data=mean(data);

new_data= data -mean_data;
K= model_order+1;

%for row we take i and column we take i_c

for i = 1:1:model_order
    j(i,1) = sum(new_data(K : N) .* new_data(K - i : N - i));
    for i_c = 1 : 1 : model_order
    k(i,i_c) = sum(new_data(K - i :N - i) .* new_data (K - i_c : N - i_c));
    end
end

%inverse of k
k_inv= eye(size(k))/k;

coeffs= -(k_inv*j);
avg = mean(data(K : N));

end

%ALTERNATVIE METHOD
%[coeffs,avg] = armcov(data,model_order);
%coeffs = coeffs(2:end);
%coeffs = coeffs';
%avg = mean(coeffs);
%avg = abs(avg);
