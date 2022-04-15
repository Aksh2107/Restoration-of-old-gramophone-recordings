% You know enough now that I don't have to give you the starter line
function [restored2, Ak2, Au2, ik2] = interpolateAR(block, detected_missing, coeffs)
B = length(block);
coeff = flip(coeffs);
L = length(coeffs);
N = not(detected_missing);
B = length(block);
K = B - L;



R = [coeff, 1, zeros(1, B - (L+1))];
C = [coeff(1), zeros(1,K-1)];
A = toeplitz(C,R)


%size(block)
Au2 = A .* detected_missing;
Au2(:, all(Au2 == 0,1)) = [];

Ak2 = A .* N;
Ak2(:, all(Ak2 == 0,1)) = [];

yk = nonzeros(block);
Avg = mean(yk);
yk2 = yk - Avg;
ik2 = yk2

yu = -inv(Au2'*Au2)*Au2'*Ak2*ik2;

yu = yu + Avg

F = find(detected_missing);

for t = 1 : length(yu)
    k = F(t)
    block(k) = yu(t)
end

restored2 = block ;

end