%% Make Y(2-dim) and X(3-dim)

function [Y0, YLm] = makeYX(Y,p)

% Data information: Sample Szie(T), Number of Variables(k)
T = rows(Y);
k = cols(Y);

% Y(t): Dependent Variables
Y0 = Y(p+1:T,:);

% 2-dimension, (T-p) by p*k
% Y(t-1), Y(t-2), ..., Y(t-p) : Lagged Dependent Variables
YL = zeros(T-p, p*k);
for i = 1:p
    YL(:, k*(i-1)+1:k*i) = Y(p+1-i:T-i, :);
end

% 3-dimension, k by p*k by (T-p)
YLm = zeros(k, k*k*p, T-p);
for t = 1: (T-p)
    % I(k) (X) [Y(t-1), Y(t-2), ..., Y(t-p)]
    xt = kron(eye(k), YL(t,:));
    YLm(:,:,t) = xt;

end

end