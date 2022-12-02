%% Generate Impulse Response Function

function ImpulseRespm = Gen_ImRes(Omega, F, mlag, n0, ImpulseRespm, iter)

% # of variables
k = rows(Omega);

% Cholesky Decompositions (Lower triangular matrix)
Binv = chol(Omega)';

% Initial values
FF = eye(rows(F));


for j = 1:(mlag+1)
    psi_j = FF(1:k, 1:k); % First black of F^j Matrix, k by k
    theta = psi_j*Binv;   % Impulse Response, k by k
    theta = vec(theta);   % k^2 by 1

    % shock 1 => 1, 2, ..., k
    % shock 2 => 1, 2, ..., k
    % ...
    % shock k => 1, 2, ..., k
    for i = 1:k^2
        ImpulseRespm(iter-n0, j, i) = theta(i);
    end

    % Update, use F^j in nest lag
    FF = FF*F;

end

end