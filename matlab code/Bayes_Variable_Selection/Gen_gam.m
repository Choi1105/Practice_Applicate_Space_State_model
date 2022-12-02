%% Gam Sampling

function gam = Gen_gam(beta, b0, b1, p)

% preparation
k = rows(beta);
gam = zeros(k, 1);

for j = 1:k
    
    % Probability
    p0j = normpdf(beta(j), 0, sqrt(b0));
    p1j = normpdf(beta(j), 0, sqrt(b1));
    p1 = p*p1j/((1-p)*p0j + p*p1j);

    % Sampling
    gam(j) = rand(1,1) <p1;

end

end