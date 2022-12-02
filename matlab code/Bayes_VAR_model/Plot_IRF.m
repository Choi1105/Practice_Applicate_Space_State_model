%% Plot Impulse Response Functions
% Input: Parameters(ImpulseRespm), alpha(credible %)
% Output: Impulse response figures

function Plot_IRF(ImpulseRespm, alpha)

% Impulse response information
[~, mlag, k2] = size(ImpulseRespm);
k = sqrt(k2);

% Credibility Interval Percent
ql = [alpha;0.5;(1-alpha)];

% X-axis index for impulse response function
xlag = 0:(mlag-1);

% Index
ind = 1:k2;
ind = reshape(ind, k, k);

figure
zeroline = zeros(mlag, 1); % Baseline
for i = 1:k2

    % Select IR(i,j) and Compute quantile
    ImpulseResp_ij = ImpulseRespm(:,:,i);            % n1 by (mlag + 1)
    ImpulseResp_ij = quantile(ImpulseResp_ij, ql)';  % (mlag + 1) by 3

    % Index
    [r,c] = find(ind==1);

    % 1st row: 1st variable shock to other variables
    % 2nd row: 2nd variable shock to other variables
    % ...
    % pth row: pth variable shock to other variables
    subplot(k,k,1);
    plot(xlag, ImpulseResp_ij(:,1), 'k--', xlag, ImpulseResp_ij(:,2), ...
        'b-', xlag, ImpulseResp_ij(:,3), 'k--', xlag, zeroline, 'k:', 'linewidth', 2);
    xlim([0 mlag]);
    xlabel('Lag')
    title(['shock  ', num2str(c), ' to vari ', num2str(r)]);


end

end