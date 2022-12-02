function [bhat, Yhat, ehat, sig2hat, varbhat, stde, t_val, TSS, RSS, R2, R2_, SC, AIC] = OLS_OUT(Y,X,printi)

T = rows(X);
k = cols(X);

% OLS estimator
bhat = inv(X'*X)*X'*Y;
Yhat = X*bhat;
ehat = Y - Yhat;
sig2hat = ehat'*ehat/ (T-k);
varbhat = sig2hat*inv(X'*X);
stde = sqrt(diag(varbhat));

% Model comparison
RSS = ehat'*ehat;
TSS = (Y - mean(Y))'*(Y - mean(Y));
R2 = 1 - RSS/TSS;
R2_ = 1 - (RSS*(T-1))/(TSS*(T-k));
SC = log(RSS/T) + k*log(T)/T;
AIC = log(RSS/T) + 2*k/T;

% Simple Hypothesis Testing
b0 = zeros(k,1);
t_val = (bhat - b0)./ stde;
p_val = 2*(1 - cdf('t',abs(t_val),T-k));

% Confidence Interval 90% = 1.64, 95% = 1.96, 99% = 2.57
alpha = 0.05;
UB = bhat + tinv(1-alpha/2, T-k)*stde;
LB = bhat - tinv(1-alpha/2, T-k)*stde;

if printi > 0
    
    % Table
    disp('======================================================================================')
    disp(['    Estimates   ',    '     S.E.   ', '  t value  ',   ' p value  ', '        95% CI'])
    disp('======================================================================================')
    disp([bhat stde t_val p_val LB UB]);
    disp('======================================================================================')
    disp(['S.E. of regression is      ',   num2str(sqrt(sig2hat))]);
    disp(['R2 is   ', num2str(R2)]);
    disp(['adjusted R2 is   ', num2str(R2_)]);
    disp(['SC is    ', num2str(SC)]);
    disp(['AIC is   ', num2str(AIC)])'
    disp('======================================================================================')

    if printi == 2
        % Figure
        subplot(2,1,1);
        plot([Y Yhat]);
        legend('Y', 'Yhat');
        title('Actual and Fitted');
        subplot(2,1,2);
        plot([ehat zeros(T,1)]);
        legend('ehat');
        title('Residual');
    end
    
end

end

