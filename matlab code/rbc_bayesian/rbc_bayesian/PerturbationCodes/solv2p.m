% Copyright Andrew Binning 2015
% From:
% Binning, A. (2013). Solving second and third-order approximations to dsge models: a recur-
% sive sylvester equation solution. Norges Bank Working Paper 2013/18. 
% URL http://www.norges-bank.no/no/om/publisert/publikasjoner/working-papers/2013/18/.
% Based on: 
% Kamenik, O. (2005). Solving sdge models: A new algorithm for the sylvester equation.
% Working Papers 2005/10, Czech National Bank, Research Department.
% Martin, C. D. M. & Van Loan, C. F. (2006). Shifted kronecker product systems. SIAM J. Matrix
% Analysis Applications, 29 (1), 184{198.

function y = solv2p(alpha,beta_sq,F,K,d,k)

beta = sqrt(beta_sq);
n = length(K);
m = length(F);
y = zeros(n*m^k,1);

if k == 0
    M = eye(n) + 2*alpha*K + (alpha^2 + beta^2)*K^2; % M is n x n
    y = M\d; % y is n x 1, and d is n x 1

else
    jj = 1;
    while (jj <= m)
        if jj == m || abs(F(jj,jj+1)) == 0

            ind = 1 + (jj-1)*n*m^(k-1):jj*n*m^(k-1);
            yj = solv2p(F(jj,jj)*alpha,F(jj,jj)^2*beta^2,F,K,d(ind),k-1);
            y(ind) = yj;
            z = KPvecprod(F,K,k-1,yj);
            z = 2*alpha*z;
            w = KPvecprod(F^2,K^2,k-1,yj);
            w = (alpha^2 + beta^2)*w;
            
            g = F^2;

            for kk = jj+1:m 
                ind = 1 + (kk-1)*n*m^(k-1):kk*n*m^(k-1);
                d(ind) = d(ind) - F(kk,jj)*z - g(kk,jj)*w;
            end
             jj = jj + 1;
        else

            gamma = F(jj,jj);
            delta1 = F(jj,jj+1);
            delta2 = F(jj+1,jj);
            E = [gamma, -delta1; -delta2, gamma];
            
            ind1 = 1 + (jj-1)*n*m^(k-1):jj*n*m^(k-1);
            ind2 = 1 + jj*n*m^(k-1):(jj+1)*n*m^(k-1);

            d_hat1 = KPvecprod2(E^2,F^2,K^2,k-1,[d(ind1);d(ind2)]);
            d_hat2 = KPvecprod2(E,F,K,k-1,[d(ind1);d(ind2)]);
            d_hat =  [d(ind1);d(ind2)] + 2*alpha*d_hat2 + (alpha^2 + beta^2)*d_hat1;
            
            delta = sqrt(-delta1*delta2);
            
            a1 = alpha*gamma - beta*delta;
            b1 = alpha*delta + gamma*beta;
            a2 = alpha*gamma + beta*delta;
            b2 = alpha*delta - gamma*beta;
            
            y1_temp = solv2p(a2,b2^2,F,K,d_hat(1:n*m^(k-1)),k-1);
            y2_temp = solv2p(a2,b2^2,F,K,d_hat(n*m^(k-1)+1:2*n*m^(k-1)),k-1);
            
            y1 = solv2p(a1,b1^2,F,K,y1_temp,k-1);
            y2 = solv2p(a1,b1^2,F,K,y2_temp,k-1);
            y(ind1) = y1;
            y(ind2) = y2;
            
            z1 = KPvecprod(F,K,k-1,y1);
            z1 = 2*alpha*z1;
            
            z2 = KPvecprod(F,K,k-1,y2);
            z2 = 2*alpha*z2;

            w1 = KPvecprod(F^2,K^2,k-1,y1);
            w1 = (alpha^2+beta^2)*w1;
            
            w2 = KPvecprod(F^2,K^2,k-1,y2);
            w2 = (alpha^2+beta^2)*w2;
            
            g = F^2;

            for kk = jj + 2:m 
                ind = 1 + (kk-1)*n*m^(k-1):kk*n*m^(k-1);
                d(ind) = d(ind) - F(kk,jj)*z1 - F(kk,jj+1)*z2 - g(kk,jj)*w1 - g(kk,jj+1)*w2;
            end
             jj = jj + 2;

        end
    end
end