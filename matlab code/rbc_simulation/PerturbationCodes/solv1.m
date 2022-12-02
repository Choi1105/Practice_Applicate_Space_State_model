% Copyright Andrew Binning 2015
% From:
% Binning, A. (2013). Solving second and third-order approximations to dsge models: a recur-
% sive sylvester equation solution. Norges Bank Working Paper 2013/18. 
% URL http://www.norges-bank.no/no/om/publisert/publikasjoner/working-papers/2013/18/.
% Based on: 
% Kamenik, O. (2005). Solving sdge models: A new algorithm for the sylvester equation.
% Working Papers 2005/10, Czech National Bank, Research Department.
% Martin, C. D. M. & Van Loan, C. F. (2006). Shifted kronecker product systems. SIAM J. Matrix
% Analysis Applications, 29 (1), 184-198.

function y = solv1(r,F,K,d,k)

m = length(F);
n = length(K);
y = zeros(n*m^k,1);

if k == 0
    
    M = [eye(n) + r*K];
    y = M\d;
    
else
    jj = 1;
    while (jj <= m)
        if jj == m || abs(F(jj,jj+1)) == 0
            ind = 1+(jj-1)*n*m^(k-1):jj*n*m^(k-1);
            yj = solv1(r*F(jj,jj),F,K,d(ind),k-1);
            y(ind) = yj;
            z = KPvecprod(F,K,k-1,yj);
            z = r*z;
            for kk = jj + 1:m
                ind = 1 + (kk-1)*n*m^(k-1):kk*n*m^(k-1);
                d(ind) = d(ind) - F(kk,jj)*z;
            end
            jj = jj + 1;
        else
            alpha = F(jj,jj);
            beta1 = F(jj,jj+1);
            beta2 = -F(jj+1,jj);
            ind = 1+(jj-1)*n*m^(k-1):(jj+1)*n*m^(k-1);
            yj = solv2(r*alpha,r*beta1,r*beta2,F,K,d(ind),k-1);
            y(ind) = yj;
            z1 = KPvecprod(F,K,k-1,yj(1:n*m^(k-1)));
            z1 = r*z1;
            z2 = KPvecprod(F,K,k-1,yj(1+n*m^(k-1):2*n*m^(k-1)));
            z2 = r*z2;
            for kk = jj+2:m
                ind = 1 + (kk-1)*n*m^(k-1):kk*n*m^(k-1);
               d(ind) = d(ind) - F(kk,jj)*z1 - F(kk,jj+1)*z2; 
            end
            jj = jj + 2;
        end
    end
end