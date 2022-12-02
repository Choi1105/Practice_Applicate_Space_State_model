% Copyright Andrew Binning 2015
% From:
% Binning, A. (2013). Solving second and third-order approximations to dsge models: a recur-
% sive sylvester equation solution. Norges Bank Working Paper 2013/18. 
% URL http://www.norges-bank.no/no/om/publisert/publikasjoner/working-papers/2013/13/.
% Based on:
% Martin, C. D. M. & Van Loan, C. F. (2006). Shifted kronecker product systems. SIAM J. Matrix
% Analysis Applications, 29 (1), 184-198.

function y = KPvecprod2(E,F,K,k,x)

n = length(K);
m = length(F);
p = length(E);

N = p*n*m^k;

x = (K*reshape(x,n,[])).';

if k > 0
    
    for i=1:k
        x = (F*reshape(x,m,[])).';
    end
    
end

x = (E*reshape(x,p,[])).';

y = reshape(x,N,1);
