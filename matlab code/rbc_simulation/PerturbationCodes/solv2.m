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

function y = solv2(alpha,beta1,beta2,F,K,d,k)

m = length(F);
n = length(K);

E = [alpha, -beta1; beta2, alpha];

dhat = KPvecprod2(E,F,K,k,d);

dhat = d + dhat;

y1 = solv2p(alpha,beta1*beta2,F,K,dhat(1:n*m^k),k);
y2 = solv2p(alpha,beta1*beta2,F,K,dhat(n*m^k+1:2*n*m^k),k);

y = [y1;y2];