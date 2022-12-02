function corrtable = crosscorrGDP(data_d,numlags)

T = size(data_d,1);
n = size(data_d,2);

%Create business cycle statistics table
%Column 1: (Std dev(x))/(Std dev(GDP))
%Other columns: corr(GDP(t),x(tau)) for tau = t - nlag, ..., t + nlag

corrtable      = zeros(n,1 + 2*numlags + 1);
corrtable(:,1) = 100*(var(data_d).^0.5)';
corrtable(:,1) = corrtable(:,1)/corrtable(1,1);
for i = 1:n
    for j = 1:2*numlags + 1
        lag                = circshift(data_d(:,i),numlags - j + 1);
        corrtable(i,j + 1) = corr(data_d(1 + numlags:T - numlags,1),lag(1 + numlags:T - numlags));
    end
end