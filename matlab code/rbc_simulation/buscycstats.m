function buscyctable = buscycstats(data_d)

T = size(data_d,1);
n = size(data_d,2);

%Create business cycle statistics table
%Column 1: Std dev(x)
%Column 2: corr(x,GDP)
%Column 3: corr(x,Trade Balance)
%Column 4: corr(x,x_1)

buscyctable      = zeros(n,4);
buscyctable(:,1) = 100*(var(data_d).^0.5)';
for i = 1:n
    buscyctable(i,2) = corr(data_d(:,i),data_d(:,1));
    buscyctable(i,3) = corr(data_d(:,i),data_d(:,4));
    buscyctable(i,4) = corr(data_d(1:T-1,i),data_d(2:T,i));
end