function [yd]=ldtrend(y)
%TREND extract a linear trend
%	y-yd is the detrended series
[T,k]=size(y);
tdum = 1:T;
const=ones(T,1);
X=[const tdum'];
bet=inv(X'*X)*X'*y;
yd=y-X*bet;