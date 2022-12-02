function [parameter]=transformdomain(x)
% parameter([1 3 5 7])=exp(x([1 3 5 7]));      % to make them positive    
parameter([1 3 5 6 7])=exp(x([1 3 5 6 7]));      % to make them positive    
parameter([2])= 1/(1+exp(x(2))); % to make them lie between 0 and 1
parameter([4])= 1/(1+exp(x(4)));
% parameter(6)  = x(6);
parameter=parameter';
