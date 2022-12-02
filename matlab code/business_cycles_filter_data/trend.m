function [yd]=trend(y,phi)
%TREND uses the Hodrick-Prescott method of detrending to find yd where
%	y-yd is the detrended series
[T,k]=size(y);
A=[1+phi 	-2*phi 	phi 	zeros(1,T-3);
   -2*phi	1+5*phi -4*phi	phi	zeros(1,T-4);
		zeros(T-4,T);
   zeros(1,T-4)	phi	-4*phi	1+5*phi	-2*phi;
   zeros(1,T-3) phi -2*phi	1+phi];

for i=3:T-2;
	A(i,i-2)=phi;
	A(i,i-1)=-4*phi;
	A(i,i)=1+6*phi;
	A(i,i+1)=-4*phi;
	A(i,i+2)=phi;
end;
yd=pentle(A,y);
