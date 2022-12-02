function generateMfuncNum3rd_file(nameSteadyState,strucOfArrays,nameOfFunction)
% This function generates an M function which evaluates all second order 
% derivatives to a DSGE model at the steady state.  
% The results are saved in the file nameOfFunction, for instance 
% nameOfFunction = 'numDerivF_3rd.m';
% In this version of generateMfuncNum3rd_file we write the output directly
% to the file and not to the screen and then to the file, as in
% generateMfuncNum3rd. 

% We start deleting the old version of the file - if it exists
if exist(nameOfFunction,'file') > 0
    delete(nameOfFunction)
end

text = ['function [fxxx,fxxxp,fxxy,fxxyp,fxxpxp,fxxpy,fxxpyp,fxyy,fxyyp,fxypyp,fxpxpxp,fxpxpy,fxpxpyp,fxpyy,fxpyyp,fxpypyp,fyyy,fyyyp,fyypyp,fypypyp,errorMes] = ',nameOfFunction(1:end-2),'(params)'];
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '%Initializing the flag for errors';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'errorMes = 0;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = ' ';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '%The steady state of the model';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = [nameSteadyState,';'];
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = ' ';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = '%Setting dimension for the matrices';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'n      = ny+nx;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '%For fx:';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxxx   = zeros(n,nx,nx,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxxxp  = zeros(n,nx,nx,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxxy   = zeros(n,nx,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxxyp  = zeros(n,nx,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fxxpxp = zeros(n,nx,nx,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxxpy  = zeros(n,nx,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxxpyp = zeros(n,nx,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fxyy   = zeros(n,nx,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxyyp  = zeros(n,nx,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fxypyp = zeros(n,nx,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '%For fxp:';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxpxpxp= zeros(n,nx,nx,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxpxpy = zeros(n,nx,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxpxpyp= zeros(n,nx,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fxpyy  = zeros(n,nx,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxpyyp = zeros(n,nx,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fxpypyp= zeros(n,nx,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '%For fy:';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fyyy   = zeros(n,ny,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fyyyp  = zeros(n,ny,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fyypyp = zeros(n,ny,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '%For fyp:';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fypypyp = zeros(n,ny,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = ' ';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

names = fieldnames(strucOfArrays);
nn = size(names,1);
for i=1:nn
    DispSymMatrixMatlab_file(strucOfArrays.(names{i}),names{i},nameOfFunction)
    text = ['if any(any(any(any(isnan(',names{i},'),1),2),3),4) ~= 0 || any(any(any(any(isinf(',names{i},'),1),2),3),4) ~= 0 || any(any(any(any(imag(',names{i},'),1),2),3),4) ~= 0'];
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = '   errorMes = max(errorMes,1);';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = 'end';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = ' ';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
end

% Some testing of the derivatives
text = 'if errorMes == 1';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxx   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxxp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxy   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxyp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxpxp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxpy  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxpyp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxyy   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxyyp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxypyp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpxpxp= NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpxpy = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpxpyp= NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpyy  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpyyp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpypyp= NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyyy   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyyyp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyypyp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fypypyp= NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'end';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'end';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
end





