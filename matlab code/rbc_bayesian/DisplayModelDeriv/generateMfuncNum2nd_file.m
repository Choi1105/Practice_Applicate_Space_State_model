function generateMfuncNum2nd_file(nameSteadyState,strucOfArrays,nameOfFunction)
% This function generates an M function which evaluates all second order 
% derivatives to a DSGE model at the steady state.  
% The results are saved in the file nameOfFunction, for instance 
% nameOfFunction = 'numDerivF_2nd.m';
% In this version of generateMfuncNum2nd_file we write the output directly
% to the file and not to the screen and then to the file, as in
% generateMfuncNum2nd. 

% We start deleting the old version of the file - if it exists
if exist(nameOfFunction,'file') > 0
    delete(nameOfFunction)
end


text = ['function [fxx,fxxp,fxy,fxyp,fxpx,fxpxp,fxpy,fxpyp,fyx,fyxp,fyy,fyyp,fypx,fypxp,fypy,fypyp,errorMes] = ',nameOfFunction(1:end-2),'(params)'];
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = ' ';
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
text = 'n     = ny+nx;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxx   = zeros(n,nx,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxxp  = zeros(n,nx,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxy   = zeros(n,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxyp  = zeros(n,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fxpxp = zeros(n,nx,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxpy  = zeros(n,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxpyp = zeros(n,nx,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fyy   = zeros(n,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fyyp  = zeros(n,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'fypyp = zeros(n,ny,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = ' ';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

names = fieldnames(strucOfArrays);
nn = size(names,1);
for i=1:nn
    DispSymMatrixMatlab_file(strucOfArrays.(names{i}),names{i},nameOfFunction)
    text = ['if any(any(any(isnan(',names{i},'),1),2),3) ~= 0 || any(any(any(isinf(',names{i},'),1),2),3) ~= 0 || any(any(any(imag(',names{i},'),1),2),3) ~= 0'];
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = '   errorMes = max(errorMes,1);';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = 'end';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = ' ';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
end

% The remaining derivatives by symmetry
text = 'fxpx  = permute(fxxp,[1,3,2]);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fyx   = permute(fxy,[1,3,2]);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fyxp  = permute(fxpy,[1,3,2]);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fypx  = permute(fxyp,[1,3,2]);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fypxp = permute(fxpyp,[1,3,2]);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fypy  = permute(fyyp,[1,3,2]);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

% Some testing of the derivatives
text = ' ';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'if errorMes == 1';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxx   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxxp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxy   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxyp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpx  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpxp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpy  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxpyp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyx   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyxp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyy   = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyyp  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fypx  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fypxp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fypy  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fypyp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'end';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'end';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
end





