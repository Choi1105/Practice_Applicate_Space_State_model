function generateMfuncNum1st_file(nameSteadyState,strucOfArrays,nameOfFunction)
% This function generates an M function which evaluates all first order 
% derivatives to a DSGE model at the steady state.  
% The results are saved in the file nameOfFunction, for instance 
% nameOfFunction = 'numDerivF_1st.m';
% In this version of generateMfuncNum1st_file we write the output directly
% to the file and not to the screen and then to the file, as in
% generateMfuncNum1st.

% We start deleting the old version of the file - if it exists
if exist(nameOfFunction,'file') > 0
    delete(nameOfFunction)
end

text = ['function [f,fx,fxp,fy,fyp,eta,errorMes] = ',nameOfFunction(1:end-2),'(params)'];
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
text = 'n   = ny+nx;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'f   = zeros(n,1);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fx  = zeros(n,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fxp = zeros(n,nx);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fy  = zeros(n,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'fyp = zeros(n,ny);';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = ' ';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

names = fieldnames(strucOfArrays);
nn = size(names,1);
for i=1:nn
    DispSymMatrixMatlab_file(strucOfArrays.(names{i}),names{i},nameOfFunction)
    text = ['if any(any(isnan(',names{i},'),1),2) ~= 0 || any(any(isinf(',names{i},'),1),2) ~= 0 || any(any(imag(',names{i},'),1),2) ~= 0'];
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = '   errorMes = max(errorMes,1);';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = 'end';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
    text = ' ';
    dlmwrite(nameOfFunction,text,'-append','delimiter','');
end

text = 'if errorMes == 1';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fx  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fxp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fy  = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    fyp = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = '    eta = NaN;';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
text = 'end';
dlmwrite(nameOfFunction,text,'-append','delimiter','');

text = 'end';
dlmwrite(nameOfFunction,text,'-append','delimiter','');
end





