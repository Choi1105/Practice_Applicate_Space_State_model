addpath(genpath(pwd))
clc
% **************************** USER SETTING *******************************
order_app = 3;         % The order of analytical derivatives which we compute
logApprox = 1;                        
% *************************************************************************

%% Model specific part
% The DSGE model
unitFree = 0;
[f,x,xp,y,yp,logTransY,logTransX] = growth_model(unitFree,logApprox);

% A string with the output and name of the function computing the steady
% state
nameSteadyState = '[GAMA, DELTA, ALFA, PSSI, OMEGA, SIGMAA, RHOA, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, YY, la, lap, a, ap, tby, tbyp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, nx, ny, ne, sig, eta]=growth_model_ss(params);';

%% The analytical derivatives of the model
[fx,fxp,fy,fyp,...                                                 %first order derivatives
    fxx,fxxp,fxy,fxyp, fxpxp,fxpy,fxpyp, fyy,fyyp, fypyp,...                %second order derivatives
    fxxx,fxxxp,fxxy,fxxyp , fxxpxp,fxxpy,fxxpyp, fxyy,fxyyp, fxypyp, ...    %third order derivatives related to: fxx,fxxp,fxy,fxyp   
    fxpxpxp,fxpxpy,fxpxpyp, fxpyy,fxpyyp, fxpypyp, ...                      %third order derivatives related to: fxpxp,fxpy,fxpyp
    fyyy,fyyyp, fyypyp, ...                                                 %third order derivatives related to: fyy,fyyp
    fypypyp...                                                              %third order derivatives related to: fypyp
    ] = Anal_derivatives(f,x,xp,y,yp,order_app);

% The first order derivatives
strucOfArraysF1st = struct('f',f,'fx',fx,'fxp',fxp,'fy',fy,'fyp',fyp);
if logApprox == 1
    nameOfFunction = 'numDerivF_1st.m';
end
generateMfuncNum1st_file(nameSteadyState,strucOfArraysF1st,nameOfFunction)

% The second order derivatives
strucOfArraysF2nd = struct('fxx',fxx,'fxxp',fxxp,'fxy',fxy,'fxyp',fxyp,...
    'fxpxp',fxpxp,'fxpy',fxpy,'fxpyp',fxpyp,...
    'fyy',fyy,'fyyp',fyyp,'fypyp',fypyp);
if logApprox == 1
    nameOfFunction = 'numDerivF_2nd.m';
end
generateMfuncNum2nd_file(nameSteadyState,strucOfArraysF2nd,nameOfFunction);

% The third order derivatives
strucOfArraysF3rd = struct('fxxx',fxxx,'fxxxp',fxxxp,'fxxy',fxxy,'fxxyp',fxxyp,...
    'fxxpxp',fxxpxp,'fxxpy',fxxpy,'fxxpyp',fxxpyp,...
    'fxyy',fxyy,'fxyyp',fxyyp,...
    'fxypyp',fxypyp,...
    'fxpxpxp',fxpxpxp,'fxpxpy',fxpxpy,'fxpxpyp',fxpxpyp,...
    'fxpyy',fxpyy,'fxpyyp',fxpyyp,...
    'fxpypyp',fxpypyp,...
    'fyyy',fyyy,'fyyyp',fyyyp,'fyypyp',fyypyp,'fypypyp',fypypyp);
if logApprox == 1
    nameOfFunction = 'numDerivF_3rd.m';
end
generateMfuncNum3rd_file(nameSteadyState,strucOfArraysF3rd,nameOfFunction);








