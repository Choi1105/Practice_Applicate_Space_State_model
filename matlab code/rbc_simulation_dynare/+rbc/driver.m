%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'rbc';
M_.dynare_version = '5.1';
oo_.dynare_version = '5.1';
options_.dynare_version = '5.1';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps_a'};
M_.exo_names_tex(1) = {'eps\_a'};
M_.exo_names_long(1) = {'eps_a'};
M_.endo_names = cell(14,1);
M_.endo_names_tex = cell(14,1);
M_.endo_names_long = cell(14,1);
M_.endo_names(1) = {'c'};
M_.endo_names_tex(1) = {'c'};
M_.endo_names_long(1) = {'c'};
M_.endo_names(2) = {'k'};
M_.endo_names_tex(2) = {'k'};
M_.endo_names_long(2) = {'k'};
M_.endo_names(3) = {'a'};
M_.endo_names_tex(3) = {'a'};
M_.endo_names_long(3) = {'a'};
M_.endo_names(4) = {'h'};
M_.endo_names_tex(4) = {'h'};
M_.endo_names_long(4) = {'h'};
M_.endo_names(5) = {'d'};
M_.endo_names_tex(5) = {'d'};
M_.endo_names_long(5) = {'d'};
M_.endo_names(6) = {'y'};
M_.endo_names_tex(6) = {'y'};
M_.endo_names_long(6) = {'y'};
M_.endo_names(7) = {'invest'};
M_.endo_names_tex(7) = {'invest'};
M_.endo_names_long(7) = {'invest'};
M_.endo_names(8) = {'tb'};
M_.endo_names_tex(8) = {'tb'};
M_.endo_names_long(8) = {'tb'};
M_.endo_names(9) = {'mu_c'};
M_.endo_names_tex(9) = {'mu\_c'};
M_.endo_names_long(9) = {'mu_c'};
M_.endo_names(10) = {'tb_y'};
M_.endo_names_tex(10) = {'tb\_y'};
M_.endo_names_long(10) = {'tb_y'};
M_.endo_names(11) = {'g_y'};
M_.endo_names_tex(11) = {'g\_y'};
M_.endo_names_long(11) = {'g_y'};
M_.endo_names(12) = {'g_c'};
M_.endo_names_tex(12) = {'g\_c'};
M_.endo_names_long(12) = {'g_c'};
M_.endo_names(13) = {'g_invest'};
M_.endo_names_tex(13) = {'g\_invest'};
M_.endo_names_long(13) = {'g_invest'};
M_.endo_names(14) = {'r'};
M_.endo_names_tex(14) = {'r'};
M_.endo_names_long(14) = {'r'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
M_.param_names(1) = {'gamma'};
M_.param_names_tex(1) = {'gamma'};
M_.param_names_long(1) = {'gamma'};
M_.param_names(2) = {'delta'};
M_.param_names_tex(2) = {'delta'};
M_.param_names_long(2) = {'delta'};
M_.param_names(3) = {'alpha'};
M_.param_names_tex(3) = {'alpha'};
M_.param_names_long(3) = {'alpha'};
M_.param_names(4) = {'psi'};
M_.param_names_tex(4) = {'psi'};
M_.param_names_long(4) = {'psi'};
M_.param_names(5) = {'omega'};
M_.param_names_tex(5) = {'omega'};
M_.param_names_long(5) = {'omega'};
M_.param_names(6) = {'theta'};
M_.param_names_tex(6) = {'theta'};
M_.param_names_long(6) = {'theta'};
M_.param_names(7) = {'rho_a'};
M_.param_names_tex(7) = {'rho\_a'};
M_.param_names_long(7) = {'rho_a'};
M_.param_names(8) = {'RSTAR'};
M_.param_names_tex(8) = {'RSTAR'};
M_.param_names_long(8) = {'RSTAR'};
M_.param_names(9) = {'beta'};
M_.param_names_tex(9) = {'beta'};
M_.param_names_long(9) = {'beta'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 14;
M_.param_nbr = 9;
M_.orig_endo_nbr = 14;
M_.aux_vars = [];
M_.predetermined_variables = [ 2 5 ];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 14;
M_.eq_nbr = 14;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 7 0;
 2 8 0;
 3 9 21;
 0 10 22;
 4 11 0;
 5 12 0;
 6 13 0;
 0 14 0;
 0 15 23;
 0 16 0;
 0 17 0;
 0 18 0;
 0 19 0;
 0 20 0;]';
M_.nstatic = 6;
M_.nfwrd   = 2;
M_.npred   = 5;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 6;
M_.ndynamic   = 8;
M_.dynamic_tmp_nbr = [9; 3; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'r' ;
  2 , 'name' , 'mu_c' ;
  3 , 'name' , 'y' ;
  4 , 'name' , '4' ;
  5 , 'name' , '5' ;
  6 , 'name' , 'invest' ;
  7 , 'name' , '7' ;
  8 , 'name' , '8' ;
  9 , 'name' , '9' ;
  10 , 'name' , '10' ;
  11 , 'name' , 'g_y' ;
  12 , 'name' , 'g_c' ;
  13 , 'name' , 'g_invest' ;
  14 , 'name' , '14' ;
};
M_.mapping.c.eqidx = [2 3 12 ];
M_.mapping.k.eqidx = [5 6 8 9 ];
M_.mapping.a.eqidx = [5 8 9 14 ];
M_.mapping.h.eqidx = [2 5 8 9 ];
M_.mapping.d.eqidx = [1 4 ];
M_.mapping.y.eqidx = [3 5 10 11 ];
M_.mapping.invest.eqidx = [3 6 13 ];
M_.mapping.tb.eqidx = [3 4 10 ];
M_.mapping.mu_c.eqidx = [2 7 9 ];
M_.mapping.tb_y.eqidx = [10 ];
M_.mapping.g_y.eqidx = [11 ];
M_.mapping.g_c.eqidx = [12 ];
M_.mapping.g_invest.eqidx = [13 ];
M_.mapping.r.eqidx = [1 4 7 ];
M_.mapping.eps_a.eqidx = [14 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 2 3 5 6 7 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(14, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(9, 1);
M_.endo_trends = struct('deflator', cell(14, 1), 'log_deflator', cell(14, 1), 'growth_factor', cell(14, 1), 'log_growth_factor', cell(14, 1));
M_.NNZDerivatives = [46; -1; -1; ];
M_.static_tmp_nbr = [9; 3; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(7) = 0.1;
rho_a = M_.params(7);
M_.params(1) = 2;
gamma = M_.params(1);
M_.params(2) = 0.1255088099999999;
delta = M_.params(2);
M_.params(3) = 0.32;
alpha = M_.params(3);
M_.params(5) = 1.6;
omega = M_.params(5);
M_.params(6) = 1.4*M_.params(5);
theta = M_.params(6);
M_.params(4) = 0.001;
psi = M_.params(4);
M_.params(8) = 1.1;
RSTAR = M_.params(8);
M_.params(9) = 1/M_.params(8);
beta = M_.params(9);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.05)^2;
options_.irf = 20;
options_.loglinear = true;
options_.order = 1;
var_list_ = {'g_y';'g_c';'g_invest';'tb_y'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'rbc_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbc_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbc_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbc_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbc_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbc_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbc_results.mat'], 'oo_recursive_', '-append');
end
disp('Note: 2 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
