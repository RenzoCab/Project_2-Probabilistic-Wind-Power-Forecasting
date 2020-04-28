close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet       = 'MTLOG';
% epsilon can be 0.035 (B), 0.020 (A) or 0.028 (C).
epsilon       = 0.020;
% dataKind can be classic or comparable.
dataKind      = 'comparable';
[Ta_Tra_Comp] = load_data_eps(epsilon,dataSet,dataKind);

Date              = Ta_Tra_Comp.Date;
Time              = Ta_Tra_Comp.Time;
Forecast          = Ta_Tra_Comp.Forecast;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
Real_ADME         = Ta_Tra_Comp.Real_ADME;
Error             = Ta_Tra_Comp.Error;
Error_Transitions = Ta_Tra_Comp.Error_Transitions;
Lamparti_Data     = Ta_Tra_Comp.Error_Lamp;
Lamparti_Tran     = Ta_Tra_Comp.Error_Lamp_Transitions;

dt         = Time(1,2);
[M, N_ini] = size(Forecast);
N          = N_ini - 1; % We have N_ini measurements but N samples.

num_days   = height(Ta_Tra_Comp); % Maximum 127 for MTLOG.
[~, batch] = new_batch_fixed(Ta_Tra_Comp,num_days,N);

%% Initial parameters:

est           = initial_guess(Real_ADME, M, N, dt);
gamma         = 0.35; 
samples_gamma = create_samples_minus_eps(Forecast, Error, gamma);
theta_ini     = mean_regression_eps(samples_gamma, dt);
alpha_ini     = est/theta_ini;

%% Optimization:

theta_0 = theta_ini;
alpha   = alpha_ini;
    
x0  = [theta_0, alpha, 0.05];
fun = @(x) -likelihood_optimization_all(batch, x(1), x(2), x(3), dt, N);

x = fminsearch(fun, x0);
% Initial Delta: 0.05.
% Theta = 1.4617, Alpha = 0.066422, eps = 0.020, prod = 0.0971, Delta = 0.084315.

%% Many initial guesses for Delta:

theta_0   = theta_ini;
alpha     = alpha_ini;
Delta_vec = [0:0.05:1];
fun       = @(x) -likelihood_optimization_all(batch, x(1), x(2), x(3), dt, N);

parfor i = 1:length(Delta_vec)

    x0           = [theta_0, alpha, Delta_vec(i)];
    x            = fminsearch(fun, x0);
    results(i,:) = x;

end
%     Theta_0   Alpha     Delta     Product
%     2.2210    0.0436    0.0856    0.0969
%     1.4617    0.0664    0.0844    0.0971
%     1.3647    0.0712    0.0839    0.0971
%     1.4607    0.0665    0.0844    0.0971
%     1.4422    0.0673    0.0842    0.0971
%     1.4653    0.0663    0.0842    0.0971
%     1.4585    0.0665    0.0843    0.0970
%     1.2390    0.0784    0.0837    0.0972
%     2.6163    0.0370    0.0901    0.0968
%     2.7354    0.0353    0.0877    0.0964
%     2.7330    0.0354    0.0910    0.0967
%     2.7365    0.0354    0.0922    0.0968
%     2.7901    0.0347    0.0913    0.0967
%     2.7938    0.0350    0.0923    0.0977
%     2.7364    0.0354    0.0909    0.0968
%     2.3509    0.0412    0.0883    0.0968
%     2.3509    0.0412    0.0880    0.0969
%     2.0711    0.0466    0.7088    0.0965
%     2.3536    0.0411    0.0882    0.0968
%     2.4210    0.0400    0.0888    0.0969
%     2.3508    0.0412    0.0894    0.0968