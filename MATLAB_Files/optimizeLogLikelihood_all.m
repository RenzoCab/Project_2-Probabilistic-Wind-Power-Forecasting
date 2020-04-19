close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet = 'UTEP5';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon       = 0.028;
[Ta_Tra_Comp] = load_data_eps(epsilon,dataSet);

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
% (A) Theta = 1.2627, Alpha = 0.0779, prod = 0.0983.
% (B) Theta = 0.7398, Alpha = 0.1383, prod = 0.1023.
% (C) Theta = 1.4666, Alpha = 0.0717, prod = 0.1051.

%% Optimization:

theta_0 = theta_ini;
alpha   = alpha_ini;
    
x0    = [theta_0, alpha, 0.2];
fun   = @(x) -likelihood_optimization_all(batch, x(1), x(2), x(3), dt, N);

x = fminsearch(fun, x0);
% fminsearch:
% (A) Theta = 0.95907, Alpha = 0.086308, eps = 0.018, prod = 0.0828.
% (B) Theta = 0.76586, Alpha = 0.11728,  eps = 0.035, prod = 0.0898.
% (C) Theta = 1.347,   Alpha = 0.073388, eps = 0.028, prod = 0.0989.

