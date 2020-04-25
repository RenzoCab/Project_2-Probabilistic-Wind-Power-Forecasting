close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon       = 0.018;
approx        = 1;
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

gamma         = 0; 
samples_gamma = create_samples_minus_eps(Forecast, Error, gamma);
rang          = [-6:0.01:6];
est           = initial_guess(Real_ADME, M, N, dt);
theta_ini     = mean_regression_M1(batch, dt, rang, gamma);
theta_ini     = rang(theta_ini == min(theta_ini(:)));
alpha_ini     = est/theta_ini;

%% Optimization:

theta_0 = theta_ini;
alpha   = alpha_ini;
    
x0    = [theta_0, alpha];
fun   = @(x) -likelihood_optimization_MODEL1(batch, x(1), x(2), dt, approx);

[x,f] = fminsearch(fun, x0);
% fminsearch:
% (A) Theta_0 = -1.2201, Alpha = -0.075142, eps = 0.018, Product = 0.091681, f = -52342.
% (B) Theta = 0.02584, Alpha = 3.7390, eps = 0.035, prod = 0.0966, f = -30347.
% (C) Theta = 0.02479, Alpha = 3.9996, eps = 0.028, prod = 0.0992, f = -50455.