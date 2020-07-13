close all;
clear all;
clc;

likelihood    = 'normal'; % 'normal' or 'lamperti'.
% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet       = 'MTLOG';
% epsilon can be 0.035 (B), 0.020 (A) or 0.028 (C).
epsilon       = 0.020;
% dataKind can be classic or comparable.
dataKind      = 'comparable';
dataNumber    = 'testing';
[Ta_Tra_Comp] = load_data_eps(epsilon,dataSet,dataKind,dataNumber);

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

x0    = [theta_0, alpha];
fun   = @(x) fixed_point_function(batch, x, dt);

[x,f] = fminsearch(fun, x0);

%     Thet_0 Alpha  Eval   Prod
% V1: 2.8150 0.1276 0.0487 0.3592
% V2: 1.8708 0.0431 0.0667 0.0806

%     Thet_0 Alpha  Eval   Prod (testing data - 08/06/2020)
% V1: 2.5656 0.1276 0.0324 0.0831
% V2: 3.9657 0.0431 0.0371 0.1471

% The values for V2 are also fixed points of V1 (see level sets plot).
