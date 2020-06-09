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
% (A) Theta = 1.5441, Alpha = 0.0717, prod = 0.1107.
% (B) Theta = 1.0657, Alpha = 0.1039, prod = 0.1107.
% (C) Theta = 1.1625, Alpha = 0.0953, prod = 0.1108.

%% Optimization:

theta_0 = theta_ini;
alpha   = alpha_ini;

if strcmp(likelihood,'normal')
    
    x0    = [theta_0, alpha];
    fun   = @(x) -likelihood_optimization(batch, x(1), x(2), dt);
    
    [x,f] = fminunc(fun, x0);
    
    % fminsearch -> Theta_0 = 1.278, Alpha = 0.0760, prod = 0.0971 (24/04/2020).
    % fmincon    -> Theta_0 = 1.577, Alpha = 0.0615, prod = 0.0970 (24/04/2020).
    % fminunc    -> Theta_0 = 1.544, Alpha = 0.0628, prod = 0.0970 (24/04/2020).
    % eval       -> Theta_0 = 1.930, Alpha = 0.0500, prod = 0.0965 (20/05/2020).
    
    % For the testing data:
    % Initial    -> Theta_0 = 1.963, Alpha = 0.0527, prod = 0.1036 (09/06/2020).
    % fminsearch -> Theta_0 = 1.642, Alpha = 0.0542, prod = 0.0890 (08/06/2020).
    % fmincon    -> Theta_0 = 1.632, Alpha = 0.0545, prod = 0.0890 (08/06/2020).
    % fminunc    -> Theta_0 = 1.963, Alpha = 0.0454, prod = 0.0891 (08/06/2020).
    % eval       -> Theta_0 = 1.590, Alpha = 0.0560, prod = 0.0890 (08/06/2020).
    
elseif strcmp(likelihood,'lamperti')
    
    x0 = [0.1 0.1];
    x1 = [1.14, 0.073]; % Output from 'normal' optimization.
    
    while norm(x1-x0)/norm(x1) > 0.01
        
        x0        = x1;
        [x1,fval] = likelihood_nested_optimization_L(batch, x0, dt);
    
        % Theta_0 = 1.57, Alpha = 0.034 (04/04/2020).
    
    end
    
else
    
	disp('Please, choose normal or lamperti.');
    
end