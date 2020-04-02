close all;
clear all;
clc;

likelihood    = 'normal'; % 'normal' or 'lamperti'.
epsilon       = 0.018;
[Ta_Tra_Comp] = load_data_eps(epsilon);

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

num_days   = 127; % Maximum 127.
[~, batch] = new_batch_fixed(Ta_Tra_Comp,num_days,N);

%% Initial parameters:

est           = initial_guess(Real_ADME, M, N, dt);
gamma         = 0.35; 
samples_gamma = create_samples_minus_eps(Forecast, Error, gamma);
theta_ini     = mean_regression_eps(samples_gamma, dt);
alpha_ini     = est/theta_ini;

%% Optimization:

theta   = theta_ini;
alpha   = alpha_ini;

if strcmp(likelihood,'normal')
    
    x0      = [theta, alpha];
    fun     = @(x) -likelihood_optimization(batch, x(1), x(2), dt);
    options = optimset('PlotFcns',@optimplotfval);

    x = fminsearch(fun, x0, options);
    % Theta = 3.912, Alpha = 0.019308 (20/03/2020).

elseif strcmp(likelihood,'lamperti')
    
    x0      = [theta, alpha];
    fun     = @(x) -likelihood_optimization_L(batch, x(1), x(2), dt);
    options = optimset('PlotFcns',@optimplotfval);

    x = fminsearch(fun, x0, options);
    % Theta = 3.912, Alpha = 0.019308 (20/03/2020).
    
else
    
	disp('Please, choose normal or lamperti.');
    
end