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

gamma_vec = [0:0.05:0.46];

for i = 1:length(gamma_vec)

    disp(num2str(i/length(gamma_vec)));
    gamma         = gamma_vec(i); 
    samples_gamma = create_samples_minus_eps(Forecast, Error, gamma);
    rang          = [-4:0.1:4];
    est           = initial_guess(Real_ADME, M, N, dt);
    theta_ini     = mean_regression_M1(batch, dt, rang, gamma_vec(i));
    theta_ini     = rang(theta_ini == min(theta_ini(:)));
    alpha_ini     = est/theta_ini;
    theta_vec(i)  = theta_ini;
    alpha_vec(i)  = alpha_ini;
    
end

plot(gamma_vec,theta_vec);
hold on;
plot(gamma_vec,alpha_vec);
grid minor; xlabel('Gamma');
sup_lim = max([max(theta_vec) max(alpha_vec)]);
% ylim([0 sup_lim]);
xlim([min(gamma_vec) max(gamma_vec)]);
legend('Theta_0','Alpha');