close all;
clear all;
clc;

% dataSet can be 1 (Gaussian Proxy), 2 (Shoji-Ozaki) or 3 (Beta Proxy).
approx        = 3;
% dataSet can be AWSTP, MTLOG or UTEP5.
dataSet       = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon       = 0.018;
[Ta_Tra_Comp] = load_data_eps(epsilon,dataSet);

Table_Training = Ta_Tra_Comp; % We copy it so we can modify it.

Date              = Ta_Tra_Comp.Date;
Time              = Ta_Tra_Comp.Time;
Forecast          = Ta_Tra_Comp.Forecast;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
Real_ADME         = Ta_Tra_Comp.Real_ADME;
Error             = Ta_Tra_Comp.Error;
Error_Transitions = Ta_Tra_Comp.Error_Transitions;
Lamparti_Data     = Ta_Tra_Comp.Error_Lamp;
Lamparti_Tran     = Ta_Tra_Comp.Error_Lamp_Transitions;

dt             = Time(1,2);
[M, N_ini]     = size(Forecast);
N              = N_ini - 1; % We have N_ini measurements but N samples.

num_days = height(Ta_Tra_Comp); % Maximum 127 for MTLOG.
[Table_Training, batch] = new_batch_fixed(Table_Training,num_days,N);

theta_0_vec = -5:0.2:5;

for k = 1:length(theta_0_vec)

    for i = 1:length(batch(1,:))/2

        j = i*2;
        p1 = batch(1,j-1); p2 = batch(1,j);
        v1 = batch(3,j-1); v2 = batch(3,j);
        n  = 100;
        x1 = v1 + p1;
        x2 = v2 + p2;

        alpha   = 0.093 / theta_0_vec(k);
        theta_0 = theta_0_vec(k);

        for j = 1:200
            sim_path_f(j) = sde_FE_MODEL1(x1,alpha,theta_0,dt,p1);
        end
        e_x1(i) = (mean(sim_path_f)-x2)^2;

    end
    
    toPlot(k) = mean(e_x1);

end

plot(theta_0_vec,toPlot);
grid on;
