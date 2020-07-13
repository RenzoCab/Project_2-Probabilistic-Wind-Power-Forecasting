close all;
clear all;
clc;

% Lamperti has no sense here. I wanted to plot it as a function of the
% parameters. Still, as the Lamperti transform also depends on the 
% parameters, the data set changed in each evaluation of the Likelihood,
% which has no sense. See plotLamperti.m.

likelihood    = 'lamperti'; % 'normal' or 'lamperti'.
% dataSet can be AWSTP, MTLOG or UTEP5.
dataSet = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon       = 0.018;
% dataKind can be classic or comparable.
dataKind   = 'comparable';
dataNumber = 'training';
[Ta_Tra_Comp] = load_data_eps(epsilon,dataSet,dataKind,dataNumber);

% AWSTP) (B)
% T_0 = 1.42, Alpha = 0.069, prod = 0.0980, f = -36753.
% AIC = -73502, BIA = -73487.
% MTLOG) (A)
% T_0 = 1.93, Alpha = 0.050, prod = 0.0965, f = -36852.
% AIC = -73700, BIA = -73685.
% UTEP5) (B)
% T_0 = 1.38, Alpha = 0.078, prod = 0.1076, f = -36261.
% AIC = -72518, BIA = -72503.

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

%% Initial parameters:

est           = initial_guess(Real_ADME, M, N, dt);
gamma         = 0.35; 
samples_gamma = create_samples_minus_eps(Forecast, Error, gamma);
theta_ini     = mean_regression_eps(samples_gamma, dt);
alpha_ini     = est/theta_ini;
% AWSTP) (B)
% T_0 = 1.0657, Alpha = 0.1039, prod = 0.1108.
% MTLOG) (A)
% T_0 = 1.5441, Alpha = 0.0717, prod = 0.1108.
% UTEP5) (B)
% T_0 = 1.1625, Alpha = 0.0953, prod = 0.1108.

%% Log-Likelihood plot:

val = [];

vec_theta = 0.1:0.01:5;
vec_alpha = 0.01:0.001:1;
num_divs  = 10000;

for i = 1:length(vec_theta)

    vec = ones(1,length(vec_alpha));

    parfor j = 1:length(vec_alpha)
        
        theta_0 = vec_theta(i);
        alpha   = vec_alpha(j);
        
        if strcmp(likelihood,'normal')
            
            batch_complete = batch_with_theta(batch, alpha, theta_0);
            vec(j) = -log_LH_evaluation(batch_complete, alpha, theta_0, dt);
            
        elseif strcmp(likelihood,'lamperti')
            
            batch_complete = batch_with_theta_L(batch, alpha, theta_0);
            vec(j) = -log_LH_evaluation_L(batch_complete, theta_0, alpha, dt);
            
        end
        
    end

    val(i,:) = vec(:);

    disp(num2str(i/length(vec_theta)*100));

end

theta_0_opt_1 = 1.278;
alpha_opt_1   = 0.076;
theta_0_opt_2 = 1.577;
alpha_opt_2   = 0.062;
theta_0_opt_3 = 1.544;
alpha_opt_3   = 0.063;

% To save the matrix with the evaluations, we use:
% save([pwd '/Results/likelihood/',likelihood,'/Log-Likelihood.mat'],'val');
% To load the matrix with the evaluations, we use:
% load([pwd '/Results/likelihood/',likelihood,'/Log-Likelihood.mat'],'val');

% To save the matrix with the evaluations, we use:
% save([pwd '/Results/likelihood/',likelihood,'/Log-Likelihood_lam.mat'],'val');
% To load the matrix with the evaluations, we use:
% load([pwd '/Results/likelihood/',likelihood,'/Log-Likelihood_lam.mat'],'val');

val_norm = (val-min(min(val)))/max(max((val-min(min(val)))));

minMatrix = min(val(:));
[row,col] = find(val==minMatrix);
theta_0_opt_4 = vec_theta(row);
alpha_opt_4   = vec_alpha(col);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on;
[X,Y] = meshgrid(vec_theta,vec_alpha);
% contourf(X,Y,val',num_divs); % colorbar;
% contourf(X,Y,val',[-14800 -14500 0 10000]); colorbar;
imagesc(vec_theta,vec_alpha,val');
xlim([min(vec_theta) max(vec_theta)]);
ylim([min(vec_alpha) max(vec_alpha)]);
% contourf(X,Y,val'); colorbar;
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Negative Log-Likelihoog for ',num2str(num_days),' days']);
if strcmp(likelihood,'normal')
    plot(theta_ini,alpha_ini,'-p','MarkerFaceColor','red','MarkerSize',12);
    plot(theta_0_opt_1,alpha_opt_1,'-p','MarkerFaceColor','cyan','MarkerSize',12);
    plot(theta_0_opt_2,alpha_opt_2,'-p','MarkerFaceColor','green','MarkerSize',12);
    plot(theta_0_opt_3,alpha_opt_3,'-p','MarkerFaceColor','yellow','MarkerSize',12);
    plot(theta_0_opt_4,alpha_opt_4,'-p','MarkerFaceColor','white','MarkerSize',12);
    legend('Level sets','Initial guess','Optimal value (fminsearch)',...
        'Optimal value (fmincon)','Optimal value (fminunc)','Optimal value (evaluations)');
elseif strcmp(likelihood,'lamperti')
    legend('Level sets');
end
pause(0.1);