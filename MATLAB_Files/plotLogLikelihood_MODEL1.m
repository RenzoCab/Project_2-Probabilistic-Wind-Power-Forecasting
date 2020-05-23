close all;
clear all;
clc;

% Lamperti has no sense here. I wanted to plot it as a function of the
% parameters. Still, as the Lamperti transform also depends on the 
% parameters, the data set changed in each evaluation of the Likelihood,
% which has no sense. See plotLamperti.m.

% dataSet can be 1 (Gaussian Proxy), 2 (Shoji-Ozaki) or 3 (Beta Proxy).
approx        = 3;
% dataSet can be AWSTP, MTLOG or UTEP5.
dataSet       = 'MTLOG';
% epsilon can be 0.035 (B), 0.020 (A) or 0.028 (C).
epsilon       = 0.020;
% dataKind can be classic or comparable.
dataKind = 'comparable';
[Ta_Tra_Comp] = load_data_eps(epsilon,dataSet,dataKind);

% AWSTP) (B)
% 1) T_0 = 0.0400, Alpha = 2.6300, prod = 0.1052, f = -29115.
% AIC = -58226, BIA = -58211.
% 2) T_0 = 0.0400, Alpha = 2.6000, prod = 0.1040, f = -29115.
% AIC = -58226, BIA = -58211.
% 3) T_0 = 0.0400, Alpha = 2.6000, prod = 0.1040, f = -29146.
% AIC = -58288, BIA = -58273.
% MTLOG) (A)
% 1) T_0 = 0.0400, Alpha = 2.6300, prod = 0.1052, f = -29115.
% AIC = -58226, BIA = -58211.
% 2) T_0 = 0.0400, Alpha = 2.6000, prod = 0.1040, f = -29115.
% AIC = -58226, BIA = -58211.
% 3) T_0 = 0.0400, Alpha = 2.6000, prod = 0.1040, f = -29145.
% AIC = -58286, BIA = -58271.
% UTEP5) (C)
% 1) T_0 = 0.0400, Alpha = 2.6300, prod = 0.1052, f = -29115.
% AIC = -58226, BIA = -58211.
% 2) T_0 = 0.0400, Alpha = 2.6000, prod = 0.1040, f = -29115.
% AIC = -58226, BIA = -58211.
% 3) T_0 = 0.0400, Alpha = 2.6000, prod = 0.1040, f = -29145.
% AIC = -58286, BIA = -58271.

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

%% Log-Likelihood plot:

val = [];

vec_theta = 0.01:0.01:3;
vec_alpha = 0.01:0.01:3;
num_divs  = 10000;

for i = 1:length(vec_theta)

    vec = ones(1,length(vec_alpha));

    parfor j = 1:length(vec_alpha)
        
        theta_0        = vec_theta(i);
        alpha          = vec_alpha(j);
        batch_complete = batch_with_theta(batch, alpha, theta_0);
        vec(j)         = -log_LH_evaluation_MODEL1(batch_complete, alpha, theta_0, dt, approx);

    end

    val(i,:) = vec(:);

    disp(num2str(i/length(vec_theta)*100));

end

% To save the matrix with the evaluations, we use:
% save([pwd '/Results/likelihood/model1/Log-Likelihood.mat'],'val');
% To load the matrix with the evaluations, we use:
% load([pwd '/Results/likelihood/model1/Log-Likelihood.mat'],'val');

minMatrix = min(val(:));
[row,col] = find(val==minMatrix);
theta_0_opt_4 = vec_theta(row);
alpha_opt_4   = vec_alpha(col);
alpha_plot    = theta_0_opt_4*alpha_opt_4 ./ vec_theta;

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on;
[X,Y] = meshgrid(vec_theta,vec_alpha);
contourf(X,Y,val',num_divs); % colorbar;
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Negative Log-Likelihoog for ',num2str(num_days),' days']);
plot(vec_theta,alpha_plot,'LineWidth',2);
% plot(theta_ini,alpha_ini,'-p','MarkerFaceColor','red','MarkerSize',12);
plot(theta_0_opt_4,alpha_opt_4,'-p','MarkerFaceColor','white','MarkerSize',12);
legend('Level sets','Curve $\theta_0\alpha=C$','Optimal value (evaluations)',...
    'interpreter','Latex');
ylim([min(vec_alpha) max(vec_alpha)]);
pause(0.1);
% To save the plot and change the fonts size, we use:
set(gca,'FontSize',16);
saveas(gcf,[pwd '/Results/likelihood/model1/Log-Likelihood_',num2str(approx),'_',dataSet],'epsc');

%% Plot with the important line:

% alpha_1 = 0.0728./vec_theta;
% alpha_2 = 0.0828./vec_theta;
% alpha_3 = 0.0928./vec_theta;
% alpha_4 = 0.1028./vec_theta;
% alpha_5 = 0.1128./vec_theta;
% 
% parfor i = 1:length(alpha_2)
%     
%     disp(num2str(i/length(alpha_2)*100));
%         
%     batch_complete = batch_with_theta(batch, alpha_2(i), vec_theta(i));
%     
%     vec_1(i) = -log_LH_evaluation_MODEL1(batch_complete, alpha_1(i), vec_theta(i), dt);
%     vec_2(i) = -log_LH_evaluation_MODEL1(batch_complete, alpha_2(i), vec_theta(i), dt);
%     vec_3(i) = -log_LH_evaluation_MODEL1(batch_complete, alpha_3(i), vec_theta(i), dt);    
%     vec_4(i) = -log_LH_evaluation_MODEL1(batch_complete, alpha_4(i), vec_theta(i), dt);    
%     vec_5(i) = -log_LH_evaluation_MODEL1(batch_complete, alpha_5(i), vec_theta(i), dt);    
%     
% end
% 
% figure;
% plot(vec_theta,vec_1,'LineWidth',2);
% hold on;
% plot(vec_theta,vec_2,'LineWidth',2);
% plot(vec_theta,vec_3,'LineWidth',2);
% plot(vec_theta,vec_4,'LineWidth',2);
% plot(vec_theta,vec_5,'LineWidth',2);
% title(['Negative Log-Likelihoog for ',num2str(num_days),' days']);
% grid minor;
% xlabel('$\theta_0$','interpreter','latex');
% ylabel('Negative Log-Likelihoog evaluation');
% legend('Curve $\theta_0\alpha=0.0728$','Curve $\theta_0\alpha=0.0828$',...
%     'Curve $\theta_0\alpha=0.0928$','Curve $\theta_0\alpha=0.1028$',...
%     'Curve $\theta_0\alpha=0.1128$','Interpreter','Latex');
% xlim([min(vec_theta) max(vec_theta)]);
% saveas(gcf,[pwd '/Results/likelihood/model1/Log-Likelihood_Evaluations'],'epsc');