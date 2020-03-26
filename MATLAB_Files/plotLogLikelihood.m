close all;
clear all;
clc;

epsilon       = 0.018;
[Ta_Tra_Comp] = load_data_eps(epsilon);

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

num_days = 127;
[Table_Training, batch] = new_batch_fixed(Table_Training,num_days,N);

%% Initial parameters:

est           = initial_guess(Real_ADME, M, N, dt);
gamma         = 0.35; 
samples_gamma = create_samples_minus_eps(Forecast, Error, gamma);
theta_ini     = mean_regression_eps(samples_gamma, dt);
alpha_ini     = est/theta_ini;

%% Log-Likelihood plot:

val = [];
ini_theta = 10;
ini_alpha = 10;
len_theta = 500;
len_alpha = 300;
div_theta = 100; 
div_alpha = 300;

for i = ini_theta:1:len_theta

    vec = ones(1,len_alpha);

    parfor j = ini_alpha:len_alpha
        theta          = i/div_theta;
        alpha          = (j-1+ini_alpha)/div_alpha;
        batch_complete = batch_with_theta(batch, alpha, theta);
        vec(j)         = -log_LH_evaluation(batch_complete, alpha, dt);
    end

    val(i+1-ini_theta,:) = vec(ini_alpha:len_alpha);

    disp(num2str(i/len_theta*100));

end

figure;
hold on;
vec_theta = [ini_theta:len_theta] / div_theta;
vec_alpha = [ini_alpha:len_alpha] / div_alpha;
[X,Y] = meshgrid(vec_theta,vec_alpha);
contourf(X,Y,val',20); colorbar;
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Negative Log-Likelihoog for ',num2str(num_days),' days']);
plot(theta_ini,alpha_ini,'-p','MarkerFaceColor','red','MarkerSize',15);
legend('Level sets','Initial guess');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/Log-Likelihood'],'epsc');

%% Log-Likelihood plot (more refined):

val = [];
ini_theta = 10;
ini_alpha = 10;
len_theta = 300;
len_alpha = 200;
div_theta = 100; 
div_alpha = 1000;

for i = ini_theta:1:len_theta

    vec = ones(1,len_alpha);

    parfor j = ini_alpha:len_alpha
        theta          = i/div_theta;
        alpha          = (j-1+ini_alpha)/div_alpha;
        batch_complete = batch_with_theta(batch, alpha, theta);
        vec(j)         = -log_LH_evaluation(batch_complete, alpha, dt);
    end

    val(i+1-ini_theta,:) = vec(ini_alpha:len_alpha);

    disp(num2str(i/len_theta*100));

end

figure;
hold on;
vec_theta = [ini_theta:len_theta] / div_theta;
vec_alpha = [ini_alpha:len_alpha] / div_alpha;
[X,Y] = meshgrid(vec_theta,vec_alpha);
contourf(X,Y,val',30); colorbar;
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Negative Log-Likelihoog for ',num2str(num_days),' days']);
plot(theta_ini,alpha_ini,'-p','MarkerFaceColor','red','MarkerSize',25);
legend('Level sets','Initial guess');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/Log-Likelihood_refined'],'epsc');

%% Log-Likelihood plot (more more refined):

val = [];
ini_theta = 150;
ini_alpha = 10;
len_theta = 450;
len_alpha = 200;
div_theta = 100; 
div_alpha = 1000;

for i = ini_theta:1:len_theta

    vec = ones(1,len_alpha);

    parfor j = ini_alpha:len_alpha
        theta          = i/div_theta;
        alpha          = (j-1+ini_alpha)/div_alpha;
        batch_complete = batch_with_theta(batch, alpha, theta);
        vec(j)         = -log_LH_evaluation(batch_complete, alpha, dt);
    end

    val(i+1-ini_theta,:) = vec(ini_alpha:len_alpha);

    disp(num2str(i/len_theta*100));

end

figure;
hold on;
vec_theta = [ini_theta:len_theta] / div_theta;
vec_alpha = [ini_alpha:len_alpha] / div_alpha;
[X,Y] = meshgrid(vec_theta,vec_alpha);
contourf(X,Y,val',30); colorbar;
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Negative Log-Likelihoog for ',num2str(num_days),' days']);
legend('Level sets','Initial guess');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/Log-Likelihood_more_refined'],'epsc');
