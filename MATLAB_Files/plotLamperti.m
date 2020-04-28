close all;
clear all;
clc;

likelihood    = 'lamperti';
% dataSet can be AWSTP, MTLOG or UTEP5.
dataSet = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon       = 0.018;
% dataKind can be classic or comparable.
dataKind = 'comparable';
[Ta_Tra_Comp] = load_data_eps(epsilon,dataSet,dataKind);

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

dt         = Time(1,2);
[M, N_ini] = size(Forecast);
N          = N_ini - 1; % We have N_ini measurements but N samples.

num_days                = height(Ta_Tra_Comp); % Maximum 127 for MTLOG.
[Table_Training, batch] = new_batch_fixed(Table_Training,num_days,N);

%% Log-Likelihood plot:

val_theta_0 = [];
val_alpha   = [];
val_fval    = [];
val_fval_2  = [];

vec_theta_0 = 0.5:0.005:3;
vec_alpha   = 0.01:0.001:0.1;
num_divs    = 20;
new_alpha   = 0.0965./vec_theta_0;

for i = 1:length(vec_theta_0)

    aux_theta_0 = ones(1,length(vec_alpha));
    aux_alpha   = ones(1,length(vec_alpha));
    aux_fval    = ones(1,length(vec_alpha));
    aux_fval_2  = ones(1,length(vec_alpha));
    theta_0     = vec_theta_0(i);
    
    parfor j = 1:length(vec_alpha)
        
        alpha          = vec_alpha(j);
        x0             = [theta_0,alpha];
        [x1,fval]      = likelihood_nested_optimization_L(batch, x0, dt);
        aux_theta_0(j) = x1(1);
        aux_alpha(j)   = x1(2);
        aux_fval(j)    = fval;
        
        batch_L        = batch_with_theta_L(batch, alpha, theta_0);
        aux_fval_2(j)  = -likelihood_optimization_L(batch_L, theta_0, alpha, dt);
        
    end

    val_theta_0(i,:) = aux_theta_0(:);
    val_alpha(i,:)   = aux_alpha(:);
    val_fval(i,:)    = aux_fval(:);
    val_fval_2(i,:)  = aux_fval_2(:);

    disp(num2str(i/length(vec_theta_0)*100));

end

% load([pwd '/Results/likelihood/',likelihood,'/val_theta_0.mat'],'val_theta_0');
% load([pwd '/Results/likelihood/',likelihood,'/val_alpha.mat'],'val_alpha');
% save([pwd '/Results/likelihood/',likelihood,'/val_theta_0.mat'],'val_theta_0');
% save([pwd '/Results/likelihood/',likelihood,'/val_alpha.mat'],'val_alpha');
% save([pwd '/Results/likelihood/',likelihood,'/val_fval.mat'],'val_fval');
% save([pwd '/Results/likelihood/',likelihood,'/val_fval_2.mat'],'val_fval_2');

for i = 1:length(vec_theta_0)
    for j = 1:length(vec_alpha)
        rel_change_theta_0(i,j) = abs(vec_theta_0(i)-val_theta_0(i,j)) / abs(vec_theta_0(i));
        rel_change_alpha(i,j)   = abs(vec_alpha(j)-val_alpha(i,j)) / abs(vec_alpha(j));
    end
end

total_change = rel_change_theta_0 + rel_change_alpha;

figure;
[X,Y] = meshgrid(vec_theta_0,vec_alpha);
hold on;
contourf(X,Y,total_change',[0:0.05:1]); colorbar;
% colormap summer;
plot(vec_theta_0,new_alpha,'LineWidth',2,'Color','r');
ylim([min(vec_alpha) max(vec_alpha)]);
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Fixed-Point Lamperti ',num2str(num_days),' days']);
legend('Level sets','Curve $\theta_0\alpha=0.097$','Interpreter','Latex');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/',likelihood,'/Log-Likelihood'],'epsc');

figure;
[X,Y] = meshgrid(vec_theta_0,vec_alpha);
hold on;
contourf(X,Y,val_fval'); colorbar;
% colormap summer;
% plot(vec_theta_0,new_alpha,'LineWidth',2,'Color','r');
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Value in optimal point ',num2str(num_days),' days']);
% legend('Level sets','Curve $\theta_0\alpha=0.097$','Interpreter','Latex');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/',likelihood,'/Log-Likelihood_f'],'epsc');

figure;
[X,Y] = meshgrid(vec_theta_0,vec_alpha);
hold on;
contourf(X,Y,val_fval_2'); colorbar;
% colormap summer;
% plot(vec_theta_0,new_alpha,'LineWidth',2,'Color','r');
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Value in point ',num2str(num_days),' days']);
% legend('Level sets','Curve $\theta_0\alpha=0.097$','Interpreter','Latex');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/',likelihood,'/Log-Likelihood_f_2'],'epsc');

figure;
[X,Y] = meshgrid(vec_theta_0,vec_alpha);
hold on;
contourf(X,Y,val_fval_2'-val_fval'); colorbar;
% colormap summer;
% plot(vec_theta_0,new_alpha,'LineWidth',2,'Color','r');
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Difference ',num2str(num_days),' days']);
% legend('Level sets','Curve $\theta_0\alpha=0.097$','Interpreter','Latex');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/',likelihood,'/Log-Likelihood_f_extra'],'epsc');

%% Plot with zoom:

% figure;
% hold on;
% [X,Y] = meshgrid(vec_theta_0(1:61),vec_alpha);
% total_change = total_change(1:61,:);
% contourf(X,Y,total_change',num_divs, 'LineStyle','none'); colorbar;
% xlabel('$\theta_0$','interpreter','latex');
% ylabel('$\alpha$','interpreter','latex');
% title(['Fixed-Point Lamperti ',num2str(num_days),' days']);
% legend('Level sets');
% pause(0.1);
% saveas(gcf,[pwd '/Results/likelihood/',likelihood,'/Log-Likelihood_Zoom'],'epsc');