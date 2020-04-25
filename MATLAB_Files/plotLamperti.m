close all;
clear all;
clc;

likelihood    = 'lamperti';
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

dt         = Time(1,2);
[M, N_ini] = size(Forecast);
N          = N_ini - 1; % We have N_ini measurements but N samples.

num_days                = 127;
[Table_Training, batch] = new_batch_fixed(Table_Training,num_days,N);

%% Log-Likelihood plot:

val_theta_0 = [];
val_alpha   = [];

% vec_theta_0 = 0.5:0.1:3;
% vec_alpha   = 0.01:0.005:0.1;
% num_divs    = 20;

vec_theta_0 = 0.5:0.025:3;
vec_alpha   = 0.01:0.001:0.1;
num_divs    = 20;
new_alpha   = 0.083./vec_theta_0;

for i = 1:length(vec_theta_0)

    aux_theta_0 = ones(1,length(vec_alpha));
    aux_alpha   = ones(1,length(vec_alpha));
    theta_0     = vec_theta_0(i);
    
    parfor j = 1:length(vec_alpha)
        
        alpha          = vec_alpha(j);
        x0             = [theta_0,alpha];
        [x1,fval]      = likelihood_nested_optimization_L(batch, x0, dt);
        aux_theta_0(j) = x1(1);
        aux_alpha(j)   = x1(2);
        
    end

    val_theta_0(i,:) = aux_theta_0(:);
    val_alpha(i,:)   = aux_alpha(:);

    disp(num2str(i/length(vec_theta_0)*100));

end

% load([pwd '/Results/likelihood/',likelihood,'/val_theta_0.mat'],'val_theta_0');
% load([pwd '/Results/likelihood/',likelihood,'/val_alpha.mat'],'val_alpha');
% save([pwd '/Results/likelihood/',likelihood,'/val_theta_0.mat'],'val_theta_0');
% save([pwd '/Results/likelihood/',likelihood,'/val_alpha.mat'],'val_alpha');

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
colormap summer
plot(vec_theta_0,new_alpha,'LineWidth',2,'Color','r');
ylim([min(vec_alpha) max(vec_alpha)]);
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Fixed-Point Lamperti ',num2str(num_days),' days']);
legend('Level sets','Curve $\theta_0\alpha=0.083$','Interpreter','Latex');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/',likelihood,'/Log-Likelihood'],'epsc');

%% Plot with zoom:

figure;
hold on;
[X,Y] = meshgrid(vec_theta_0(1:61),vec_alpha);
total_change = total_change(1:61,:);
contourf(X,Y,total_change',num_divs, 'LineStyle','none'); colorbar;
xlabel('$\theta_0$','interpreter','latex');
ylabel('$\alpha$','interpreter','latex');
title(['Fixed-Point Lamperti ',num2str(num_days),' days']);
legend('Level sets');
pause(0.1);
saveas(gcf,[pwd '/Results/likelihood/',likelihood,'/Log-Likelihood_Zoom'],'epsc');
