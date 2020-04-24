close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet                = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon                = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet);

t  = Table_Testing_Complete.Time(1,:);
dt = t(2);

theta_0 = 1.2656;
alpha   = 0.07052;
n       = 100;
p_end   = 1/2;
p_lims  = 0.03;

p = p_end-p_lims:0.001:p_end+p_lims;
x = -0.4:0.001:0.2;
% Here we are being confused. The vector x is the error vector. The idea is
% that for each p, we compute the measurement as p+x, which would imply
% that this vector x is the error vector v, but we called it x.

for i = 1:length(p)
    for j = 1:length(x)
                        
%         moment_1_M1(x,theta_0,p1,p2,dt,n)
        Exp_X   = moment_1_M1(p(i)+x(j),theta_0,p(i),p_end,dt,n);
%         moment_2_M1(x,m1,theta_0,alpha,p1,p2,dt,n)
        Exp_X_2 = moment_2_M1(p(i)+x(j),Exp_X,theta_0,alpha,p(i),p_end,dt,n);
                
        for k = 1:100
%             sde_FE_MODEL1(sim_path_i,alpha,theta_0,dt,P)
            sim_path_x(k)   = sde_FE_MODEL1(p(i)+x(j),alpha,theta_0,dt,p(i));
%             sde_FE_MODEL1_2(sim_path_i,alpha,theta_0,dt,P) 
            sim_path_x_2(k) = sde_FE_MODEL1_2((p(i)+x(j))^2,alpha,theta_0,dt,p(i));
        end
        
        Exp_X_emp   = mean(sim_path_x);
        Exp_X_emp_2 = mean(sim_path_x_2);

        exp_x_ode(i,j)   = Exp_X(end);
        exp_x_2_ode(i,j) = Exp_X_2(end);
        exp_x_emp(i,j)   = Exp_X_emp;
        exp_x_emp_2(i,j) = Exp_X_emp_2;
        
    end
end

error_v   = abs(exp_x_ode-exp_x_emp);
error_v_2 = abs(exp_x_2_ode-exp_x_emp_2);

figure;
[V,Y] = meshgrid(p-p_end,x);
surf(V,Y,exp_x_emp');
title('$E[V]$ Empirical','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(x) max(x)]);
saveas(gcf,[pwd '/Results/moments/model1/1'],'epsc');

figure;
surf(V,Y,exp_x_ode');
title('$E[V]$ ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(x) max(x)]);
saveas(gcf,[pwd '/Results/moments/model1/2'],'epsc');

figure;
contourf(V,Y,error_v'); colorbar;
title('$E[V]$ Abs-error Empirical-ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(x) max(x)]);
saveas(gcf,[pwd '/Results/moments/model1/3'],'epsc');

figure;
surf(V,Y,exp_x_emp_2');
title('$E[V^2]$ Empirical','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(x) max(x)]);
saveas(gcf,[pwd '/Results/moments/model1/4'],'epsc');

figure;
surf(V,Y,exp_x_2_ode');
title('$E[V^2]$ ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(x) max(x)]);
saveas(gcf,[pwd '/Results/moments/model1/5'],'epsc');

figure;
contourf(V,Y,error_v_2',10); colorbar;
title('$E[V^2]$ Abs-error Empirical-ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(x) max(x)]);
saveas(gcf,[pwd '/Results/moments/model1/6'],'epsc');

%% Now we plot w.r.t. p_end: STILL NOT DONE

error_v = []; error_v_2 = [];
p_end = 0.2:0.1:0.8;

for l = 1:length(p_end)
    
    p = p_end(l)-p_lims:0.001:p_end(l)+p_lims;
    disp(['Completed: ',num2str(l/length(p_end)*100),'%']);
    num_trunc = 0;

    for i = 1:length(p)
        for j = 1:length(x)
            
            x_to_use = p(i)+x(j);
            
            if 0 > x_to_use
                x_to_use   = 0;
                num_trunc = num_trunc + 1;
                flag      = 0;
            elseif x_to_use > 1
                x_to_use   = 1;
                num_trunc = num_trunc + 1;
                flag      = 0;
            else
                flag    = 1;
            end
            
            Exp_X   = moment_1_M1(x_to_use,theta_0,p(i),p_end(l),dt,n);
            Exp_X_2 = moment_2_M1(x_to_use,Exp_X,theta_0,alpha,p(i),p_end(l),dt,n);

            for k = 1:100
                sim_path_x(k)   = sde_FE_MODEL1(x_to_use,alpha,theta_0,dt,p(i));
                sim_path_x_2(k) = sde_FE_MODEL1_2(x_to_use^2,alpha,theta_0,dt,p(i));
            end

            Exp_X_emp   = mean(sim_path_x);
            Exp_X_emp_2 = mean(sim_path_x_2);

            exp_x_ode(i,j)   = Exp_X(end)*flag;
            exp_x_2_ode(i,j) = Exp_X_2(end)*flag;
            exp_x_emp(i,j)   = Exp_X_emp*flag;
            exp_x_emp_2(i,j) = Exp_X_emp_2*flag;

        end
    end

    error_v(l)   = sum(sum(abs(exp_x_ode-exp_x_emp))) / (length(p)*length(x) - num_trunc);
    error_v_2(l) = sum(sum(abs(exp_x_2_ode-exp_x_emp_2))) / (length(p)*length(x) - num_trunc);
    
end

figure;
plot(p_end,error_v);
hold on;
plot(p_end,error_v_2);
grid minor;
xlabel('Next Forecast Value');
ylabel('Error');
title('Moments mean absolute error');
legend('Error first moment (empirical-ODE)','Error second moment (empirical-ODE)');
xlim([min(p_end) max(p_end)]);
saveas(gcf,[pwd '/Results/moments/model1/7'],'epsc');