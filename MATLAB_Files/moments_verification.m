close all;
clear all;
clc;

epsilon                = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon);

t  = Table_Testing_Complete.Time(1,:);
dt = t(2);

theta_0 = 1.2656;
alpha   = 0.07052;
n       = 100;
p_end   = 1/2;
p_lims  = 0.03;

p = p_end-p_lims:0.001:p_end+p_lims;
v = -0.3:0.001:0.3;

for i = 1:length(p)
    for j = 1:length(v)
                
        P_dot        = (p_end-p(i)) / dt;
        Theta_t      = theta_t(theta_0, alpha, p(i), P_dot);
        Theta_t_next = theta_t(theta_0, alpha, p_end, 0);
        
        Exp_V   = moment_1(v(j),theta_0,alpha,Theta_t,Theta_t_next,p(i),p_end,dt,n);
        Exp_V_2 = moment_2(v(j),theta_0,Theta_t,Theta_t_next,p(i),p_end,alpha,Exp_V,dt,n);
                
        for k = 1:100
            sim_path_v(k)   = sde_FE_error(v(j),alpha,theta_0,Theta_t,dt,p(i));
            sim_path_v_2(k) = sde_FE_error_2(v(j)^2,alpha,theta_0,Theta_t,dt,p(i));
        end
        
        Exp_V_emp   = mean(sim_path_v);
        Exp_V_emp_2 = mean(sim_path_v_2);

        exp_v_ode(i,j)   = Exp_V(end);
        exp_v_2_ode(i,j) = Exp_V_2(end);
        exp_v_emp(i,j)   = Exp_V_emp;
        exp_v_emp_2(i,j) = Exp_V_emp_2;
        
    end
end

error_v   = abs(exp_v_ode-exp_v_emp);
error_v_2 = abs(exp_v_2_ode-exp_v_emp_2);

figure;
[V,Y] = meshgrid(p-p_end,v);
surf(V,Y,exp_v_emp');
title('$E[V]$ Empirical','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/classic/1'],'epsc');

figure;
surf(V,Y,exp_v_ode');
title('$E[V]$ ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/classic/2'],'epsc');

figure;
contourf(V,Y,error_v'); colorbar;
title('$E[V]$ Abs-error Empirical-ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/classic/3'],'epsc');

figure;
surf(V,Y,exp_v_emp_2');
title('$E[V^2]$ Empirical','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/classic/4'],'epsc');

figure;
surf(V,Y,exp_v_2_ode');
title('$E[V^2]$ ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/classic/5'],'epsc');

figure;
contourf(V,Y,error_v_2',10); colorbar;
title('$E[V^2]$ Abs-error Empirical-ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/classic/6'],'epsc');

%% Now we plot w.r.t. p_end:

error_v = []; error_v_2 = [];
p_end = 0.2:0.1:0.8;

for l = 1:length(p_end)
    
    p = p_end(l)-p_lims:0.001:p_end(l)+p_lims;
    disp(['Completed: ',num2str(l/length(p_end)*100),'%']);
    num_trunc = 0;

    for i = 1:length(p)
        for j = 1:length(v)

            if v(j) < -p(i)
                v_toUse   = -p(i);
                num_trunc = num_trunc + 1;
                flag      = 0;
            elseif v(j) > 1-p(i)
                v_toUse   = 1-p(i);
                num_trunc = num_trunc + 1;
                flag      = 0;
            else
                v_toUse = v(j);
                flag    = 1;
            end
            
            P_dot        = (p_end(l)-p(i)) / dt;
            Theta_t      = theta_t(theta_0, alpha, p(i), P_dot);
            Theta_t_next = theta_t(theta_0, alpha, p_end(l), 0);

            Exp_V   = moment_1(v_toUse,theta_0,alpha,Theta_t,Theta_t_next,p(i),p_end(l),dt,n);
            Exp_V_2 = moment_2(v_toUse,theta_0,Theta_t,Theta_t_next,p(i),p_end(l),alpha,Exp_V,dt,n);

            for k = 1:100
                sim_path_v(k)   = sde_FE_error(v_toUse,alpha,theta_0,Theta_t,dt,p(i));
                sim_path_v_2(k) = sde_FE_error_2(v_toUse^2,alpha,theta_0,Theta_t,dt,p(i));
            end

            Exp_V_emp   = mean(sim_path_v);
            Exp_V_emp_2 = mean(sim_path_v_2);

            exp_v_ode(i,j)   = Exp_V(end)*flag;
            exp_v_2_ode(i,j) = Exp_V_2(end)*flag;
            exp_v_emp(i,j)   = Exp_V_emp*flag;
            exp_v_emp_2(i,j) = Exp_V_emp_2*flag;

        end
    end

    error_v(l)   = sum(sum(abs(exp_v_ode-exp_v_emp))) / (length(p)*length(v) - num_trunc);
    error_v_2(l) = sum(sum(abs(exp_v_2_ode-exp_v_emp_2))) / (length(p)*length(v) - num_trunc);
    
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
saveas(gcf,[pwd '/Results/moments/classic/7'],'epsc');