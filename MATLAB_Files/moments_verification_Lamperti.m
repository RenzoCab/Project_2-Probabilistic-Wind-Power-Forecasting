close all;
clear all;
clc;

epsilon                = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon);

t  = Table_Testing_Complete.Time(1,:);
dt = t(2);

theta_0 = 3.912;
alpha   = 0.019308;
p_end   = 1/2;
n       = 100;

lims = 0.1;

p = lims:0.01:1-lims;
x = lims:0.01:1-lims;

for i = 1:length(p)
    for j = 1:length(x)
                
        v            = x(j) - p(i);
        P_dot        = (p_end-p(i)) / dt;
        Theta_t      = theta_t(theta_0, alpha, p(i), P_dot);
        Theta_t_next = theta_t(theta_0, alpha, p_end, 0);
        
        Exp_V   = moment_1(v,Theta_t,Theta_t_next,dt,n);
        Exp_V_2 = moment_2(v,Theta_t,Theta_t_next,p(i),p_end,alpha,Exp_V,dt,n);
        
        for k = 1:100
            sim_path_x(k)   = sde_FE(x(j),alpha,theta_0,Theta_t,dt,p(i),P_dot);
            sim_path_x_2(k) = sim_path_x(k)^2;
            
            sim_path_v(k)   = sde_FE_error(v,alpha,theta_0,Theta_t,dt,p(i));
            sim_path_v_2(k) = sim_path_v(k)^2;
        end
        Exp_X   = mean(sim_path_x);
        Exp_X_2 = mean(sim_path_x_2);
        Exp_V_emp   = mean(sim_path_v);
        Exp_V_emp_2 = mean(sim_path_v_2);
        
        error_x_1(i,j) = (Exp_X -  + p_end);
        error_x_2(i,j) = (Exp_X_2 - Exp_V_2(end) -2*p_end*Exp_X - p_end^2);      
        error_V_1(i,j) = (Exp_V_emp - Exp_V(end));
        error_V_2(i,j) = (Exp_V_emp_2 - Exp_V_2(end));       
        exp_x(i,j)   = Exp_X;
        exp_x_2(i,j) = Exp_X_2;
        exp_v(i,j)   = Exp_V(end);
        exp_v_2(i,j) = Exp_V_2(end);
        exp_v_emp(i,j)   = Exp_V_emp;
        exp_v_emp_2(i,j) = Exp_V_emp_2;
        
    end
end

[X,Y] = meshgrid(p,x);
surf(X,Y,error_x_1');
title('Error in X first moment');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,error_x_2');
title('Error in X second moment');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

[X,Y] = meshgrid(p,x);
surf(X,Y,error_V_1');
title('Error in V first moment');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,error_V_2');
title('Error in V second moment');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,exp_x');
title('$E[X]$ empirical','interpreter','latex');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,exp_x_2');
title('$E[X^2]$ empirical','interpreter','latex');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,exp_v_emp');
title('$E[V]$ empirical','interpreter','latex');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,exp_v_emp_2');
title('$E[V^2]$ empirical','interpreter','latex');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,exp_v');
title('$E[V]$ ODE','interpreter','latex');
xlabel('Forecast (p)'); ylabel('Real Production (X)');

figure;
[X,Y] = meshgrid(p,x);
surf(X,Y,exp_v_2');
title('$E[V^2]$ ODE','interpreter','latex');
xlabel('Forecast (p)'); ylabel('Real Production (X)');