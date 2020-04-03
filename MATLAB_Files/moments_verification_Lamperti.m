close all;
clear all;
clc;

epsilon                = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon);

t  = Table_Testing_Complete.Time(1,:);
dt = t(2);

theta_0 = 1.2656;
alpha   = 0.07052;
n       = 500;
p_end   = 1/2;
p_lims  = 0.03;

p = p_end-p_lims:0.001:p_end+p_lims;
v = -0.3:0.01:0.3;

for i = 1:length(p)
    
    for j = 1:length(v)
        
        disp([num2str(i/length(p)),' and ',num2str(j/length(v))]);
                
        z = lamperti_transform(theta_0,alpha,v(j),p(i));
        
        P_dot        = (p_end-p(i)) / dt;
        Theta_t      = theta_t(theta_0, alpha, p(i), P_dot);
        
        Exp_Z   = moment_1_L(z,theta_0,alpha,p(i),p_end,dt,n);
        Exp_Z_2 = moment_2_L(Exp_Z,theta_0,alpha,p(i),p_end,dt,n) + Exp_Z(end)^2;
                
        parfor k = 1:100
            sim_path_z(k)   = sde_Lamperti_FE(z,alpha,theta_0,Theta_t,dt,p(i),P_dot);
            sim_path_z_2(k) = sde_Lamperti_FE_2(z^2,alpha,theta_0,Theta_t,dt,p(i),P_dot);
        end
        
        Exp_Z_emp   = mean(sim_path_z);
        Exp_Z_emp_2 = mean(sim_path_z_2);

        exp_z_ode(i,j)   = Exp_Z(end);
        exp_z_2_ode(i,j) = Exp_Z_2(end);
        exp_z_emp(i,j)   = Exp_Z_emp;
        exp_z_emp_2(i,j) = Exp_Z_emp_2;
        
    end
end

error_z   = abs(exp_z_ode-exp_z_emp);
error_z_2 = abs(exp_z_2_ode-exp_z_emp_2);

figure;
[V,Y] = meshgrid(p-p_end,v);
surf(V,Y,exp_z_emp');
title('$E[Z]$ Empirical','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/lamperti/1'],'epsc');

figure;
surf(V,Y,exp_z_ode');
title('$E[Z]$ ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/lamperti/2'],'epsc');

figure;
contourf(V,Y,error_z'); colorbar;
title('$E[Z]$ Abs-error Empirical-ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/lamperti/3'],'epsc');

figure;
surf(V,Y,exp_z_emp_2');
title('$E[Z^2]$ Empirical','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/lamperti/4'],'epsc');

figure;
surf(V,Y,exp_z_2_ode');
title('$E[Z^2]$ ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/lamperti/5'],'epsc');

figure;
contourf(V,Y,error_z_2',[0:0.1:1]); colorbar;
title('$E[Z^2]$ Abs-error Empirical-ODE','interpreter','latex');
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
xlim([min(p-p_end) max(p-p_end)]); ylim([min(v) max(v)]);
saveas(gcf,[pwd '/Results/moments/lamperti/6'],'epsc');

%% Now we plot w.r.t. p_end:

error_z = []; error_z_2 = [];
p_end = 0.2:0.1:0.8;

Z_inf = lamperti_transform(theta_0,alpha,0,0) ;
Z_sup = lamperti_transform(theta_0,alpha,1,0) ;
Z_100 = (Z_sup-Z_inf)/100;

for l = 1:length(p_end)
    
    p = p_end(l)-p_lims:0.001:p_end(l)+p_lims;
    disp(['Completed: ',num2str(l/length(p_end)*100),'%']);
    num_trunc = 0;

    for i = 1:length(p)
        disp(num2str(i/length(p)));
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
            
            z = lamperti_transform(theta_0,alpha,v_toUse,p(i));
            
            if z >= Z_sup - Z_100 * 5
                z = Z_sup - Z_100 * 5;
                if flag == 1
                    flag = 0;
                    num_trunc = num_trunc + 1;
                end
            elseif z <= Z_inf + Z_100 * 5
                z = Z_inf + Z_100 * 5;
                if flag == 1
                    flag = 0;
                    num_trunc = num_trunc + 1;
                end
            end
        
            P_dot        = (p_end(l)-p(i)) / dt;
            Theta_t      = theta_t(theta_0, alpha, p(i), P_dot);

            Exp_Z   = moment_1_L(z,theta_0,alpha,p(i),p_end(l),dt,n);
            Var_Z   = moment_2_L(Exp_Z,theta_0,alpha,p(i),p_end(l),dt,n);
            Exp_Z_2 = Var_Z + Exp_Z(end)^2;

            for k = 1:100
                sim_path_z(k)   = sde_Lamperti_FE(z,alpha,theta_0,Theta_t,dt,p(i),P_dot);
                sim_path_z_2(k) = sde_Lamperti_FE_2(z^2,alpha,theta_0,Theta_t,dt,p(i),P_dot);
            end

            Exp_Z_emp   = mean(sim_path_z);
            Exp_Z_emp_2 = mean(sim_path_z_2);

            z_plot(i,j)      = z;
            exp_z_ode(i,j)   = Exp_Z(end)*flag;
            exp_z_2_ode(i,j) = Exp_Z_2(end)*flag;
            var_z_ode(i,j)   = Var_Z(end)*flag;
            exp_z_emp(i,j)   = Exp_Z_emp*flag;
            exp_z_emp_2(i,j) = Exp_Z_emp_2*flag;
            flat_plot(i,j)   = flag;

        end
    end
    
    dif_exp_z   = abs(exp_z_ode-exp_z_emp);
    dif_exp_z_2 = abs(exp_z_2_ode-exp_z_emp_2);
    
    figure;
    surf(V,Y,(exp_z_ode)'); 
    title(['First moment ODE P(t+1)=',num2str(p_end(l))]);
    saveas(gcf,[pwd '/Results/moments/lamperti/errors/fm_ODE_',num2str(l)],'epsc');
    figure;
    surf(V,Y,(exp_z_emp)'); 
    title(['First moment MC P(t+1)=',num2str(p_end(l))]);
    saveas(gcf,[pwd '/Results/moments/lamperti/errors/fm_MC_',num2str(l)],'epsc');
    figure;
    surf(V,Y,(exp_z_2_ode)');
    title(['Second moment ODE P(t+1)=',num2str(p_end(l))]);
    saveas(gcf,[pwd '/Results/moments/lamperti/errors/sm_ODE_',num2str(l)],'epsc');
    figure;
    surf(V,Y,(var_z_ode)');
    title(['Variance ODE P(t+1)=',num2str(p_end(l))]);
    saveas(gcf,[pwd '/Results/moments/lamperti/errors/var_ODE_',num2str(l)],'epsc');
    figure;
    surf(V,Y,(exp_z_emp_2)'); 
    title(['Second moment MC P(t+1)=',num2str(p_end(l))]);
    saveas(gcf,[pwd '/Results/moments/lamperti/errors/sm_MC_',num2str(l)],'epsc');
    figure;
    contourf(V,Y,(dif_exp_z)',[0:0.1:1]); colorbar;
    title(['Error first moment when P(t+1)=',num2str(p_end(l))]);
    saveas(gcf,[pwd '/Results/moments/lamperti/errors/fm_',num2str(l)],'epsc');
    figure;
    contourf(V,Y,(dif_exp_z_2)',[0:0.1:1]); colorbar;
    title(['Error second moment when P(t+1)=',num2str(p_end(l))]);
    saveas(gcf,[pwd '/Results/moments/lamperti/errors/sm_',num2str(l)],'epsc');
    
    figure;
    surf(V,Y,(flat_plot)'); 
    title(['Flag for P(t+1)=',num2str(p_end(l))]);
    figure;
    surf(V,Y,(z_plot)'); 
    title(['Z(v,p) for P(t+1)=',num2str(p_end(l))]);
    
    error_z(l)   = sum(sum(abs(dif_exp_z))) / (length(p)*length(v) - num_trunc);
    error_z_2(l) = sum(sum(abs(dif_exp_z_2))) / (length(p)*length(v) - num_trunc);
    
    pause(0.5);
    
end

figure;
plot(p_end,error_z);
hold on;
plot(p_end,error_z_2);
grid minor;
xlabel('Next Forecast Value');
ylabel('Error');
title('Lamperti moments mean absolute error');
legend('Error first moment (empirical-ODE)','Error second moment (empirical-ODE)');
xlim([min(p_end) max(p_end)]);
saveas(gcf,[pwd '/Results/moments/lamperti/7'],'epsc');