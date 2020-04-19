close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet = 'UTEP5';
[Ta_Tra_Comp, Ta_Test_Compl, Ta_Comp] = load_data(dataSet);

figure(1); figure(2); 

Time              = Ta_Tra_Comp.Time;
Forecast          = Ta_Tra_Comp.Forecast;
Error             = Ta_Tra_Comp.Error;
Real_ADME         = Ta_Tra_Comp.Real_ADME;

dt             = Time(1,2);
[M, N_ini]     = size(Forecast);
N              = N_ini - 1; % We have N_ini measurements but N samples.
addLimitation  = 1;
maxChange      = 1.01;
maxIter        = 100;

est               = initial_guess(Real_ADME, M, N, dt);
samples_minus_eps = create_samples_minus_eps(Forecast, Error, 0.2);
[theta_0,~]       = mean_regression_eps(samples_minus_eps, dt);
max_plot_theta    = 0;

for n = 1:500

    epsilon_new = 0.01 + rand(1)*0.09;
    step        = 1;
    plot_eps    = epsilon_new;
    plot_theta  = [];
    counter     = 0;

    while (step < maxIter)

        epsilon     = epsilon_new;
        samples     = create_samples_eps(Forecast, Real_ADME, epsilon);
        theta_t     = mean_regression_eps(samples, dt);
        if theta_t > theta_0
            epsilon_new = est / theta_t;
            if addLimitation
                if epsilon_new > epsilon*maxChange
                    epsilon_new = epsilon*maxChange;
                elseif epsilon_new < epsilon/maxChange
                    epsilon_new = epsilon/maxChange;
                end
            end
              disp(['========== Step: ',num2str(step),' ==========']);
            disp(['Epsilon = ',num2str(epsilon_new)]);
            disp(['Theta(t) = ',num2str(theta_t)]);
            step = step + 1;
            plot_eps(end+1)   = epsilon_new;
            plot_theta(end+1) = theta_t;
            counter           = counter + 1;
        else
        	epsilon_new = epsilon_new * 0.999;
        end

    end

    max_plot_theta = max(max(plot_theta(1:end-1)),max_plot_theta);
    
    Forecast_vec = Forecast(:);
    disp('==============================');
    disp(['Minimum forecast: ',num2str(min(Forecast_vec))]);
    disp(['Maximum forecast: ',num2str(max(Forecast_vec))]);

    set(0,'CurrentFigure',1);
    plot([0:step-2],plot_eps(1:end-1)); 
    hold on;
    title('$\epsilon$ over iterations','Interpreter','Latex');
    xlabel('Number iteration');
    xlim([0 step-2]);
    ylim([0 0.1]);

    set(0,'CurrentFigure',2);
    plot(plot_theta(1:end-1)); 
    hold on;
    plot(theta_0*ones(1,length(plot_theta(1:end-1))));
    title('k over iterations','Interpreter','Latex');
    xlabel('Number iteration');
    xlim([1 step-2]);
    ylim([theta_0*0.9 max_plot_theta*1.1]);
    
    pause(0.1)

end

set(0,'CurrentFigure',1);
grid minor;
saveas(gcf,[pwd '/Results/epsilon/epsilon'],'epsc');
set(0,'CurrentFigure',2);
grid minor;
saveas(gcf,[pwd '/Results/epsilon/theta_t'],'epsc');

%% Plot QMM formula as a function of delta:

epsilon           = [0.01:0.00001:0.05];
samples_minus_eps = create_samples_minus_eps(Forecast, Error, 0.1);
therta_sweep      = [0:0.01:10];

for i = 1:length(epsilon)

    samples                   = create_samples_eps(Forecast, Real_ADME, epsilon(i));
    theta_t(i)                = mean_regression_eps(samples, dt);
    [min_value_fun_real(i),~] = value_function_QMM_eps(samples, theta_t(i), dt);

    for j = 1:length(therta_sweep)
        min_val(j) = value_function_QMM_eps(samples, therta_sweep(j), dt);
    end
    min_value_fun_hard(i) = min(min_val);
    
    if theta_0 > theta_t(i)
        theta_t(i) = theta_0;
    end
    [min_value_fun(i),accum(i)] = value_function_QMM_eps(samples, theta_t(i), dt);
    
end

figure;
plot(epsilon, min_value_fun);
hold on;
plot(epsilon, min_value_fun_real);
plot(epsilon, min_value_fun_hard);
legend('Using fixed k', 'Using analytic minimum', 'Using samples minimum', 'interpreter', 'latex');
grid minor;
xlim([min(epsilon) max(epsilon)]);
title('LSM over $\epsilon$','interpreter','latex');
xlabel('$\epsilon$','interpreter','latex');
saveas(gcf,[pwd '/Results/epsilon/LSM'],'epsc');

figure;
plot(epsilon, theta_t); grid minor;
xlim([min(epsilon) max(epsilon)]);
title('k over $\epsilon$','interpreter','latex');
xlabel('$\epsilon$','interpreter','latex');
saveas(gcf,[pwd '/Results/epsilon/theta_over_eps'],'epsc');

figure;
plot(epsilon, accum); grid minor;
xlim([min(epsilon) max(epsilon)]);
title('Number of samples over $\epsilon$','interpreter','latex');
xlabel('$\epsilon$','interpreter','latex');
saveas(gcf,[pwd '/Results/epsilon/num_over_eps'],'epsc');

%% Plor theta_0:

epsilon = [0.01:0.001:0.49];
accum   = [];

for i = 1:length(epsilon)

    samples               = create_samples_minus_eps(Forecast, Error, epsilon(i));
    [theta_0(i),accum(i)] = mean_regression_eps(samples, dt);
    
end

figure;
plot(epsilon, theta_0);
title('$\theta_0$ approximation over $\gamma$','interpreter','latex');
xlim([min(epsilon) max(epsilon)]);
xlabel('$\gamma$','interpreter','latex');
grid minor;
saveas(gcf,[pwd '/Results/epsilon/theta_0'],'epsc');

figure;
plot(epsilon, accum); grid minor;
xlim([min(epsilon) max(epsilon)]);
title('Number of samples over $\gamma$','interpreter','latex');
xlabel('$\gamma$','interpreter','latex');
saveas(gcf,[pwd '/Results/epsilon/num_over_eps_t0'],'epsc');
