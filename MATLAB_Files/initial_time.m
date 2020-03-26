close all;
clear all;
clc;

epsilon               = 0.018;
eta                   = 0.1;
we_remove_bad_samples = 1;
Ta_Tra_Comp           = load_data_eps(epsilon);
save_plots            = 0;

Time     = Ta_Tra_Comp.Time;
Forecast = Ta_Tra_Comp.Forecast;

dt         = Time(1,2);
[~, N_ini] = size(Forecast);
N          = N_ini - 1; % We have N_ini measurements but N samples.

num_days = 127;

if we_remove_bad_samples == 0
    [~, batch] = new_batch_fixed(Ta_Tra_Comp,num_days,N);
elseif we_remove_bad_samples == 1
    [~, batch] = new_batch_fixed_removed_samples(Ta_Tra_Comp,eta,num_days,N);
end

theta_0 = 1.6290;
alpha   = 0.06;

batch_complete = batch_with_theta(batch, alpha, theta_0);

norm  = 1000;
times = 1;

for delta = 1:1:norm*times

    disp(num2str(delta/(norm*times)));
    
    [val(delta),m2(delta),xi1(delta),sig2(delta)] = ...
        first_log_LH_evaluation(batch_complete, theta_0, alpha, delta/norm, dt, N, norm, times);

end

title('Initial error Histogram and 1-element Likelihoods');
grid minor;
if save_plots == 1
saveas(gcf,[pwd ['/Results/delta/allLH_',num2str(we_remove_bad_samples)]],'epsc');
end
% The likelihoods plotted corresponds to the final value of $\delta$.
% This value is 'times' in days.

figure;
plot([1:norm*times]/norm,val);
xlabel('$\delta$','interpreter','latex');
title('Negative log-Likelihood'); grid minor;
if save_plots == 1
saveas(gcf,[pwd ['/Results/delta/LL_',num2str(we_remove_bad_samples)]],'epsc');
end

for i = 1:length(val)-1
    fd(i) = (val(i+1)-val(i))/(1/norm);
end

figure;
plot([1:norm*times-1]/norm,fd);
xlabel('$\delta$','interpreter','latex');
title('F-D for Negative log-Likelihood'); grid minor;
if save_plots == 1
saveas(gcf,[pwd ['/Results/delta/dif_LL_',num2str(we_remove_bad_samples)]],'epsc');
end

figure;
plot([1:norm*times]/norm,xi1);
xlabel('$\delta$','interpreter','latex');
title('Mean $\xi_1$','interpreter','latex'); grid minor;
if save_plots == 1
saveas(gcf,[pwd ['/Results/delta/Xi_',num2str(we_remove_bad_samples)]],'epsc');
end

figure;
plot([1:norm*times]/norm,m2);
xlabel('$\delta$','interpreter','latex');
title('Mean $m_2(t_0)$','interpreter','latex'); grid minor;
if save_plots == 1
saveas(gcf,[pwd ['/Results/delta/M2_',num2str(we_remove_bad_samples)]],'epsc');
end

figure;
plot([1:norm*times]/norm,sig2);
xlabel('$\delta$','interpreter','latex');
title('Mean $\sigma_2$','interpreter','latex'); grid minor;
if save_plots == 1
saveas(gcf,[pwd ['/Results/delta/Sig_',num2str(we_remove_bad_samples)]],'epsc');
end

%% Delta w.r.t. Eta:

close all;
clear all;
clc;

initial     = 0;
save_plots  = 1;

if initial == 1
    theta_0     = 1.6290;
    alpha       = 0.06;
else
    theta_0     = 3.912;
    alpha       = 0.019308;
end

epsilon     = 0.018;
eta         = 0.005:0.005:0.2;
num_days    = 127;
Ta_Tra_Comp = load_data_eps(epsilon);
Time        = Ta_Tra_Comp.Time;
Forecast    = Ta_Tra_Comp.Forecast;
dt          = Time(1,2);
[M, N_ini]  = size(Forecast);
N           = N_ini - 1; % We have N_ini measurements but N samples.
norm        = 1000;
times       = 1;

for j = 1:length(eta)

    [~, batch]     = new_batch_fixed_removed_samples(Ta_Tra_Comp,eta(j),num_days,N);
    batch_complete = batch_with_theta(batch, alpha, theta_0);

    for delta = 1:1:norm*times
        disp(num2str(delta/(norm*times)));

        [val(delta),m2(delta),xi1(delta),sig2(delta)] = ...
            first_log_LH_evaluation(batch_complete, theta_0, alpha, delta/norm, dt, N, norm, times);
    end

    [min_val,min_index(j)] = min(val);

end

figure;
plot(eta,min_index/norm);
xlabel('$\eta$','interpreter','latex');
ylabel('$\delta$','interpreter','latex');
title('$\delta$ over $\eta$','interpreter','latex'); grid minor;
if save_plots == 1
    if initial == 1
        saveas(gcf,[pwd '/Results/delta/eta_ini'],'epsc');
    else
        saveas(gcf,[pwd '/Results/delta/eta_opt'],'epsc');
    end
end

figure;
plot(eta,min_index/norm*24*60);
xlabel('$\eta$','interpreter','latex');
ylabel('$\delta$ in minutes','interpreter','latex');
title('$\delta$ over $\eta$','interpreter','latex'); grid minor;
if save_plots == 1
    if initial == 1
        saveas(gcf,[pwd '/Results/delta/eta_min_ini'],'epsc');
    else
        saveas(gcf,[pwd '/Results/delta/eta_min_opt'],'epsc');
    end
end