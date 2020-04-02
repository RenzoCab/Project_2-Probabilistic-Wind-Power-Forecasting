close all;
clear all;
clc;

epsilon                = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon);

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
quantil  = 0;
save     = 0;
delta    = 22; % The time is delta*10 minutes.
xlimit   = 1; % If this in 1, the plots start at time 0. Otherwise, at -delta.

if quantil
    numPaths = 5000;
else
    numPaths = 5;
end

theta_0 = 1.6290;
alpha   = 0.06;

d  = Table_Testing_Complete.Date;
p  = Table_Testing_Complete.Forecast;
t  = Table_Testing_Complete.Time(1,:);
x  = Table_Testing_Complete.Real_ADME;
v  = Table_Testing_Complete.Error;
dt = t(2); % This is the default dt. I will choose a smaller one.
% If the dt is not small enough, we always scape from [0 1].
dt = dt/1;

figure(10);

minus_t = linspace(-delta*dt,-dt,delta);
exten_t = [minus_t,t];

for i = 1 : height(Table_Testing_Complete)

    X = x(i,:);
    P = p(i,:);
    % All of them show similar results:
    P_exten_linear = interp1([0 6*dt],[P(1) P(6)],minus_t,'linear','extrap');
    P_exten_pchip  = interp1([0 6*dt],[P(1) P(6)],minus_t,'pchip','extrap');
    P_exten_spline = interp1([0 6*dt],[P(1) P(6)],minus_t,'spline','extrap');
    P_exten_makima = interp1([0 6*dt],[P(1) P(6)],minus_t,'makima','extrap');
    
    P        = [P_exten_linear,P]; % We extend the forecast.
    for j = 1:length(P)
        P(j) = max(P(j),epsilon);
        P(j) = min(P(j),1-epsilon);
    end
    P_dot    = (-P(1:end-1) + P(2:end)) / dt; % This starts at 1 end ends at (end - 1).
    for j = 1:length(P)-1
        Theta_t(j) = theta_t(theta_0, alpha, P(j), P_dot(j));
    end
    
    % We assume that the forecast is perfect one hour before it starts.
    sim_path(1:numPaths,1) = P(1);

    set(0,'CurrentFigure',10); clf(10);
    J = plot(exten_t,P); J.LineWidth = 2;
    hold on;
    J = plot(t,X); J.LineWidth = 2;
    xlim([min(exten_t) max(exten_t)]); ylim([0 1]); grid minor;
    xlabel('Time'), ylabel('Power');
    date = d(i);
    title(date{1});
    
    for k = 1:numPaths
        for j = 1 : length(exten_t)-1
            sim_path(k,j+1) = sde_FE(sim_path(k,j),alpha,Theta_t(j),dt,P(j),P_dot(j));
            if sim_path(k,j+1) > 1
                sim_path(k,j+1) = 1;
            elseif sim_path(k,j+1) < 0
                sim_path(k,j+1) = 0;
            end
        end
    end
    
    if quantil == 1
        
        for k = 1:length(sim_path(1,:))
            quant(k,:) = quantile(sim_path(:,k),[0.01 0.10 0.25 0.75 0.90 0.99]);
            men(k)     = mean(sim_path(:,k));
        end
        legend('Forecast','Real Production');
        set(0,'CurrentFigure',10); clf(10); hold on;
        ciplot(quant(:,1),quant(:,6),exten_t,[153 255 153]/255)
        ciplot(quant(:,2),quant(:,5),exten_t,[153 255 255]/255)
        ciplot(quant(:,3),quant(:,4),exten_t,[153 153 255]/255)
        J = plot(exten_t,P,'k'); J.LineWidth = 2;
        J = plot(t,X,'b');       J.LineWidth = 2;
        legend({'99% confidence interval', '90% confidence interval', '50% confidence interval',...
            'Forecast', 'Real production'});
        if xlimit
            xlim([0 max(exten_t)]);
        else
            xlim([min(exten_t) max(exten_t)]);
        end
        ylim([0 1]); grid minor;
        xlabel('Time'), ylabel('Power'); box;
        date = d(i);
        title(date{1});
        if save == 1
            saveas(gcf,[pwd '/Results/bands_testing_days/initial_guess/',num2str(i)],'epsc');
        end
        
    else
        
        plot(exten_t,sim_path,'r');
        legend('Forecast','Real production','Simulations','HandleVisibility','off');
        if xlimit
            xlim([0 max(exten_t)]);
        else
            xlim([min(exten_t) max(exten_t)]);
        end
        box;
        pause(0.1); box;
        if save == 1
            saveas(gcf,[pwd '/Results/paths_testing_days/initial_guess/',num2str(i)],'epsc');
        end
        
    end
    
    pause(0.1);
    
end