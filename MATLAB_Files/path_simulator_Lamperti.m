close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet);
whatToDo               = 'lamperti_clasic_optimal';
% whatToDo               = 'lamperti_clasic_IG';
% whatToDo               = 'lamperti_optimal';

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
quantil  = 1;
save     = 1;
delta    = 21; % The time is delta*10 minutes.
xlimit   = 1; % If this in 1, the plots start at time 0. Otherwise, at -delta.

if quantil
    numPaths = 5000;
else
    numPaths = 5;
end

if  strcmp(whatToDo,'lamperti_clasic_optimal')
    theta_0 = 1.180;
    alpha   = 0.070;
elseif  strcmp(whatToDo,'lamperti_clasic_IG')
    theta_0 = 1.6290;
    alpha   = 0.06;
elseif  strcmp(whatToDo,'lamperti_optimal')
    theta_0 = 2.200;
    alpha   = 0.038;
else
    error('Wrong whatToDo script!');
end

Z_inf = lamperti_transform(theta_0,alpha,0,0,2);
Z_sup = lamperti_transform(theta_0,alpha,1,0,2);
Z_100 = (Z_sup-Z_inf)/100;

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

    X   = x(i,:);
    P   = p(i,:);
    Z   = lamperti_transform(theta_0,alpha,X-P,P,2);
    
    P_exten_linear = interp1([0 6*dt],[P(1) P(6)],minus_t,'linear','extrap');
    P              = [P_exten_linear,P]; % We extend the forecast.
    
    for j = 1:length(P)
        P(j) = max(P(j),epsilon);
        P(j) = min(P(j),1-epsilon);
    end
    P_dot   = (-P(1:end-1) + P(2:end)) / dt; % This starts at 1 end ends at (end - 1).
    for j = 1:length(P)-1
        Theta_t(j) = theta_t(theta_0, alpha, P(j), P_dot(j));
    end
    Z_P = lamperti_transform(theta_0,alpha,P,0,2);
    % Z_P is the Lamperti transform of the forecast.
    
    % We assume that the forecast is perfect one hour before it starts.
    sim_path(1:numPaths,1) = Z_P(1);

    set(0,'CurrentFigure',10); clf(10);
    J = plot(exten_t,Z_P); J.LineWidth = 2;
    hold on;
    J = plot(t,Z); J.LineWidth = 2;
    ylim([Z_inf Z_sup]);
    xlabel('Time'), ylabel('Lamperti Power');
    date = d(i);
    title(date{1});
    
    for k = 1:numPaths
        for j = 1 : length(exten_t)-1
            
            sim_path(k,j+1) = sde_Lamperti_FE(sim_path(k,j),alpha,theta_0,Theta_t(j),dt,P(j),P_dot(j),2);
            
%             if sim_path(k,j+1) >= Z_sup - Z_100 * 5
%                 sim_path(k,j+1) = sim_path(k,j);
%             elseif sim_path(k,j+1) <= Z_inf + Z_100 * 5
%                 sim_path(k,j+1) = sim_path(k,j);
%             end
        end
        
    end
    
    if quantil == 1
        
        for k = 1:length(sim_path(1,:))
            quant(k,:) = quantile(sim_path(:,k),[0.01 0.05 0.25 0.75 0.95 0.99]);
            men(k)     = mean(sim_path(:,k));
        end
        set(0,'CurrentFigure',10); clf(10); hold on;
        ciplot(quant(:,1),quant(:,6),exten_t,[153 255 153]/255)
        ciplot(quant(:,2),quant(:,5),exten_t,[153 255 255]/255)
        ciplot(quant(:,3),quant(:,4),exten_t,[153 153 255]/255)
        J = plot(exten_t,Z_P,'k'); J.LineWidth = 2;
        J = plot(t,Z,'b'); J.LineWidth = 2;
        legend({'99% confidence interval', '90% confidence interval', '50% confidence interval',...
            'Forecast', 'Real production'});
        if xlimit
            xlim([0 max(exten_t)]);
        else
            xlim([min(exten_t) max(exten_t)]);
        end
        ylim([0 1]); grid minor;
        xlabel('Time'), ylabel('Lamperti Power'); box;
        date = d(i);
        ylim([Z_inf Z_sup]);
        title(date{1});
        if save == 1
            saveas(gcf,[pwd '/Results/bands_testing_days/',whatToDo,'/',num2str(i)],'epsc');
        end
        
    elseif quantil == 0
        
        plot(exten_t,sim_path,'r');
        legend('Forecast','Real Production','Simulations','HandleVisibility','off');
        if xlimit
            xlim([0 max(exten_t)]);
        else
            xlim([min(exten_t) max(exten_t)]);
        end
        box; grid minor;
        pause(0.1); box;
        if save == 1
            saveas(gcf,[pwd '/Results/paths_testing_days/',whatToDo,'/',num2str(i)],'epsc');
        end
        
%         figure; plot(exten_t,drift,'r'); grid minor;
        
    end
    
    pause(0.1);
    
end