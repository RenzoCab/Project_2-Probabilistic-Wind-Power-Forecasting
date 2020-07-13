% Supplementary Material for the Wind-Power Paper.
% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email: Renzo.CaballeroRosas@kaust.edu.sa CaballeroRenzo@hotmail.com
% Website: https://orcid.org/0000-0003-3220-0923.
% July 2020; Last revision: 13/07/2020.

close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet  = 'MTLOG';
% epsilon can be 0.035 (B), 0.020 (A) or 0.028 (C).
epsilon  = 0.020;
% dataKind can be classic or comparable.
dataKind = 'comparable';
set_type = 'testing'; % 'testing' or 'training'.
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet,dataKind,set_type);
whatToDo = 'optimal_value';
% whatToDo = 'Optimal_Lamperti';

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
quantil  = 0;
save     = 0;
delta    = 13; % The time is delta*10 minutes.
xlimit   = 1; % If this in 1, the plots start at time 0. Otherwise, at -delta.
norm_lim = 0; % If 1, x=[0,1], if 0, x=[1PM 1PM].
% Also, when norm_lim = 0, we check the real accuracy for the bands.
new_time = [13-(1/6)*delta:1/6:37]/24;
counting   = zeros(4,1);
counting_2 = zeros(4,1);

if quantil
    numPaths = 5000;
else
    numPaths = 5;
end

if  strcmp(whatToDo,'optimal_value')
    
    if  strcmp(set_type,'testing') % We do inference with training.
        theta_0 = 1.93;
        alpha   = 0.05;
    elseif strcmp(set_type,'training') % We do inference with testing.
        theta_0 = 1.59;
        alpha   = 0.056;
    end
    
elseif  strcmp(whatToDo,'Optimal_Lamperti')
    
%     theta_0 = 2.200;
%     alpha   = 0.038;
% V2:
    theta_0 = 1.8708;
    alpha   = 0.0431;    
    
else
    error('Wrong whatToDo script!');
end

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
if norm_lim == 1
    exten_t = [minus_t,t];
elseif norm_lim == 0
    exten_t = new_time;
    t       = [13:1/6:37]/24;
    t_ticks = [13:37]/24;
%     t_ticks = [13:2:37]/24;
end

for i = 1 : height(Table_Testing_Complete)

    X = x(i,:);
    P = p(i,:);
    % All of them show similar results:
    P_exten_linear = interp1([0 6*dt],[P(1) P(6)],minus_t,'linear','extrap');
    P_exten_pchip  = interp1([0 6*dt],[P(1) P(6)],minus_t,'pchip','extrap');
    P_exten_spline = interp1([0 6*dt],[P(1) P(6)],minus_t,'spline','extrap');
    P_exten_makima = interp1([0 6*dt],[P(1) P(6)],minus_t,'makima','extrap');
    
    P       = [P_exten_linear,P]; % We extend the forecast.
    for j = 1:length(P)
        P(j) = max(P(j),epsilon);
        P(j) = min(P(j),1-epsilon);
    end
    P_dot   = (-P(1:end-1) + P(2:end)) / dt; % This starts at 1 end ends at (end - 1).
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
    date_format = datetime(date{1},'InputFormat','yyyyMMdd');
    date_format_next = date_format + days(1);
    title([datestr(date_format),' and ',datestr(date_format_next)]);
    
    for k = 1:numPaths
        for j = 1 : length(exten_t)-1
            sim_path(k,j+1) = sde_FE(sim_path(k,j),alpha,theta_0,Theta_t(j),dt,P(j),P_dot(j));
            if sim_path(k,j+1) > 1
                sim_path(k,j+1) = 1;
            elseif sim_path(k,j+1) < 0
                sim_path(k,j+1) = 0;
            end
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
        J = plot(exten_t,P,'k'); J.LineWidth = 2;
        J = plot(t,X,'b'); J.LineWidth = 2;
        legend({'99% confidence interval', '90% confidence interval', '50% confidence interval',...
            'Forecast', 'Real production'});

        xticks(t_ticks);
        datetick('x','HHPM','keepticks');
        xtickangle(90);

        if xlimit
            xlim([min(t_ticks) max(t_ticks)]);
        else
            xlim([min(new_time) max(new_time)]);
        end

        ylim([0 1]); grid minor;
        xlabel('Time'), ylabel('Power'); box;
        date = d(i);
        date_format = datetime(date{1},'InputFormat','yyyyMMdd');
        date_format_next = date_format + days(1);
        title([datestr(date_format),' and ',datestr(date_format_next)]);
        if save == 1
            saveas(gcf,[pwd '/Results/bands_testing_days/',whatToDo,'/withDate/',num2str(i)],'epsc');
        end
                    
    else
        
        plot(exten_t,sim_path,'r');
        legend('Forecast','Real Production','Simulations','HandleVisibility','off');

        xticks(t_ticks);
        datetick('x','HHPM','keepticks');
        xtickangle(90);

        pause(0.1); box;
        if xlimit
            xlim([min(t_ticks) max(t_ticks)]);
        else
            xlim([min(new_time) max(new_time)]);
        end

        box;
        if save == 1
            saveas(gcf,[pwd '/Results/paths_testing_days/',whatToDo,'/withDate/',num2str(i)],'epsc');
        end
        
    end
    
    pause(0.1);
    
end