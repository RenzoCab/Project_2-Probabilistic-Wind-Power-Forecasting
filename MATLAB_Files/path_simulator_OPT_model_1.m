close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet  = 'MTLOG';
% epsilon can be 0.035 (B), 0.020 (A) or 0.028 (C).
epsilon  = 0.020;
dataKind = 'comparable';
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet,dataKind);

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
quantil  = 1;
save     = 0;
delta    = 22; % The time is delta*10 minutes.
xlimit   = 1; % If this in 1, the plots start at time 0. Otherwise, at -delta.
norm_lim = 0; % If 1, x=[0,1], if 0, x=[1PM 1PM].
new_time = [13-(1/6)*delta:1/6:37]/24;
new_ti   = [1:1/6:25];

if quantil
    numPaths = 5000;
else
    numPaths = 5;
end

% (A) Theta = 0.02081, Alpha = 4.4421, eps = 0.018, prod = 0.0924, f = -52305.
% (B) Theta = 0.02584, Alpha = 3.7390, eps = 0.035, prod = 0.0966, f = -30347.
% (C) Theta = 0.02479, Alpha = 3.9996, eps = 0.028, prod = 0.0992, f = -50455.

theta = 1.2201/1.12201*1.93;
alpha = 0.075142/1.93*1.12201;

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
        Theta_t(j) = theta_t(theta, alpha, P(j), P_dot(j));
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
            sim_path(k,j+1) = sim_path(k,j) + (-theta*(sim_path(k,j)-P(j)))*dt +...
                sqrt(dt)*sqrt(2*theta*alpha*sim_path(k,j)*(1-sim_path(k,j)))*randn(1);
            if sim_path(k,j+1) > 1
                sim_path(k,j+1) = 1;
            elseif sim_path(k,j+1) < 0
                sim_path(k,j+1) = 0;
            end
        end
    end
    
    if quantil == 1
        
        if norm_lim == 1
        
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
            
            if xlimit == 1
                xlim([0 max(exten_t)]);
            elseif xlimit == 0
                xlim([min(exten_t) max(exten_t)]);
            end
            
            ylim([0 1]); grid minor;
            xlabel('Time'), ylabel('Power'); box;
            date_format = datetime(date{1},'InputFormat','yyyyMMdd');
            date_format_next = date_format + days(1);
            title([datestr(date_format),' and ',datestr(date_format_next)]);
            if save == 1
                saveas(gcf,[pwd '/Results/bands_testing_days/model_1/',num2str(i)],'epsc');
            end
        
        elseif norm_lim == 0
            
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
%             legend({'99% confidence interval', '90% confidence interval', '50% confidence interval',...
%             	'Forecast', 'Real production'},'location','southeast');

            
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
                saveas(gcf,[pwd '/Results/bands_testing_days/model_1/withDate/',num2str(i)],'epsc');
            end
            
        end
        
    elseif quantil == 0
        
        if norm_lim == 1
            
            plot(exten_t,sim_path,'r');
            legend('Forecast','Real Production','Simulations','HandleVisibility','off');
            if xlimit
                xlim([0 max(exten_t)]);
            else
                xlim([min(exten_t) max(exten_t)]);
            end
            box;
            pause(0.1); box;
            if save == 1
                saveas(gcf,[pwd '/Results/paths_testing_days/model_1/',num2str(i)],'epsc');
            end
        
        elseif norm_lim == 0
            
            plot(new_time,sim_path,'r');
            legend('Forecast','Real Production','Simulations','HandleVisibility','off');
            box;
            
            xticks(t_ticks);
            datetick('x','HHPM','keepticks');
            xtickangle(90);
            
            pause(0.1); box;
            if xlimit
                xlim([min(t_ticks) max(t_ticks)]);
            else
                xlim([min(new_time) max(new_time)]);
            end

            if save == 1
                saveas(gcf,[pwd '/Results/paths_testing_days/model_1/withDate/',num2str(i)],'epsc');
            end
            
        end
        
    end
    
    pause(0.1);
    
end