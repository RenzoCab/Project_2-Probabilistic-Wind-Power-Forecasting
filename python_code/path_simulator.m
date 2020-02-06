close all;
clear all;
clc;

load('data/cleansed/MTLOG_0100_and_Real_24h_Testing_Data.mat');

theta = 25;
alpha = 10;

d  = Table_Testing_Complete.Date;
p  = Table_Testing_Complete.Forecast;
t  = Table_Testing_Complete.Time(1,:);
x  = Table_Testing_Complete.Real_ADME;
v  = Table_Testing_Complete.Error;
dt = t(2); % This is the default dt. I will choose a smaller one.
% If the dt is not small enough, we always scape from [0 1].
dt = dt;

figure(10);

minus_t = linspace(-6*dt,-dt,6); % Time from -01:00 (or 23:00) to 00:00 (not including the 00:00).
exten_t = [minus_t,t]; % We extend two hours backward.

for i = 1 : height(Table_Testing_Complete)

    X = x(i,:);
    P = p(i,:);
    % All of them show similar results:
    P_exten_linear = interp1([0 6*dt],[P(1) P(6)],minus_t,'linear','extrap');
    P_exten_pchip  = interp1([0 6*dt],[P(1) P(6)],minus_t,'pchip','extrap');
    P_exten_spline = interp1([0 6*dt],[P(1) P(6)],minus_t,'spline','extrap');
    P_exten_makima = interp1([0 6*dt],[P(1) P(6)],minus_t,'makima','extrap');
    
    P       = [P_exten_linear,P]; % We extend the forecast.
    P_dot   = (P(1:end-1) - P(2:end)) / dt; % This starts at 1 end ends at (end - 1).
    theta_t = max(theta, abs(P_dot)./min(P(1:end-1),1-P(1:end-1)));
    
    % We assume that the forecast is perfect one hour before it starts.
    numPaths = 1000;
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

            sim_path(k,j+1) = sim_path(k,j) + (P_dot(j)-theta_t(j)*(sim_path(k,j)-P(j)))*dt +...
                sqrt(dt)*sqrt(2*theta_t(j)*alpha*sim_path(k,j)*(1-sim_path(k,j)))*randn(1);
            if sim_path(k,j+1) > 1
                sim_path(k,j+1) = 1;
            elseif sim_path(k,j+1) < 0
                sim_path(k,j+1) = 0;
            end
            
        end
%         plot(exten_t,sim_path(k,:));
    end
    
    for k = 1:length(sim_path(1,:))
        quant(k,:) = quantile(sim_path(:,k),[0.05 0.25 0.75 0.95]);
    end
    plot(exten_t,quant(:,:),'LineWidth',3);
    
    legend('Forecast','Real Production','HandleVisibility','off');
    saveas(gcf,[pwd '/plots/matlab/',num2str(i)],'epsc');
    pause(0.01);
    
end