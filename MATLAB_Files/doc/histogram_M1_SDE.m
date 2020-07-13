close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon       = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet);
% whatToDo               = 'Optimal'; 
% whatToDo               = 'Lamperti_Optimal'; 
% whatToDo               = 'Initial_Guess'; 
whatToDo               = 'M1'; 

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
save     = 0;
delta    = 21; % The time is delta*10 minutes.
numPaths = 1;

if  strcmp(whatToDo,'Optimal')
    theta_0 = 1.180;
    alpha   = 0.070;
elseif  strcmp(whatToDo,'Lamperti_Optimal')
    theta_0 = 2.200;
    alpha   = 0.038;
elseif  strcmp(whatToDo,'Initial_Guess')
    theta_0 = 1.6290;
    alpha   = 0.06;
elseif  strcmp(whatToDo,'M1')
    theta_0 = 0.06;
    alpha   = 1.5100;
%     theta_0 = -1.2201;
%     alpha   = -0.075142;
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
    
    for k = 1:numPaths
        for j = 1 : length(exten_t)-1
            sim_path(k,j+1) = sde_FE_MODEL1(sim_path(k,j),alpha,theta_0,dt,P(j));
            if sim_path(k,j+1) > 1
                sim_path(k,j+1) = 1;
            elseif sim_path(k,j+1) < 0
                sim_path(k,j+1) = 0;
            end
        end
    end
    
    x_trans(i,:) = X(2:end) - X(1:end-1);
    s_trans(i,:) = sim_path(delta+2:end) - sim_path(delta+1:end-1);
    
end

figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
h1 = histogram(x_trans);
hold on; grid minor; title('Error (Data and Simulated) Transitions');
h2 = histogram(s_trans);
h1.Normalization = 'pdf';
h2.Normalization = 'pdf';
% h1.FaceAlpha = 1;
h2.FaceAlpha = 0;
h2.FaceColor = [0 0 0];
h2.LineWidth = 2;
h2.EdgeColor = 'r';
xlabel('Value of Transition');
ylabel('Probability');
xlim([-0.1 0.1]);
legend('Error Transitions Histogram','Simulated Error Transitions Histogram');
set(gca,'FontSize',18);
saveas(gcf,[pwd '/Results/histograms/model1/',whatToDo],'epsc');