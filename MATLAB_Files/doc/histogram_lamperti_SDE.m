close all;
clear all;
clc;

epsilon                = 0.018;
dataKind               = 'comparable';
dataSet                = 'MTLOG';
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet,dataKind);
% whatToDo               = 'Optimal'; 
whatToDo               = 'Optimal_Lamperti'; 

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
save     = 0;
delta    = 21; % The time is delta*10 minutes.
numPaths = 1;

if  strcmp(whatToDo,'Optimal')
    theta_0 = 1.18;
    alpha   = 0.07;
elseif  strcmp(whatToDo,'Optimal_Lamperti')
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

minus_t = linspace(-delta*dt,-dt,delta);
exten_t = [minus_t,t];

for i = 1 : height(Table_Testing_Complete)

    X = x(i,:);
    P = p(i,:);
    Z = lamperti_transform(theta_0,alpha,X-P,P,2);
    
    P_exten_linear = interp1([0 6*dt],[P(1) P(6)],minus_t,'linear','extrap');
    P              = [P_exten_linear,P]; % We extend the forecast.
    for j = 1:length(P)
        P(j) = max(P(j),epsilon);
        P(j) = min(P(j),1-epsilon);
    end
    Z_P = lamperti_transform(theta_0,alpha,P,0,2);
%     Z_P is the Lamperti transform of the forecast.

    P_dot   = (-P(1:end-1) + P(2:end)) / dt; % This starts at 1 end ends at (end - 1).
    for j = 1:length(P)-1
        Theta_t(j) = theta_t(theta_0, alpha, P(j), P_dot(j));
    end
    
    % We assume that the forecast is perfect one hour before it starts.
    sim_path(1:numPaths,1) = Z_P(1);
    
    for k = 1:numPaths
        for j = 1 : length(exten_t)-1
            
            sim_path(k,j+1) = sde_Lamperti_FE(sim_path(k,j),alpha,theta_0,Theta_t(j),dt,P(j),P_dot(j),2);            
            if sim_path(k,j+1) >= Z_sup - Z_100 * 5
                sim_path(k,j+1) = sim_path(k,j);
            elseif sim_path(k,j+1) <= Z_inf + Z_100 * 5
                sim_path(k,j+1) = sim_path(k,j);
            end
            
        end
    end
    
    z_trans(i,:) = Z(2:end) - Z(1:end-1);
    s_trans(i,:) = sim_path(delta+2:end) - sim_path(delta+1:end-1);
    
    % Next transformation:
    for n = 2:length(Z)-1
    
        % The first transition is from 1 to 2 because we use information
        % in time 0 to construct it.
        Delta_new_Transform(n-1) = (Z(n+1)-2*Z(n)+Z(n-1)) / sqrt(dt);
        
    end
    for n = 2:length(Z)
        new_Transform(n-1) = (Z(n)-Z(n-1)) / sqrt(dt);
    end
    
    new_Transform_trans(i,:)   = Delta_new_Transform;
    new_Transform_allData(i,:) = new_Transform;
        
end

% figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
figure;
h1 = histogram(z_trans);
hold on; grid minor;
title('Lamperti (Data and Simulated) Transitions');
h2 = histogram(s_trans);
h1.Normalization = 'pdf';
h2.Normalization = 'pdf';
h2.FaceAlpha = 0;
h2.FaceColor = [0 0 0];
h2.LineWidth = 2;
h2.EdgeColor = 'r';
xlim([-0.5 0.5]);
xlabel('Value of Transition');
ylabel('Probability');
% xlim([-0.4 0.4]);
% xlim([-1 1]);
set(gca,'FontSize',12);
legend('Lamperti Transitions Histogram','Simulated Lamperti Transitions Histogram');
saveas(gcf,[pwd '/Results/histograms/lamperti/',whatToDo],'epsc');

% figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
figure;
h3 = histogram(new_Transform_allData);
h3.Normalization = 'pdf';
set(gca,'FontSize',12);
disp(['Mean of the second transformation: ',num2str(mean(new_Transform_trans(:)))]);
disp(['Variance of the second transformation: ',num2str(var(new_Transform_trans(:)))]);