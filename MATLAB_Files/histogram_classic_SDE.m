close all;
clear all;
clc;

epsilon                = 0.020;
dataKind               = 'comparable';
dataSet                = 'MTLOG';
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet,dataKind);
whatToDo               = 'Optimal'; 
% whatToDo               = 'Lamperti_Optimal'; 
% whatToDo               = 'Initial_Guess'; 

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
save     = 0;
delta    = 13; % The time is delta*10 minutes.

if  strcmp(whatToDo,'Optimal')
    theta_0 = 1.930;
    alpha   = 0.050;
elseif  strcmp(whatToDo,'Lamperti_Optimal')
    theta_0 = 2.450;
    alpha   = 0.097/2.450;
    % V1: 
    theta_0 = 2.8150; alpha = 0.1276;
    % V2: 
    theta_0 = 1.8708; alpha = 0.0431;
elseif  strcmp(whatToDo,'Initial_Guess')
    theta_0 = 1.6290;
    alpha   = 0.06;
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

total_x_s_trans = [];
total_v_s_trans = [];

for l = 1:1
    for i = 1 : height(Table_Testing_Complete)

        X = x(i,:);
        P = p(i,:);
        % All of them show similar results:
        P_exten_linear = interp1([0 6*dt],[P(1) P(6)],minus_t,'linear','extrap');
        P_exten_pchip  = interp1([0 6*dt],[P(1) P(6)],minus_t,'pchip','extrap');
        P_exten_spline = interp1([0 6*dt],[P(1) P(6)],minus_t,'spline','extrap');
        P_exten_makima = interp1([0 6*dt],[P(1) P(6)],minus_t,'makima','extrap');

        initial_error_V     = X(1)-P(1);
        initial_error_slope = initial_error_V / (delta*dt);
        for j = 1:delta
            X_exten(j) = initial_error_slope*(j-1)*dt;
        end

        P       = [P_exten_linear,P]; % We extend the forecast.
        for j = 1:length(P)
            P(j) = max(P(j),epsilon);
            P(j) = min(P(j),1-epsilon);
        end
        P_dot   = (-P(1:end-1) + P(2:end)) / dt; % This starts at 1 end ends at (end - 1).
        for j = 1:length(P)-1
            Theta_t(j) = theta_t(theta_0, alpha, P(j), P_dot(j));
        end

        % We assume that the forecast is perfect *delta* hours before it starts.
        sim_path(1,1) = P(1);

        for j = 1 : length(exten_t)-1
            sim_path(j+1) = sde_FE(sim_path(j),alpha,theta_0,Theta_t(j),dt,P(j),P_dot(j));
            if sim_path(j+1) > 1
                sim_path(j+1) = 1;
            elseif sim_path(j+1) < 0
                sim_path(j+1) = 0;
            end
        end

        X     = [X_exten,X];
        V     = X - P;
        V_sim = sim_path - P;
        
        x_trans(i,:)   = X(delta+2:end) - X(delta+1:end-1);
        x_s_trans(i,:) = sim_path(delta+2:end) - sim_path(delta+1:end-1);
        v_trans(i,:)   = V(delta+2:end) - V(delta+1:end-1);
        v_s_trans(i,:) = V_sim(delta+2:end) - V_sim(delta+1:end-1);

    end
    
    total_x_s_trans = [total_x_s_trans,x_s_trans];
    total_v_s_trans = [total_v_s_trans,v_s_trans];

end

% figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
figure;
h1 = histogram(v_trans);
hold on; grid minor; title('Error (Data and Simulated) Transitions');
h2 = histogram(total_v_s_trans);
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
set(gca,'FontSize',12);
saveas(gcf,[pwd '/Results/histograms/classic/',whatToDo],'epsc');

%% With no Delta:

close all;
clear all;
clc;

epsilon                = 0.020;
dataKind               = 'comparable';
dataSet                = 'MTLOG';
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet,dataKind);
whatToDo               = 'Optimal'; 
whatToDo               = 'Lamperti_Optimal'; 
% whatToDo               = 'Initial_Guess'; 

% PARAMETERS:
% set(0,'defaultAxesFontSize',18);
save     = 0;
numPaths = 1;

if  strcmp(whatToDo,'Optimal')
    theta_0 = 1.930;
    alpha   = 0.050;
elseif  strcmp(whatToDo,'Lamperti_Optimal')

% V1: 
%     theta_0 = 2.8150; alpha = 0.1276;
% V2: 
%     theta_0 = 1.8708; alpha = 0.0431;
elseif  strcmp(whatToDo,'Initial_Guess')
    theta_0 = 1.6290;
    alpha   = 0.06;
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

total_x_s_trans = [];
total_v_s_trans = [];

for l = 1:5
    for i = 1 : height(Table_Testing_Complete)

        X = x(i,:);
        P = p(i,:);

        for j = 1:length(P)
            P(j) = max(P(j),epsilon);
            P(j) = min(P(j),1-epsilon);
        end
        P_dot   = (-P(1:end-1) + P(2:end)) / dt; % This starts at 1 end ends at (end - 1).
        for j = 1:length(P)-1
            Theta_t(j) = theta_t(theta_0, alpha, P(j), P_dot(j));
        end

        % We assume that the forecast is perfect *delta* hours before it starts.
        sim_path(1,1) = X(1);

        for j = 1 : length(t)-1
            sim_path(j+1) = sde_FE(sim_path(j),alpha,theta_0,Theta_t(j),dt,P(j),P_dot(j));
            if sim_path(j+1) > 1
                sim_path(j+1) = 1;
            elseif sim_path(j+1) < 0
                sim_path(j+1) = 0;
            end
        end

        V     = X - P;
        V_sim = sim_path - P;
        
        x_trans(i,:)   = X(2:end) - X(1:end-1);
        x_s_trans(i,:) = sim_path(2:end) - sim_path(1:end-1);
        v_trans(i,:)   = V(2:end) - V(1:end-1);
        v_s_trans(i,:) = V_sim(2:end) - V_sim(1:end-1);

    end
    
    total_x_s_trans = [total_x_s_trans,x_s_trans];
    total_v_s_trans = [total_v_s_trans,v_s_trans];
    
end

% figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
figure;
h1 = histogram(v_trans);
hold on; grid minor; title('Error (Data and Simulated) Transitions');
h2 = histogram(total_v_s_trans);
h1.Normalization = 'pdf';
h2.Normalization = 'pdf';
% h1.FaceAlpha = 1;
h2.FaceAlpha = 0;
h2.FaceColor = [0 0 0];
h2.LineWidth = 1;
h2.EdgeColor = 'r';
xlabel('Value of Transition');
ylabel('Probability');
xlim([-0.1 0.1]);
legend('Error Transitions Histogram','Simulated Error Transitions Histogram',...
    'location','southoutside');
set(gca,'FontSize',12);
saveas(gcf,[pwd '/Results/histograms/classic_noDelta/',whatToDo],'epsc');