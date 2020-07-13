close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet  = 'MTLOG';
% epsilon can be 0.035 (B), 0.020 (A) or 0.028 (C).
epsilon  = 0.020;
% dataKind can be classic or comparable.
dataKind = 'comparable';
Table_Testing_Complete = load_data_eps_test(epsilon,dataSet,dataKind);

delta   = 13; % The time is delta*10 minutes.
theta_0 = 2.2210;
alpha   = 0.0436;

d  = Table_Testing_Complete.Date;
p  = Table_Testing_Complete.Forecast;
t  = Table_Testing_Complete.Time(1,:);
x  = Table_Testing_Complete.Real_ADME;
dt = t(2); % This is the default dt. I will choose a smaller one.

minus_t = linspace(-delta*dt,-dt,delta);
exten_t = [minus_t,t];
NT      = 1000000;
NX      = 1000;

for i = 1 : height(Table_Testing_Complete)

    X = x(i,:);
    P = p(i,:);
    % All of them show similar results:
    P_exten_linear = interp1([0 6*dt],[P(1) P(6)],minus_t,'linear','extrap');
    P              = [P_exten_linear,P]; % We extend the forecast.
    dP             = (P(2:end) - P(1:end-1)) / dt; dP(end+1) = dP(end);
    
    pde_time = linspace(exten_t(1),exten_t(end),NT);
    short_t  = linspace(t(1),t(end),NT);
    P        = interp1(exten_t,P,pde_time);
    dP       = interp1(exten_t,dP,pde_time);
    X        = interp1(t,X,short_t);
    dt       = pde_time(2) - pde_time(1);
    
    pde_space    = linspace(0,1,NX);
    pde_val(1,:) = normpdf(pde_space,P(1),0.1);
    d_g          = zeros(1,NX);
    dd_h         = zeros(1,NX);
    
    for j = 1:length(pde_time)-1
        
        ddp = dP(j);
        pp  = P(j);
        tt  = theta_t(theta_0, alpha, pp, ddp);
        
        g = (ddp-tt*(pde_space+pp)).*pde_val(j,:);
        h = alpha*theta_0*pde_space.*(1-pde_space).*pde_val(j,:);
        
        for k = 2:NX-1
            
            if g(k) < 0
                d_g(k) = (g(k+1)-g(k))/dt;
            elseif g(k) >= 0
                d_g(k) = (g(k)-g(k-1))/dt;
            end
            dd_h(k) = (h(k+1)-2*h(k)-h(k-1))/dt^2;
            
        end
        
        plot(pde_val(j,:));
        hold on;
        
        pde_val(j+1,:) = pde_val(j,:) + dt*(-d_g+dd_h);
        
    end
    
end