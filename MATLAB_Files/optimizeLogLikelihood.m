close all;
clear all;
clc;

likelihood    = 'normal'; % 'normal' or 'lamperti'.
epsilon       = 0.018;
[Ta_Tra_Comp] = load_data_eps(epsilon);

Date              = Ta_Tra_Comp.Date;
Time              = Ta_Tra_Comp.Time;
Forecast          = Ta_Tra_Comp.Forecast;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
Real_ADME         = Ta_Tra_Comp.Real_ADME;
Error             = Ta_Tra_Comp.Error;
Error_Transitions = Ta_Tra_Comp.Error_Transitions;
Lamparti_Data     = Ta_Tra_Comp.Error_Lamp;
Lamparti_Tran     = Ta_Tra_Comp.Error_Lamp_Transitions;

dt         = Time(1,2);
[M, N_ini] = size(Forecast);
N          = N_ini - 1; % We have N_ini measurements but N samples.

num_days   = 127; % Maximum 127.
[~, batch] = new_batch_fixed(Ta_Tra_Comp,num_days,N);

%% Initial parameters:

est           = initial_guess(Real_ADME, M, N, dt);
gamma         = 0.35; 
samples_gamma = create_samples_minus_eps(Forecast, Error, gamma);
theta_ini     = mean_regression_eps(samples_gamma, dt);
alpha_ini     = est/theta_ini;

%% Optimization:

theta_0 = theta_ini;
alpha   = alpha_ini;

if strcmp(likelihood,'normal')
    
    x0    = [theta_0, alpha];
    fun   = @(x) -likelihood_optimization(batch, x(1), x(2), dt);
    
    x = fminunc(fun, x0);
    
    % fminsearch -> Theta_0 = 1.135, Alpha = 0.073 (05/04/2020).
    % fmincon    -> Theta_0 = 1.366, Alpha = 0.061 (05/04/2020).
    % fminunc    -> Theta_0 = 1.628, Alpha = 0.051 (05/04/2020).
    
%     D_theta_0 = 4;
%     D_alpha   = 1;
%     div_fact  = 10;
%     x_end = [theta_0, alpha];  
%     opt   = x_end;
%     circ_fact = 1;
%        
%     for k = 1:50
% 
%         disp('=============================');
%         disp(x_end);
%         disp(k);
%         disp('=============================');
%         x0 = [0 0];
%         
%         while not(isequal(x0,x_end))
%             
%             x0   = x_end;
%             eval = fun(x0);
%             
%             x1 = x0 + [0 D_alpha/(div_fact*sqrt(k))];
%             
%             x2 = x0 + [D_theta_0*circ_fact/(div_fact*sqrt(k)) D_alpha*circ_fact/(div_fact*sqrt(k))];
%             
%             x3 = x0 + [D_theta_0/(div_fact*sqrt(k)) 0];
%             
%             x4 = x0 + [D_theta_0*circ_fact/(div_fact*sqrt(k)) -D_alpha*circ_fact/(div_fact*sqrt(k))];
%             
%             x5 = x0 + [0 -D_alpha/(div_fact*sqrt(k))];
%             
%             x6 = x0 + [-D_theta_0*circ_fact/(div_fact*sqrt(k)) -D_alpha*circ_fact/(div_fact*sqrt(k))];
%             
%             x7 = x0 + [-D_theta_0/(div_fact*sqrt(k)) 0];
%             
%             x8 = x0 + [-D_theta_0*circ_fact/(div_fact*sqrt(k)) D_alpha*circ_fact/(div_fact*sqrt(k))];
%             
%             eval1 = fun(x1);
%             if eval1 < eval
%                 eval = eval1;
%                 x_end   = x1;
%             end
%             
%             eval2 = fun(x2);
%             if eval2 < eval
%                 eval = eval2;
%                 x_end   = x2;
%             end
%             
%             eval3 = fun(x3);
%             if eval3 < eval
%                 eval = eval3;
%                 x_end   = x3;
%             end
%             
%             eval4 = fun(x4);
%             if eval4 < eval
%                 eval = eval4;
%                 x_end   = x4;
%             end
%             
%             eval5 = fun(x5);
%             if eval5 < eval
%                 eval = eval5;
%                 x_end   = x5;
%             end
%             
%             eval6 = fun(x6);
%             if eval6 < eval
%                 eval = eval6;
%                 x_end   = x6;
%             end
%             
%             eval7 = fun(x7);
%             if eval7 < eval
%                 eval = eval7;
%                 x_end   = x7;
%             end
%             
%             eval8 = fun(x8);
%             if eval8 < eval
%                 eval = eval8;
%                 x_end   = x8;
%             end
%             
%             opt(end+1,:) = x_end;
%             
%         end
%         
%     end
    
elseif strcmp(likelihood,'lamperti')
    
    x0 = [0.1 0.1];
    x1 = [1.14, 0.073]; % Output from 'normal' optimization.
    
    while norm(x1-x0)/norm(x1) > 0.01
        
        x0        = x1;
        [x1,fval] = likelihood_nested_optimization_L(batch, x0, dt);
    
        % Theta_0 = 1.57, Alpha = 0.034 (04/04/2020).
    
    end
    
else
    
	disp('Please, choose normal or lamperti.');
    
end