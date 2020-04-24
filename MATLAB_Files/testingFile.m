close all;
clear all;
clc;

% TESTING:
% v     = 0.2
% th1   = 1
% th2   = 2
% dt    = 0.05
% n     = 10
% p1    = 0.5
% p2    = 0.6
% alpha = 0.1
% m1    = moment_1(v,th1,th2,dt,n)
% m2    = moment_2(v,th1,th2,p1,p2,alpha,m1,dt,n)

epsilon = 0.018;
Ta_Tra_Comp = load_data_eps(epsilon);

Table_Training = Ta_Tra_Comp; % We copy it so we can modify it.

Date              = Ta_Tra_Comp.Date;
Time              = Ta_Tra_Comp.Time;
Forecast          = Ta_Tra_Comp.Forecast;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
Real_ADME         = Ta_Tra_Comp.Real_ADME;
Error             = Ta_Tra_Comp.Error;
Error_Transitions = Ta_Tra_Comp.Error_Transitions;
Lamparti_Data     = Ta_Tra_Comp.Error_Lamp;
Lamparti_Tran     = Ta_Tra_Comp.Error_Lamp_Transitions;

dt             = Time(1,2);
[M, N_ini]     = size(Forecast);
N              = N_ini - 1; % We have N_ini measurements but N samples.

num_days = 127;

% we_remove_badSamples = 0;
% if we_remove_badSamples == 0
%     [Table_Training, batch] = new_batch_fixed(Table_Training,num_days,N);
% elseif we_remove_badSamples == 1
%     [Table_Training, batch] = new_batch_fixed_2(Table_Training,num_days,N);
% end

ini_theta = 10;
ini_alpha = 10;
len_theta = 500;
len_alpha = 300;
div_theta = 100; 
div_alpha = 300;

% Initial parameters:
if 1 == 2
    est = initial_guess(Real_ADME, M, N, dt);
    limit_theta = 2000;
    div_theta   = 1000;
    for i = 1:limit_theta
        val(i) = initial_theta(Error, Forecast, M, N, dt, i/div_theta);
        val2(i) = initial_theta_0(Error, M, N, dt, i/div_theta);
        disp(num2str(i/limit_theta*100));
    end
    theta_ini = find(min(val)==val) / div_theta;
    alpha_ini = est/theta_ini;
end

%% Testing Epsilon:

est = initial_guess(Real_ADME, M, N, dt);

epsilon     = 0.2;
epsilon_new = 0.1;
step        = 1;
plot_eps    = epsilon_new;
while abs(epsilon-epsilon_new)/epsilon > 0.05
    epsilon     = epsilon_new;
    samples     = create_samples_eps(Forecast, Real_ADME, epsilon);
    theta_t     = mean_regression_eps(samples, dt);
    epsilon_new = est / theta_t;
    disp(['========== Step: ',num2str(step),' ==========']);
    disp(['Epsilon = ',num2str(epsilon_new)]);
    disp(['Theta(t) = ',num2str(theta_t)]);
    step = step + 1;
    plot_eps(end+1) = epsilon_new;
end

Forecast_vec = Forecast(:);
disp('==============================');
disp(['Minimum forecast: ',num2str(min(Forecast_vec))]);
disp(['Maximum forecast: ',num2str(max(Forecast_vec))]);
plot(plot_eps); grid minor;

%%

if 1 == 2

    plot([1/div_theta:1/div_theta:2],val); grid minor; hold on;
    title('Value function');
    xlabel('$\theta_0$','interpreter','latex');
    plot(theta_ini,val(theta_ini*div_theta),'*r');
    disp(['Initial \theta_0 = ',num2str(theta_ini)]);
    disp(['Initial \alpha = ',num2str(alpha_ini)]);
    saveas(gcf,[pwd '/Results/initial_guess/initial_theta_1'],'epsc');

    figure;
    plot([1/div_theta:1/div_theta:2],val2); grid minor; hold on;
    title('Value function');
    xlabel('$\theta_0$','interpreter','latex');
    theta_ini_2 = find(min(val2)==val2) / div_theta;
    plot(theta_ini_2,val2(theta_ini_2*div_theta),'*r');
    disp(['Initial \theta_0 using initial_theta_0 is ',num2str(theta_ini_2)]);
    saveas(gcf,[pwd '/log_likelihood_2'],'epsc');
    saveas(gcf,[pwd '/Results/initial_guess/initial_theta_2'],'epsc');

end

if 1 == 1
    val = [];
    % TO PLOT THE LIKELIHOOD:
    
    for i = 0:5720/2-1
        n(i+1) = batch(3,2+i*2) - batch(3,1+i*2);
    end
%     plot(n);
%     grid minor;
%     xlim([1 5720/2]);
    lev1 = []; lev2 = []; lev3 = []; lev4 = []; lev5 = [];
    for j = 1:5720/2
        
        if     0.000 <= abs(n(j)) && abs(n(j)) < 0.010
            lev1(end+1) = j;
        elseif 0.010 <= abs(n(j)) && abs(n(j)) < 0.025
            lev2(end+1) = j;
        elseif 0.025 <= abs(n(j)) && abs(n(j)) < 0.050
            lev3(end+1) = j;
        elseif 0.050 <= abs(n(j)) && abs(n(j)) < 0.100
            lev4(end+1) = j;
        elseif 0.100 <= abs(n(j)) && abs(n(j)) < 0.150
            lev5(end+1) = j;
        end

    end
    
    for z = 0:length(batch(1,:))-1
        
        if 1==1%ismember((z+1) / 2, lev5)
        
            auxBatch = batch(:,1+2*z:2+2*z);

            for i = ini_theta:1:len_theta

                vec = ones(1,len_alpha);

%                 parfor j = ini_alpha:len_alpha
                for j = ini_alpha:len_alpha
                    theta          = i/div_theta;
                    alpha          = (j-1+ini_alpha)/div_alpha;
                    batch_complete = batch_with_theta(auxBatch, theta);
                    vec(j)         = Inf;
                    vec(j)         = -log_LH_evaluation(batch_complete, alpha, dt);
                end

                val(i+1-ini_theta,:) = vec(ini_alpha:len_alpha);

                disp(num2str(i/len_theta*100));

            end

            pause(0.1);
            figure;
            hold on;
            vec_theta = [ini_theta:len_theta] / div_theta;
            vec_alpha = [ini_alpha:len_alpha] / div_alpha;
            [X,Y] = meshgrid(vec_theta,vec_alpha);
            contourf(X,Y,val',20); colorbar;
            xlabel('$\theta_0$','interpreter','latex');
            ylabel('$\alpha$','interpreter','latex');
            title(['Negative Log-likelihoog for ',num2str(num_days),' days']);
            title(['Log-likelihoog 1 Transition: ',num2str(n(z+1))]);
            plot(1.145,0.0823,'-p','MarkerFaceColor','red','MarkerSize',15);
            pause(0.1);
            saveas(gcf,[pwd '/Results/likelihood/level_sets/',num2str(z+1)],'epsc');

        end

    end

    vec_theta = [ini_theta:len_theta] / div_theta;
    vec_alpha = [ini_alpha:len_alpha] / div_alpha;
    [X,Y] = meshgrid(vec_theta,vec_alpha);
    
%     s = surf(X,Y,val'); s.EdgeColor = 'none';
    
    contourf(X,Y,val',20); colorbar;
    
    xlabel('$\theta_0$','interpreter','latex');
    ylabel('$\alpha$','interpreter','latex');
    title(['Negative Log-likelihoog for ',num2str(num_days),' days']);
    view(-45,45); box on; grid minor;
    if 1 == 0
        saveas(gcf,[pwd '/Results/likelihood/neg_log_likelihood'],'epsc');
    end

end

if 1 == 2
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>
    
    theta_0 = 3;
    alpha = 0.03;

    batch_complete = batch_with_theta(batch, alpha, theta_0);
    
%     % >>>
%     
%     norm = 100;
%     times = 1;
%     
%     for delta = 1:1:norm*times
% 
%         pause(0.1);
%         disp(num2str(delta/(norm*times)));
%         [val(delta),m1(delta),m2(delta),xi1(delta),xi2(delta),sig2(delta)] = ...
%             first_log_LH_evaluation(batch_complete, theta_0, alpha, delta/norm, dt, N);
% 
%     end
%     title('Initial error Histogram and All 1-element Likelihoods');
%     grid minor;
%     %>>>>>>>>>>>>>>>>>>>>>>>>>
    
%     figure;
%     plot([1:norm*times]/norm,val);
%     title('-log-likelihoog')
%     figure;
%     plot([1:norm*times]/norm,xi1);
%     title('xi1')
%     figure;
%     plot([1:norm*times]/norm,m2);
%     title('m2')
%     figure;
%     plot([1:norm*times]/norm,sig2);
%     title('sig2')
%     
%     [vall,ind] = min(val);
%     
    % >>>
    
end

theta = 1.145;
alpha = 0.0823;
x0 = [theta, alpha];
num_days = 127;
[Table_Training, batch] = new_batch(Table_Training,num_days,N);
fun = @(x) -likelihood_optimization(batch, x(1), x(2), dt);
options = optimset('PlotFcns',@optimplotfval);

x = fminsearch(fun, x0, options);

%% TESTING:

close all; clear all; clc;

theta = 1.145;
alpha = 0.0823;

dt = 0.0069;
[Ta_Tra_Comp, Ta_Test_Compl, Ta_Comp] = load_data();
Table_Training = Ta_Tra_Comp; % We copy it so we can modify it.
num_days = 127; N = 144;
[Table_Training, batch] = new_batch(Table_Training,num_days,N);
Out = likelihood_optimization(batch, theta, alpha, dt);
disp(num2str(Out));

%% HISTOGRAMS:

close all;
clear all;
clc;

epsilon = 0.018;
Ta_Tra_Comp = load_data_eps(epsilon);

Table_Training = Ta_Tra_Comp; % We copy it so we can modify it.

Date              = Ta_Tra_Comp.Date;
Time              = Ta_Tra_Comp.Time;
Forecast          = Ta_Tra_Comp.Forecast;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
Real_ADME         = Ta_Tra_Comp.Real_ADME;
Error             = Ta_Tra_Comp.Error;
Error_Transitions = Ta_Tra_Comp.Error_Transitions;
Lamparti_Data     = Ta_Tra_Comp.Error_Lamp;
Lamparti_Tran     = Ta_Tra_Comp.Error_Lamp_Transitions;

Forecast_Dot = Forecast_Dot(:)*Time(1,2);

histogram(Forecast_Dot);

%% Testing Z^2 computation:

close all;
clear all;
clc;

epsilon                = 0.018;
Table_Testing_Complete = load_data_eps_test(epsilon);

t  = Table_Testing_Complete.Time(1,:);
dt = t(2);

theta_0 = 1.2656;
alpha   = 0.07052;

z1 = lamperti_transform(theta_0,alpha,-0.07,0.2-0.03);
z2 = lamperti_transform(theta_0,alpha,0.02,0.2-0.03);

P_dot        = (-0.03) / dt;
Theta_t      = theta_t(theta_0, alpha, 0.2-0.03, P_dot);

n = -0.04:0.01:-0.03;

Z_inf = lamperti_transform(theta_0,alpha,0,0);
Z_sup = lamperti_transform(theta_0,alpha,1,0);
Z_inf_2 = lamperti_transform(theta_0,alpha,0,0)^2;
Z_sup_2 = lamperti_transform(theta_0,alpha,1,0)^2;
Z_100 = (Z_sup-Z_inf)/100;

for i = 1:length(n)
    
    figure(i)
    figure(i*1000)
    zz = lamperti_transform(theta_0,alpha,n(i),0.2-0.03);
    disp(num2str(zz*sqrt(2*alpha*theta_0)/(pi/2)));
    disp(['i=',num2str(i),', n = ',num2str(n(i)),', z = ',num2str(zz),', and zz^2 = ',num2str(zz^2)]);
    for k = 1:100
        sim_path_z_2(k) = sde_Lamperti_FE_2_TEXSTING(zz^2,alpha,theta_0,Theta_t,dt,0.2-0.03,P_dot,i,i*1000);
    end

end

%% Testing Plot Likelihood for Model 1 and only 1 transition:

xx = 0:0.01:1;

for i = 1:length(xx)
    
    dist_to_plot(i) = log_dist_L(xx(i),m1(end),sqrt(m2(end)-m1(end)^2));
    
end

dist_to_plot_val = log_dist_L(x2,m1(end),sqrt(m2(end)-m1(end)^2));

figure;
plot(xx,dist_to_plot);
hold on;
plot(x2,dist_to_plot_val,'x');



























