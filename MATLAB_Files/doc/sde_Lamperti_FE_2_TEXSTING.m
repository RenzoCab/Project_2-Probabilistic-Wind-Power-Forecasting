function sim_path_f = sde_Lamperti_FE_2_TEXSTING(sim_path_i,alpha,theta_0,Theta_t,dt,P,P_dot,m,mm) % 02/04/2020 20:40

%     b = sde_Lamperti_drift(sim_path_i,alpha,theta_0,Theta_t,P,P_dot);
% 
%     sim_path_f = sim_path_i + dt * (2*b*sim_path_i + 1) + ...
%         2 * sim_path_i * sqrt(dt) * randn(1);
    
    n  = 200;
    dt = dt / n;
    fv = sim_path_i;
    
    for i = 1:n
        
        forplot(i) = sim_path_i;
        b = sde_Lamperti_drift(sim_path_i,alpha,theta_0,Theta_t,P,P_dot);
        bforplot(i) = b;
        
        sim_path_i = sim_path_i + dt * (2*b*sim_path_i + 1) + ...
        2 * sim_path_i * sqrt(dt) * randn(1);
        
    end
    
    sim_path_f = sim_path_i;
    
    set(0,'CurrentFigure',m);
    hold on;
    plot(forplot);
    set(0,'CurrentFigure',mm);
    hold on;
    plot(bforplot);
    
    
end