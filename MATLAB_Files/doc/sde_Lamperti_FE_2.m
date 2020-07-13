function sim_path_f = sde_Lamperti_FE_2(sim_path_i,alpha,theta_0,Theta_t,dt,P,P_dot)

    % 02/04/2020 20:40
    
    n  = 50;
    dt = dt / n;
    fv = sim_path_i;
    
    for i = 1:n

        z = -sqrt(sim_path_i);
        b = sde_Lamperti_drift_cand2(z,alpha,theta_0,Theta_t,P,P_dot);

        sim_path_i = sim_path_i + dt * (2*b*z + 1) + ...
        2 * z * sqrt(dt) * randn(1);
        sim_path_i = max(0,sim_path_i); % ARTIFICIAL.
        
    end
    
    sim_path_f = sim_path_i;
    
end