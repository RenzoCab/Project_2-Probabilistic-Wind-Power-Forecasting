function sim_path_f = sde_Lamperti_FE(sim_path_i,alpha,theta_0,Theta_t,dt,P,P_dot,candidate) 

    % 03/04/2020 20:37

    if candidate == 1
    
        sim_path_f = sim_path_i + ...
            dt * sde_Lamperti_drift_cand1(sim_path_i,alpha,theta_0,Theta_t,P,P_dot) + ...
            sqrt(dt)*randn(1);
    
    elseif candidate == 2
        
        sim_path_f = sim_path_i + ...
            dt * sde_Lamperti_drift_cand2(sim_path_i,alpha,theta_0,Theta_t,P,P_dot) + ...
            sqrt(dt)*randn(1);
    
    else
        
        error('Choose a correct CANDIDATE.');
        
    end
    
end