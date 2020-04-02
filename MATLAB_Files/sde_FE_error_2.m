function sim_path_f = sde_FE_error_2(sim_path_i,alpha,theta_0,Theta_t,dt,P) % 01/04/2020 19:53

    % SDE for (V_t)^2.

    sim_path_f = sim_path_i + (2*alpha*theta_0*(1-P-sim_path_i)*(P+sim_path_i) - ...
        sim_path_i^2*2*Theta_t) * dt + 2*sim_path_i*sqrt(2)*sqrt(alpha*theta_0 * ...
        (1-P-sim_path_i)*(P+sim_path_i)) * sqrt(dt) * randn(1);
    
end