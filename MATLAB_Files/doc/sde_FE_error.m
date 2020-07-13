function sim_path_f = sde_FE_error(sim_path_i,alpha,theta_0,Theta_t,dt,P) % 01/04/2020 19:53

    % SDE for V_t.

    sim_path_f = sim_path_i - Theta_t*sim_path_i*dt +...
                sqrt(dt)*sqrt(2*theta_0*alpha*(sim_path_i+P)*(1-sim_path_i-P))*randn(1);
    
end