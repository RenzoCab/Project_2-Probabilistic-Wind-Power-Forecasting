function sim_path_f = sde_FE_MODEL1(sim_path_i,alpha,theta_0,dt,P) % 02/02/2020 18:28

    % SDE for X_t in MODEL 1.

    sim_path_f = sim_path_i -theta_0*(sim_path_i-P)*dt +...
                sqrt(dt)*sqrt(2*theta_0*alpha*sim_path_i*(1-sim_path_i))*randn(1);
    
end