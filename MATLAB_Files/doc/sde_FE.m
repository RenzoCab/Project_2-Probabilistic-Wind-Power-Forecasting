function sim_path_f = sde_FE(sim_path_i,alpha,theta_0,Theta_t,dt,P,P_dot) % 02/02/2020 18:28

    % SDE for X_t.

    sim_path_f = sim_path_i + (P_dot-Theta_t*(sim_path_i-P))*dt +...
                sqrt(dt)*sqrt(2*theta_0*alpha*sim_path_i*(1-sim_path_i))*randn(1);
    
end