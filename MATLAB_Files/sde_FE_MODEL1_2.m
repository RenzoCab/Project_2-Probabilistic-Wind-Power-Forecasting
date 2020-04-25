function sim_path_f = sde_FE_MODEL1_2(sim_path_i,alpha,theta_0,dt,P) % 01/04/2020 19:53

    % SDE for X_t^2 in MODEL 1.

    sim_path_f = sim_path_i + 2*theta_0*(sqrt(sim_path_i)*(alpha+P) - ...
        sim_path_i*(1+alpha)) * dt + 2*sqrt(sim_path_i)*sqrt(2*alpha*theta_0 * ...
        sqrt(sim_path_i)*(1-sqrt(sim_path_i))) * sqrt(dt) * randn(1);
    
end