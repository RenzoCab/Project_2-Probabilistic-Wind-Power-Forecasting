function drift = sde_Lamperti_drift_cand1(sim_path_i,alpha,theta_0,Theta_t,P,P_dot)

    drift = ((alpha*theta_0-Theta_t)*sin(sqrt(2*alpha*theta_0)*sim_path_i) - ...
                Theta_t*(1-2*P)+2*P_dot) / (sqrt(2*alpha*theta_0)*cos(sqrt(2*alpha*theta_0)*sim_path_i));
    
end