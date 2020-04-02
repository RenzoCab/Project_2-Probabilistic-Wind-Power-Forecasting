function sim_path_f = sde_Lamperti_drift(sim_path_i,alpha,theta_0,Theta_t,P,P_dot)

    sim_path_f = ((alpha*theta_0-Theta_t)*sin(sqrt(2*alpha*theta_0)*sim_path_i) - ...
                Theta_t*(1-2*P)+2*P_dot) / (sqrt(2*alpha*theta_0)*cos(sqrt(2*alpha*theta_0)*sim_path_i));
            
%     sim_path_f = ((alpha*theta_0-Theta_t)*cos(sqrt(2*alpha*theta_0)*sim_path_i) - ...
%         Theta_t*(1-2*P)+2*P_dot) / (sqrt(alpha*theta_0)*(1-cos(2*sqrt(2*alpha*theta_0)*sim_path_i)));
    
end