function drift_prime = sde_Lamperti_drift_prime_cand2(z,alpha,theta_0,theta_t,p,p_dot)

    drift_prime = (alpha*theta_0 - theta_t + cos(sqrt(2*alpha*theta_0)*z) * ...
        (2*p_dot + theta_t*(2*p - 1))) / (sin(sqrt(2*theta_0*alpha)*z))^2;
    
end