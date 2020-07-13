function m2 = moment_2_L(m1,theta_0,alpha,p1,p2,dt,n)
    % 29/03/2020 16:50

    % This function returns the approximated variance.

    p_dot = (p2-p1) / dt;
    m2(1) = 0;
    p_t   = @(i) p1 + (p2-p1) * i/n;
    ds    = dt/n;
    
    for i = 2:n
        Theta_t = theta_t(theta_0, alpha, p_t(i), p_dot);
        b_prime = sde_Lamperti_drift_prime_cand2(m1(i-1),alpha,theta_0,Theta_t,p_t(i-1),p_dot);
        m2(i)   = m2(i-1) + ds * ...
            (1+2*m2(i-1)*b_prime);
    end
    
end