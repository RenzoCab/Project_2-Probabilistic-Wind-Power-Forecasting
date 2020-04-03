function m1 = moment_1_L(z,theta_0,alpha,p1,p2,dt,n)
    % 29/03/2020 16:28
    p_dot   = (p2-p1) / dt;
    m1(1)   = z;
    p_t     = @(i) p1 + (p2-p1) * i/n;
    ds      = dt/n;
    
    for i = 2:n
        Theta_t = theta_t(theta_0, alpha, p_t(i), p_dot);
        m1(i)   = m1(i-1) + ds * ...
            sde_Lamperti_drift_cand2(m1(i-1),alpha,theta_0,Theta_t,p_t(i-1),p_dot);
    end
    
end