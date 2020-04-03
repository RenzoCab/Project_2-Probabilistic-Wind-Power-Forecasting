function m1 = moment_1_L(z,theta_0,alpha,p1,p2,dt,n)
    % 29/03/2020 16:28
    p_dot   = (p2-p1) / dt;
    m1(1)   = z;
    p_t     = @(i) p1 + (p2-p1) * i/n;
    ds      = dt/n;
    
    for i = 2:n
        Theta_t = theta_t(theta_0, alpha, p_t(i), p_dot);
        m1(i)   = m1(i-1) + ds * ...
            ((alpha*theta_0-Theta_t) * sin(m1(i-1)*sqrt(2*alpha*theta_0)) - ...
            Theta_t*(1-2*p_t(i)) + 2*p_dot) / ...
            (sqrt(2*alpha*theta_0) * cos(m1(i-1)*sqrt(2*alpha*theta_0)));
    end
    
end