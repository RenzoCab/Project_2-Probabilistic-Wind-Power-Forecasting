function m2 = moment_2(v,theta_0,th1,th2,p1,p2,alpha,m1,dt,n) % 02/04/2020 15:13

    p      = @(i) p1 + (p2-p1) * i/n; % We interpolate p(t).
    pe_dot = (p2-p1) / dt;
    m2(1)  = v^2;
    ds     = dt/n;
    
    if th1 == th2 % In this case, m_1 has exact solution and theta_t = theta_0.
        theta = @(i) th2;
        m1    = @(i) v*exp(-th1*dt*(i/n));
    else
        theta  = @(i) theta_t(theta_0, alpha, p(i), pe_dot);
    end
    
    for i = 2:n
        m2(i) = m2(i-1) + (-2*m2(i-1)*(theta(i-1)+alpha*theta_0) + ...
            2*alpha*theta_0*m1(i-1)*(1-2*p(i-1)) + ...
            2*alpha*theta_0*p(i-1)*(1-p(i-1))) * ds;
    end
    
end