function m1 = moment_1(v,theta_0,alpha,th1,th2,p1,p2,dt,n) 
    % 02/04/2020 15:16

    if th1 == th2 % We have the exact solution.
        m1 = v*exp(-th1*dt);
    else % Otherwise, we compute F-E.
        p      = @(i) p1 + (p2-p1) * i/n; % We interpolate p(t).
        pe_dot = (p2-p1) / dt;
        theta  = @(i) theta_t(theta_0, alpha, p(i), pe_dot);
        m1(1) = v;
        
        ds    = dt/n;
        for i = 2:n
            m1(i) = m1(i-1) * (1 - theta(i-1)*ds);
        end
    end
    
end