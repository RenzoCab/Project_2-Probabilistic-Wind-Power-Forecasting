function m2 = moment_2_M1_gamma(x,m1,theta_0,gamma,p1,p2,dt,n) % 19/04/2020 22:00
    p     = @(i) p1 + (p2-p1) * i/n; % We interpolate p(t).
    m2(1) = x^2;
    ds    = dt/n;
    for i = 2:n
        m2(i) = m2(i-1) + theta_0*(-2*m2(i-1)+2*p(i-1)*m1(i-1))*ds + ...
            gamma*(-2*m2(i-1) + 2*m1(i-1))*ds;
    end

end