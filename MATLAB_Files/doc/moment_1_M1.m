function m1 = moment_1_M1(x,theta_0,p1,p2,dt,n) % 19/04/2020 21:56
    p     = @(i) p1 + (p2-p1) * i/n; % We interpolate p(t).
    m1(1) = x;
    ds    = dt/n;
    for i = 2:n
        m1(i) = m1(i-1) - theta_0*(m1(i-1)-p(i-1))*ds;
    end
end