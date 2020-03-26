function m1 = moment_1(v,th1,th2,dt,n) % 02/02/2020 18:28

    if th1 == th2 % We have the exact solution.
        m1 = v*exp(-th1*dt);
    else % Otherwise, we compute F-E.
        m1(1) = v;
        theta = @(i) th1 + (th2-th1) * i/n;
        ds    = dt/n;
        for i = 2:n
            m1(i) = m1(i-1) * (1 - theta(i-1)*ds);
        end
    end
    
end