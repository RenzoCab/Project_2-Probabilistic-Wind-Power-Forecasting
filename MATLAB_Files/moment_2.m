function m2 = moment_2(v,th1,th2,p1,p2,alpha,m1,dt,n) % 02/02/2020 18:28

    if th1 == th2
        theta = @(i) th2;
        m1    = @(i) v*exp(-th1*dt*(i/n));
    else
        theta = @(i) th1 + (th2-th1) * i/n;
    end
    p     = @(i) p1  + (p2-p1)   * i/n;
    m2(1) = v^2;
    ds    = dt/n;
    
    for i = 2:n
        m2(i) = m2(i-1) + 2*theta(i-1)*ds * (-(1+alpha)*m2(i-1) + ...
            alpha*m1(i-1)*(1-2*p(i-1)) + alpha*p(i-1)*(1-p(i-1))); 
    end
    
end