function m1 = moment_1_M1_SO(x,theta_0,p1,p2,dt) % 21/04/2020 14:27

    L = -theta_0;
    M = theta_0 * (p2-p1) / dt;
    f = -theta_0 * (x - p1);
    
    m1 = 1 + f/(x*L) * (exp(L*dt)-1) + M/(x*L^2) * (exp(L*dt) - 1 - L*dt);

end