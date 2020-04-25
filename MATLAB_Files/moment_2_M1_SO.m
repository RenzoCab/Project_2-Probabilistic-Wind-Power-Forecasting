function m2 = moment_2_M1_SO(x,theta_0,alpha,dt) % 21/04/2020 14:27
    
    L = -theta_0;
    g = sqrt(2*theta_0*alpha*x*(1-x));
    
    m2 = g * sqrt((exp(2*L*dt)-1) / (2*L));
    
end