function m2 = moment_2_ini(theta,p0,p1,alpha,delta) % 18/03/2020 17:21

    a  = p0;
    b  = (p1-p0)/delta;
    k1 = 2*theta*(1+alpha);
    k2 = 2*theta*alpha;
    
    t  = delta;
    
    m2 = (-exp(-k1*t)*k2/k1^3) * ((a-1)*a*(exp(k1*t)-1)*k1^2 + (2*a-1)*b*k1*(1+(k1*t-1)*exp(k1*t))...
        + b^2*(-2+exp(k1*t)*(2-2*k1*t+k1^2*t^2)));
    
end