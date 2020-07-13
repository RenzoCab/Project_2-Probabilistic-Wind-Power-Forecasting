function [val] = mean_regression_M1_fe(x1,theta_0,dt,p1,p2)
    
    n   = 50;
    p   = @(i) p1 + (p2-p1) * i/n;
    val = x1;
    
    for i = 2:n
        val = val - theta_0*(val - p(i))*dt;
    end

end
