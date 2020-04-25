function [val] = mean_regression_M1_eq(x1,theta_0,dt,m,n)

    val = x1*exp(-theta_0*dt) + theta_0*exp(-theta_0*dt) * ...
        (m + (exp(theta_0*dt)-1)*n*theta_0 + exp(theta_0*dt)*m*(theta_0*dt-1)) / ...
        (theta_0^2);

end
