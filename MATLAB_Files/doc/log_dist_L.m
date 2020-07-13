function [val] = log_dist_L(z,mu,sig) % 29/03/2020 17:49

    val = -(1/2)*((z-mu)/sig)^2 - log(sig*sqrt(2*pi));
    
end