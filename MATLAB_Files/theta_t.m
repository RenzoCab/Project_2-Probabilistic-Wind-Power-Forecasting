function [theta_t] = theta_t(theta_0, alpha, pe, pe_dot) % 13/03/2020 19:09

    theta_t = max(theta_0, (theta_0*alpha + 2*abs(pe_dot))/(2*min(pe,1-pe)));

end