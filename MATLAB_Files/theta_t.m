function [theta_t] = theta_t(theta_0, alpha, pe, pe_dot) % 13/03/2020 19:09

    theta_t = max(theta_0, (theta_0*alpha + abs(pe_dot))/(min(pe,1-pe)));

end