function [xi1,xi2] = moments_matching_MODEL1(m1,m2) 

    % 19/04/2020 15:10

    mu   = m1;
    sig2 = m2 - m1^2;
    xi1  = (-sig2*mu + mu^2 - mu^3) / (sig2);
    xi2  = (mu-1) * (sig2 - mu + mu^2) / (sig2);
    
end