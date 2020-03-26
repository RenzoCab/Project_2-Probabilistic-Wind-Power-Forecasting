function [val] = log_dist(v,xi1,xi2) % 03/02/2020 11:20
    a   = -1;
    b   =  1;
    val = -betaln(xi1,xi2) + (xi1-1)*log((v-a)/(b-a)) + ...
        (xi2-1)*log((b-v)/(b-a));
end