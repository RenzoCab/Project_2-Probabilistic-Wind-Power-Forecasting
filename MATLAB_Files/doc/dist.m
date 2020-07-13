function [val] = dist(v,xi1,xi2) % 03/02/2020 11:20
    a   = -1;
    b   =  1;
    val = (1/abs(b-a)) * (1/beta(xi1,xi2)) * ((v-a)/(b-a))^(xi1-1)...
        * (1-(v-a)/(b-a))^(xi2-1);
end