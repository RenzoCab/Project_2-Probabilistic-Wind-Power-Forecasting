function [val] = lamperti_transform(theta_0,alpha,v,p) 

    % 29/03/2020 17:49

    val = (1/sqrt(2*alpha*theta_0)) * asin(2*(v+p)-1);
    
end