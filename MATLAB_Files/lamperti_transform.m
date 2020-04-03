function [val] = lamperti_transform(theta_0,alpha,v,p,candidate) 

    % 29/03/2020 17:49

    if candidate == 1
        val = (1/sqrt(2*alpha*theta_0)) * asin(2*(v+p)-1);
    elseif candidate == 2
        val = -sqrt(2/(alpha*theta_0)) * asin(sqrt(1-v-p));
    end
    
end