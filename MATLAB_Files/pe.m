function [pe] = pe(p, eps) % 18/03/2020 09:51
    for i = 1:length(p)    
        if p(i) < eps
                pe(i) = eps;
            elseif p(i) >= eps && p(i) < 1-eps
                pe(i) = p(i);
            elseif p(i) >= 1-eps
                pe(i) = 1-eps;
        end
    end
end