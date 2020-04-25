function [val] = mean_regression_M1(batch_complete, dt, rang, gamma) % 16/03/2020 11:08

    num_of_samples = length(batch_complete(1,:))/2;

    theta_0   = rang;
    max_error = 0;
    
    for k = 1:length(theta_0)
        
        value     = [];
        value_2   = [];
        index     = [];
        count_1   = 0;
        count_2   = 0;
        val_1     = 0;
        val_2     = 0;
        extra_div = 1;
        dt        = dt * extra_div;
        theta_pos = 0;
        theta_neg = 0;
        
        for i = 1:num_of_samples/extra_div

            j = i*2*extra_div; % This is the real index (parfor must go one-by-one).
            % Recall that: j is t_n and j-1 is t_{n-1}.
            p1 = batch_complete(1,j-1); p2 = batch_complete(1,j);
            v1 = batch_complete(3,j-1); v2 = batch_complete(3,j);
            x1 = v1 + p1; x2 = v2 + p2;

            n = p1;
            m = (p2 - p1) / dt;

            value(i)   = mean_regression_M1_eq(x1,theta_0(k),dt,m,n);
            value_2(i) = mean_regression_M1_eq(x1,theta_0(k),dt,p1,p2);
            aux_1(i)   = value(i);
            aux_2(i)   = value_2(i);
            value(i)   = (value(i) - x2)^2;
            aux_3(i)   = x1;
            aux_4(i)   = x2;
            aux_5(i)   = p1;
            aux_6(i)   = p2;
            
            if not(p1 > gamma && p1 < 1-gamma)
                index(end+1) = i;
            end
            
            if abs(x2-p2) < abs(x1-p1)
                count_2 = count_2 + 1;
                val_2   = val_2 + abs(x2-p2);
            else
                count_1 = count_1 + 1;
                val_1   = val_1 + abs(x1-p1);
            end
            
            if x1 < p1 && x2 < p2 && x2 < x1
                theta_neg = theta_neg + 1;
            elseif x1 < p1 && x2 < p2 && x2 > x1
                theta_pos = theta_pos + 1;
            elseif x1 > p1 && x2 > p2 && x2 < x1
                theta_pos = theta_pos + 1;
            elseif x1 > p1 && x2 > p2 && x2 > x1
                theta_neg = theta_neg + 1;
            end
            % This shows that it is more likely to obtain a negative Theta_0.
            
        end
        
        for i = length(index):-1:1
            value(index(i)) = [];
        end
        
        % Error Forward-Euler and Exact solution:
        max_error = max(max_error,max(abs(aux_1-aux_2)));
        
        val(k) = mean(value);

    end

end
