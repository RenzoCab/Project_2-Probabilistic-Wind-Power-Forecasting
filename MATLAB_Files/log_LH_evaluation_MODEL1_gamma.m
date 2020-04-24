function [value] = log_LH_evaluation_MODEL1_gamma(batch_complete, gamma, theta_0, dt)

    % 19/04/2020 22:09

    for i = 1:length(batch_complete(1,:))/2

        j = i*2; % This is the real index (parfor must go one-by-one).
        % Recall that: j is t_n and j-1 is t_{n-1}.
        p1 = batch_complete(1,j-1); p2 = batch_complete(1,j);
        v1 = batch_complete(3,j-1); v2 = batch_complete(3,j);
        n  = 100;
        x1 = v1 + p1;
        x2 = v2 + p2;

        m1 = moment_1_M1(x1,theta_0,p1,p2,dt,n);
        m2 = moment_2_M1_gamma(x1,m1,theta_0,gamma,p1,p2,dt,n);
                
        val(i) = log_dist_L(x2,m1(end),sqrt(m2(end)-m1(end)^2));
        
    end
    
    value = sum(val);
    
end

