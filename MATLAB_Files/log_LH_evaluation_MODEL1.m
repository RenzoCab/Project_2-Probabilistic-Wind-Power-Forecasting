function [value] = log_LH_evaluation_MODEL1(batch_complete, alpha, theta_0, dt, approx)

    % 19/04/2020 22:09
    
    % This makes (t_0,alpha) to be (t_0,gamma).
%     theta_0 = theta_0;
%     alpha   = alpha / theta_0;

    for i = 1:length(batch_complete(1,:))/2

        j = i*2; % This is the real index (parfor must go one-by-one).
        % Recall that: j is t_n and j-1 is t_{n-1}.
        p1 = batch_complete(1,j-1); p2 = batch_complete(1,j);
        v1 = batch_complete(3,j-1); v2 = batch_complete(3,j);
        n  = 100;
        x1 = v1 + p1;
        x2 = v2 + p2;

        if approx == 1
        
            m1 = moment_1_M1(x1,theta_0,p1,p2,dt,n);
            m2 = moment_2_M1(x1,m1,theta_0,alpha,p1,p2,dt,n);
            val(i) = log_dist_L(x2,m1(end),sqrt(m2(end)-m1(end)^2));
        
        elseif approx == 2

            m1 = moment_1_M1_SO(x1,theta_0,p1,p2,dt);
            m2 = moment_2_M1_SO(x1,theta_0,alpha,dt);
            A = m1; B = m2;
            val(i) = log_dist_L(x2,A*x1,B);
        
        elseif approx == 3

            m1 = moment_1_M1(x1,theta_0,p1,p2,dt,n);
            m2 = moment_2_M1(x1,m1,theta_0,alpha,p1,p2,dt,n);
            [xi1,xi2] = moments_matching_MODEL1(m1(end),m2(end));
            val(i) = log_dist_extra(x2,xi1,xi2);
        
        else
            
            error('Choose a number between 1, 2 or 3.');
            
        end
        
    end
    
    value = sum(val);
    
end

