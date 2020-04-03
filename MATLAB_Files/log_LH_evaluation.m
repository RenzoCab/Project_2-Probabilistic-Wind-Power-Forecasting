function [value] = log_LH_evaluation(batch_complete, alpha, theta_0, dt) % 03/02/2020 19:42

    for i = 1:length(batch_complete(1,:))/2

        j = i*2; % This is the real index (parfor must go one-by-one).
        % Recall that: j is t_n and j-1 is t_{n-1}.
        p1 = batch_complete(1,j-1); p2 = batch_complete(1,j);
        v1 = batch_complete(3,j-1); v2 = batch_complete(3,j);
        th1 = batch_complete(4,j-1); th2 = batch_complete(4,j);
        n = 100; % 50 discretizations for the ODEs.

        m1        = moment_1(v1,theta_0,alpha,th1,th2,p1,p2,dt,n);
        m2        = moment_2(v1,theta_0,th1,th2,p1,p2,alpha,m1,dt,n);
        [xi1,xi2] = moments_matching(m1(end),m2(end));
        if xi1<0 || xi2<0
            val(i)    = 0;
        else
            val(i)    = log_dist(v2,xi1,xi2);
        end
        
    end
    
    value = sum(val);
    
end
