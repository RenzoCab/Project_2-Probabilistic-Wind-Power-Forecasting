function [value] = log_LH_evaluation_L(batch_complete, theta0, alpha, dt) % 03/02/2020 19:42

    for i = 1:length(batch_complete(1,:))/2

        j = i*2; % This is the real index (parfor must go one-by-one).
        % Recall that: j is t_n and j-1 is t_{n-1}.
        p1 = batch_complete(1,j-1); p2 = batch_complete(1,j);
        z1 = batch_complete(5,j-1); z2 = batch_complete(5,j);
        n = 10; % 10 discretizations for the ODEs.
        
        mu        = moment_1_L(z1,theta0,alpha,p1,p2,dt,n);
        sig       = moment_2_L(mu,theta0,alpha,p1,p2,dt,n);
        val(i)    = log_dist_L(z2,mu(end),sig(end));
        
    end
    
    value = sum(val);
    
end
