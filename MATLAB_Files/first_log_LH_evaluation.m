function [value,m2,xi1,sig2] = first_log_LH_evaluation(batch_complete, theta_0, alpha, delta, dt, N, norm, times) 

    % 09/02/2020 18:35

    for i = 0:length(batch_complete(1,:))/((N-1)*2)-1
    
        p1  = batch_complete(1,1 + i*(N-1)*2);
        p2  = batch_complete(1,12 + i*(N-1)*2); % This is the measurement at 01:10.
        p0  = p1 - delta*(p2-p1)/(6*dt);
        p0  = min(1-0.018,p0); p0 = max(0.018,p0);
        n   = 50;
        th1 = batch_complete(4,1  + i*(N-1)*2);
        th0 = theta_t(theta_0, alpha, p0, (p1-p0)/(delta));
        v0  = 0; v1 = batch_complete(3,1  + i*(N-1)*2);

        v_ini(i+1) = v1;
        % (v,theta_0,alpha,th1,th2,p1,p2,dt,n) 
        m1   = moment_1(v0,theta_0,alpha,th0,th1,p0,p1,delta,n);
        % (v,theta_0,th1,th2,p1,p2,alpha,m1,dt,n)
        m2   = moment_2(v0,theta_0,th0,th1,p0,p1,alpha,m1,delta,n);
        mom1(i+1)  = m1(end);
        mom2(i+1)  = m2(end);
        
        [xi1(i+1),xi2(i+1)] = moments_matching(m1(end),m2(end));
        sig2(i+1) = 1/(4*(xi1(i+1)*2+1));
        val(i+1)   = -log_dist(v1,xi1(i+1),xi2(i+1));
    
    end

    if delta == 1/norm
        figure
        histogram((v_ini+1)/2,100,'Normalization','pdf')
        hold on
    end
    if delta == times
        X = 0:.01:1;
        for i = 1:length(xi1)
            y1 = betapdf(X,xi1(i),xi1(i));
            plot(X,y1);
        end
        histogram((v_ini+1)/2,100,'Normalization','pdf')
    end
    
    value = sum(val);
    m2    = mean(mom2);
    xi1   = mean(xi1);
    sig2  = mean(sig2);
                
end