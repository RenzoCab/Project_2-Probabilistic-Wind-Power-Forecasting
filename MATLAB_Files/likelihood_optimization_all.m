function [val] = likelihood_optimization_all(batch, theta_0, alpha, delta, dt, N)
    % 10/02/2020 11:20
    disp(['Theta = ',num2str(theta_0),', Alpha = ',num2str(alpha),', Delta = ',num2str(delta)]);
    batch_complete = batch_with_theta(batch, alpha, theta_0);
    if theta_0 < 0 || alpha < 0 || delta < 0
        val = - Inf;
    else
        val = log_LH_evaluation(batch_complete, alpha, theta_0, dt) - ...
            first_log_LH_evaluation(batch_complete, theta_0, alpha, delta, dt, N, 1e6, 1e6);
    end

end