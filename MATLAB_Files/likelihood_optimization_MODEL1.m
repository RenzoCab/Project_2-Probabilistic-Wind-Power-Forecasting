function [val] = likelihood_optimization_MODEL1(batch, theta_0, alpha, dt, approx)
    % 19/04/2020 22:08
    
    disp(['Theta_0 = ',num2str(theta_0),', Alpha = ',num2str(alpha),', Product = ',num2str(alpha*theta_0)]);
    batch_complete = batch_with_theta(batch, alpha, theta_0);
    
    val = log_LH_evaluation_MODEL1(batch_complete, alpha, theta_0, dt, approx);

end