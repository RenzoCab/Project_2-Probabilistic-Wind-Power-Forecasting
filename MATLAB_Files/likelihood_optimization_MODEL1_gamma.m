function [val] = likelihood_optimization_MODEL1_gamma(batch, theta_0, gamma, dt)
    % 19/04/2020 22:08
    
    disp(['Theta_0 = ',num2str(theta_0),', Gamma = ',num2str(gamma)]);
    batch_complete = batch_with_theta(batch, gamma/theta_0, theta_0);
    
    val = log_LH_evaluation_MODEL1_gamma(batch_complete, gamma, theta_0, dt);

end