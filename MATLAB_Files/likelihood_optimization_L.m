function [val] = likelihood_optimization_L(batch, theta0, alpha, dt)

    % 29/03/2020 11:20
    
    disp(['Theta_0 = ',num2str(theta0),', Alpha = ',num2str(alpha)]);
    
%     batch = batch_with_theta_L(batch, alpha, theta0);
    val   = log_LH_evaluation_L(batch, theta0, alpha, dt);

end