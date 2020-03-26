function [val] = likelihood_optimization(batch, theta0, alpha, dt)
    % 10/02/2020 11:20
    disp(['Theta = ',num2str(theta0),', Alpha = ',num2str(alpha)]);
    batch_complete = batch_with_theta(batch, alpha, theta0);
    val            = log_LH_evaluation(batch_complete, alpha, dt);

end