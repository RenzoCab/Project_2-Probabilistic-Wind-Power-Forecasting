function [batch_theta] = batch_with_theta(batch, alpha, theta_0) % 09/02/2020 18:51
    batch(4,1)   = theta_t(theta_0, alpha, batch(1,1), batch(2,1));
    batch(4,end) = theta_t(theta_0, alpha, batch(1,end), batch(2,end));
    for i = 2:2:length(batch(1,:))-1
        batch(4,i:i+1) = theta_t(theta_0, alpha, batch(1,i), batch(2,i));
    end
    batch_theta = batch;
end