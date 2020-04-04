function [batch_theta] = batch_with_theta_L(batch, alpha, theta_0) % 29/03/2020 19:43

    batch(4,1)   = theta_t(theta_0, alpha, batch(1,1), batch(2,1));
    batch(4,end) = theta_t(theta_0, alpha, batch(1,end), batch(2,end));
    batch(5,1)   = lamperti_transform(theta_0,alpha,batch(3,1),batch(1,1),2);
    batch(5,end) = lamperti_transform(theta_0,alpha,batch(3,end),batch(1,end),2);
    
    for i = 2:2:length(batch(1,:))-1
        batch(4,i:i+1) = theta_t(theta_0, alpha, batch(1,i), batch(2,i));
        batch(5,i:i+1) = lamperti_transform(theta_0, alpha, batch(3,i), batch(1,i),2);
    end
    
    batch_theta = batch;
    
end