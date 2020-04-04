function [fval] = likelihood_nested_optimization_L(batch, x0, dt)

    disp('================================================');
    disp(['Initial values: Theta_0 = ',num2str(x0(1)),' and alpha = ',num2str(x0(2))]);
    disp('================================================');

    batch   = batch_with_theta_L(batch, x0(2), x0(1));
    fun     = @(x) -likelihood_optimization_L(batch, x(1), x(2), dt);

    [x0,fval] = fminsearch(fun, x0);

end