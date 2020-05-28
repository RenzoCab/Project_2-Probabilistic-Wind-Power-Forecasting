function [fval] = fixed_point_function(batch, x0, dt)

    [x1,~] = likelihood_nested_optimization_L(batch, x0, dt);

    fval = abs(x1(1)-x0(1))/x0(1) + abs(x1(2)-x0(2))/x0(2);
%     fval = abs(x1(1)*x1(2)-x0(1)*x0(2)) / (x0(1)*x0(2));

end