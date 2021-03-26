function [x_dot, x_double_dot] = get_derivatives(x, h)
%Uses simple gradient method to find derivatives. 

x_dot = gradient(x, h);
x_double_dot = gradient(x_dot, h);

end

