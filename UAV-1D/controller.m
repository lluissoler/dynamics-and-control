function [u] = controller(~, x, x_des, params)

%   x: 2x1 vector containing the current state [z; z_dot]
%   x_des: 2x1 vector containing desired state
%   using a PD controller


kp = 10; % proportional gain
kd = 9; % derivative gain
z_dot_dot = 0; % desired acceleration
e = x_des-x; % error

u = params.mass*(z_dot_dot + kp*e(1) + kd*e(2) + params.gravity);

end

